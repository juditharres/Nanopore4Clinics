#!/bin/bash
set -euo pipefail

# Function to display help
usage() {
    echo "Usage: $0 -s SAMPLE_ID -R RUN_ID -i IN_DIR -o OUT_DIR -t THREADS -r REFERENCE -m CLAIR3_MODEL -b BED_FILE"
    echo "Options:"
    echo "  -s SAMPLE_ID     Sample ID"
    echo "  -R RUN_ID        Run ID"
    echo "  -i IN_DIR        Input directory with bam_pass basecalled data"
    echo "  -o OUT_DIR       Output directory"
    echo "  -t THREADS       Number of threads"
    echo "  -r REFERENCE     Reference genome file"
    echo "  -m CLAIR3_MODEL  Clair3 model file"
    echo "  -b BED_FILE      BED file for STR calling"
    echo "  -h               Display this help message"
    exit 1
}

# Parse command-line arguments
while getopts "s:i:o:t:r:m:R:b:h" opt; do
    case ${opt} in
        s ) SAMPLE_ID="${OPTARG}" ;;
        i ) IN_DIR="${OPTARG}" ;;
        o ) OUT_DIR="${OPTARG}" ;;
        t ) THREADS="${OPTARG}" ;;
        r ) REFERENCE="${OPTARG}" ;;
        m ) CLAIR3_MODEL="${OPTARG}" ;;
        R ) RUN_ID="${OPTARG}" ;;
        b ) BED_FILE="${OPTARG}" ;;
        h ) usage ;;
        * ) usage ;;
    esac
done

# Check if all required arguments are provided
if [ -z "${SAMPLE_ID:-}" ] || [ -z "${IN_DIR:-}" ] || [ -z "${OUT_DIR:-}" ] || \
   [ -z "${THREADS:-}" ] || [ -z "${REFERENCE:-}" ] || [ -z "${CLAIR3_MODEL:-}" ] || \
   [ -z "${RUN_ID:-}" ] || [ -z "${BED_FILE:-}" ]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

# Find UBAM pass files
bam_files=$(find "${IN_DIR}" -name "*.bam")

if [ -z "$bam_files" ]; then
    echo "Error: No BAM files found in ${IN_DIR}. Exiting script."
    exit 1
fi

# If BAM files are found, proceed with the script
echo "BAM files found. Continuing with the analysis."
echo "$bam_files" > "${OUT_DIR}/result.txt"

# Concatenate multiple uBAM files into one uBAM for alignment
samtools cat -o "${OUT_DIR}/${SAMPLE_ID}.bam" --threads "${THREADS}" -b "${OUT_DIR}/result.txt"

# Check if the BAM file has a good EOF block
if ! samtools quickcheck -vvv "${OUT_DIR}/${SAMPLE_ID}.bam" 2>&1 | grep -q "good EOF block"; then
    echo "Error: BAM file ${OUT_DIR}/${SAMPLE_ID}.bam failed samtools quickcheck (EOF block missing or corrupted)."
    exit 1
fi

echo "BAM file ${OUT_DIR}/${SAMPLE_ID}.bam passed samtools quickcheck with a good EOF block."

# Alignment
samtools fastq -T Mm,Ml,MM,ML,mm,ml "${OUT_DIR}/${SAMPLE_ID}.bam" --reference "${REFERENCE}" --threads "${THREADS}" | \
    sentieon minimap2 -a -y -x map-ont -t "${THREADS}" --secondary no \
    -R "@RG\tID:${RUN_ID}\tSM:${SAMPLE_ID}\tPL:ONT" "${REFERENCE}.mmi" - | \
    samtools sort --reference "${REFERENCE}" --threads "${THREADS}" -O CRAM -o "${OUT_DIR}/${SAMPLE_ID}.sorted.cram"

# Alignment QC
alfred qc -r "${REFERENCE}" -j "${OUT_DIR}/sample_alfred.json.gz" -o "${OUT_DIR}/sample_alfred.tsv.gz" "${OUT_DIR}/${SAMPLE_ID}.sorted.cram"
mosdepth -b 1000 -n -t "${THREADS}" -f "${REFERENCE}" "${OUT_DIR}/mosdepth/" "${OUT_DIR}/${SAMPLE_ID}.sorted.cram"

# Variant Calling
run_clair3.sh -b "${OUT_DIR}/${SAMPLE_ID}.sorted.cram" \
              -f "${REFERENCE}" \
              --ctg_name="$(seq 1 22 | sed 's/^/chr/' | tr '\n' ',' | sed 's/,$//')chrX,chrY,chrM" \
              --gvcf --enable_phasing --longphase_for_phasing \
              -t "${THREADS}" -p ont --remove_intermediate_dir \
              -m "${CLAIR3_MODEL}" --sample_name="${SAMPLE_ID}" \
              -o "${OUT_DIR}/"

# Structural variant calling
sniffles --input "${OUT_DIR}/${SAMPLE_ID}.sorted.cram" \
         --allow-overwrite \
         --vcf "${OUT_DIR}/${SAMPLE_ID}_sniffles2.vcf.gz" \
         --reference "${REFERENCE}" \
         -t "${THREADS}" \
         --snf "${OUT_DIR}/${SAMPLE_ID}_sniffles2.snf"

# Structural variant calling QC
SURVIVOR stats "${OUT_DIR}/${SAMPLE_ID}_sniffles2.vcf.gz" -1 -1 -1 "${OUT_DIR}/${SAMPLE_ID}_survivor.tsv"

# CNV Calling
python spectre.py CNVCaller --bin-size 1000 --coverage "${OUT_DIR}/mosdepth/" \
                            --sample-id "${SAMPLE_ID}" --output-dir "${OUT_DIR}/spectre/" \
                            --reference "${REFERENCE}"

# Extract gender using goleft indexcov
goleft indexcov "${OUT_DIR}/${SAMPLE_ID}.sorted.cram.crai" --extranormalize \
                --directory "${OUT_DIR}/goleft_straglr/" --sex chrX,chrY --fai "${REFERENCE}.fai"

awk -F'\t' 'NR>1{if($5==2) print $2"\tfemale"; if($5==1) print $2"\tmale"}' \
    "${OUT_DIR}/goleft_straglr/"*.ped > "${OUT_DIR}/gender.txt"

# STR Calling
straglr-genotype "${OUT_DIR}/${SAMPLE_ID}.sorted.cram" "${REFERENCE}" \
                 --loci "${BED_FILE}" --sample "${SAMPLE_ID}" \
                 --tsv "${OUT_DIR}/goleft_straglr/${SAMPLE_ID}.straglr.tsv" \
                 --vcf "${OUT_DIR}/goleft_straglr/${SAMPLE_ID}.straglr.vcf" \
                 --sex "$(cat "${OUT_DIR}/gender.txt")" --debug

echo "Analysis completed successfully."