#!/bin/bash
#SBATCH --job-name=parent_finder
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=08:00:00
#SBATCH --array=1-2040
#SBATCH --output=logs/parent_finder_%A_%a.out
#SBATCH --error=logs/parent_finder_%A_%a.err

# NOTE: Submit this with: sbatch --array=1-N%100 
# where N = number of progeny files

# Create directories
mkdir -p logs
mkdir -p results

# Configuration
PROGENY_FOLDER="./analysis/preprocessing"
PARENT_FOLDER="./analysis/parent_classification"

# Get the progeny file for this array task
PROGENY_FILE=$(ls ${PROGENY_FOLDER}/*.txt | grep -E '[0-9]+_[0-9]+_[0-9]+\.txt' | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Skip if no file found
if [ -z "$PROGENY_FILE" ]; then
    echo "No progeny file found for task ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# Extract progeny ID
PROGENY_ID=$(basename "$PROGENY_FILE" .txt)

# Output file
OUTPUT_FILE="results/${PROGENY_ID}_top20.txt"

# Run the script
echo "Processing: $PROGENY_ID"
python ./code/intersecion.sh "$PROGENY_FILE" "$PARENT_FOLDER" "$OUTPUT_FILE"

echo "Done: $PROGENY_ID
