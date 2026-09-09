#!/bin/bash

# Define the target Nextflow work directory path
TARGET_DIR="work/ArcSimSyn_offtarget_unique"
PATTERN="${1:-*_dge.rds}"

# Check if the directory exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Nextflow work directory '$TARGET_DIR' does not exist."
    exit 1
fi

echo "Searching Nextflow work tree '$TARGET_DIR'..."
echo "Searching Nextflow work tree for subdirectories containing '$PATTERN'..."

# Find unique directories containing the pattern
# We use a temporary file to store the list securely
DIR_LIST=$(mktemp)

# Search deeply nested Nextflow subdirectories
find "$TARGET_DIR" -type f -name "$PATTERN" -exec dirname {} \; | sort -u > "$DIR_LIST"

# Check if any directories were found
if [ ! -s "$DIR_LIST" ]; then
    echo "No directories containing '$PATTERN' were found in $TARGET_DIR."
    rm "$DIR_LIST"
    exit 0
fi

# Count the number of unique directories found
DIR_COUNT=$(wc -l < "$DIR_LIST")
echo "Found $DIR_COUNT Nextflow task directory(ies) containing matching files."

# Ask for user confirmation
read -p "Are you sure you want to delete these $DIR_COUNT directories and all their contents? (y/N): " CONFIRM
if [[ "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Deleting directories..."
    while IFS= read -r dir; do
        if [ -d "$dir" ]; then
            rm -rf "$dir"
            echo "Deleted: $dir"
        fi
    done < "$DIR_LIST"
    echo "Cleanup complete."
else
    echo "Operation canceled. No directories were deleted."
fi

# Clean up temporary file
rm "$DIR_LIST"
