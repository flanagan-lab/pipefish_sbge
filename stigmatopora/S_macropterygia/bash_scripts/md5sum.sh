#!/bin/bash

# Loop through all files in the current directory
for file in *; do
    # Only process regular files (skip directories)
    if [ -f "$file" ]; then
        # Calculate and display MD5 sum
        md5sum "$file"
    fi
done
