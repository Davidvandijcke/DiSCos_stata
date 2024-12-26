#!/bin/bash

# Create or overwrite the output file
output_file="combined_output.txt"
touch "$output_file"

# Loop through all files in the current directory
for file in *; do
    # Skip the output file itself if it exists
    if [ "$file" != "$output_file" ]; then
        echo "=== Content of $file ===" >> "$output_file"
        
        # Check if file is binary using 'file' command
        if file "$file" | grep -q "text"; then
            # If it's a text file, directly cat it
            cat "$file" >> "$output_file"
        else
            # For binary files, try to extract text using strings
            strings "$file" >> "$output_file"
        fi
        
        echo -e "\n\n" >> "$output_file"
    fi
done

echo "All content has been extracted to $output_file"