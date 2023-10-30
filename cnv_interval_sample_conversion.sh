#!/bin/bash

input_file="interval_sample_results.txt"
output_file="output.txt"

# Use awk to process the input file and convert columns 3 and 4 to integers
awk -F'\t' -v OFS='\t' 'NR==1 {print; next} { $3=int($3); $4=int($4); print }' "$input_file" > "$output_file"

rm -rf "$input_file"
mv "$output_file" "$input_file"

echo "Conversion complete. Results are in $output_file"
