#!/bin/bash

# Set variables
assembly_file=$1
output_file=$2

# Function to calculate average transcript length
calculate_avg_length() {
    # Extract first sequence (assuming it's representative)
    first_seq=$(head -n1 $1 | sed 's/>//')
    
    # Calculate length
    avg_length=$(echo "scale=2; $(wc -c < $assembly_file) / $(echo "$first_seq" | wc -m)" | bc)
    
    echo "$avg_length"
}

# Function to calculate total assembly size
calculate_total_size() {
    num_sequences=$(wc -l < $1)
    avg_length=$(calculate_avg_length)
    total_size=$(echo "scale=2; $num_sequences * $avg_length" | bc)
    
    echo "$total_size"
}

# Main script
echo "Assembly Size Estimation Report" > $2
echo "" >> $2

echo "Number of sequences: $(wc -l < $assembly_file)" >> $2
echo "Average transcript length: $(calculate_avg_length)" >> $2
echo "Total assembly size (MB): $(($(calculate_total_size) / 1000000))" >> $2

echo "" >> $2
echo "Estimated Assembly Size Breakdown:" >> $2
echo "-----------------------------------" >> $2
echo "Sequences: $(wc -l < $assembly_file)" >> $2
echo "Size per sequence (bases): $(echo "$(calculate_total_size) / $(wc -l < $1)" | bc)" >> $2
echo "Total bases: $(calculate_total_size)" >> $2
echo "Total MB: $(($(calculate_total_size) / 1000000))" >> $2

echo "" >> $2
echo "Note: This is an approximation based on the first sequence in the file." >> $2
