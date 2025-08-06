#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p results
mkdir -p results-err

# Loop through all files in the 'files' directory
for file in ../../GILS-RVND-TSP/leitor-instancias/BB-instances/*; do
    filename=$(basename "$file")           # Extract filename without path
    ./solve.out "$file" > "results/${filename}.txt" 2> "results-err/${filename}.txt" # Run command and redirect output
done
