#!/bin/bash

# Function to process each alignment file
run_iqtree2() {
    alignment_file="$1"
    # Count the number of sequences in the alignment file
    num_sequences=$(grep -c "^>" "$alignment_file")
    # Only process the alignment file if it has more than 5 sequences
    if [ "$num_sequences" -gt 5 ]; then
        # Extract the HG_name from the alignment file name
        HG_name=$(basename "$alignment_file" | cut -d '.' -f 1)
        # Define the tree folder path
        tree_folder="../results/molecular_evolution_analyses/phylogenetic_trees/$HG_name"
        # Create the tree folder if it doesn't exist
        mkdir -p "$tree_folder"
        # Run iqtree2 command
        ~/mauricio_PROGRAMAS/iqtree-2.2.2.6-Linux/bin/iqtree2 -s "$alignment_file" -m TEST --threads-max 32 --alrt 1000 --ufboot 1000 --prefix "../results/molecular_evolution_analyses/phylogenetic_trees/$HG_name/$HG_name"
    else
        echo "Skipping alignment file $alignment_file: Number of sequences is less than or equal to 5"
    fi
}

# Function to export run_iqtree2 so it's available in parallel
export -f run_iqtree2

# Use find to get a list of alignment files and pipe them to parallel
find ../results/molecular_evolution_analyses/alignments_gblocks/ -name "*gblo" | parallel -j 32 run_iqtree2

