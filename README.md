# FindConservedSeq
An easy-to-use tool for finding conserved sequences of DNA stands of the same strain

## How to use:
Still a work in-progress
1) Ignore the main file
2) To use:
- Call parse_nucleotides function from the parse class and input your fasta file. This function will return the array of instances for each sequence
- To get the areas that are most conserved call the construct_top_row function from the parse class and input the instances recieved from parse_nucleotides
- To make your instances aligned with their mutations call the shift_mutations function from the mutations class and input your instances
- To get the linker strand call the create_linker_strand function from the conserved_features class and input your target strand
- To get the incumbent strand call the create_incumbent_strand function from the conserved_features class and input your target strand
- To get the longest conserved strand call the get_longest_conserved function from the conserved features class
- To get the most conserved sequence call the calculate_most_conserved_seq function for the conserved features class
- Note: you need to do the first two steps in order for any of the other functions to work since the data is based off of your fasta file
