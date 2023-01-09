# Outliers Analysis

This script reads in a SLAC result file in JSON format and performs an analysis to identify problematic stretches of sequences. It provides a number of options that can be specified through keyword arguments to customize the behavior of the analysis.

## Example invocation 
`hyphy LIBPATH=/path/to/hyphy/res/ find-outliers-slac.bf  --slac test.SLAC.json`

## Keyword Arguments

The following keyword arguments are supported by the script:

| Keyword Argument |	Description	| Default Value |
| ---------------- | ------------ | ------------- |
| code | 	Which genetic code should be used	| "Universal"
| slac |	A SLAC result file (JSON) |	N/A (required)
| window-width | The width of a sliding window to use for filtering |	5 |
| window-fraction	| The fraction of multiple-hits in a window used for filtering	| 0.4 |
| hit-threshold	| Number of substitutions on a branch to count towards 'filtering' (e.g. 2 for 2 or more changes) |	1.5
| filter-partials	| Filter partially resolved codons | 	"No"
| output |	Output FASTA file with suspect sites masked |	N/A (required)
| outlier-coord-output |	Output coordinates that were masked |	N/A (required)


## Discovering sites to mask

The script parses a SLAC result file (which stands for Single-Likelihood Ancestor Counting) in JSON format, which contains data about the evolutionary history of the sequences in a phylogenetic tree. 

The algorithm works by iterating over the sequences in the tree and examining each one in a "sliding window" of a specified width. If the proportion of multiple hits in the current window (i.e., the number of residues in "outliers.not_gaps" that are not identical) is greater than or equal to the specified "window-fraction" of the total number of residues in the window, the script adds the current column to a list of "suspect" columns. The "window-fraction" argument is set to a default value of 0.4 but can be modified by the user. The script then moves the sliding window to the next set of columns and repeats the process until it has iterated over the entire MSA. The script stores the suspect column ranges for each sequence in the MSA. 

## Hit Threshold
In this script, the hit-threshold argument is used to specify the number of substitutions that must occur on a branch of the phylogenetic tree in order for the sequences in that branch to be considered problematic and subject to filtering. The hit threshold is used as a criterion for identifying stretches of sequence that may have unusually high rates of evolution and could potentially be problematic or difficult to interpret.

For example, if the hit threshold is set to 1.5, this means that any sequence in the tree that has more than 1.5 substitutions on its branch will be considered problematic and subject to filtering. The script will then use this threshold to identify stretches of sequence that have a high number of substitutions, and may perform further analysis or filtering on these stretches.

If an internal node is determined to exceed the hit threshold, then all of its descendents are marked for further outlier analysis.

## Filter Partials

In this script, the filter-partials argument is a keyword argument that specifies whether or not to filter partially resolved codons. A codon is a sequence of three nucleotides that specifies a particular amino acid in a protein. A partially resolved codon is one that has not been fully determined due to gaps or missing data in the sequence.

The filter-partials argument allows the user to specify whether partially resolved codons should be filtered out of the analysis or not. If the argument is set to "Yes", the script will treat partially resolved codons as problematic and will filter them out of the analysis. If the argument is set to "No", partially resolved codons will be treated as normal and will not be filtered out.

## Output
A FASTA file containing sequences with problematic stretches masked will be written to a path specied by `--output`. 
A JSON file listing each site masked will be written to a path specified by `--outlier-coord-output`.

## Summary
Overall, the script is designed to identify stretches of sequence-level residues that may be problematic in the context of the SLAC analysis. It does this by sliding a window across the multiple sequence alignment and identifying regions where a specified fraction of the residues in the window are multiple hits. The script also allows the user to filter out partially resolved codons and to specify a threshold for the number of substitutions on a branch that should be counted towards the "filtering" process.
