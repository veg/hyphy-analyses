# Outliers Analysis

This script reads in a SLAC result file in JSON format and performs an analysis to identify problematic stretches of sequences. It provides a number of options that can be specified through keyword arguments to customize the behavior of the analysis.

Keyword Arguments

The following keyword arguments are supported by the script:

| Keyword Argument |	Description	| Default Value |
| ---------------- | ------------ | ------------- |
| code | 	Which genetic code should be used	| "Universal"
| slac |	A SLAC result file (JSON) |	N/A (required)
| window-width | The width of a sliding window to use for filtering |	5 |
| window-fraction	| The fraction of multiple-hits in a window used to filtering	| 0.4 |
| hit-threshold	| Number of substitution on a branch to count towards 'filtering' (e.g. 2 for 2 or more changes) |	1.5
| filter-partials	| Filter partially resolved codons | 	"No"
| output |	Output FASTA file with suspect sites masked |	N/A (required)
| outlier-coord-output |	Output coordinates that were masked |	N/A (required)

If the number of multiple hits in the current window (i.e., the number of residues in "outliers.not_gaps" that are not identical) is greater than or equal to the specified "window-fraction" of the total number of residues in the window, the script adds the current column to a list of "suspect" columns. The "window-fraction" argument is set to a default value of 0.4 but can be modified by the user. The script then moves the sliding window to the next set of columns and repeats the process until it has iterated over the entire MSA. The script stores the suspect column ranges for each sequence in the MSA. 

Overall, the script is designed to identify stretches of sequence-level residues that may be problematic in the context of the SLAC analysis. It does this by sliding a window across the multiple sequence alignment and identifying regions where a specified fraction of the residues in the window are multiple hits. The script also allows the user to filter out partially resolved codons and to specify a threshold for the number of substitutions on a branch that should be counted towards the "filtering" process.
