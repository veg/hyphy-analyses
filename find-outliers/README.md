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
