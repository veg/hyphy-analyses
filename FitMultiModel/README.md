## Fit Multiple Hit Models

This analysis fits three codon models to a codon alignment and compares the fits of all three to determine if a model that allows multiple simultaneous nucleotide substitutions fits better. The models it compares are:
-   [Muse Gaut](https://www.ncbi.nlm.nih.gov/pubmed/7968485) + REV which does NOT permit simultaneous substitutions.
-   [MG+REV+MH](https://www.nature.com/articles/s41559-018-0584-5) which allows for two simultaneous nucleotide substitutions.
-   MG+REV+TRIP which is introduced here and allows for three simultaneous nucleotide substitutions.

## Invocation

This analysis has one **required** argument

- `--alignment` the alignment file and tree (in FASTA, PHYLIP, MEGA or NEXUS formats)

HyPhy will write Markdown output to the screen and a JSON file with detailed fit results.
See example at the end of the document

### Complete options list

```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

code
        Which genetic code should be used
        default value: Universal

alignment [required]
        An in-frame codon alignment in one of the formats supported by HyPhy

tree [conditionally required]
        A phylogenetic tree (optionally annotated with {})
        applies to: Please select a tree file for the data:

rates
        The number omega rate classes to include in the model [2-10, default 3]
        default value: 3 [computed at run time]

triple-islands
        Use a separate rate parameter for synonymous triple-hit substitutions
        default value: No

output
        Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'FITTER.json')
        default value: fitter.codon_data_info[terms.json.json] [computed at run time]
        
save-fit
        Write model fit files (HyPhy NEXUS) to this file path with extensions .MODEL_NAME.bf; default is NOT to save, or 'dev/null'
        default value: /dev/null
```

## Example run

```
HYPHYMP FitMultiModel.bf --alignment p51.nex --triple-islands Yes
```

The following data are output to the screen.

Analysis Description
--------------------
Examine whether or not a codon alignment is better fit by models which
permit multiple instantaneous substitutions. v0.2 adds a separate rate
for codon-island triple-hit rates

- __Requirements__: in-frame codon alignment and a phylogenetic tree

- __Written by__: Sergei L Kosakovsky Pond, Sadie Wisotsky and Alexander Lucaci

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.2

>code –> Universal
>Loaded a multiple sequence alignment with **8** sequences, **440** codons, and **1** partitions from `/home/swisotsky/hyphy-analyses/FitMultiModel/p51.nex`
The number of omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer): 
>rates –> 3

>triple-islands –> Yes


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -3320.50, AIC-c =  6683.09 (21 estimated parameters)

### Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases
* Log(L) = -3178.60, AIC-c =  6413.66 (28 estimated parameters)
* non-synonymous/synonymous rate ratio for *test* =   0.2459

### Fitting Standard MG94
* Log(L) = -3121.04, AIC-c =  6308.73 (33 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.2967
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.186               |    86.040     |                                   |
|               4.466               |    13.674     |                                   |
|              80.143               |     0.286     |                                   |


### Fitting MG94 with double instantaneous substitutions
* Log(L) = -3121.03, AIC-c =  6310.74 (34 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.2946
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0127
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.186               |    85.954     |                                   |
|               4.448               |    13.757     |                                   |
|              79.054               |     0.289     |                                   |

### Fitting MG94 with double and triple instantaneous substitutions
* Log(L) = -3121.03, AIC-c =  6314.82 (36 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.2947
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0127
* rate at which 3 nucleotides are changed instantly within a single codon =   0.0000
* rate at which 3 nucleotides are changed instantly within a single codon between synonymous codon islands =   0.0000
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.186               |    85.952     |                                   |
|               4.444               |    13.756     |                                   |
|              78.643               |     0.291     |                                   |


### Fitting MG94 with double and triple instantaneous substitutions [only synonymous islands]
* Log(L) = -3121.03, AIC-c =  6312.78 (35 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.2946
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0126
* rate at which 3 nucleotides are changed instantly within a single codon =   0.0000
* rate at which 3 nucleotides are changed instantly within a single codon between synonymous codon islands =   0.0000
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.186               |    85.959     |                                   |
|               4.451               |    13.751     |                                   |
|              78.760               |     0.290     |                                   |


### Summary of rate estimates and significance testing
|                Model                 |   Log-L    |   omega    | 2-hit rate |              p-value               | 3-hit rate |              p-value               |
|:------------------------------------:|:----------:|:----------:|:----------:|:----------------------------------:|:----------:|:----------------------------------:|
|            Standard MG94             |  -3121.04  |    0.2967  |    N/A     |                N/A                 |    N/A     |                N/A                 |
|        Standard MG94 + 2 hits        |  -3121.03  |    0.2946  |    0.0127  |       0.8622 (2-hit rate = 0)      |    N/A     |                N/A                 |
|   Standard MG94 + 3 hits (islands)   |  -3121.03  |    0.2946  |    0.0000  |    0.9960 (3-hit island vs 2-hit)  |    0.0000  |          1.0000 (3-hit = 0)        |
|     Standard MG94 + 2 or 3 hits      |  -3121.03  |    0.2947  |    0.0127  |      0.9986 (2&3-hit rates = 0)    |    0.0000  |      1.0000 (3-hit rate(s) = 0)    |

### No individual sites showed sufficiently strong preference for multiple-hit models

### Writing detailed analysis report to `~/FitMultiModel/p51.nex.FITTER.json'
