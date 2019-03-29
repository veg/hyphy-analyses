## Fit Multiple Hit Models

This analysis fits three codon models to a codon alignment and compares the fits of all three to determine if a model that allows multiple simultaneous nucleotide substitutions fits better. The models it compares are:
-   [Muse Gaut](https://www.ncbi.nlm.nih.gov/pubmed/7968485) + REV which does NOT permit simultaneous substitutions.
-   [MG+REV+MH](https://www.nature.com/articles/s41559-018-0584-5) which allows for two simultaneous nucleotide substitutions.
-   MG+REV+TRIP which is introduced here and allows for three simultaneous nucleotide substitutions.

## Invokation

This analysis has one **required** argument

- `--alignment` the alignment file and tree (in FASTA, PHYLIP, MEGA or NEXUS formats)

HyPhy will write Markdown output to the screen and a JSON file with detailed fit results.
See example at the end of the document

### Complete options list

```
Analysis options description
----------------------------
code
	Which genetic code should be used
	default value: Universal

alignment [required]
	An in-frame codon alignment in one of the formats supported by HyPhy

tree
	A phylogenetic tree
	default value: null [computed at run time]

rates
	The number omega rate classes to include in the model
	default value: 3 [2 - 10 accepted]

output
	Write the resulting JSON to this file
	default: same path as the alignment file + 'FITTER.json'
```


## Example run

```
HYPHYMP FitMultiModel.bf --alignment p51.nex
```

The following data are output to the screen.
Analysis Description
--------------------
Examine whether or not a codon alignment is better fit by models which
permit multiple instantaneous substitutions

- __Requirements__: in-frame codon alignment and a phylogenetic tree

- __Written by__: Sergei L Kosakovsky Pond, Sadie Wisotsky and Alexander Lucaci

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1


>code –> Universal

>alignment –> p51.nex
>Loaded a multiple sequence alignment with **8** sequences, **440** codons, and **1** partitions from `/Users/Sadie/Dropbox/hyphy-analyses/FitMultiModel/p51.nex`
The number of omega rate classes to include in the model (permissible range = [2,10], default value = 3, integer):
>rates –> 3


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -3320.50, AIC-c =  6683.09 (21 estimated parameters)

>output –> /Users/Sadie/Dropbox/hyphy-analyses/FitMultiModel/p51.nex.FITTER.json

### Branches to test for selection in the fitter analysis
* Selected 13 branches to test in the fitter analysis: `D_CD_83_ELI, D_CD_83_NDK, Node3, D_UG_94_94UG114, Node2, D_CD_84_84ZR085, Node1, B_US_83_RF, B_FR_83_HXB2, B_US_86_JRFL, Node10, B_US_90_WEAU160, Node9`

### Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases
* Log(L) = -3178.59, AIC-c =  6413.65 (28 estimated parameters)
* non-synonymous/synonymous rate ratio for *test* =   0.2455

### Fitting Standard MG94
* Log(L) = -3129.78, AIC-c =  6326.21 (33 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.2680
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.000               |     3.973     |                                   |
|               0.314               |    86.314     |                                   |
|               7.502               |     9.713     |                                   |


### Fitting MG94 with double instantaneous substitutions
* Log(L) = -3127.18, AIC-c =  6323.03 (34 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.2442
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0639
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.000               |     6.534     |                                   |
|               0.323               |    83.420     |                                   |
|               7.268               |    10.046     |                                   |


### Fitting MG94 with double and triple instantaneous substitutions
* Log(L) = -3121.09, AIC-c =  6312.90 (35 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.3051
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0098
* rate at which 3 nucleotides are changed instantly within a single codon =   0.0000
* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.173               |    85.538     |                                   |
|               4.204               |    14.194     |                                   |
|              95.390               |     0.268     |                                   |


### Summary of rate estimates and significance testing
|               Model               |   Log-L    |   omega    | 2-hit rate |             p-value             | 3-hit rate |           p-value            |
|:---------------------------------:|:----------:|:----------:|:----------:|:-------------------------------:|:----------:|:----------------------------:|
|           Standard MG94           |  -3129.78  |    0.2680  |    N/A     |               N/A               |    N/A     |             N/A              |
|      Standard MG94 + 2 hits       |  -3127.18  |    0.2442  |    0.0639  |      0.0223 (2-hit rate = 0)    |    N/A     |             N/A              |
|    Standard MG94 + 2 or 3 hits    |  -3121.09  |    0.3051  |    0.0098  |    0.0002 (2&3-hit rates = 0)   |    0.0000  |    0.0005 (3-hit rate = 0)   |

### 1 individual site which showed sufficiently strong preference for multiple-hit models
|   Site   | Evidence Ratio (2-hit)  | Evidence Ratio (3-hit)  |                       Substitutions                        |
|:--------:|:-----------------------:|:-----------------------:|:----------------------------------------------------------:|
|   245    |           2.0845        |         298.5884        |         AAA->AAC(1)ACA(1)CAT(1), GTG->AAA(1)AAG(1)         |

### Writing detailed analysis report to `/Users/Sadie/Dropbox/hyphy-analyses/FitMultiModel/p51.nex.FITTER.json'
