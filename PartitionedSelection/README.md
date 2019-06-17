## Partitioned analysis of selection.

> Analysis inspired by Marina Escalera Zamudio marinaescalera@gmail.com

The purpose of this analysis is to compare **selective pressures** between partitions in an alignment of coding sequences. The input is comprised of a `NEXUS` file with at least two defined partitions (see `data/Measles_Fusion_Complete.nex` for an example) as `CHARSET` in the `ASSUMPTION` block. 

The analysis fits a series of codon-substitution models (jointly to all partitions) and performs statistical comparisons between them to determine whether or not `dN/dS` distributions differ between partitions. The workflow is as follows, and is controlled by several options, described below. 

1. Fit the GTR nucleotide model to obtain nucleotide substitution biases and (per partition) branch lengths. These estimates are used as initial values for subsequent optimizations. Nucleotide substitution biases and equilibirum frequency estimates are shared by all partitions. 
2. (unless turned off by the `--mg No` option) Fit the [Muse-Gaut 94](https://www.ncbi.nlm.nih.gov/pubmed/7968485) model with a single `dN/dS` ratio per partition. Report these estimates and corresponding confidence intervals. 
	* 	Perform a likelihood ratio test for equality of partition `dN/dS` estimates. If the test is singificant, then **mean selective pressures differ between partitions**.
3. (unless turned off by the `--bs-rel No` option) Fit the [Branch-site REL](https://academic.oup.com/mbe/article/32/5/1365/1134918) model with `K` rates (controlled by the `--rates` argument, default **2**) to partitions. Each partition is therefore described by its own K-rate distribution of `dN/dS` values, which vary across both **sites and branches**
	* 	If any of the partitions have estimates of `dN/dS > 1`, perform the test for positive selection. This test compares the null of having all `dN/dS` values constrained to be ≤1, versus the alternative that `dN/dS > 1` somewhere along the tree in at least one of the partitions. If the test is singificant, then **there is evidence for positive selection on some of the partitions**.
	* 	Compare the alternative model (different distributions on each of the partitions) to the null model, where a single K-rate `dN/dS` distribution is If the test is singificant, then **there is evidence for differences in selective pressure between some of the partitions**.

### Important options

#### Branch subset

The test can be run on all branches in the tree or a subset of the branches in the tree (e.g. internal branches only). In the case where only **some** of the branches are selected, the "background" branches get their own (per partition) nuisance parameters (`dN/dS`) which are estimated but not used for testing. An important use-case is to test for selection on internal branches of pathogen trees to reduce `dN/dS` inflation issues due to "within-host" variation at the tips (see [Pond et al](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020062) and [Pybus et al](https://academic.oup.com/mbe/article/24/3/845/1246056)

#### Branch lengths

Branch lengths between partitions can be either compeletely independent or proporitional (default). The latter scenario assumes that relative lineage specific evolutionary rates (other than due to selective pressures) do not vary between partitions. The former scenario is more general, but has many more parameters, and assumes general heterotachy.


## Invokation

This analysis has one **required** argument

- `--alignment` the alignment file (NEXUS with at least two partitions) and tree (either one shared by all partitions, or one per partition)

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
	defaut value: Universal

alignment [required]
	An in-frame codon alignment in one of the formats supported by HyPhy

branches 
	The set of branches to test [All, Internal, Leaves]
	defaut value: All

rates
	The number omega rate classes to include in the model [1-10, default 2]
	defaut value: partitioned.rate_classes [computed at run time]

branch-lengths
	Proportional branch lengths across partitions [Proportion or Unlinked]
	defaut value: Proportional

bs-rel
	Run BS-REL analyses [Yes, No]
	defaut value: Yes

mg
	Run MG94 analyses [Yes, No]
	defaut value: Yes

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'partitioned.json')
	defaut value: partitioned.codon_data_info[terms.json.json] [computed at run time]

```
 

##Example run 


```
hyphy PartitionedAnalysis.bf --alignment data/Measles_Fusion_Complete.nex --branches Internal
```

---

Analysis Description
--------------------
Fit several models of partition-level dN/dS variation to _partitioned_
data and perform inference on whether or not there is evidence of
different selective pressures among partitions

- __Requirements__: in-frame codon alignment with at least two partitions and a phylogenetic
tree

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1


>code –> Universal
>Loaded a multiple sequence alignment with **56** sequences, **550** codons, and **4** partitions from `/Users/sergei/Development/hyphy-analyses/PartitionedSelection/data/Measles_Fusion_Complete.nex`

>branches –> Internal
The number omega rate classes to include in the model (permissible range = [1,10], default value = 2, integer): 
>rates –> 2

>branch-lengths –> Proportional

>bs-rel –> Yes

>mg –> Yes


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -4971.94, AIC-c = 10836.17 (444 estimated parameters)

### Fitting MG94xREV with separate dN/dS ratios for each partition
* Log(L) = -4862.24, AIC-c =  9975.51 (125 estimated parameters)
* non-synonymous/synonymous rate ratio for *G1|background* =   0.2249
* non-synonymous/synonymous rate ratio for *G1|test* =   0.0838 (95% profile CI   0.0500-  0.1304)
* non-synonymous/synonymous rate ratio for *G2|background* =   0.1992
* non-synonymous/synonymous rate ratio for *G2|test* =   0.0351 (95% profile CI   0.0126-  0.0756)
* non-synonymous/synonymous rate ratio for *G3|background* =   0.3343
* non-synonymous/synonymous rate ratio for *G3|test* =   0.1020 (95% profile CI   0.0253-  0.2651)
* non-synonymous/synonymous rate ratio for *G4|background* =   0.1416
* non-synonymous/synonymous rate ratio for *G4|test* =   0.1203 (95% profile CI   0.0491-  0.2377)

### Checking for difference in partition-level dN/dS using a likelihood ratio test
* Log(L) = -4864.73, AIC-c =  9974.44 (122 estimated parameters)
* Test for dN/dS equality across partitions:  p-value =   0.1730, LRT =     4.98

### Fitting partition-specific BS-REL models with 2 rate classes
* Log(L) = -4845.65, AIC-c = 10015.01 (161 estimated parameters)
* Inferred rate distribution for `G1`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.001     |   98.788    |                                   |
|      Diversifying selection       |     9.314     |    1.212    |                                   |

* _background_ rate distribution for `G1`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   97.737    |                                   |
|      Diversifying selection       |    13.622     |    2.263    |                                   |

* Inferred rate distribution for `G2`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.034     |   100.000   |                                   |
|      Diversifying selection       |     1.001     |    0.000    |       Not supported by data       |

* _background_ rate distribution for `G2`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   98.907    |                                   |
|      Diversifying selection       |    34.938     |    1.093    |                                   |

* Inferred rate distribution for `G3`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.121     |   100.000   |                                   |
|      Diversifying selection       |     1.001     |    0.000    |       Not supported by data       |

* _background_ rate distribution for `G3`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |    1.545    |                                   |
|        Negative selection         |     0.318     |   98.455    |                                   |

* Inferred rate distribution for `G4`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.141     |   100.000   |                                   |
|      Diversifying selection       |     1.001     |    0.000    |       Not supported by data       |

* _background_ rate distribution for `G4`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.087     |   28.956    |                                   |
|        Negative selection         |     0.133     |   71.044    |                                   |


### Testing for positive selection on any of the partitions
* Log(L) = -4846.39, AIC-c = 10008.41 (157 estimated parameters)
* Inferred rate distribution for `G1`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   91.471    |                                   |
|         Neutral evolution         |     1.000     |    8.529    |                                   |

* _background_ rate distribution for `G1`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   97.749    |                                   |
|      Diversifying selection       |    13.820     |    2.251    |                                   |

* Inferred rate distribution for `G2`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.034     |   100.000   |                                   |
|         Neutral evolution         |     1.000     |    0.000    |       Not supported by data       |

* _background_ rate distribution for `G2`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   98.957    |                                   |
|      Diversifying selection       |    37.388     |    1.043    |                                   |

* Inferred rate distribution for `G3`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.121     |   100.000   |                                   |
|         Neutral evolution         |     1.000     |    0.000    |       Not supported by data       |

* _background_ rate distribution for `G3`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |    6.244    |                                   |
|        Negative selection         |     0.333     |   93.756    |                                   |

* Inferred rate distribution for `G4`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.141     |   100.000   |                                   |
|         Neutral evolution         |     1.000     |    0.000    |       Not supported by data       |

* _background_ rate distribution for `G4`

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.093     |   24.374    |                                   |
|        Negative selection         |     0.128     |   75.626    |                                   |

* Test for positive selection on any of the partitions:  p-value =   0.8299, LRT =     1.48

### Testing for equality of rate distributions across partitions
* Log(L) = -4860.30, AIC-c = 10026.11 (152 estimated parameters)
* Test for selective differences between partitions:  p-value =   0.0006, LRT =    29.28
* Inferred joint rate distribution

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   98.437    |                                   |
|      Diversifying selection       |    11.486     |    1.563    |                                   |

----
## Partitioned analyses test summary (corrected p-values)
* Likelihood ratio test for Test for selective differences between partitions, **p =   0.0017**.
* Likelihood ratio test for Test for dN/dS equality across partitions, **p =   0.3460**.
* Likelihood ratio test for Test for positive selection on any of the partitions, **p =   0.8299**.

--- 


