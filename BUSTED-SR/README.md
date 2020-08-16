## A test for episodic conservation

This analysis is a modification of [BUSTED[S]](https://academic.oup.com/mbe/article/37/8/2430/5739973) where **test** and **background** branches have independent synonymous rate distributions. Unless there are two branch classes (test and background) and the `SRV` option is set to `Yes` (also default), the analysis will be identical to BUSTED[S]. Can be used to examine, e.g., the effects of heterotachy (not selective pressure) affecting subtrees. For instance, on the example dataset, there is no benefit of adding a separate distribution of synonymous rates (higher c-AIC), and the p-value for the selection test is largely unaffected.


#### Performing the full (dN/dS > 1 allowed) branch-site model fit
* Log(L) = -3416.83, AIC-c =  6919.64 (42 estimated parameters)

...

* The following rate distribution for site-to-site **synonymous** rate variation was inferred along **background** branches

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.497               |    62.677     |                                   |
|               1.845               |    37.323     |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.466               |    72.967     |                                   |
|               2.441               |    27.033     |                                   |

...

Likelihood ratio test for episodic diversifying positive selection, **p =   0.3763**.

**whereas starndard BUSTED[S]** has

## Performing the full (dN/dS > 1 allowed) branch-site model fit
* Log(L) = -3416.44, AIC-c =  6912.59 (39 estimated parameters)

...

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.417               |    44.192     |                                   |
|               1.462               |    55.808     |                                   |

...

Likelihood ratio test for episodic diversifying positive selection, **p =   0.2835**.

---

## Invokation

This analysis has one **required** argument

- `--alignment` the alignment file and tree
HyPhy will write Markdown output to the screen and a JSON file with detailed fit results.
See example at the end of the document.

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

branches
	Branches to test
	default value: All

srv
	Include synonymous rate variation in the model
	default value: Yes

rates
	The number omega rate classes to include in the model [1-10, default 3]
	default value: busted.rate_classes [computed at run time]

syn-rates
	The number alpha rate classes to include in the model [1-10, default 3]
	default value: busted.synonymous_rate_classes [computed at run time]

grid-size
	The number of points in the initial distributional guess for likelihood fitting
	default value: 250 [computed at run time]

starting-points
	The number of initial random guesses to seed rate values optimization
	default value: 1 [computed at run time]

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'BUSTED.json')
	default value: busted.codon_data_info[terms.json.json] [computed at run time]

save-fit
	Save BUSTED model fit to this file (default is not to save)
	default value: /dev/null

```

##Example run

Test a subset of branches, use two rate classes, turn off synonymous rate variation, set options to improve convergence, save model fit to file.

```
hyphy BUSTED-SR.bf --alignment data/example.nex --branches TEST --srv No --rates 2 --syn-rates 2 --starting-points 20 --grid-size 500 --save-fit data/fit.lf
```

---

Analysis Description
--------------------
BUSTED (branch-site unrestricted statistical test of episodic
diversification) uses a random effects branch-site model fitted jointly
to all or a subset of tree branches in order to test for alignment-wide
evidence of episodic diversifying selection. Assuming there is evidence
of positive selection (i.e. there is an omega > 1), BUSTED will also
perform a quick evidence-ratio style analysis to explore which
individual sites may have been subject to selection. v2.0 adds support
for synonymous rate variation, and relaxes the test statistic to 0.5
(chi^2_0 + chi^2_2). Version 2.1 adds a grid search for the initial
starting point. Version 2.2 changes the grid search to LHC, and adds an
initial search phase to use adaptive Nedler-Mead. Version 3.0 implements
the option for branch-site variation in synonymous substitution rates.
Version 3.1 adds HMM auto-correlation option for SRV, and binds SRV
distributions for multiple branch sets. This version allows test and
background branches to have independent synonymous rate distributions,
if selected.

- __Requirements__: in-frame codon alignment and a phylogenetic tree (optionally annotated
with {})

- __Citation__: *Gene-wide identification of episodic selection*, Mol Biol Evol.
32(5):1365-71

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 3.1


>code –> Universal
>Loaded a multiple sequence alignment with **10** sequences, **187** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/BUSTED-SR/data/example.nex`

>branches –> TEST

>srv –> Yes
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer):
>rates –> 2
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer):
>syn-rates –> 2
The number of points in the initial distributional guess for likelihood fitting (permissible range = [1,10000], default value = 250, integer):
>grid-size –> 500
The number of initial random guesses to 'seed' rate values optimization (permissible range = [1,50], default value = 1, integer):
>starting-points –> 20


### Branches to test for selection in the BUSTED analysis
* Selected 7 branches to test in the BUSTED analysis: `RHMONKEY, BABOON, Node11, HUMAN, CHIMP, Node14, Node10`


### Performing the full (dN/dS > 1 allowed) branch-site model fit
* Log(L) = -3416.83, AIC-c =  6919.64 (42 estimated parameters)
* For *test* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   53.670    |                                   |
|      Diversifying selection       |     1.831     |   46.330    |                                   |

* For *background* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.216     |   55.342    |                                   |
|      Diversifying selection       |     2.693     |   44.658    |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred along **background** branches

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.497               |    62.677     |                                   |
|               1.845               |    37.323     |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.466               |    72.967     |                                   |
|               2.441               |    27.033     |                                   |


### Performing the constrained (dN/dS > 1 not allowed) model fit
* Log(L) = -3417.11, AIC-c =  6918.11 (41 estimated parameters)
* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   29.004    |                                   |
|         Neutral evolution         |     1.000     |   70.996    |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred along **background** branches

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.499               |    63.103     |                                   |
|               1.857               |    36.897     |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.466               |    76.533     |                                   |
|               2.741               |    23.467     |                                   |

----
## Branch-site unrestricted statistical test of episodic diversification [BUSTED]
Likelihood ratio test for episodic diversifying positive selection, **p =   0.3763**.
---
