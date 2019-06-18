## A test for episodic conservation

> Analysis suggested by Sammy Cheng (scheng25@u.rochester.edu)

This analysis is a modification of [BUSTED](https://academic.oup.com/mbe/article/32/5/1365/1134918) to look for **purifying** instead of **diversifying** selection. Recall that BUSTED fits a branch site random effects model to an alignment, where the &omega; value at each branch/site is draw from a K-bin general discrete distribution (K=3 by default), with dN/dS ratio &omega;<sub>i</sub> and corresponding weights p<sub>i</sub>.  &omega;<sub>i</sub> ≤ 1 for `i` < K, and  &omega;<sub>K</sub> ≥ 1. 

The modification of BUSTED, which we call BUSTEC (Branch-site unrestricted statistical test of episodic conservation) compares the fit of the above model with the model where **all** &omega;<sub>i</sub> ≥ 1 (in the original BUSTED, &omega;<sub>K</sub> ≤ 1).

Significance is assessed by a LRT (conservative here) with 2K-2 degrees of freedom.

## Invokation

This analysis has one **required** argument

- `--alignment` the alignment file and tree 
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

tree [conditionally required]
	A phylogenetic tree (optionally annotated with {})
	applies to: Please select a tree file for the data:

branches
	Branches to test
	defaut value: All

srv
	Include synonymous rate variation in the model
	defaut value: Yes

rates
	The number omega rate classes to include in the model [1-10, default 3]
	defaut value: busted.rate_classes [computed at run time]

syn-rates
	The number alpha rate classes to include in the model [1-10, default 3]
	defaut value: busted.synonymous_rate_classes [computed at run time]

grid-size
	The number of points in the initial distributional guess for likelihood fitting
	defaut value: 250 [computed at run time]

starting-points
	The number of initial random guesses to seed rate values optimization
	defaut value: 1 [computed at run time]

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'BUSTED.json')
	defaut value: busted.codon_data_info[terms.json.json] [computed at run time]

save-fit
	Save BUSTED model fit to this file (default is not to save)
	defaut value: /dev/null

```
 

##Example run 

Test a subset of branches, use two rate classes, turn off synonymous rate variation, set options to improve convergence, save model fit to file.

```
hyphy BUSTEC.bf --alignment data/example.fas --branches TEST --srv No --starting-points 20 --grid-size 500 --save-fit data/fit.lf 
```

---

Analysis Description
--------------------
BUSTEC (branch-site unrestricted statistical test of episodic
conservation) uses a random effects branch-site model fitted jointly to
all or a subset of tree branches in order to test for alignment-wide
evidence of conservation (negative selection). Assuming there is
evidence of negative selection (i.e. there is an omega < 1), BUSTED will
also perform a quick evidence-ratio style analysis to explore which
individual sites may have been subject to selection. v2.0 adds support
for synonymous rate variation, and relaxes the test statistic to 0.5
(chi^2_0 + chi^2_2). Version 2.1 adds a grid search for the initial
starting point. Version 2.2 changes the grid search to LHC, and adds an
initial search phase to use adaptive Nedler-Mead. 

- __Requirements__: in-frame codon alignment and a phylogenetic tree (optionally annotated
with {})

- __Citation__: *Gene-wide identification of episodic selection*, Mol Biol Evol.
32(5):1365-71

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 2.2


>code –> Universal
>Loaded a multiple sequence alignment with **16** sequences, **267** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/BUSTED-conservation/data/data.fas.txt`

>branches –> TEST

>srv –> No
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer): 
>rates –> 3
The number of points in the initial distributional guess for likelihood fitting (permissible range = [1,10000], default value = 250, integer): 
>grid-size –> 500
The number of initial random guesses to 'seed' rate values optimization (permissible range = [1,50], default value = 1, integer): 
>starting-points –> 20


### Branches to test for selection in the BUSTED analysis
* Selected 5 branches to test in the BUSTED analysis: `XP_026276666_1_uncharacterized_protein_LOC113205312_Frankliniella_occidentalis, XP_026285291_1_uncharacterized_protein_LOC113211199_Frankliniella_occidentalis, XP_026285289_1_uncharacterized_protein_LOC113211197_Frankliniella_occidentalis, Node21, Node19`


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -8873.79, AIC-c = 17821.80 (37 estimated parameters)

### Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases
* Log(L) = -7991.99, AIC-c = 16074.95 (45 estimated parameters)
* non-synonymous/synonymous rate ratio for *background* =   0.0944
* non-synonymous/synonymous rate ratio for *test* =   0.1337

### Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model
* Log(L) = -7855.57, AIC-c = 15802.11 (45 estimated parameters)
* non-synonymous/synonymous rate ratio for *background* =   0.0287
* non-synonymous/synonymous rate ratio for *test* =   0.1403

### Performing the full (dN/dS < 1 and dN/dS > 1 allowed) branch-site model fit
* Log(L) = -7610.06, AIC-c = 15327.48 (53 estimated parameters)
* For *test* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.002     |   68.390    |                                   |
|        Negative selection         |     0.100     |   30.463    |                                   |
|      Diversifying selection       |    236.673    |    1.146    |                                   |

* For *background* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   74.917    |                                   |
|        Negative selection         |     0.025     |   18.861    |                                   |
|        Negative selection         |     0.862     |    6.222    |                                   |

### Performing the constrained (dN/dS < 1 not allowed) model fit
* Log(L) = -7627.61, AIC-c = 15362.59 (53 estimated parameters)
* For *test* branches under the null (no dN/dS < 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|         Neutral evolution         |     1.000     |   75.397    |                                   |
|      Diversifying selection       |    25.610     |   24.149    |                                   |
|      Diversifying selection       |9999999171.5...|    0.454    |                                   |

----
## Branch-site unrestricted statistical test of episodic conservation [BUSTEC]
Likelihood ratio test for episodic purifying selection, **p =   0.0000**.




