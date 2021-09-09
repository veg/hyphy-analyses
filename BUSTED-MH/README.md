## BUSTED[S] with Multiple Hits

> Analysis suggested by Ray Cui <rcui@age.mpg.de>

This modification of BUSTED[S], which we call which we call BUSTED[SMSH], allows for a gene-wide branch-site test for episodic selection that also allows for site-to-site synonymous rate variation and multiple simultaneous hits. We combine the framework of BUSTED[S] with the MG94 + REV + TRIP codon model that accounts for single, double and triple nucleotide changes.

The core test of BUSTED remains the same, with the analysis fitting an unconstrained model with `K` rate categories (`K` = 3 by default) such that &omega;<sub>i</sub> ≤ 1 for `i` < `K`, and  &omega;<sub>K</sub> ≥ 1 and compares the fit to a constrained model (all &omega; ≤ 1) using a likelihood ratio test.

## Invokation

This analysis has one **required** argument

- `--alignment` the alignment file and tree (in FASTA, PHYLIP, MEGA or NEXUS formats)

HyPhy will write Markdown output to the screen and a JSON file with detailed fit results.
See example at the end of the document

### Complete options list

```
Available analysis command line options
--------------------------------------
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


## Example run

```
HYPHYMP  BUSTED-MH.bf --alignment CD2.nex --srv Yes --branches PR
```

The following data are output to the screen.
Analysis Description
--------------------
BUSTED-SMSH (branch-site unrestricted statistical test of episodic
diversification) uses a random effects branch-site model fitted jointly
to all or a subset of tree branches in order to test for alignment-wide
evidence of episodic diversifying selection and. Assuming there is
evidence of positive selection (i.e. there is an omega > 1), BUSTED will
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

- __Written by__: Sergei L Kosakovsky Pond, Sadie Wisotsky

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 2.2


>code –> Universal
>Loaded a multiple sequence alignment with **10** sequences, **187** codons, and **1** partitions from `/Users/Sadie/bin/hyphy-analyses/BUSTED-MH/CD2.nex`

>branches –> All

>srv –> Yes
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer):
>rates –> 3
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer):
>syn-rates –> 3
The number of points in the initial distributional guess for likelihood fitting (permissible range = [1,10000], default value = 250, integer):
>grid-size –> 250
The number of initial random guesses to 'seed' rate values optimization (permissible range = [1,25], default value = 1, integer):
>starting-points –> 1


### Branches to test for selection in the BUSTED analysis
* Selected 16 branches to test in the BUSTED analysis: `PIG, COW, Node3, HORSE, CAT, Node2, RHMONKEY, BABOON, Node9, HUMAN, CHIMP, Node12, Node8, Node1, RAT, MOUSE`


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -3532.32, AIC-c =  7112.86 (24 estimated parameters)

### Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model
* Log(L) = -3467.01, AIC-c =  6997.09 (31 estimated parameters)
* non-synonymous/synonymous rate ratio for *test* =   0.9975

### Fitting MG94 with double instantaneous substitutions
* Log(L) = -3467.01, AIC-c =  6997.09 (31 estimated parameters)
* non-synonymous/synonymous rate ratio for *test* =   0.9888

### Fitting MG94 with double and triple instantaneous substitutions
* Log(L) = -3467.01, AIC-c =  7030.49 (47 estimated parameters)
* non-synonymous/synonymous rate ratio for *test* =   0.9896

### Performing the full (dN/dS > 1 allowed) branch-site model fit
* Log(L) = -3414.49, AIC-c =  6914.96 (42 estimated parameters)
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0000
* rate at which 3 nucleotides are changed instantly within a single codon =   0.0000
* For *test* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |    9.177    |                                   |
|        Negative selection         |     0.000     |   29.502    |       Collapsed rate class        |
|      Diversifying selection       |     1.878     |   61.322    |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.275               |    26.571     |                                   |
|               0.975               |    63.188     |                                   |
|               3.038               |    10.242     |                                   |


### Performing the constrained (dN/dS > 1 not allowed) model fit
* Log(L) = -3418.83, AIC-c =  6921.54 (41 estimated parameters)
* rate at which 2 nucleotides are changed instantly within a single codon =   0.0925
* rate at which 3 nucleotides are changed instantly within a single codon =   0.0000
* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |    2.348    |                                   |
|        Negative selection         |     0.000     |   13.307    |       Collapsed rate class        |
|         Neutral evolution         |     1.000     |   84.345    |                                   |

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.327               |    33.285     |                                   |
|               1.056               |    58.091     |                                   |
|               3.221               |     8.624     |                                   |

----
## Branch-site unrestricted statistical test of episodic diversification [BUSTED]
Likelihood ratio test for episodic diversifying positive selection, **p =   0.0065**.

Check messages.log for diagnostic messages.
