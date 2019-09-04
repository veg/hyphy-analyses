## RELAX-scan

> Modifications suggested by Dr. Vinny Lynch (see [this issue](https://github.com/veg/hyphy/issues/1012))

This analysis allows modifies the [RELAX](https://academic.oup.com/mbe/article/32/3/820/981440) method to _scan_ branches for those *individual* branches, where there is evidence of relaxation/intensification relative to the gene background. 

The test uses the _General Descriptive_ model, where variation in &omega; across sites and branches is handled via a **k**-bin distribution.

&omega;<sub>1</sub> ≤ 1, &omega;<sub>2</sub> ≤ 1, &omega;<sub>3</sub> ≥ 1.

Each branch gets an additional parameter, `k`, which _scales_ the _gene_ distribution, like so:

&omega;<sub>1</sub><sup>k</sup>, &omega;<sub>2</sub><sup>k</sup>, &omega;<sub>3</sub><sup>k</sup>, 

If k < 1, we infer that selection is relaxed on branch **b** relative to gene average (shrunk towards &omega; = 1), otherwise it is intensified (pulled away from &omega; = 1). 

Upon fitting the general exploratory model, **RELAX-scan** proceeds to fit a series of null hypotheses (one per branch), where `k=1` is enforced at a given branch, and the significance of `k≠1` is tested via a one degree of freedom likelihood ratio test; p-values are corrected for multiple testing using the Holm-Bonferroni procedure. The tests can be carried out using one of two approaches

1. **Fast** (default). All global model parameters are fixed at their MLEs from the general exploratory during branch-level tests. This will generally **overestimate** significance. 

2. **Slow**. All global model parameters are refitted during each branch-level test. This will accurately estimate significance, but will be **much** slower. 

A sensible procedure might be to first use the **fast** mode to see if there is any significant results, and then run the **slow** mode to confirm that the significance is still there. 

> The branch-level series of tests are `MPI enabled` and can be run in parallel on a cluster. 

## Invokation

This analysis has three **required** arguments

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
	defaut value: Universal

alignment [required]
	An in-frame codon alignment in one of the formats supported by HyPhy

tree [conditionally required]
	A phylogenetic tree (optionally annotated with {})
	applies to: Please select a tree file for the data:

rates
	The number omega rate classes to include in the model [2-10, default 3]
	defaut value: relax.rate_classes [computed at run time]

fast-test
	Use a faster test (could have some false positives)
	defaut value: Yes

```
 

##Example run 1 (fast test)

**2** (`EU912347PIPISTRELLUSABRAMUS_OPN1LW` and `Node13`) branches are found to be **relaxed** relative to the gene background.

> The dataset is Bat M/LWS1 in echolocating bats from the RELAX paper (figure 5H)


```
hyphy RELAX.bf --alignment data/Fig4H.nex
```

### Output [partial]

### Fitting the alternative model to test K != 1
* L### Fitting the general descriptive (separate k per branch) model

### * Log(L) = -1616.21, AIC-c =  3374.42 (69 estimated parameters)
* The following baseline rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.274     |    0.000    |       Not supported by data       |
|        Negative selection         |     0.507     |   100.000   |                                   |
|      Diversifying selection       |     1.001     |    0.000    |       Not supported by data       |


##### Because some of the rate classes were collapsed to 0, the model is likely overparameterized. RELAX will reduce the number of site rate classes by one and repeat the fit now.
----

Fitting the general descriptive (separate k per branch) model
* Log(L) = -1616.21, AIC-c =  3370.19 (67 estimated parameters)
* The following baseline rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.017     |   100.000   |                                   |
|      Diversifying selection       |     1.001     |    0.000    |       Not supported by data       |

* Branch-level relaxation or intensification parameter distribution has mean  3.50, median  0.70, and 95% of the weight in  0.00 -  9.04

### Testing individual branches for relaxation/intensification of selection

|              Branch name               |   k MLE    |    LRT     |Uncorrected p-value|
|:--------------------------------------:|:----------:|:----------:|:-----------------:|
|EU912339_HARPYIONYCTERIS_CELEBENSIS_O...|      8.84  |      0.73  |       0.3941      |
|   EU912342_CYNOPTERUS_SPHINX_OPN1LW    |      8.80  |      0.69  |       0.4045      |
|  EU912341_CHAEREPHON_PLICATUS_OPN1LW   |      0.82  |      0.37  |       0.5439      |
|  EU912338_ACERODON_CELEBENSIS_OPN1LW   |      0.70  |      0.93  |       0.3337      |
|  EU912340_ARTIBEUS_JAMAICENSIS_OPN1LW  |      0.79  |      0.79  |       0.3748      |
|  EU912343_HIPPOSIDEROS_ARMIGER_OPN1LW  |      0.70  |      0.85  |       0.3572      |
|EU912344_MINIOPTERUS_SCHREIBERSII_OPN...|      0.66  |      3.41  |       0.0649      |
|    EU912345_MYOTIS_RICKETTI_OPN1LW     |      0.60  |      1.50  |       0.2202      |
|    EU912346_MEGADERMA_SPASMA_OPN1LW    |      0.69  |      2.90  |       0.0884      |
|   EU912348_PTEROPUS_GIGANTEUS_OPN1LW   |      7.76  |      0.11  |       0.7374      |
|EU912349_RHINOLOPHUS_FERRUMEQUINUM_OP...|      0.46  |      8.74  |       0.0031      |
|EU912350_ROUSETTUS_LESCHENAULTII_OPN1...|      0.70  |      0.93  |       0.3355      |
|EU912351_TAPHOZOUS_MELANOPOGON_OPN1LW...|      8.50  |      1.18  |       0.2764      |
|  EU912347_PIPISTRELLUS_ABRAMUS_OPN1LW  |      0.28  |     21.39  |       0.0000      |
|                 Node11                 |      8.59  |      0.21  |       0.6452      |
|                 Node18                 |      9.04  |      0.14  |       0.7044      |
|                 Node13                 |      0.00  |     10.17  |       0.0014      |
|                 Node16                 |     10.82  |      0.55  |       0.4586      |
|                 Node19                 |      0.38  |      2.51  |       0.1132      |
|                 Node2                  |      0.49  |      2.57  |       0.1089      |
|                 Node24                 |      7.91  |      0.62  |       0.4324      |
|                 Node22                 |      0.27  |      4.75  |       0.0293      |
|                 Node5                  |      0.58  |      2.68  |       0.1015      |
|                 Node9                  |      8.79  |      0.70  |       0.4024      |
|                 Node7                  |      0.39  |     11.62  |       0.0007      |

###Tests of individual branches for relaxation/intensification of selection

|              Branch name               |   k MLE    |    LRT     | Corrected p-value |
|:--------------------------------------:|:----------:|:----------:|:-----------------:|
|  EU912338_ACERODON_CELEBENSIS_OPN1LW   |     50.00  |      3.98  |         0.88      |
|EU912339_HARPYIONYCTERIS_CELEBENSIS_O...|      2.03  |      0.33  |         1.00      |
|   EU912342_CYNOPTERUS_SPHINX_OPN1LW    |      2.03  |      0.31  |         1.00      |
|  EU912341_CHAEREPHON_PLICATUS_OPN1LW   |     15.49  |     -0.02  |         1.00      |
|  EU912340_ARTIBEUS_JAMAICENSIS_OPN1LW  |     12.44  |     -0.60  |         1.00      |
|  EU912343_HIPPOSIDEROS_ARMIGER_OPN1LW  |      0.50  |      0.85  |         1.00      |
|   EU912348_PTEROPUS_GIGANTEUS_OPN1LW   |      1.98  |      0.05  |         1.00      |
|EU912344_MINIOPTERUS_SCHREIBERSII_OPN...|      0.44  |      3.62  |         0.97      |
|    EU912346_MEGADERMA_SPASMA_OPN1LW    |     45.99  |      3.78  |         0.94      |
|    EU912345_MYOTIS_RICKETTI_OPN1LW     |      0.31  |      1.94  |         1.00      |
|  EU912347_PIPISTRELLUS_ABRAMUS_OPN1LW  |      0.10  |     21.52  |       0.00 (*)    |
|EU912350_ROUSETTUS_LESCHENAULTII_OPN1...|     50.00  |      4.03  |         0.89      |
|                 Node11                 |      1.99  |      0.09  |         1.00      |
|EU912351_TAPHOZOUS_MELANOPOGON_OPN1LW...|      2.06  |      0.54  |         1.00      |
|EU912349_RHINOLOPHUS_FERRUMEQUINUM_OP...|      0.22  |      8.41  |         0.09      |
|                 Node13                 |      0.00  |      9.76  |       0.04 (*)    |
|                 Node18                 |      1.99  |      0.07  |         1.00      |
|                 Node16                 |      2.05  |      0.42  |         1.00      |
|                 Node24                 |      2.02  |      0.28  |         1.00      |
|                 Node2                  |      0.25  |      3.54  |         0.96      |
|                 Node22                 |      0.08  |      6.17  |         0.27      |
|                 Node19                 |      0.16  |      2.36  |         1.00      |
|                 Node5                  |      0.34  |      2.80  |         1.00      |
|                 Node9                  |      2.03  |      0.31  |         1.00      |
|                 Node7                  |      0.21  |      7.70  |         0.12      |


## RELAX scan result
**2** branches had significant relaxation/intensification of selection at significance level of 0.05

##Example run 2 (slow test)

**0** branches are found to be **different** relative to the gene background.

###Tests of individual branches for relaxation/intensification of selection

|              Branch name               |   k MLE    |    LRT     | Corrected p-value |
|:--------------------------------------:|:----------:|:----------:|:-----------------:|
|  EU912340_ARTIBEUS_JAMAICENSIS_OPN1LW  |      0.58  |      0.00  |         0.97      |
|  EU912341_CHAEREPHON_PLICATUS_OPN1LW   |      0.61  |      0.00  |         1.00      |
|   EU912342_CYNOPTERUS_SPHINX_OPN1LW    |      5.95  |      0.06  |         1.00      |
|  EU912338_ACERODON_CELEBENSIS_OPN1LW   |      0.52  |      0.01  |         1.00      |
|EU912339_HARPYIONYCTERIS_CELEBENSIS_O...|      6.02  |      0.06  |         1.00      |
|EU912344_MINIOPTERUS_SCHREIBERSII_OPN...|      0.49  |      0.00  |         1.00      |
|   EU912348_PTEROPUS_GIGANTEUS_OPN1LW   |      8.83  |      0.03  |         1.00      |
|    EU912346_MEGADERMA_SPASMA_OPN1LW    |      0.51  |      0.00  |         1.00      |
|  EU912343_HIPPOSIDEROS_ARMIGER_OPN1LW  |      0.52  |      0.01  |         1.00      |
|  EU912347_PIPISTRELLUS_ABRAMUS_OPN1LW  |      0.21  |      0.01  |         1.00      |
|                 Node11                 |      5.89  |      0.05  |         1.00      |
|    EU912345_MYOTIS_RICKETTI_OPN1LW     |      0.44  |      0.01  |         1.00      |
|EU912349_RHINOLOPHUS_FERRUMEQUINUM_OP...|      0.34  |      0.00  |         1.00      |
|                 Node18                 |      6.80  |      0.04  |         1.00      |
|EU912350_ROUSETTUS_LESCHENAULTII_OPN1...|      0.52  |      0.01  |         1.00      |
|EU912351_TAPHOZOUS_MELANOPOGON_OPN1LW...|      6.29  |      0.06  |         1.00      |
|                 Node16                 |      6.12  |      0.06  |         1.00      |
|                 Node2                  |      0.36  |      0.01  |         1.00      |
|                 Node24                 |      5.89  |      0.07  |         1.00      |
|                 Node19                 |      0.28  |      0.04  |         1.00      |
|                 Node5                  |      0.43  |      0.01  |         1.00      |
|                 Node22                 |      0.20  |      0.05  |         1.00      |
|                 Node7                  |      0.28  |      0.01  |         1.00      |
|                 Node9                  |      5.96  |      0.07  |         1.00      |
|                 Node13                 |      0.00  |      0.58  |         1.00      |


```
hyphy RELAX.bf --alignment data/Fig4H.nex --fast-test No
```


