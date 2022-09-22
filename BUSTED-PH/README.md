## Is there evidence of selection associated with a trait/phenotype?

> In collaboration with the [Nathan Clark lab](https://clark.genetics.utah.edu)


The `BUSTED-PH` analysis builds upon [BUSTED](https://academic.oup.com/mbe/article/32/5/1365/1134918) which is a method to test for evidence of epidosic diversifying selection. The analysis seeks to answer the question: 

>Is a specific feature/phenotype/trait associated with positive selection?

The analysis requires that a coding alignment with a tree that includes a subset of branches designated as **test** set (i.e. those that have feature/phenotype/trait); and a subset of branches in the **background** set (i.e. those that lack feature/phenotype/trait; see an example tree in `examples/positive-test.fas`). For details on how to generate annotations interactively, see [http://www.hyphy.org/tutorials/phylotree/](http://www.hyphy.org/tutorials/phylotree/). 

The default options is that all branchs **not** in the test set are in the background set. If there are two or more branch sets labeled in the tree, there is also an option to designate only a subset of them as **background** and treat the rest as **nuisance**. The comparisons will **exclude** any nuisance branches (see an example below).


The BUSTED-PH analysis proceeds as follows

1). An **unrestricted** branch-site REL model is fitted; this model endows **test** and **background** sets with independent &omega; distributions

2). A **constrained** model is fitted, where &omega; ≤ 1 is enforced on the **test** branches. A likelihood ratio test of this model and the unrestricted model is used to determine whether or not **test** branches are subject to episodic diversifying selection.

3). A **constrained** model is fitted, where &omega; ≤ 1 is enforced on the **background** branches. A likelihood ratio test of this model and the unrestricted model is used to determine whether or not **background** branches are subject to episodic diversifying selection.

4). A **constrained** model is fitted, where the &omega; distribution is the same for **test** and **background** branches. A likelihood ratio test of this model and the unrestricted model is used to determine whether or selective regimes **differ** between **test** and **background** branches. 

The following decision tree is then applied to provide guidance on whether or not episodic selection is associated with the presence of a trait/phenotype.

```
├── Test in (2) is significant at p≤0.05?
    ├── No
    	|==> No evidence of selection association with phenotype/trait
    ├── Yes
    	|-- Test in (3) is significant at p≤0.05?
    		├── Yes
    			|-- Test in (4) is significant at p≤0.05?
    				├── Yes
    					|==> Selection is occuring both 
    					 	  with / without phenotype/trait
    					 	  but with different selective regimes
					├── No
    					|==> Selection is occuring both 
    					 	  with / without phenotype/trait
    					 	  with similar selective regimes        		├── No
    			|-- Test in (4) is significant at p≤0.05?
    				├── Yes
    					|==> Strong evidence of selection 
    						  being associated with phenotype/trait
					├── No
    					|==> Evidence of selection 
    						  being associated with phenotype/trait,
    						  but no statistically significant 
    						  differences between test and background
    						  selective regimes		    	


```

### Example

Command line call

``` 
hyphy BUSTED-PH.bf 
--alignment examples/positive-test.fas 
--srv No 
--branches Primates
```

#### Output 


Analysis Description
--------------------
BUSTED-PH (BUSTED-PHenotype) uses several random effects branch-site
model fitted to subsets of tree branches that represent tree
partitioning by discrete trait to identify evidence of positive
diversifying selection associated with the presence of the discrete
trait. This is accomplished by fitting a BUSTED model to with separate
rate distributions for PH+/PH- branches, testing that (i) PH+ subset is
under episodic diversifying selection, (ii) PH- is NOT under
diversifying selection and (iii) PH+ and PH- have different
distributions 

- __Requirements__: in-frame codon alignment and a phylogenetic tree with exactly one
non-trivial branch subset ({} format)

- __Citation__: *TBD*

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1


>code –> Universal
>Loaded a multiple sequence alignment with **10** sequences, **500** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/BUSTED-PH/examples/positive.replicate.1`

>srv –> No
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer): 
>rates –> 3
The number of points in the initial distributional guess for likelihood fitting (permissible range = [1,10000], default value = 250, integer): 
>grid-size –> 250
The number of initial random guesses to 'seed' rate values optimization (permissible range = [1,25], default value = 1, integer): 
>starting-points –> 1


### Branches to test for selection in the BUSTED analysis
* Selected 7 branches to test in the BUSTED analysis: `Baboon, Chimp, Human, Node12, Node8, Node9, RhMonkey`


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model

>kill-zero-lengths –> Yes
* Log(L) = -8696.59, AIC-c = 17441.25 (24 estimated parameters)

### Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases
* Log(L) = -8013.31, AIC-c = 16091.05 (32 estimated parameters)
* non-synonymous/synonymous rate ratio for *background* =   0.1369
* non-synonymous/synonymous rate ratio for *test* =   0.3531

### Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model
* Log(L) = -7931.73, AIC-c = 15927.88 (32 estimated parameters)
* non-synonymous/synonymous rate ratio for *background* =   0.1073
* non-synonymous/synonymous rate ratio for *test* =   2.7179

### Performing the full (dN/dS > 1 allowed) branch-site model fit
* Log(L) = -7931.66, AIC-c = 15943.99 (40 estimated parameters)
* For *test* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.614     |    1.266    |                                   |
|        Negative selection         |     0.971     |   19.929    |                                   |
|      Diversifying selection       |     3.342     |   78.805    |                                   |

* For *background* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.107     |   100.000   |                                   |
|        Negative selection         |     0.121     |    0.000    |       Not supported by data       |
|        Negative selection         |     0.400     |    0.000    |       Not supported by data       |


### Performing the constrained (dN/dS > 1 not allowed) model fit
* Log(L) = -7938.59, AIC-c = 15955.80 (39 estimated parameters)
* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|         Neutral evolution         |     1.000     |    1.057    |                                   |
|         Neutral evolution         |     1.000     |    3.034    |       Collapsed rate class        |
|         Neutral evolution         |     1.000     |   95.909    |       Collapsed rate class        |

* For *background* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.109     |   100.000   |                                   |
|        Negative selection         |     0.115     |    0.000    |       Not supported by data       |
|        Negative selection         |     0.374     |    0.000    |       Not supported by data       |


### No evidence for episodic diversifying positive selection on background branches under the unconstrained model, skipping constrained model fitting

### Performing the shared distribution (same on test and background brances) model fit
* Log(L) = -8034.94, AIC-c = 16140.38 (35 estimated parameters)
* For the shared rates model (same between test and background), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.146     |    0.110    |                                   |
|        Negative selection         |     0.147     |   99.448    |       Collapsed rate class        |
|      Diversifying selection       |    22.787     |    0.441    |                                   |

----
## Branch-site unrestricted statistical test of episodic diversification and association with phenotype/trait [BUSTED-PH]
Likelihood ratio test for episodic diversifying positive selection on test branches , **p =   0.0005**.
Likelihood ratio test for episodic diversifying positive selection on background branches , **p =   0.5000**.
Likelihood ratio test for differences in distributions between **test** and **background** , **p =   0.0000**.


## Analysis summary (p = 0.05)
Selection is associated with the phenotype / trait

### Example (with nuisance)

Command line call

``` 
hyphy BUSTED-PH.bf 
--alignment examples/test-contrast-nuisance.fas 
--tree examples/test-contrast-nuisance.nwk 
--srv No 
--branches TEST 
--comparison REFERENCE 
```

Output 
----

Analysis Description
--------------------
BUSTED-PH (BUSTED-PHenotype) uses several random effects branch-site
model fitted to subsets of tree branches that represent tree
partitioning by discrete trait to identify evidence of positive
diversifying selection associated with the presence of the discrete
trait. This is accomplished by fitting a BUSTED model to with separate
rate distributions for PH+/PH- branches, testing that (i) PH+ subset is
under episodic diversifying selection, (ii) PH- is NOT under
diversifying selection and (iii) PH+ and PH- have different
distributions. Version 0.2 adds support for >1 labeled classes. 

- __Requirements__: in-frame codon alignment and a phylogenetic tree with exactly one
non-trivial branch subset ({} format)

- __Citation__: *TBD*

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.2


>code –> Universal
/Users/sergei/Development/hyphy-analyses/BUSTED-PH/examples/test-contrast-nuisance.nwk
examples/test-contrast-nuisance.nwk

>Loaded a multiple sequence alignment with **57** sequences, **250** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/BUSTED-PH/examples/test-contrast-nuisance.fas`

>branches –> TEST

>comparison –> REFERENCE

>srv –> No
The number omega rate classes to include in the model (permissible range = [1,10], default value = 3, integer): 
>rates –> 3
The number of points in the initial distributional guess for likelihood fitting (permissible range = [1,10000], default value = 250, integer): 
>grid-size –> 250
The number of initial random guesses to 'seed' rate values optimization (permissible range = [1,25], default value = 1, integer): 
>starting-points –> 1


### Branches to test for selection in the BUSTED analysis
* Selected 45 branches to test in the BUSTED analysis: `Amissidens_hainesi_RB0374, Arius_aff_leptaspis_RB0493, Arius_augustus_RB0475, Arius_berneyi_RB0526, Arius_dioctes_RB0463, Arius_leptaspis_RB0394, Arius_mastersi_RB0482, Arius_midgleyi_RB0543, Arius_proximus_RB0375, Brustiarius_nox_RB0369, Brustiarius_solidus_RB0355, Cinetodus_froggatti_RB0449, Cochlefelis_danielsi_RB0425, Cochlefelis_spatula_RB0458, Doiichthys_novaeguineae_RB0474, Nedystoma_dayi_RB0466, Nemapteryx_aff_armiger_RB0457, Nemapteryx_armiger_RB0383, Neoarius_aff_graeffei_RB0347, Neoarius_graeffei_RB0537, Node102, Node103, Node106, Node107, Node110, Node68, Node70, Node71, Node72, Node75, Node76, Node79, Node81, Node83, Node85, Node88, Node89, Node91, Node92, Node95, Node98, Node99, Potamosilurus_coatesi_RB0362, Potamosilurus_latirostris_RB0346, Potamosilurus_macrorhynchus_RB0528`


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model

>kill-zero-lengths –> Yes

### Deleted 17 zero-length internal branches: `Node11, Node16, Node25, Node29, Node3, Node31, Node43, Node49, Node51, Node63, Node66, Node7, Node70, Node71, Node75, Node91, Node99`
* Log(L) = -5428.95, AIC-c = 11096.57 (119 estimated parameters)
* 1 partition. Total tree length by partition (subs/site)  1.108

### Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases
* Log(L) = -5392.49, AIC-c = 11008.74 (111 estimated parameters)
* 1 partition. Total tree length by partition (subs/site)  1.107
* non-synonymous/synonymous rate ratio for *background* =   0.6206
* non-synonymous/synonymous rate ratio for *nuisance* =   0.6645
* non-synonymous/synonymous rate ratio for *test* =   0.5960

### Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model
* Log(L) = -5389.85, AIC-c = 11003.46 (111 estimated parameters)
* non-synonymous/synonymous rate ratio for *background* =   0.5563
* non-synonymous/synonymous rate ratio for *nuisance* =   0.6865
* non-synonymous/synonymous rate ratio for *test* =   0.5003

### Performing the full (dN/dS > 1 allowed) branch-site model fit
|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.341     |   99.228    |                                   |
|         Neutral evolution         |     1.000     |    0.236    |                                   |
|      Diversifying selection       |    403.342    |    0.536    |                                   |

* For *background* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.463     |   99.155    |                                   |
|         Neutral evolution         |     1.000     |    0.071    |                                   |
|      Diversifying selection       |   10000.000   |    0.774    |                                   |

* For *nuisance* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.405     |   92.161    |                                   |
|        Negative selection         |     0.514     |    5.579    |                                   |
|      Diversifying selection       |    265.362    |    2.261    |                                   |


### Performing the constrained (dN/dS > 1 not allowed) model fit
* Log(L) = -5088.27, AIC-c = 10422.66 (122 estimated parameters)
* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   64.737    |                                   |
|         Neutral evolution         |     1.000     |   10.252    |                                   |
|         Neutral evolution         |     1.000     |   25.012    |       Collapsed rate class        |

* For *background* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.454     |   99.155    |                                   |
|         Neutral evolution         |     1.000     |    0.085    |                                   |
|      Diversifying selection       |   10000.000   |    0.759    |                                   |

* For *nuisance* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.404     |   92.161    |                                   |
|        Negative selection         |     0.405     |    5.588    |       Collapsed rate class        |
|      Diversifying selection       |    240.849    |    2.252    |                                   |


### Performing the constrained background (dN/dS > 1 not allowed on background branches) model fit
* Log(L) = -5062.92, AIC-c = 10371.97 (122 estimated parameters)
* For *test* branches under the null (no dN/dS > 1 on background branches model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.343     |   99.227    |                                   |
|         Neutral evolution         |     1.000     |    0.278    |                                   |
|      Diversifying selection       |    807.685    |    0.495    |                                   |

* For *background* branches under the null (no dN/dS > 1 on background branches model), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   59.705    |                                   |
|         Neutral evolution         |     1.000     |    3.231    |                                   |
|         Neutral evolution         |     1.000     |   37.064    |       Collapsed rate class        |

* For *nuisance* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.399     |    5.656    |                                   |
|        Negative selection         |     0.412     |   92.240    |       Collapsed rate class        |
|      Diversifying selection       |    265.362    |    2.104    |                                   |


### Performing the shared distribution (same on test and background brances) model fit
* Log(L) = -5021.97, AIC-c = 10281.94 (118 estimated parameters)
* For the shared rates model (same between test and background), the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.379     |   99.217    |                                   |
|         Neutral evolution         |     1.000     |    0.235    |                                   |
|      Diversifying selection       |    807.685    |    0.548    |                                   |

* For *nuisance* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.411     |   92.161    |                                   |
|        Negative selection         |     0.435     |    5.584    |                                   |
|      Diversifying selection       |    265.362    |    2.255    |                                   |

----
## Branch-site unrestricted statistical test of episodic diversification and association with phenotype/trait [BUSTED-PH]
Likelihood ratio test for episodic diversifying positive selection on test branches , **p =   0.0000**.
Likelihood ratio test for episodic diversifying positive selection on background branches , **p =   0.0000**.
Likelihood ratio test for differences in distributions between **test** and **background** , **p =   0.4780**.


## Analysis summary (p = 0.05)
Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.
There is no significant difference between test and background branches in terms of selective pressure