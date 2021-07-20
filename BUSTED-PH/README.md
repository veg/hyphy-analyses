## Is there evidence of selection associated with a trait/phenotype?

> In collaboration with the [Nathan Clark lab](https://clark.genetics.utah.edu)


The `BUSTED-PH` analysis builds upon [BUSTED](https://academic.oup.com/mbe/article/32/5/1365/1134918) which is a method to test for evidence of epidosic diversifying selection. The analysis seeks to answer the question: 

>Is a specific feature/phenotype/trait associated with positive selection?

The analysis requires that a coding alignment with a tree that includes a subset of branches designated as **test** set (i.e. those that have feature/phenotype/trait); the rest of the branches are placed in the **background** set (i.e. those that lack feature/phenotype/trait; see an example tree in `examples/positive-test.fas`). For details on how to generate annotations interactively, see [http://www.hyphy.org/tutorials/phylotree/](http://www.hyphy.org/tutorials/phylotree/)


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
hyphy BUSTED-PH.bf --alignment examples/positive-test.fas --srv No
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