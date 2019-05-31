## Modifying the RELAX test

> Modifications suggested by Dr. Rongfeng (Ray) Cui

This analysis allows for _limited_ call-back based modification of the testing procedure described in the [RELAX](https://academic.oup.com/mbe/article/32/3/820/981440) paper. In the original RELAX method, there are (at least) two sets of branches in a phylogenetic tree. One is designated as **R**eference, and the other as **T**est. Each partition of branches is endowed with a 3-bin branch-site REL (for details see the original paper) model to account for variation in selective pressures (&omega;) across sites in branches. 
&omega;<sub>1</sub> ≤ 1, &omega;<sub>2</sub> ≤ 1, &omega;<sub>3</sub> ≥ 1.

In RELAX, &omega;<sub>T</sub> := &omega;<sub>R</sub> <sup>K</sup>. A likelihood ratio test is performed to check if K ≠ 1. If K < 1, we infer that selection is relaxed on **T** relative to **R** (shrunk towards &omega; = 1), otherwise it is intensified (pulled away from &omega; = 1). 

This analysis, `RELAX-mod` allows for sligthly different parametrizations of the test. To modify the test, it is necessary to create a file in the `constraint-definition` directory. This file must contain several following HyPhy Batch Language functions that will be called by **RELAX-mod**. 

**<span style = 'color:red'>The modification apply ONLY to the MINIMAL set of model runs; All models will continue using the standard RELAX parameterization </span>**


We discuss the example of `increase-only` parameterization 
> to define your own, make a copy of this file, rename it, and edit the functions as needed, then supply the name of the file as the `--parameterization` argument

&omega;<sub>T</sub> := &omega;<sub>R</sub> <sup>K</sup>, if &omega;<sub>R</sub> > 1.

&omega;<sub>T</sub> := &omega;<sub>R</sub> <sup>1/K</sup>, if &omega;<sub>R</sub> ≤ 1.

K ≥ 1. 

In this parameterization, negative selection on **T**est is relaxed (moved closer to 1), while positive selection is intensified (moved away from 1), assuming K > 1. When K = 1, there is no difference between the two sets of branches (this is the null hypothesis).


```
relax.skip_convergence_checking = TRUE;
/** this variable must be set because some convergence testing heuristics 
in RELAX assume the standard w^K parameterization **/

/**
 * @name relax.declare_parameter_with_ranges
 
 * Define and constrain the range of the "K" parameter 
 
 * @param group_index {Number}  the index of the branch group (1 in standard testing, 1 : N-1, for testing on N groups of branches)
 * @param parameter_name {String} the name of the K parameter (internal variable name)
 * @return nothing
 */

lfunction    relax.declare_parameter_with_ranges (group_index, parameter_name) {
    parameters.DeclareGlobalWithRanges (parameter_name, 1, 1, 50);
    /** declare the variable set the initial value to 1, and the permitted range of values to [1,50] **/
}


/**
 * @name relax.impose_distributional_constraints_alternative 
 
 * Define the functional constraints on omega
 
 * @param rate_index {Number}  the index of the rate group (i.e. the 'i' in omega_i)
 * @param test_parameter {String} the name of the test omega parameter for group i
 * @param reference_parameter {String} the name of the reference omega parameter for group i
 * @param relax_parameter {String} the name of the reference rela (K) parameter for group
 * @return nothing
 */


lfunction    relax.impose_distributional_constraints_alternative (rate_index, test_parameter, reference_parameter, relax_parameter) {
    if (rate_index < ^"relax.rate_classes") { // only the last omega ratio is allowed to be >=1
        parameters.SetConstraint (test_parameter, reference_parameter + "^(1/" + relax_parameter + ")", ^"terms.global");      
        /**
            if omega < 1, then omega_T := omega_R ^ (1/K)
        */  
    } else {
        parameters.SetConstraint (test_parameter, reference_parameter + "^(" + relax_parameter + ")", ^"terms.global");
       /**
            if omega >= 1, then omega_T := omega_R ^ K
        */  
    }
}

/**
 * @name relax.impose_distributional_constraints_null
 
 * Apply the necessary constraints to define the null hypothesis
 
 * @param rate_index {Number}  the index of the rate group (i.e. the 'i' in omega_i)
 * @param test_parameter {String} the name of the test omega parameter for group i
 * @param reference_parameter {String} the name of the reference omega parameter for group i
 * @param relax_parameter {String} the name of the reference rela (K) parameter for group
 * @return nothing
 */
 
lfunction    relax.impose_distributional_constraints_null (rate_index, test_parameter, reference_parameter, relax_parameter) {
    parameters.SetConstraint (relax_parameter, ^"terms.parameters.one", ^"terms.global");
        // set K := 1
}


```


## Invokation

This analysis has three **required** arguments

- `--alignment` the alignment file and tree (in FASTA, PHYLIP, MEGA or NEXUS formats)
- `--test ` the name of the test branch set

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

mode
	Run mode
	defaut value: Classic mode
	applies to: Group test mode

reference-group [conditionally required]
	The set of branches to use as reference
	applies to: Select the set of branches to use as reference

test [conditionally required]
	Branches to use as the test set
	applies to: Choose the set of branches to use as the _test_ set

reference [conditionally required]
	Branches to use as the reference test
	applies to: Choose the set of branches to use as the _reference_ set

rates
	The number omega rate classes to include in the model [2-10, default 3]
	defaut value: relax.rate_classes [computed at run time]

models
	Which version of the test to run (All or Minimal)
	defaut value: All

parameterization
	Which parameterization should the model use
	defaut value: standard

```
 

##Example run 1

Relax **and** intensify parameterization 
> The dataset is SWS1 for HDC vr LDC bats from the RELAX paper (figure 5E)


```
hyphy RELAX.bf --alignment data/Fig4E.nex --test T --models Minimal --parameterization increase-only 
```

### Output [partial]

### Fitting the alternative model to test K != 1
* Log(L) = -4982.41, AIC-c = 10132.32 (83 estimated parameters)
* Relaxation/intensification parameter (K) =     4.58
* The following rate distribution was inferred for **test** branches

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.580     |   88.724    |                                   |
|        Negative selection         |     0.663     |    9.804    |                                   |
|      Diversifying selection       |    11.280     |    1.471    |                                   |

* The following rate distribution was inferred for **reference** branches

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.082     |   88.724    |                                   |
|        Negative selection         |     0.152     |    9.804    |                                   |
|      Diversifying selection       |     1.697     |    1.471    |                                   |

### Fitting the null (K := 1) model
* Log(L) = -5044.43, AIC-c = 10254.31 (82 estimated parameters)
* The following rate distribution for test/reference branches was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |    2.480    |                                   |
|        Negative selection         |     0.077     |   81.718    |                                   |
|      Diversifying selection       |     1.066     |   15.802    |                                   |

----
### Test for relaxation (or intensification) of selection [RELAX]
Likelihood ratio test **p =   0.0000**.

>Evidence for *intensification of selection* among **test** branches _relative_ to the **reference** branches at P<=0.05



##Example run 2

Standard parameterization 


```
hyphy RELAX.bf --alignment data/Fig4E.nex --test T --models Minimal 
```

### Output [partial]

### Fitting the alternative model to test K != 1
* Log(L) = -4984.87, AIC-c = 10137.23 (83 estimated parameters)
* Relaxation/intensification parameter (K) =     0.17
* The following rate distribution was inferred for **test** branches

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.641     |   72.126    |                                   |
|        Negative selection         |     0.754     |   26.683    |                                   |
|         Neutral evolution         |     1.000     |    1.191    |                                   |

* The following rate distribution was inferred for **reference** branches

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.068     |   72.126    |                                   |
|        Negative selection         |     0.182     |   26.683    |                                   |
|         Neutral evolution         |     1.000     |    1.191    |                                   |

### Fitting the null (K := 1) model
* Log(L) = -5044.42, AIC-c = 10254.29 (82 estimated parameters)
* The following rate distribution for test/reference branches was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.042     |   22.080    |                                   |
|        Negative selection         |     0.075     |   60.576    |                                   |
|      Diversifying selection       |     1.017     |   17.344    |                                   |

----
### Test for relaxation (or intensification) of selection [RELAX]
Likelihood ratio test **p =   0.0000**.
>Evidence for *relaxation of selection* among **test** branches _relative_ to the **reference** branches at P<=0.05
]


