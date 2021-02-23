## Are dN/dS ratios different between genes?


This analysis builds upon [BUSTED](https://academic.oup.com/mbe/article/32/5/1365/1134918) which is a method to test for evidence of epidosic diversifying selection. BUSTED infers a distribution of &omega; (dN/dS) ratios for each gene it analyzes. If we BUSTED to two or more genes, we could ask whether or not these distributions differ between the genes in a statistically significant way. Here's an example of two such distributions (example files included in `examples`).

#### Gene 1

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   32.430    |                                   |
|        Negative selection         |     0.000     |    5.170    |       Collapsed rate class        |
|      Diversifying selection       |     1.723     |   62.400    |                                   |


#### Gene 2

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.307     |   68.513    |                                   |
|        Negative selection         |     0.820     |    0.000    |       Not supported by data       |
|      Diversifying selection       |     3.369     |   31.487    |                                   |

They **look** different, e.g. in gene 2 a smaller fraction of the alignment (31.5% vs 62.4%) appears to be under stronger selection (dN/dS ~ 3.4 vs dN/dS ~ 1.7) compared to gene 1. But **are they different**? (turns out, no, see below). A statistically principled test for **any** distributional difference would compare

1. Alternative: the distributions of &omega; are different.
2. Null: the distributions of &omega; are exactly the same.

Signififcance can be assessed by a likelihood ratio test with 2N-1 degrees of freedom, where **N** is the number of dN/dS rate categories (3 here).

This is what `BUSTED-compare.bf` does.

## Invokation

This analysis requires that you **first** perform two independent BUSTED analysies **and** save the resulting model fits. For example

```
hyphy busted --alignment examples/model1_replicate1.fas \
			 --save-fit  examples/model1_replicate1.fit
```

This analysis has three **required** argument

- `--fit1` the file with the first BUSTED fit
- `--fit2` the file with the second BUSTED fit
- `--output` save a `JSON` file with analysis results here.


HyPhy will write Markdown output to the screen and a JSON file with detailed fit results. 

See example at the end of the document

### Complete options list 

```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

fit1 [required]
	First fit file

fit2 [required]
	Second fit file

output [required]
	Write comparison JSON to

```
 

##Example positive control run 

Compare fits from two datasets simulated under different models (expected positive result)

```
hyphy BUSTED-compare.bf --fit1 examples/model1-replicate1.fit --fit2 examples/model2-replicate1.fit --output results/model12.json
```

---

Analysis Description
--------------------
 Read two fit files from BUSTED analyses and conduct an LRT test to
determine if the inferred dN/dS rate distributions are different between
them 

- __Requirements__: Two BUSTED fit files

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1



### Loaded distributions report

>dN/dS distribution for /Users/sergei/Development/hyphy-analyses/BUSTED-compare/examples/model1-replicate1.fit

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   32.430    |                                   |
|        Negative selection         |     0.000     |    5.170    |       Collapsed rate class        |
|      Diversifying selection       |     1.723     |   62.400    |                                   |


>dN/dS distribution for /Users/sergei/Development/hyphy-analyses/BUSTED-compare/examples/model2-replicate1.fit

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.000     |   40.834    |                                   |
|        Negative selection         |     0.706     |   22.480    |                                   |
|      Diversifying selection       |    15.896     |   36.686    |                                   |


### Fitting the constrained model

- Independent model likelihood  -7907.308
- Degrees of freedom = 5
Current Max: -7934.301      (99 % done) LF Evals/Sec: 59.51   CPU Load: 2.559   
### Model fit results

- Constrained model likelihood  -7934.301
- p-value 2.110283059408857e-10

Inferred shared dN/dS distribution

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.018     |   30.485    |                                   |
|        Negative selection         |     0.044     |   21.955    |                                   |
|      Diversifying selection       |     4.069     |   47.560    |                                   |

### Example `JSON` output

```
{
 "distributions":{
   "0":    [
[0, 0.3243033731484445],
    [0, 0.05169799679720603],
    [1.72311268758951, 0.6239986300543495] 
    ],
   "1":    [
[0.306610934700062, 0.6851311276395798],
    [0.8203718498689812, 0],
    [3.369098635615948, 0.3148688723604202] 
    ],
   "shared":    [
[0.1168235325171251, 0.249566428272481],
    [0.1263383582567932, 0.2512010138826047],
    [2.147400118054277, 0.4992325578449143] 
    ]
  },
 "files":  [
["/Users/sergei/Development/hyphy-analyses/BUSTED-compare/examples/model1-replicate1.fit", "/Users/sergei/Development/hyphy-analyses/BUSTED-compare/examples/model1-replicate2.fit"] 
  ],
 "fits":{
   "constrained":-6963.464144056938,
   "unconstrained":-6961.841970096701
  },
 "test":{
   "LRT":3.244347920473956,
   "p-value":0.6623721626107559
  }
}

```

