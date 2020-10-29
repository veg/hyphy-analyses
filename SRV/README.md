### Modeling synonymous rate variation in selection tests (BUSTED and RELAX)

Starting with `2.5.21` release of HyPhy, BUSTED and RELAX analyses support several options for modeling site-to-site variation in synonymous substitutions rates (**SRV**), which is controlled by the `--srv option` command line argument. The following values for `option` are supported (they are case sensitive). Setting different options alters the model being fitted, screen output, and `JSON` files.

1. `No` Synonymous rates are assumed to be constant across sites in the alignment, so the rate at site **s** is &alpha;<sub>s</sub> = 1. This is the default setting for both [BUSTED](https://academic.oup.com/mbe/article/32/5/1365/1134918) and [RELAX](https://academic.oup.com/mbe/article/32/3/820/981440), and corresponds to the original implementations of the methods.

2. `Yes`  Synonymous rates are assumed to be constant across branches at each site but vary from site to site in an indepdentent manner according to a K rate (K is set via the `--syn-rates` argument, default is 3) unit-mean general discrete distribution, &alpha;<sub>s</sub> ~ GDD (weights, rates). This is the model implemeted in [BUSTED[S]](https://academic.oup.com/mbe/article/37/8/2430/5739973). Biological intepretation of this model is that there is some selection operating on synonymous changes, but it does not vary in time, and has no site-to-site dependance. 

3. `HMM`  Synonymous rates are assumed to be constant across branches at each site, vary from site to site in an indepdentent manner according to a K rate (K is set via the `--syn-rates` argument, default is 3) unit-mean general discrete distribution, &alpha;<sub>s</sub> ~ GDD (weights, rates), and are spatially autocorrelated with the [Felsenstein-Churchill Hidden Markov Model](https://academic.oup.com/mbe/article/13/1/93/1055515); with the switching rate given by &lambda;. Biological intepretation of this model is that there is some selection operating on synonymous changes, but it does not vary in time, and adjacent sites tend to have the same synonymous rate (depending on &lambda;).

4. `Branch-site` Synonymous rates are assumed to be vary across branches **and** sites using the same unrestricted branch site model as the one used for &omega; ratios.  At each branch and site, the synonymous rate is drawn (independently from other sites and branches) from a K rate (K is set via the `--syn-rates` argument, default is 3) unit-mean general discrete distribution, &alpha;<sub>b,s</sub> ~ GDD (weights, rates). Biological intepretation of this model is that there is some selection operating on synonymous changes, and it varies both between sites and branches; this is also the model most amenable to handling heterotachy. 

#### Differences in screen output. 

For all options except `No`, a report will be written to the screen identifying the GDD distribution inferred for synonymous rates. For example:

---

* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.421               |    52.820     |                                   |
|               1.269               |    43.803     |                                   |
|               6.574               |     3.377     |                                   |

---

For the `HMM` option, there will be two additional outputs: the switching rate &lambda;, and the collection of rate switch-points inferred using the Viterbi algorithm - these are sites where &alpha;<sub>s</sub> switch from one rate-category to another, in the most likely rate _path_. For example:

---

* HMM switching rate =   0.1734 (95% profile CI   0.0992-  0.2652)
* The following rate distribution for site-to-site **synonymous** rate variation was inferred

|               Rate                | Proportion, % |               Notes               |
|-----------------------------------|---------------|-----------------------------------|
|               0.582               |    67.827     |                                   |
|               1.106               |    17.064     |                                   |
|               2.758               |    15.109     |                                   |

* The following switch points for synonymous rates were inferred

|               Site                |  From class   |   To class    |
|-----------------------------------|---------------|---------------|
|                 6                 |       2       |       1       |
|                12                 |       1       |       2       |
|                16                 |       2       |       1       |
|                18                 |       1       |       2       |
|                24                 |       2       |       0       |
|                34                 |       0       |       2       |
|                35                 |       2       |       0       |
|                41                 |       0       |       1       |

...

---

#### Differences in `JSON` output. 

For settings other than `No`, the output `JSON` file will contain the following records **in addition** to the standard fields for BUSTED or RELAX.

1. `:"Synonymous site-posteriors"` For each site, records the estimated posterior probabilities (estimated via an Empricial Bayes procedure) for each synonymous rate class as an K x N (N = number of sites) matrix. `(i,j)` &rightarrow; posterior probability that &alpha;<sub>j</sub> == i-th rate.

2. `:fits:*:"Rate Distributions":"Synonymous site-to-site rates"` The inferred rate distribution for ,e.g.,
```
{
 "0":{
   "proportion":0.5281294386831287,
   "rate":0.4209048453725299
  },
 "1":{
   "proportion":0.4382271939443426,
   "rate":1.270180968926857
  },
 "2":{
   "proportion":0.03364336737252874,
   "rate":6.571277958359039
  }
}
```
3. For the `HMM` model, there is also `:"Viterbi synonymous rate path"`, which reports the most likely _path_ of rate assignements for each site, and `:fits:*:"Rate Distributions":"HMM rate switching parameter"` for the estimated value of &lambda; and the associated 95% profile confidence interval, like so

```
"HMM rate switching parameter":{
 "LB":0.09923831400908656,
 "MLE":0.1733505229099538,
 "UB":0.265216407206666
}
```

#### An example with intepretation [BUSTED]

Consider BUSTED with different SRV options applied to an alignment of **17** &beta;-globin sequences (originally analyzed in [Yang et al](https://www.genetics.org/content/155/1/431)). 

```
hyphy busted --alignment tests/bglobin.nex --starting-points 10 --srv [option] \
--output results/bglobin.[option].busted.json 
```

##### At-a-glance model comparison 

* HMM switching rate =   0.1734 (95% profile CI   0.0992-  0.2652)

| SRV option |  Log L   |   # parameters    |   c-AIC    | p-value | &omega;<sub>3</sub> (%) | &alpha;<sub>1</sub> (%) | &alpha;<sub>2</sub> (%)| &alpha;<sub>3</sub> (%) |
|-------|---------------|---------------|---------------|---------------|---------------|---------------|--------------|--------------|
| `No`  |  -3676.68     | 50            | 7455.49       | <0.0001      | 11.0 (4.4)    |  N/A | N/A |N/A |
| `Yes`  |  -3651.50     | 55            | 7415.58     | <0.0001      | 8.8 (4.4)    |  0.351 (43.682) | 1.269              (53.250) |  9.162  (3.069) |
| `HMM`  |  -3651.73    | 56            | 7418.13     | <0.0001      | 7.3 (3.2)    |  0.582 (67.827) | 1.106               (17.064) |  2.758  (15.109) |
| `Branch-site`  |  -3672.05    | 55            | 7456.68     | 0.0004     | 7.6 (7.1)    |  0.299 (65.609) | 0.301               (5.291) |  2.707  (29.100) |

##### Observations.

1. All models find strong evidence for episodic positive diversifying selection in this alignment at ~10% of sites with &omega; ~ 10, p < 0.001.
2. There's strong evidence of SRV (`Yes` vs `No` model, based on c-AIC). 
3. The standard SRV (`Yes`) model provides the best fit to the data using c-AIC.
4. There's some evidence of spatial rate auto-correlation (&lambda; ≠ 1/3), although there's no explicit test for it conducted here, but the reported confidence interval (which is too narrow, generally) does not include 1/3. A proper test would conduct a likelihood ratio test.
6. The estimates of the SRV distribution (&alpha;) are strongly influenced by which model is being selected.

The following are the lists of the sites with `Constrained Test Statistic` > 10 which provide an indication of which sites support models that permit positive selection.

1. `No`: 7, 10, 11, 42, 48, 50, 54, 67, 85, 110, 123, 124
2. `Yes`: 42, 48, 54, 110
3. `HMM`: 7, 10, 11, 42, 48, 50, 54, 67, 85, 110, 123, 124
4. `Branch-site`: 7, 11, 42, 48, 50, 54, 85, 110, 123

The standard `SRV` model has the fewest sites in this category (but it selects a subset of sites from other models), while the other three models largely agree.

#### An example with intepretation [RELAX]

Consider RELAX with different SRV options applied to an alignment of 33 *Bat* SWS1 gene sequences (originally analyzed in [Wertheim et al](https://academic.oup.com/mbe/article/32/3/820/981440)). 

```
./hyphy relax --alignment tests/SWS1.nex --test T --reference R --models Minimal \
--starting-points 10 --srv [option] --output results/SWS1.[option].relax.json
```

* HMM switching rate =   0.2995 (95% profile CI   0.2226-  0.3333)

| SRV option |  Log L   |   # parameters    |   c-AIC    | p-value | K | &alpha;<sub>1</sub> (%) | &alpha;<sub>2</sub> (%)| &alpha;<sub>3</sub> (%) |
|-------|---------------|---------------|---------------|---------------|---------------|---------------|--------------|--------------|
| `No`  |  -4982.88     | 81            | 10129.17      | <0.0001      | 0.00   |  N/A | N/A |N/A |
| `Yes`  |  -4958.51     | 86           | 10090.63      | <0.0001      | 0.18   | 0.407 (33.026)   | 1.180 (62.858)              |  3.170 (3.437) |
| `HMM`  |  -4961.15     | 87           | 10097.94     | <0.0001      | 0.17    | 0.856 (92.432)   | 1.675 (1.868)              |  3.113 (5.700) |
| `Branch-site`  |  -4982.97    | 86          | 10139.55     | <0.0001      | 0.17    |  0.908 (98.91) | 9.361 (1.09) | N/A |  

##### Observations.

1. All models agree that there is strong evidence for relaxation (K < 1, p < 0.0001). 
2. The `Yes` model (standard SRV) provides the best fit to the data (AIC-c)
3. There's no evidence of auto-correlation in SRV (HMM model **not** preferred over the standard model, and the CI for &lambda; overlaps 1/3).
4. The `Branch-site` model provides a poor fit to this alignment. 
6. The estimates of the SRV distribution (&alpha;) are strongly influenced by which model is being selected.
