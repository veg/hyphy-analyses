## NRM -- examine nucleotide data for non-reversibility

> Original code by Wayne Delport; project with Darren Matrin and Rita Sianga

Given an input (rooted) Newick tree, and a nucleotide sequence alignment, fit **4** nucleotide models, compare their fits and estimate equilibrim nucleotide frequencies.

#### Model 1: GTR + &Gamma;4

This is the standard 6 parameter general time reversible nucleotide model and a 4-bin discretized gamma distribution site-to-site rate variation. The equilibritum nucleotide frequencies (&pi;<sub>A</sub>,&pi;<sub>C</sub>,&pi;<sub>G</sub>,&pi;<sub>T</sub>) are estimated by counts from the data.

**Q** matrix

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   r<sub>AC</sub> &pi;<sub>C</sub>|   &pi;<sub>G</sub> |  r<sub>AT</sub>&pi;<sub>T</sub> |   &pi;<sub>A</sub> |
|     C    |     r<sub>AC</sub> &pi;<sub>A</sub>| *  |   r<sub>CG</sub>   &pi;<sub>G</sub> |  r<sub>CT</sub>&pi;<sub>T</sub> |   &pi;<sub>C</sub> |
|     G    |     &pi;<sub>A</sub>|   r<sub>CG</sub>   &pi;<sub>C</sub> |  * | r<sub>GT</sub>&pi;<sub>T</sub> |   &pi;<sub>C</sub> |
|     T    |     r<sub>AT</sub> &pi;<sub>A</sub>|  r<sub>CT</sub>   &pi;<sub>C</sub> |  r<sub>GT</sub>&pi;<sub>G</sub> |   * | &pi;<sub>T</sub> |

#### Model 1: stGTR + &Gamma;4

This is the "paired strand" 6 parameter non-reversible nucleotide model and a 4-bin discretized gamma distribution site-to-site rate variation. The nucleotide frequencies at the root of the tree are estimated by counts from the data.

**Q** matrix

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   r<sub>AC</sub> &pi;<sub>C</sub>|   &pi;<sub>G</sub> |  r<sub>AT</sub>&pi;<sub>T</sub> |   Induced by the model |
|     C    |     r<sub>GT</sub> &pi;<sub>A</sub>| *  |   r<sub>CG</sub>   &pi;<sub>G</sub> |  r<sub>CT</sub>&pi;<sub>T</sub> |  " |
|     G    |     r<sub>CT</sub>&pi;<sub>A</sub>|   r<sub>CG</sub>   &pi;<sub>C</sub> |  * | r<sub>GT</sub>&pi;<sub>T</sub> |   "|
|     T    |     r<sub>AT</sub> &pi;<sub>A</sub>|  &pi;<sub>C</sub> |  r<sub>AC</sub>&pi;<sub>G</sub> |   * | "|

#### Model 1: NRM + &Gamma;4

This is the 12 parameter parameter non-reversible nucleotide model and a 4-bin discretized gamma distribution site-to-site rate variation. The nucleotide frequencies at the root of the tree are estimated by counts from the data.
**Q** matrix

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   r<sub>AC</sub>|   1|  r<sub>AT</sub> |   Induced by the model |
|     C    |     r<sub>CA</sub> | *  |   r<sub>CG</sub>    |  r<sub>CT</sub>|  " |
|     G    |     r<sub>GA</sub>|   r<sub>GC</sub>    |  * | r<sub>GT</sub> |   "|
|     T    |     r<sub>TA</sub> | r<sub>TC</sub>  |  r<sub>TG</sub>|   * | "|

#### Model 1: NRMF + &Gamma;4

This is the 12+3 parameter parameter non-reversible nucleotide model and a 4-bin discretized gamma distribution site-to-site rate variation. The nucleotide frequencies at the root of the tree are estimated by maximum likelihood.
**Q** matrix

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   r<sub>AC</sub>|   1|  r<sub>AT</sub> |   Induced by the model |
|     C    |     r<sub>CA</sub> | *  |   r<sub>CG</sub>    |  r<sub>CT</sub>|  " |
|     G    |     r<sub>GA</sub>|   r<sub>GC</sub>    |  * | r<sub>GT</sub> |   "|
|     T    |     r<sub>TA</sub> | r<sub>TC</sub>  |  r<sub>TG</sub>|   * | "|

## Invokation

This analysis has one **required** argument

- `--alignment` the file containing the alignment 

Detailed analysis results written to the `--output` file (default == `alignment file.NRM.json`)

### Complete options list



```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

alignment [required]
	Sequence alignment to screen for recombination

tree [conditionally required]
	A phylogenetic tree (optionally annotated with {})
	applies to: Please select a tree file for the data:

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'NRM.json')
	default value: nrm.nuc_data_info[terms.json.json] [computed at run time]
```


## Example run

```
hyphy NRM.bf --alignment test.nex
```

Analysis Description
--------------------
 Perform a fit of GTR, strand non-reversible, and fully non-reversible
model (with Gamma rate variation) to a nucleotide alignment. Report
estimated rate matrices, and perform nested model fits. 

- __Requirements__: a nucleotide alignment and a phylogenetic tree (rooted)

- __Citation__: TBD

- __Written by__: Wayne Delport and Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1

>Loaded a multiple sequence alignment with **16** sequences, **864** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/NucleotideNonREV/test.nex`


### Fitting the GTR + G model with empirical base frequencies

>Log(L) = -2078.21, AIC-c =  4226.60 (35 estimated parameters)
1 partition. Total tree length by partition (subs/site)  0.236
_Gamma shape parameter_ =   0.0955

#### Rate and equilibrium frequency estimates.

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   0.3489 |   1.0000 |   0.0394 |   0.3951 |
|     C    |   0.3489 |     *    |   0.0569 |   0.4478 |   0.1802 |
|     G    |   1.0000 |   0.0569 |     *    |   0.4201 |   0.1808 |
|     T    |   0.0394 |   0.4478 |   0.4201 |     *    |   0.2439 |

### Fitting the strandGTR + G model with empirical base frequencies

>Log(L) = -2082.52, AIC-c =  4235.23 (35 estimated parameters)
1 partition. Total tree length by partition (subs/site)  0.233
_Gamma shape parameter_ =   0.0765

#### Rate and equilibrium frequency estimates.

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   0.5995 |   1.0000 |   0.0578 |   0.4001 |
|     C    |   0.4982 |     *    |   0.0795 |   1.2234 |   0.1795 |
|     G    |   1.2234 |   0.0795 |     *    |   0.4982 |   0.1656 |
|     T    |   0.0578 |   1.0000 |   0.5995 |     *    |   0.2547 |

### Fitting the NRM + G model with empirical root frequencies

>Log(L) = -2075.79, AIC-c =  4233.82 (41 estimated parameters)
1 partition. Total tree length by partition (subs/site)  0.238
_Gamma shape parameter_ =   0.1042

#### Rate and equilibrium frequency estimates.

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   0.3516 |   1.0000 |   0.0805 |   0.4645 |
|     C    |   0.9882 |     *    |   0.0000 |   0.8712 |   0.1465 |
|     G    |   2.8400 |   0.1423 |     *    |   0.4051 |   0.1780 |
|     T    |   0.0704 |   0.3966 |   0.6565 |     *    |   0.2110 |

### Fitting the NRM + G model with estimated root frequencies


### 
>Log(L) = -2075.73, AIC-c =  4239.76 (44 estimated parameters)
1 partition. Total tree length by partition (subs/site)  0.238
_Gamma shape parameter_ =   0.1038

#### Estimates for root frequencies
- A : 0.4001362558802573
- C : 0.1789561515326505
- G : 0.1799435082820515
- T : 0.2409640843050408

#### Rate and equilibrium frequency estimates.

| From/To  |     A    |     C    |     G    |     T    | Frequency| 
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|     A    |     *    |   0.3551 |   1.0000 |   0.0814 |   0.4653 |
|     C    |   0.9963 |     *    |   0.0000 |   0.8806 |   0.1466 |
|     G    |   2.8620 |   0.1435 |     *    |   0.4083 |   0.1773 |
|     T    |   0.0711 |   0.4007 |   0.6632 |     *    |   0.2109 |

### Model comparison results

|   Null   |   Alt.   |      LRT      | Deg. freedom  |       p       |    Delta c-AIC     |
|:--------:|:--------:|:-------------:|:-------------:|:-------------:|:------------------:|
|   GTR    |  stGTR   |     null      |     null      |     null      |        -8.6313     |
|   GTR    |   NRM    |      4.8418   |       6       |      0.5643   |        -7.2200     |
|  stGTR   |   NRM    |     13.4731   |       6       |      0.0361   |         1.4113     |
|   NRM    |   NRMF   |      0.1031   |       3       |      0.9915   |        -5.9318     |