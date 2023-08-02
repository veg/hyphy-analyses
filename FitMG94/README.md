## Standard MG94 fit

This analysis fits the [Muse Gaut](https://www.ncbi.nlm.nih.gov/pubmed/7968485) (+ GTR) model of codon substitution to an alignment and a tree and reports parameter estimates and trees scaled on expected number of synonymous and non synonymous substitutions per nucleotide site.

## Invokation

This analysis has one **required** argument

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

rooted
	Accept rooted trees
	default value: No

code
	Which genetic code should be used
	default value: Universal

alignment [required]
	An in-frame codon alignment in one of the formats supported by HyPhy

tree [conditionally required]
	A phylogenetic tree
	applies to: Please select a tree file for the data:

type
	Model type: global (single dN/dS for all branches) or local (separate dN/dS)
	default value: terms.global [computed at run time]
	applies to: Model Type

frequencies
	Equilibrium frequency estimator
	default value: CF3x4

lrt
	Perform LRT to test which for dN/dS == 1 (global model only)
	default value: No

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')
	default value: fitter.codon_data_info[terms.json.json] [computed at run time]

save-fit
	Save MG94 model fit to this file (default is not to save)
	default value: /dev/null
```

### Model type

* `global` : dN/dS ratio is shared by all branches
* `local` : synonymous and non-synonymous rates are inferred separately for each branch

### Frequencies

* `CF3x4` : corrected F3x4 estimator
* `F3x4` : standard F3x4 estimator
* `F1x4` : standard F1x4 estimator 

See [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011230) for details on CF3x4 and references to other estimators. 
 

##Example run 


```
hyphy FitMG94.bf --alignment CD2.nex --lrt Yes
```

--- 

Analysis Description
--------------------
Fit an MG94xREV model with several selectable options frequency
estimator and report the fit results including dN/dS ratios, and
synonymous and non-synonymous branch lengths

- __Requirements__: in-frame codon alignment and a phylogenetic tree

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1


>code –> Universal

>alignment –> CD2.nex
>Loaded a multiple sequence alignment with **10** sequences, **187** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/FitMG94/CD2.nex`

>type –> global

>frequencies –> CF3x4

>output –> /Users/sergei/Development/hyphy-analyses/FitMG94/CD2.nex.FITTER.json


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -3532.32, AIC-c =  7112.86 (24 estimated parameters)

### Fitting Standard MG94
* Log(L) = -3467.00, AIC-c =  6997.08 (31 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.9946

### Running the likelihood ratio tests for dN/dS=1

>Testing _non-synonymous/synonymous rate ratio_ == 1

Likelihood ratio test for _non-synonymous/synonymous rate ratio == 1_, **p =   0.9441**.

### **Synonymous tree** 
((((PIG:0.04984627036053099,COW:0.06403799551220878)Node3:0.02628206419849433,HORSE:0.05462441699415268,CAT:0.07072042159004242)Node2:0.01660606657724141,((RHMONKEY:0.0009603406563057343,BABOON:0.0004580422901830288)Node9:0.006698871580291192,(HUMAN:0,CHIMP:0.0004723658790532533)Node12:0.004621767454293861)Node8:0.02833161659368989)Node1:0.07347572784488095,RAT:0.01731155959494845,MOUSE:0.03105524853252499)

### **Non-synonymous tree** 
((((PIG:0.1429749925174539,COW:0.1836813840427375)Node3:0.07538533785866487,HORSE:0.15668024015681,CAT:0.2028487121410343)Node2:0.04763149233539052,((RHMONKEY:0.002754563122905964,BABOON:0.001313811294966029)Node9:0.0192144990415551,(HUMAN:0,CHIMP:0.001354895913669316)Node12:0.01325670230522016)Node8:0.08126410744860087)Node1:0.2107518087683148,RAT:0.04965507119491638,MOUSE:0.08907635204099817)

>save-fit –> /dev/null

### Writing detailed analysis report to `/Users/sergei/Development/hyphy-analyses/FitMG94/CD2.nex.FITTER.json'

##Local model fits 

Specifying `--type local` will fit a model where each branch has its own &alpha; (dS) and &beta; (dN) rates. Further specifying `--lrt Yes` will test whether or not dN ≠ dS along a given branch. The `JSON` output file will include such branch level parameters in the corresponding dictionary:

```
"branch attributes":{
   "0":{
     "BABOON":{
       "Confidence Intervals":{
         "LB":0,
         "MLE":0,
         "UB":0.6860409126311118
        },
       "LRT":{
         "Corrected P-value":1,
         "FDR":0.5456841416635031,
         "LRT":2.669070286591705,
         "p-value":0.1023157765619068
        },
       "Nucleotide GTR":0.001679681026555922,
       "Standard MG94":0.001815827116959503,
       "dN":1.267499201636807e-10,
       "dS":0.008603986131448501,
       "nonsynonymous":1e-10,
       "original name":"BABOON",
       "synonymous":0.001815827116959503
      },
```

### Example invokation

```
hyphy FitMG94.bf --alignment CD2.nex --type local --lrt Yes
```

Analysis Description
--------------------
Fit an MG94xREV model with several selectable options frequency
estimator and report the fit results including dN/dS ratios, and
synonymous and non-synonymous branch lengths. v0.2 adds LRT test for
dN/dS != 1. v0.3 adds LRT test support for dN/dS != 1 for local models

- __Requirements__: in-frame codon alignment and a phylogenetic tree

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.3

rooted: No

>code –> Universal
>Loaded a multiple sequence alignment with **10** sequences, **187** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/FitMG94/CD2.nex`

>type –> local

>frequencies –> CF3x4

>lrt –> Yes


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model

>kill-zero-lengths –> Yes
* Log(L) = -3532.32, AIC-c =  7112.86 (24 estimated parameters)
* 1 partition. Total tree length by partition (subs/site)  1.694

### Fitting Standard MG94
* Log(L) = -3450.61, AIC-c =  6995.58 (46 estimated parameters)

### Running the likelihood ratio tests for dN/dS=1 and estimating confidence intervals for dN/dS along each branch

|            Branch            |     Length     |     dN/dS      |Approximate dN/dS CI|LRT p-value dN != dS|
|:----------------------------:|:--------------:|:--------------:|:------------------:|:------------------:|
|             PIG              |     0.192      |     1.345      |   0.966 - 1.809    |       0.3794       |
|             COW              |     0.253      |     1.918      |   1.456 - 2.479    |       0.0516       |
|            Node5             |     0.103      |     1.482      |   0.836 - 2.304    |       0.4878       |
|            HORSE             |     0.209      |     1.245      |   0.926 - 1.635    |       0.4828       |
|             CAT              |     0.277      |     1.604      |   1.241 - 2.038    |       0.1190       |
|            Node4             |     0.066      |     0.664      |   0.275 - 1.202    |       0.5664       |
|           RHMONKEY           |     0.004      |10000000000.0...|0.000 - 10000.000...|       0.2845       |
|            BABOON            |     0.002      |     0.000      |   0.000 - 0.686    |       0.1023       |
|            Node11            |     0.026      |     0.401      |   0.139 - 0.817    |       0.1605       |
|            HUMAN             |     0.000      |     1.000      |0.000 - 10000.000...|       1.0000       |
|            CHIMP             |     0.002      |10000000000.0...|0.000 - 10000.000...|       0.4355       |
|            Node14            |     0.018      |     0.368      |   0.052 - 0.924    |       0.2652       |
|            Node10            |     0.110      |     1.915      |   1.270 - 2.728    |       0.2750       |
|            Node3             |     0.290      |     0.432      |   0.317 - 0.573    |       0.0015       |
|             RAT              |     0.066      |     1.089      |   0.610 - 1.717    |       0.9070       |
|            MOUSE             |     0.122      |     0.525      |   0.343 - 0.755    |       0.1057       |

### **Synonymous tree** 
((((PIG:0.03937408670808301,COW:0.03878451353271926)Node5:0.0195616708158665,HORSE:0.04559327378804292,CAT:0.04931224992892905)Node4:0.02284730872773409,((RHMONKEY:0,BABOON:0.001815827116959503)Node11:0.01202579316224195,(HUMAN:0,CHIMP:0)Node14:0.008700053742021586)Node10:0.01688166070650646)Node3:0.1292228327513776,RAT:0.01590572225225523,MOUSE:0.04854724584993934)

### **Non-synonymous tree** 
((((PIG:0.1521831461112254,COW:0.2138117251333718)Node5:0.08333270489251039,HORSE:0.1632014620903291,CAT:0.2273175155627533)Node4:0.04358798007485778,((RHMONKEY:0.003666894925291958,BABOON:0)Node11:0.01385522359919074,(HUMAN:0,CHIMP:0.001837936132619463)Node14:0.00920920510293554)Node10:0.09291437133059131)Node3:0.1603787381398532,RAT:0.04977984176707917,MOUSE:0.07324864703524042)
**Combined tree** 
((((PIG:0.1915572328193086,COW:0.2525962386660908)Node5:0.102894375708377,HORSE:0.2087947358783719,CAT:0.2766297654916826)Node4:0.06643528880259186,((RHMONKEY:0.003666894925291958,BABOON:0.001815827116959503)Node11:0.02588101676143268,(HUMAN:0,CHIMP:0.001837936132619463)Node14:0.01790925884495711)Node10:0.1097960320370977)Node3:0.2896015708912303,RAT:0.06568556401933438,MOUSE:0.1217958928851798)

### Writing detailed analysis report to `/Users/sergei/Development/hyphy-analyses/FitMG94/CD2.nex.FITTER.json'