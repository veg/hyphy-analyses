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
HYPHYMP FitMG94.bf --alignment CD2.nex --lrt Yes
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


