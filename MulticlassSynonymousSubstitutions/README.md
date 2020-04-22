## Codon models with multiple classes of synonymous codons.

> Analysis described by Jody Hey

This analysis extends the standard Muse-Gaut 94 codon substitution models to permit **multiple classes** of synonymous substitution rates. This is accomplished by providing a tab-separated annotation file, like this [example](codons.tsv), with three columns

```
AminoAcid	CODON	GROUP
VAL	GTT	SELECTED
ALA	GCA	SELECTED
GLY	GGA	SELECTED
THR	ACT	SELECTED
ALA	GCT	SELECTED
VAL	GTA	SELECTED
PRO	CCT	SELECTED
ARG	AGA	SELECTED
LEU	TTA	SELECTED
PRO	CCG	NEUTRAL
LEU	TTG	NEUTRAL
...
```

Each codon gets a **class** in the third column of the annotation, e.g. `SELECTED` and `NEUTRAL`. One of the **classes** is designated as the **reference** and others are estimated relative to it. The rate matrix for the modified model looks like this:

```
rate (codon_i, codon_j) = 
	0                      if more that one nucleotide changes
	pi_j*nuc(i,j)*omega    one nucleotide changes, the change is non-synonymous
	pi_j*nuc(i,j)          one nucleotide changes, the change is synonymous
	                       class(codon_i) = class(codon_j) = reference
	pi_j*nuc(i,j)*alpha_k  one nucleotide changes, the change is synonymous
	                       class(codon_i) = class(codon_j) = class k
	pi_j*nuc(i,j)*alpha_kl one nucleotide changes, the change is synonymous
	                       class(codon_i) = class_k, class(codon_j) = class_l or 
	                       class(codon_i) = class_l, class(codon_j) = class_k
		            
```

Here, `nuc(i,j)` is the nucleotide rate componentm and `pi_j` is the equilibrium frequency component

For example

```
rate (GTT, GTA) = \pi_A^3 * nuc (A,T)* alpha_SELECTED
rate (TTA, TTG) = \pi_G^3 * nuc (A,G)* alpha_NEUTRAL_SELECTED
```

## Invokation

This analysis has three **required** arguments

- `--alignment` the alignment file and tree (in FASTA, PHYLIP, MEGA or NEXUS formats)
- `--classes` the TSV file with codon annotations
- `--neutral` the label to use as neutral reference

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
	default value: Universal

alignment [required]
	An in-frame codon alignment in one of the formats supported by HyPhy

tree [conditionally required]
	A phylogenetic tree
	applies to: Please select a tree file for the data:

classes [required]
	A TSV file with three columns (AA, Codon, Class) which is used to partition synonymous substitutions into groups

neutral [required]
	Neutral reference class

type
	Model type: global (single dN/dS for all branches) or local (separate dN/dS)
	default value: terms.global [computed at run time]
	applies to: Model Type

frequencies
	Equilibrium frequency estimator
	default value: CF3x4

ci
	Compute profile confidence intervals
	default value: No

lrt
	Perform LRT to test which rates are different from the neutral rate
	default value: No

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')
	default value: fitter.codon_data_info[terms.json.json] [computed at run time]

save-fit
	Save MG94 model fit to this file (default is not to save)
	default value: /dev/null
```
 

##Example run 


```
HYPHYMP MultipleSynClasses.bf --alignment adh.nex --classes codons.tsv --neutral NEUTRAL \
--ci Yes --LRT Yes
```

---
Analysis Description
--------------------
Fit an MG94xREV model where synonymous substitutions are partitioned
into several classes and within- and between-class rates are estimated.
There are several selectable options for the frequency estimators and
report the fit results including dN/dS ratios, and synonymous and
non-synonymous branch lengths. v0.2 adds the ability to compute
confidence intervals. v0.3 adds the ability to perform LRTs.

- __Requirements__: in-frame codon alignment and a phylogenetic tree, and a TSV file with
codon class partitioning

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.3


>code –> Universal
>Loaded a multiple sequence alignment with **23** sequences, **254** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/MulticlassSynonymousSubstitutions/adh.nex`

>neutral –> NEUTRAL

>type –> global

>frequencies –> CF3x4

>ci –> Yes

>lrt –> Yes


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model

>kill-zero-lengths –> Yes
* Log(L) = -5137.00, AIC-c = 10376.30 (51 estimated parameters)

### Fitting MG94 with multiple classes of synonymous substitutions
* Log(L) = -4685.54, AIC-c =  9492.34 (60 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.0792 (95% profile CI   0.0696-  0.0894)
* synonymous rate between codon classes NEUTRAL and SELECTED =   0.8389 (95% profile CI   0.7198-  0.9716)
* synonymous rate within codon class SELECTED =   0.8284 (95% profile CI   0.7246-  0.9425)

### **Synonymous tree** 
((((YAK:0.01989871958581278,(MEL:0.008738330348007717,(SIM:0.001157539284531204,MA:0.009358742393211925)Node7:0.002829554607199886)Node5:0.01025974609139345)Node3:0.01004509523547339,ERE:0.01474045858146464)Node2:0.06615206004117313,((((PSE:0.006193346202900719,PS:0)Node14:0.002056057770263325,PER:0.001011948524134472)Node13:0.005218714281742668,MIR:0.005328680474428641)Node12:0.0212174039667197,(SUB:0.04843722966640891,AMB:0.04102008182477609)Node19:0.01305501221123294)Node11:0.03548937481000465)Node1:0.07247275764812215,(MET:0.07344089154937126,(CRA:0.1330869758937619,(NIG:0.03326708268926625,(MIM:0.0180932566199034,(ADI:0.02424174598060548,(PIC:0.02411941936502098,(((SIL:0.001212563642176521,HET:0.0006400592349326339)Node36:0.004939716366866063,(PLA:0,DIF:0.003693210307501258)Node39:0.001878529644238937)Node35:0.008049306700026453,AFF:0.01400823181115291)Node34:0.003983057682403317)Node32:0.007010584781876435)Node30:0.005654878108566886)Node28:0.01220119070150073)Node26:0.03840466078404714)Node24:0.03652063406404995)Node22:0.0860041888844728,LEB:0.111098748076788)

### **Non-synonymous tree** 
((((YAK:0.00666376047233357,(MEL:0.002926325993797264,(SIM:0.0003876412497883467,MA:0.003134091988265523)Node7:0.0009475722326983959)Node5:0.003435823604889195)Node3:0.003363940492869647,ERE:0.004936342000078596)Node2:0.02215325870414798,((((PSE:0.002074051822298759,PS:0)Node14:0.0006885406088180433,PER:0.0003388852506857493)Node13:0.0017476633005007,MIR:0.001784489206054756)Node12:0.007105366617647384,(SUB:0.01622084752983242,AMB:0.01373696426330471)Node19:0.004371913175814156)Node11:0.01188481962504523)Node1:0.02426995845908043,(MET:0.02459417089874263,(CRA:0.0445686832018864,(NIG:0.01114060981152972,(MIM:0.006059140024540202,(ADI:0.008118170013365289,(PIC:0.008077204801401233,(((SIL:0.0004060680203104981,HET:0.0002143455216453786)Node36:0.001654231395548195,(PLA:0,DIF:0.001236796606787126)Node39:0.0006290893007161771)Node35:0.00269558307940421,AFF:0.004691131056341968)Node34:0.001333861821036342)Node32:0.002347731850581946)Node30:0.00189372753624807)Node28:0.004085982113644806)Node26:0.01286110191072674)Node24:0.01223017172793807)Node22:0.02880141668773867,LEB:0.03720518010051679)

### Running likelihood ratio tests to compare all rates to the neutral rate (=1)

>Testing _non-synonymous/synonymous rate ratio_ == 1

Likelihood ratio test for _non-synonymous/synonymous rate ratio == 1_, uncorrected **p =   0.0000**.

>Testing _synonymous rate between codon classes NEUTRAL and SELECTED_ == 1

Likelihood ratio test for _synonymous rate between codon classes NEUTRAL and SELECTED == 1_, uncorrected **p =   0.3166**.

>Testing _synonymous rate within codon class SELECTED_ == 1

Likelihood ratio test for _synonymous rate within codon class SELECTED == 1_, uncorrected **p =   0.2698**.
#### Holm-Bonferroni Corrected p-values
- Likelihood ratio test for _non-synonymous/synonymous rate ratio_ == 1, corrected **p =   0.0000**.
- Likelihood ratio test for _synonymous rate between codon classes NEUTRAL and SELECTED_ == 1, corrected **p =   0.3166**.
- Likelihood ratio test for _synonymous rate within codon class SELECTED_ == 1, corrected **p =   0.5395**.

### Writing detailed analysis report to `/Users/sergei/Development/hyphy-analyses/MulticlassSynonymousSubstitutions/adh.nex.FITTER.json'

--- 


