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
Analysis options description
----------------------------
code
	Which genetic code should be used
	defaut value: Universal

alignment [required]
	An in-frame codon alignment in one of the formats supported by HyPhy

tree
	A phylogenetic tree
	defaut value: null [computed at run time]

classes [required]
	A TSV file with three columns (AA, Codon, Class) which is used to 
	partition synonymous substitutions into groups

neutral [required]
	Neutral reference class

type
	Model type: global (rates shared by all branches) or 
                local (separate rates for each branch)
	defaut value: terms.global [computed at run time]

frequencies
	Equilibrium frequency estimator
	defaut value: CF3x4

output
	Write the resulting JSON to this file (default is to save to 
	the same path as the alignment file + 'MG94.json')
	defaut value: fitter.codon_data_info[terms.json.json] [computed at run time]

save-fit
	Save MG94 model fit to this file (default is not to save)
	defaut value: /dev/null
```
 

##Example run 


```
HYPHYMP MultipleSynClasses.bf --alignment adh.nex --classes codons.tsv --neutral NEUTRAL
```

---

Analysis Description
--------------------
Fit an MG94xREV model where synonymous substitutions are partitioned
into several classes and within- and between-class rates are estimated.
There are several selectable options for the frequency estimators and
report the fit results including dN/dS ratios, and synonymous and
non-synonymous branch lengths

- __Requirements__: in-frame codon alignment and a phylogenetic tree, and a TSV file with
codon class partitioning

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1


>code –> Universal

>alignment –> adh.nex
>Loaded a multiple sequence alignment with **23** sequences, **254** codons, and **1** partitions from `/Users/sergei/Development/hyphy-analyses/MulticlassSynonymousSubstitutions/adh.nex`

>classes –> codons.tsv

>neutral –> NEUTRAL

>type –> global

>frequencies –> CF3x4

>output –> /Users/sergei/Development/hyphy-analyses/MulticlassSynonymousSubstitutions/adh.nex.FITTER.json


### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model
* Log(L) = -5137.11, AIC-c = 10376.52 (51 estimated parameters)

### Fitting MG94 with multiple classes of synonymous substitutions
* Log(L) = -4685.53, AIC-c =  9492.34 (60 estimated parameters)
* non-synonymous/synonymous rate ratio =   0.0789
* synonymous rate between codon classes NEUTRAL and SELECTED =   0.8355
* synonymous rate within codon class SELECTED =   0.8253

### **Synonymous tree** 
((((YAK:0.01991098348831301,(MEL:0.008737313343428603,(SIM:0.00115800211749067,MA:0.009363889194226845)Node7:0.002830832078527877)Node5:0.01026265900442833)Node3:0.01005223647747053,ERE:0.01475390925283582)Node2:0.06608609536439869,((((PSE:0.006194011946607071,PS:0)Node14:0.002056542136941789,PER:0.001012340881462104)Node13:0.005218962073655949,MIR:0.005329521652100825)Node12:0.02123137422148275,(SUB:0.04844573622631574,AMB:0.04104433512857298)Node19:0.01305359234869906)Node11:0.03554321874784543)Node1:0.0724923496394842,(MET:0.07340774741650495,(CRA:0.1331290565207244,(NIG:0.03328724078545436,(MIM:0.01809211497090321,(ADI:0.02423594093834776,(PIC:0.02413529390977663,(((SIL:0.001212075571806737,HET:0.0006404690372100284)Node36:0.004940060050832476,(PLA:0,DIF:0.00369271512384738)Node39:0.001876907364173814)Node35:0.008049856513402936,AFF:0.01400627884493204)Node34:0.003982039851475526)Node32:0.007013058944921452)Node30:0.005659775889049311)Node28:0.01219861408601335)Node26:0.03841424975855891)Node24:0.03650346975555927)Node22:0.08609632984879025,LEB:0.1110680623916397)

### **Non-synonymous tree** 
((((YAK:0.006666178593080257,(MEL:0.002925244310768755,(SIM:0.0003876980225958669,MA:0.00313502131781534)Node7:0.0009477599242430161)Node5:0.003435928607120699)Node3:0.003365479342495905,ERE:0.004939594976974992)Node2:0.02212556273159468,((((PSE:0.00207374938902397,PS:0)Node14:0.0006885283781735687,PER:0.0003389307774206013)Node13:0.001747303606269564,MIR:0.001784318849415227)Node12:0.007108244171866166,(SUB:0.01621958704086302,AMB:0.01374160489670956)Node19:0.004370330472564053)Node11:0.01189983629312581)Node1:0.02427037065308534,(MET:0.02457684497008286,(CRA:0.0445714832326756,(NIG:0.01114453698768573,(MIM:0.006057223119760644,(ADI:0.008114170290041252,(PIC:0.008080473759294198,(((SIL:0.0004058017643737288,HET:0.0002144284328238349)Node36:0.001653927470670641,(PLA:0,DIF:0.001236317559269934)Node39:0.0006283868247690583)Node35:0.002695084408990231,AFF:0.004689289017772657)Node34:0.001333183206659759)Node32:0.002347965555698744)Node30:0.001894887658128883)Node28:0.004084084552285427)Node26:0.01286105478215332)Node24:0.01222132742968185)Node22:0.02882497046505308,LEB:0.03718548309402383)

>save-fit –> /dev/null

### Writing detailed analysis report to `/Users/sergei/Development/hyphy-analyses/MulticlassSynonymousSubstitutions/adh.nex.FITTER.json'

--- 


