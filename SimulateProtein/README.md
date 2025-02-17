## Simulate codon data under the MG94 class of models.

This analysis generates sequence [Muse Gaut](https://www.ncbi.nlm.nih.gov/pubmed/7968485) (+ GTR) model of codon substitution and an (optionally annotated) tree with branch lengths to generate alignments of parametrically simulated codon sequences as FASTA. How synonymous and non-synonymous rates vary across branches and sites is specified via **plug-in scripts** located in the `modules` directory. This analysis has many options and is a work in progress. We begin by providing several usage examples, based on the trees supplied with this script.

## Examples

#### Generate data from a phylogeny with specified branch length using a constant omega 

Use &omega; = 0.25 for the entire alignment, generate 10 replicates with 300 sites each, write to `data/example1.replicate.N`. 

Important default settings: equal base frequencies, CF3x4 frequency estimator, HKY85 model with &kappa; = 2

```
hyphy SimulateMG94.bf --tree CD2.nwk --omega 0.25 --sites 300 --replicates 10 --output data/example1
```

#### Generate data from a phylogeny with specified branch length using a constant omega with more customization

Use &omega; = 0.25 for the entire alignment, generate 5 replicates with 200 sites each, write to `data/example2.replicate.N`. Use F3x4 frequency estimator, TrN 93 model with A-C rates (relative to AG) = 0.1 , and relative C-T rate = 2, and custom base frequencies:

| Base | Codon position 1 | Codon position 2 | Codon position 3 |
|:----:|:----------------:|:----------------:|:----------------:|
|   A  |		0.3			   |       0.2        |      0.4         |
|   C  |		0.15		   |       0.31       |      0.05        |
|   G  |		0.23			|       0.19      |      0.11         |
|   T  |    auto          |    auto          |    auto          |


(new lines added for formatting only)

```
hyphy SimulateMG94.bf 
--tree CD2.nwk 
--omega 0.25 
--sites 25 
--replicates 5 
--output data/example2 
--frequency-estimator F3x4 
--base-frequencies 0.3,0.15,0.23,0.2,0.31,0.19,0.4,0.05,0.11 
--CT 2 --AC 0.1 
```

#### Generate data from a phylogeny with grouped branches, where each branch group gets its own &omega;

Assume you have a labeled tree, using `{group}` to assign groups to branches, like in `CD2.nwk`:

```
((((Pig{Others}:0.147969,Cow{Others}:0.21343){Others}:0.085099,Horse{Others}:0.165787,Cat{Others}:0.264806){Others}:0.058611,((RhMonkey{Primates}:0.002015,Baboon{Primates}:0.003108){Primates}:0.022733,(Human{Primates}:0.004349,Chimp{Primates}:0.000799){Primates}:0.011873){Primates}:0.101856){Others}:0.340802,Rat{Others}:0.050958,Mouse{Others}:0.09795)
```

Then we can ask HyPhy to generate data where branches in different groups are evolved using different &omega;, here &omega; = 1.5 for Primates and &omega; = 0.25 for Others.


Important default settings: CF3x4 frequency estimator, HKY85 model with &kappa; = 2

(new lines added for formatting only)

```
hyphy SimulateMG94.bf 
--tree CD2.nwk 
--omega-Primates 1.5 
--omega-Others 0.25 
--sites 400 
--replicates 1 
--output data/example3 
--branch-variation partitioned 
--base-frequencies HIV
```

#### Generate data from a phylogeny with branches annotated with synonymous and non-synonymous lengths

Assume you have a labeled tree with comment blocks used to associate values with branches, like in `CD2.dnds.nwk`:

```
((((PIG[s=0.04985418650068867;ns=0.1429621849403452]:0.1928163714410339,
COW[s=0.06406199560179504;ns=0.1837045894379381]:0.2477665850397331
```

Each branch is now given two values : its synonymous length (expected synonymous susbtitution per **nucleotide** site) and its non-synonymous lengths. 

Then we can ask HyPhy to generate data where each branch has enforced synonymous and non-synonymous lengths.

Important default settings: CF3x4 frequency estimator, HKY85 model with &kappa; = 2

(new lines added for formatting only)

```
hyphy SimulateMG94.bf 
--tree CD2.dnds.nwk 
--sites 400 
--replicates 10
--output data/example4
--branch-variation dsdn 
--base-frequencies HIV
```



### Complete options list 

Note: these options will generally be extended by modules in `modules/branch-variation` and `modules/sitee-variation`, in a module-specific way (e.g. the `omega` keyword arguments in the examples above).

```
seed
	Random seed (0 to use default initialization)
	defaut value: 0

code
	Which genetic code should be used
	defaut value: Universal

tree [required]
	A phylogenetic tree with branch lengths or annotations

sites
	How many codon sites to simulate
	defaut value: 500 [computed at run time]

replicates
	How many replicates
	defaut value: 1 [computed at run time]

base-frequencies
	Base frequencies to use. 'equal' or 9 comma-separated values 
	[A in first codon position, C-1, G-1, A-2, C-2, G-2...] 
	or 12 comma-separated values 
	[A in first codon position, C-1, G-1, T-1, A-2, C-2, G-2, T-2...] 
	to specify positional nucleotide frequencies
	defaut value: equal

frequency-estimator
	Equilibrium frequency estimator
	defaut value: CF3x4

AC
	The AC substitution rate relative to the AG rate (=1)
	defaut value: 0.5

AT
	The AT substitution rate relative to the AG rate (=1)
	defaut value: 0.5

CG
	The CG substitution rate relative to the AG rate (=1)
	defaut value: 0.5

CT
	The CT substitution rate relative to the AG rate (=1)
	defaut value: 1.0

GT
	The GT substitution rate relative to the AG rate (=1)
	defaut value: 0.5

branch-variation
	The model for describing branch-to-branch variation in omega ratios
	defaut value: constant

site-variation
	The model for describing site-to-site variation in relative alpha (dS) 
	and beta (dN) rates
	defaut value: constant

output [required]
	Write simulated alignments (as FASTA) to the following prefix path, 
	using the syntax ${path}.replicate.index
```


### frequency-estimator

* `CF3x4` : corrected F3x4 estimator
* `F3x4` : standard F3x4 estimator
* `F1x4` : standard F1x4 estimator 

See [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011230) for details on CF3x4 and references to other estimators. 
 
### branch-variation

Maps to a file in `modules/branch-variation`, can be extended to allow new functionaility.

* `constant` : requires a tree with branch lengths, prompts for `omega`. &omega; does not vary across branches.
* `partitioned` : requires a tree with branch lengths and branches labeled with `{group}`, prompts for `omega` for each group. &omega; does not vary within groups.
* `dsdn` : requires a tree with branches labeled with `[s=value;ns=value]`, no prompts.


### site variation

Maps to a file in `modules/site-variation`, can be extended to allow new functionaility.

* `constant` : no site-to-site rate variation

 