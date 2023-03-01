# RELAX simulations

For convenience, we will simulate data using existing RELAX fits (as opposed to the more general simulation framework, where everything is specified as parameters).

### Generate a RELAX fit to use.

Here, we use an example from the [RELAX paper](http://bit.ly/veg_manuscripts_relax) (Bat SWS, Figure 4E).

```
hyphy relax --alignment data/Fig4E.nex 
			--save-fit fits/Fig4E.RELAX.fit 
			--test T 
			--starting-points 5
			--models Minimal
```

### Adjust model parameters to generate replicates

Weak relaxation (`K=0.95`)

```
$mkdir simulations/relaxation-0.95

$hyphy SimulateFromRELAX.bf 
	--fit fits/Fig4E.RELAX.fit 
	--output simulations/relaxation-0.95/replicate 
	--relax.K 0.95 
	--replicates 100
	--skip-tree True
```

Note, that you can check the console for which other parameters are adjustable (e.g. `relax.reference.omega1`, which you would set via a command line argument like `--relax.reference.omega1 <numerical value>`)

```
Value for relax.reference.omega3 to use during simulations (permissible range = [1,10000000000], default value = 1.384239219972456): relax.reference.omega3: 1.384239219972456
Value for relax.reference.omega2 to use during simulations (permissible range = [0,1], default value = 0.1616966251476609): relax.reference.omega2: 0.1616966251476609
Value for relax.reference.omega1 to use during simulations (permissible range = [0,1], default value = 0): relax.reference.omega1: 0
Value for relax.reference.theta_CG to use during simulations (permissible range = [0,10000], default value = 0.09452912049498727): relax.reference.theta_CG: 0.09452912049498727
Value for relax.reference.theta_CT to use during simulations (permissible range = [0,10000], default value = 0.5502496826079097): relax.reference.theta_CT: 0.5502496826079097
Value for relax.reference.theta_AC to use during simulations (permissible range = [0,10000], default value = 0.1859973662919255): relax.reference.theta_AC: 0.1859973662919255
Value for relax.reference.theta_GT to use during simulations (permissible range = [0,10000], default value = 0.1615908302349176): relax.reference.theta_GT: 0.1615908302349176
Value for relax.reference.theta_AT to use during simulations (permissible range = [0,10000], default value = 0.07925365227703088): relax.reference.theta_AT: 0.07925365227703088
Value for relax.reference.bsrel_mixture_aux_1 to use during simulations (permissible range = [0,1], default value = 0.9941025559305088): relax.reference.bsrel_mixture_aux_1: 0.9941025559305088
Value for relax.reference.bsrel_mixture_aux_0 to use during simulations (permissible range = [0,1], default value = 0.3514762875042471): relax.reference.bsrel_mixture_aux_0: 0.3514762875042471
Value for relax.K to use during simulations (permissible range = [0,50], default value = 9.918818434544907e-06): relax.K: 0.95
```

Also, let us simulate even weaker intensification (`K=1.025`)

```

$mkdir simulations/intensification-1.025

$hyphy SimulateFromRELAX.bf 
	--fit fits/Fig4E.RELAX.fit 
	--output simulations/intensification-1.025/replicate 
	--relax.K 1.025 --replicates 100 --skip-tree True
```

### Create a partitioned tree file 

This is a hack for now, because the simulator does not remember tree labeling. So, copy and paste the tree from the original fit file... (see `simulations/relaxation-0.95/tree.nwk`). 

### Fit individual replicates using RELAX and save fits
Here, I just use a simple shell one-liner.

```
for FILE in simulations/relaxation-0.95/replicate*.nex; 
do 
	hyphy relax --alignment $FILE --save-fit ${FILE}.RELAX.fit \
	--tree simulations/relaxation-0.95/tree.nwk \
	--test T \
	--models Minimal \
	--starting-points 5
done
```

and

```
for FILE in simulations/intensification-1.025/replicate*.nex; 
do 
	hyphy relax --alignment $FILE --save-fit ${FILE}.RELAX.fit \
	--tree simulations/intensification-1.025/tree.nwk \
	--test T \
	--models Minimal \
	--starting-points 5
done
```

We can check that the power to detect relaxation/intensification in individual files is low, and the point estimates of `K` are noisy

```
$jq '[.["test results"]["p-value"]] + [.["test results"]["relaxation or intensification parameter"]]' simulations/relaxation-0.95/*.RELAX.json 
...

[
  0.9187605371634595,
  0.9902332056009115
]
[
  0.7427206303141718,
  1.036488005745004
]
[
  0.3700943392636589,
  0.9208115333687038
]
[
  0.2737605257498569,
  0.9026783580337203
]
[
  0.5140370247247941,
  0.9378773097847889
]
[
  0.1214576954889909,
  0.7980680651257285
]
[
  0.3626431144293563,
  0.9238543107978173
]
[
  0.7021345474435101,
  0.9590518442474877
]
[
  0.01225085869510456,
  0.7878760758473555
]
[
  0.1178802974044687,
  0.8221322758351827
]
...

```
### Now check to see if a joint analysis does better

#### Relaxation 

Get a list of 5 random fits... 

```
ls /Users/sergei/Development/hyphy-analyses/RELAX-joint/simulations/relaxation-0.95/*fit 
	| shuf 
	| head -n 5 
	> simulations/relaxation-0.95/random-5.lst 
```

```
hyphy RELAX-joint.bf 
	--filelist simulations/relaxation-0.95/random-5.lst  
	--output simulations/relaxation-0.95/random-5.joint.json
	
```

Expect **two** results

1. No improvement in having separate `K` for each file, since they were all simulated with the same `K`
2. Better power to reject `K=1` because we are pooling information across genes.

Indeed, the analysis reports the following

```
- Joint model likelihood -22306.276

- Shared K =      0.921 (95% profile CI   0.8546-  0.9904)

### - p-value for file-level K (vs a single K) 0.546807741996161

...

- p-value for K!=1 0.1066079324352446

```

Let's try 10 random fits (should have more power)

```
ls /Users/sergei/Development/hyphy-analyses/RELAX-joint/simulations/relaxation-0.95/*fit 
	| shuf 
	| head -n 10 
	> simulations/relaxation-0.95/random-10.lst 
```

```
hyphy RELAX-joint.bf 
	--filelist simulations/relaxation-0.95/random-10.lst  
	--output simulations/relaxation-0.95/random-10.joint.json
	
```

Still no benefit from multiple `K` (as it should be), and a pretty damn good estimate of the K that we simulated with.

```
- Joint model likelihood -45700.486
- Shared K =      0.905 (95% profile CI   0.8616-  0.9500)
### - p-value for file-level K (vs a single K) 0.9951764875516653
```

And now we have better power to detect relaxation.

```
- p-value for K!=1 0.005281270097204471
```

#### Intensification

Should be similar to the relaxation example above, except even **harder** to detect, because the effect is even smaller (`K` is closer to 1) 

```
ls /Users/sergei/Development/hyphy-analyses/RELAX-joint/simulations/intensification-1.025/*fit 
	| shuf 
	| head -n 10
	> simulations/intensification-1.025/random-10.lst 
```

```
hyphy RELAX-joint.bf 
	--filelist simulations/intensification-1.025/random-10.lst  
	--output simulations/intensification-1.025/random-10.joint.json
```

No benefit from multiple `K` (as it should be), but a somewhat inflated estimate of the `K` that we simulated with (should be `1.025`).

```
- Independent model likelihood -45109.700
- Joint model likelihood -45112.398
- Shared K =      1.090 (95% profile CI   1.0426-  1.1408)
### - p-value for file-level K (vs a single K) 0.7985343417021018
```

Good power to detect weak intensification, however:

```
- p-value for K!=1 0.005533313667466233
```

#### Push-pull

Finally, let's mix things up by jointly fitting 10 datasets: 5 with `K=0.9` and 5 with `K=1.025`. 

```
ls /Users/sergei/Development/hyphy-analyses/RELAX-joint/simulations/intensification-1.025/*fit 
	| shuf 
	| head -n 5
	> simulations/mix-10.lst
	
ls /Users/sergei/Development/hyphy-analyses/RELAX-joint/simulations/relaxation-0.9/*fit 
	| shuf 
	| head -n 5
	>> simulations/mix-10.lst
```

```
hyphy RELAX-joint.bf 
	--filelist simulations/mix-10.lst  
	--output simulations/mix-10.json 
```

There now IS a benefit from having multiple `K` per file (since there really are two truly different values). 

```
- Independent model likelihood -44612.530
- Joint model likelihood -44622.325
- Shared K =      1.034 (95% profile CI   0.9886-  1.0828)
### - p-value for file-level K (vs a single K) 0.02061160029225839
```

Because of "conflicting" signals, the overall p-value for K≠1 is (as expected) is also not significant.

```
- p-value for K!=1 0.2795072384040402
```