## Examples of invokation for specific scenarios/models.

### PRIME

Simulate data under the property informed model of sequence evolution (PRIME), i.e. where non-synonymous substitution rates depend on residue properties.

Example call

```
hyphy LIBPATH=/Users/sergei/Development/hyphy/res SimulateMG94.bf 
--tree trees/C3-C4-rubisco.nwk 
--site-variation prime  
--model prime 
--output sims/prime
--property-set Atchley
--fraction-same 0.5 
--fraction-single 0.75 
--effect-variance 1 
--mean-effect 0
```

Option description

* `property-set`  [enum] : which property set to use (LCAP or Atchley)
* `fraction-same` [fraction, 0-1]: what fraction of simulated sites will have no property influence (all &lambda; values = 0)
* `fraction-single` [fraction, 0-1]: what fraction of sites where some &lambda; ≠ 0 (i.e. 1 - `fraction-same`) will have only ONE &lambda; ≠ 0. This set of sites is influenced by a **single** (randomly selected) property
* `mean-effect` [float, -5 to 5]: when &lambda; are being randomly generated, this is the mean of the **normal** distribution that will be used (default = 0)
* `effect-variance` [float, 0 to 10]: when &lambda; are being randomly generated, this is the variance of the **normal** distribution that will be used (default = 1.)

