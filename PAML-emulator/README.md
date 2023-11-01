## Partial PAML-style emulator

### Branch site-model A

`BranchSite.bf`

### Clade model C

`CMC.bf`

### Example run

```
hyphy CMC.bf Universal test/ATPalpha1.pml test/ATPalpha1.tree.Labeled.txt results/ATPalpha1-cmc.json

...

Universal

Loaded a 30 sequence alignment with 1045 codons from
/Users/sergei/Development/hyphy-analyses/PAML-emulator/test/ATPalpha1.pml

Base composition:
	A:    0.29513,   0.29111,   0.23327
	C:    0.16979,   0.23049,   0.26139
	G:    0.35218,   0.15767,   0.19073
	T:    0.18290,   0.32074,   0.31461
Obtaining nucleotide branch lengths and kappa to be used as starting values...

((Locusta_migratoria_MH450018:0.07791035966458586,(Chrotogonus_hemipterus_MK294068:0.05202822285002483,Atractomorpha_acutipennis_MK294069:0.05407311537760883)Node4:0.0422092064514993)Node2:0.08495320967372991,(Ctenocephalides_felis_S66043:0.1307528739634823,(((((((((((Chrysochus_auratus_AJ617743:0,Chrysochus_auratus_JQ771501:0)Node19:0.01590417419078815,Chrysochus_cobaltinus_MK765672:0.005785326531162721)Node18:0.03633047298049581,Chrysochus_auratus_HE956741:0.007772239186371064)Node17:0.01958153010738052,(Chrysochus_auratus_JQ771500:0.005796099516629186,Chrysochus_cobaltinus_MK765671:0.002610056939679774)Node24:0.01284750136209837)Node16:0.03666581302315633,Plagiodera_versicolora_JQ771522:0.05171790353713287)Node15:0.07125285460599827,Alticini_sp__SD_2012_HE956742:0.09198159500103707)Node14:0.03688152500816438,((Labidomera_clivicollis_HE956743:0.0009650350541321714,Labidomera_clivicollis_JQ771511:0.0003900305648625172)Node30:0.1004148157781499,Gastrophysa_viridula_HE956744:0.1147226102851616)Node29:0.03451835233971254)Node13:0.0188266627907306,((Tetraopes_tetrophthalmus_HE956745:0.0009417986395321217,Tetraopes_tetrophthalmus_JQ771526:0)Node35:0.1046540774133339,Megacyllene_robiniae_JQ771517:0.1060395193453841)Node34:0.03397685664771796)Node12:0.02380757359858859,(((Rhyssomatus_lineaticollis_HE956746:0.005023187627958644,Rhyssomatus_lineaticollis_JQ771524:0.007788241205389648)Node41:0.01325742259859575,Rhyssomatus_lineaticollis_JQ771523:0.02081674153324223)Node40:0.1243721692835709,Cyrtepistomus_castaneus_JQ771502:0.1226706703073612)Node39:0.04366360676343965)Node11:0.03881887250160536,Tribolium_castaneum_XM_969867:0.118903256404003)Node10:0.03827461406943879,(((Papilio_glaucus_JQ771498:0.08986992653504028,Bombyx_mori_LC029030:0.09508289990906124)Node49:0.04059321493529554,Apis_mellifera_XM_003250374:0.2333136912074149)Node48:0.02735697949697919,(Drosophila_melanogaster_AF044974:0.1095851287930945,Aedes_aegypti_XM_001662166:0.09852215417066099)Node53:0.04413518985752608)Node47:0.06939161756882033)Node9:0.07657150425028547)Node7:0.04003096802156561,(Acyrthosiphon_pisum_XM_001948888:0.2631892651342957,Pediculus_humanus_corporis_XM_002427669:0.1801289151673532)Node56:0.05716175325632992)
kappa=   3.087


11 foreground branch(es) set to: Labidomera_clivicollis_HE956743, Labidomera_clivicollis_JQ771511, Node30, Node35, Node40, Node41, Rhyssomatus_lineaticollis_HE956746, Rhyssomatus_lineaticollis_JQ771523, Rhyssomatus_lineaticollis_JQ771524, Tetraopes_tetrophthalmus_HE956745, Tetraopes_tetrophthalmus_JQ771526
Log(L) = -33253.49871910938
Checking for convergence by Latin Hypercube Sampling (this may take a bit of time...)

Convergence checks FAILED : a better score found. Restarting the optimization procedureLog(L) = -32801.27620575637
Checking for convergence by Latin Hypercube Sampling (this may take a bit of time...)

The estimation procedure appears to have converged.
** MODEL Clade model C **
Log (L) = -32801.27620575637

Inferred rate distribution:
	Class 0.  omega_0 = 0.005 weight = 0.806
	Class 1.  omega  := 1.000 weight = 0.007
	Class 2.  Foreground omega = 0.063 background omega := 0.076 weight = 0.186

===============================
** Fitting null model M2_rel **
===============================
** MODEL M2_rel model **
Log (L) = -32801.60804886476

Inferred rate distribution:
	Class 0.  omega_0 = 0.005 weight = 0.809
	Class 1.  omega  := 1.000 weight = 0.007
	Class 2.  Foreground omega = 0.075 background omega := 0.075 weight = 0.184

 p-value, Clade C model vs M2_rel = 0.4152615787
===============================
** Fitting null model M1a **
===============================
** MODEL M1a model **
Log (L) = -33253.47905199893

Inferred rate distribution:
	Class 0.  omega_0 = 0.013 weight = 0.971
	Class 1.  omega  := 1.000 weight = 0.029
	Class 2.  Foreground omega = 1.000 background omega := 1.000 weight = 0.000

 p-value, Clade C model vs M1a = 0.0000000000

```

```
cat results/ATPalpha1-cmc.json


{
 "fits":{
   "Clade model C":{
     "Class 1 omega":0.00489558177166652,
     "Class 1 weight":0.80638721731003,
     "Class 2 omega":1,
     "Class 2 weight":0.007161046455931649,
     "Class 3 omega BG":0.07613133580245579,
     "Class 3 omega FG":0.06323014617404525,
     "Class 3 weight":0.1864517362340383,
     "Log Likelihood":-32801.27620575637
    },
   "M1a model":{
     "Class 1 omega":0.01349144865812611,
     "Class 1 weight":0.9710570063211119,
     "Class 2 omega":1,
     "Class 2 weight":0.02894299367888808,
     "Class 3 omega BG":1,
     "Class 3 omega FG":1,
     "Class 3 weight":0,
     "Log Likelihood":-33253.47905199893
    },
   "M2_rel model":{
     "Class 1 omega":0.004867730912160508,
     "Class 1 weight":0.8088070192469138,
     "Class 2 omega":1,
     "Class 2 weight":0.007009274221718612,
     "Class 3 omega BG":0.07469013983656587,
     "Class 3 omega FG":0.07469013983656587,
     "Class 3 weight":0.1841837065313676,
     "Log Likelihood":-32801.60804886476
    }
  },
 "test results":{
   "Clade C model|M1a":{
     "LRT":904.4056924851175,
     "p-value":0
    },
   "Clade C model|M2_rel":{
     "LRT":0.6636862167651998,
     "p-value":0.4152615787013851
    }
  }
}          

```