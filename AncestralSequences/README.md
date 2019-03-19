## Ancestral sequence reconstruction

This analysis loads a previously fitted likelihood function file (e.g see [FitMG94](https://github.com/veg/hyphy-analyses/tree/master/FitMG94)), reconstructs ancestral sequences using joint maximum likelihood, writes a detailed JSON file [see description at the very end of this document] which can be further processed by custom scripts, and reports various summary data to the screen.

## Invokation

This analysis has one **required** argument

- `--fit` a previously generated fit file (see example below)

HyPhy will write Markdown output to the screen and a JSON file with detailed fit results. 
See example at the end of the document. 

### Complete options list 

```
Analysis options description
----------------------------
fit [required]
	Load a previously saved HyPhy fit (NEXUS + HBL)

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + '.ancestors.json')
	defaut value: ancestors.file+".ancestors.json" [computed at run time]
```


##Example run 

First we will create the fit file, by running `FitMG94.bf` on the provided example (subject 4059 from [this study](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1002286))

```
HYPHYMP /path/to/FitMG94.bf --alignment 4059.fas --save-fit 4059.fit
```

`4059.fit` is a NEXUS file with a custom HyPhy blocks that stores fitted parameter values and other data. 

Now we can run `AncestralSequences.bf` on the fit file as follows

```
HYPHYMP AncestralSequences.bf --fit 4059.fit --output 4059.json
```

The following data are output to the screen. You can use any Newick tree viewer, e.g. [Phylotree.js](phylotree.hyphy.org), to display tree strings.

--- 
Analysis Description
--------------------
Load a previously generated model fit file (NEXUS + HBL), reconstruct
ancestral sequences, and map substitutions

- __Requirements__: a fit file generated previously by a different HyPhy analysis

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1

fit: 4059.fit
output: 4059.json


### Loaded a previously saved likelihood function with **1** partitions from `/Users/sergei/Development/hyphy-analyses/AncestralSequences/4059.fit`

### Successfully reconstructed joint maximum likelihood ancestors
**186** sites have substitutions along at least one branch

### Phylogenetic tree annotated by the inferred number of substitutions
```
(((((JN563162_isolate_4059_C_5:0,(((JN563149_isolate_4059_C_12:0,JN562779_isolate_4059_C12_from_USA_envelope_glycoprotein_env_gene_complete_cds:3)Node8:6,(JN563153_isolate_4059_C_16:1,JN563148_isolate_4059_C_11:0)Node11:0)Node7:9,((JN563155_isolate_4059_C_18:0,JN563151_isolate_4059_C_14:2)Node15:12,(JN563161_isolate_4059_C_4:18,(JN563134_isolate_4059_P_4:30,(((((JN563136_isolate_4059_P_44:0,JN562781_isolate_4059_P28_from_USA_envelope_glycoprotein_env_gene_complete_cds:2)Node26:1,JN563124_isolate_4059_P_25:0)Node25:17,JN563135_isolate_4059_P_41:8)Node24:11,JN563137_isolate_4059_P_46:8)Node23:8,((JN562780_isolate_4059_P21_from_USA_envelope_glycoprotein_env_gene_complete_cds:2,JN563122_isolate_4059_P_21:0)Node33:12,(JN563138_isolate_4059_P_47:11,(JN563144_isolate_4059_P_8:5,((JN563119_isolate_4059_P_15:9,JN563133_isolate_4059_P_39:6)Node41:2,((((((JN563145_isolate_4059_P_9:5,(JN562784_isolate_4059_P26_from_USA_envelope_glycoprotein_env_gene_complete_cds:2,JN563125_isolate_4059_P_26:0)Node51:14)Node49:2,(JN563128_isolate_4059_P_31:10,((JN563142_isolate_4059_P_52:8,(JN563143_isolate_4059_P_7:5,JN563146_isolate_4059_C_1:9)Node59:2)Node57:2,JN563132_isolate_4059_P_36:3)Node56:1)Node54:1)Node48:3,(JN563130_isolate_4059_P_34:7,(JN563129_isolate_4059_P_33:12,JN563120_isolate_4059_P_19:6)Node65:4)Node63:2)Node47:0,(JN563139_isolate_4059_P_48:11,JN563140_isolate_4059_P_50:4)Node68:3)Node46:1,(((JN563131_isolate_4059_P_35:8,JN563127_isolate_4059_P_29:6)Node73:2,JN563121_isolate_4059_P_2:7)Node72:0,JN563141_isolate_4059_P_51:6)Node71:2)Node45:1,JN563123_isolate_4059_P_23:12)Node44:0)Node40:4)Node38:3)Node36:25)Node32:12)Node22:16)Node20:12)Node18:4)Node14:8)Node6:21)Node4:7,(JN563163_isolate_4059_C_6:0,JN562783_isolate_4059_C6_from_USA_envelope_glycoprotein_env_gene_complete_cds:2)Node79:9)Node3:10,JN563157_isolate_4059_C_2:1)Node2:14,((JN563152_isolate_4059_C_15:2,JN563164_isolate_4059_C_7:2)Node84:0,(JN563165_isolate_4059_C_9:7,JN563147_isolate_4059_C_10:1)Node87:0)Node83:1)Node1:4,(JN563160_isolate_4059_C_3:0,JN563158_isolate_4059_C_21:1)Node90:2,(JN563159_isolate_4059_C_22:15,JN562782_isolate_4059_C19_from_USA_envelope_glycoprotein_env_gene_complete_cds:2)Node93:3)Node93_1
```

|        Site        | Substitution  |                     Branches                     |
|--------------------|---------------|--------------------------------------------------|
|         2          |   GTA->GTC    |(6) JN562779_isolate_4059_C12_from_USA_envelope...|
|         2          |   GTC->GTA    |                    (1) Node2                     |
|         3          |   CCT->ACT    |(4) JN562779_isolate_4059_C12_from_USA_envelope...|
|         3          |   CCT->AAT    |                    (1) Node15                    |
|         3          |   CCT->CAT    |      (2) JN563135_isolate_4059_P_41, Node33      |
|         3          |   CAT->AAT    |(1) JN562780_isolate_4059_P21_from_USA_envelope...|
|         3          |   AAT->CCT    |                    (1) Node2                     |
|         4          |   TTA->CTA    |          (1) JN563140_isolate_4059_P_50          |
|         4          |   CTA->TTA    |                    (1) Node2                     |
|         5          |   AAT->CAA    |          (1) JN563134_isolate_4059_P_4           |
|         5          |   AAT->CAG    |                    (1) Node36                    |
|         7          |   ACT->TCT    |                    (1) Node36                    |
|         8          |   AAA->AAG    |                    (1) Node8                     |
|         8          |   AAA->AAT    |      (2) JN563134_isolate_4059_P_4, Node36       |
|         8          |   GAT->AAA    |                    (1) Node2                     |
|         9          |   GTA->GTG    |                    (1) Node36                    |
|         9          |   GTG->GTA    |                (2) Node22, Node79                |
|         9          |   AAA->GTG    |                    (1) Node2                     |
|         12         |   GTG->ACC    |          (1) JN563134_isolate_4059_P_4           |
|         13         |   AAC->AAA    |          (1) JN563134_isolate_4059_P_4           |
|         14         |   AAC->AAT    |          (1) JN563134_isolate_4059_P_4           |
|         14         |   AAC->ACC    |                    (1) Node36                    |
|         14         |   AAT->AAC    |                    (1) Node2                     |
|         15         |   AAA->AAT    |          (1) JN563134_isolate_4059_P_4           |
|         15         |   TCG->AAA    |                    (1) Node2                     |
|         16         |   ACT->ACG    |          (1) JN563137_isolate_4059_P_46          |
|         16         |   ACT->AAT    |                    (1) Node36                    |
|         16         |   ACC->ACT    |                    (1) Node2                     |
|         17         |   GAT->AAT    |          (1) JN563138_isolate_4059_P_47          |
|         17         |   AAT->GAT    |                    (1) Node22                    |
|         17         |   ATA->AAT    |                    (1) Node2                     |
|         18         |   ACT->ACC    |          (1) JN563120_isolate_4059_P_19          |
|         18         |   CCT->ACT    |                    (1) Node20                    |
|         18         |   AAG->CCT    |                    (1) Node2                     |
|         19         |   AAT->ATT    |          (1) JN563138_isolate_4059_P_47          |
|         20         |   ACC->GCC    |          (1) JN563138_isolate_4059_P_47          |
|         20         |   ACC->AAT    |          (1) JN563123_isolate_4059_P_23          |
|         21         |   ACC->GCC    |                (2) Node25, Node32                |
|         22         |   AAT->ATT    |                    (1) Node38                    |
|         23         |   AAT->GGT    |(2) JN563134_isolate_4059_P_4, JN563138_isolate...|
|         23         |   AAT->AGT    |(3) Node51, JN563143_isolate_4059_P_7, JN563132...|
|         23         |   AAT->ACT    |                    (1) Node22                    |
|         23         |   ACT->AAT    |                    (1) Node36                    |
|         23         |   GAT->AAT    |                    (1) Node2                     |
|         24         |   GGT->GGC    |          (1) JN563134_isolate_4059_P_4           |
|         24         |   GGT->GGA    |          (1) JN563133_isolate_4059_P_39          |
|         24         |   GGT->GAT    |                    (1) Node22                    |
|         24         |   GAT->GAC    |          (1) JN563138_isolate_4059_P_47          |
|         24         |   GAT->GGT    |                    (1) Node38                    |
|         24         |   GAC->GGT    |                    (1) Node2                     |
|         25         |   GGG->GGA    |          (1) JN563133_isolate_4059_P_39          |
|         27         |   GGA->GAA    |(5) Node26, JN563119_isolate_4059_P_15, JN56314...|
|         28         |   ATG->GGG    |                    (1) Node22                    |
|         30         |   AAT->AAG    |          (1) JN563135_isolate_4059_P_41          |
|         30         |   AAT->CAT    |                    (1) Node36                    |
|         30         |   GAT->AAT    |                    (1) Node2                     |
|         31         |   GAA->AAA    |                    (1) Node36                    |
|         31         |   AAA->GAA    |                    (1) Node2                     |
|         33         |   ATA->ATG    |                    (1) Node20                    |
|         33         |   ATG->ATA    |                    (1) Node3                     |
|         41         |   ACC->ACT    |(5) JN563138_isolate_4059_P_47, JN563128_isolat...|
|         44         |   GAA->GAG    |(3) JN563119_isolate_4059_P_15, JN563128_isolat...|
|         45         |   GAA->GAG    |(6) Node51, JN563128_isolate_4059_P_31, JN56314...|
|         45         |   GAA->CAA    |          (1) JN563130_isolate_4059_P_34          |
|         45         |   GAG->GAA    |                    (1) Node45                    |
|         45         |   GGA->GAG    |                    (1) Node32                    |
|         46         |   AAT->ATT    |          (1) JN563134_isolate_4059_P_4           |
|         46         |   AAT->ACT    |                    (1) Node25                    |
|         46         |   AAT->AAC    |          (1) JN563130_isolate_4059_P_34          |
|         48         |   GTG->GTT    |                    (1) Node41                    |
|         48         |   GTG->GTA    |          (1) JN563131_isolate_4059_P_35          |
|         49         |   AAA->CAA    |          (1) JN563161_isolate_4059_C_4           |
|         49         |   AAA->CAG    |                    (1) Node32                    |
|         49         |   CGG->AAA    |                    (1) Node3                     |
|         51         |   GAA->GAG    |                    (1) Node8                     |
|         51         |   GAA->AAA    |                    (1) Node32                    |
|         55         |   CTG->CTT    |      (2) JN563161_isolate_4059_C_4, Node22       |
|         55         |   CTT->CTG    |                (2) Node24, Node3                 |
|         60         |   ATA->GTA    |                    (1) Node24                    |
|         60         |   GTA->ATA    |                    (1) Node20                    |
|         63         |   ATA->AGA    |                    (1) Node8                     |
|         64         |   GAT->GGT    |                    (1) Node24                    |
|         65         |   AGT->AAT    |                    (1) Node24                    |
|         65         |   AGT->GAT    |          (1) JN563130_isolate_4059_P_34          |
|         65         |   GGT->AGT    |                    (1) Node18                    |
|         66         |   AGT->AAT    |          (1) JN563161_isolate_4059_C_4           |
|         66         |   AAT->AGT    |                    (1) Node3                     |
|         67         |   AGT->AAT    |      (2) JN563134_isolate_4059_P_4, Node24       |
|         67         |   AGT->GAT    |                    (1) Node63                    |
|         67         |   AAA->AGT    |                    (1) Node3                     |
|         68         |   AGT->AAT    |(3) JN563161_isolate_4059_C_4, Node24, JN563139...|
|         68         |   AGT->ACT    |          (1) JN563134_isolate_4059_P_4           |
|         68         |   AAT->AAC    |                    (1) Node25                    |
|         68         |   AAT->AGT    |                    (1) Node3                     |
|         76         |   AAG->AAA    |                    (1) Node8                     |
|         76         |   AAG->AAT    |                    (1) Node18                    |
|         76         |   AGT->GAT    |                    (1) Node24                    |
|         76         |   AGT->GGT    |(3) Node33, JN563139_isolate_4059_P_48, JN56312...|
|         76         |   AGT->AAT    |                    (1) Node63                    |
|         76         |   AAT->AAG    |                (2) Node65, Node3                 |
|         76         |   AAT->AGT    |                    (1) Node22                    |
|         77         |   AAT->GAT    |          (1) JN563161_isolate_4059_C_4           |
|         77         |   AAT->AAA    |                    (1) Node20                    |
|         77         |   AAA->GAA    |          (1) JN563134_isolate_4059_P_4           |
|         77         |   AAA->AGA    |                    (1) Node22                    |
|         77         |   AGA->AGT    |                    (1) Node24                    |
|         77         |   AGA->AAA    |(5) Node33, JN563146_isolate_4059_C_1, JN563120...|
|         78         |   AAC->AGT    |          (1) JN563161_isolate_4059_C_4           |
|         78         |   AAC->AGC    |                    (1) Node24                    |
|         78         |   AAC->GAA    |                    (1) Node3                     |
|         78         |   GAA->AAC    |                    (1) Node18                    |
|         79         |   ACC->AGC    |          (1) JN563161_isolate_4059_C_4           |
|         80         |   TAT->TAC    |          (1) JN563142_isolate_4059_P_52          |
|         80         |   TAT->CAT    |                    (1) Node65                    |
|         80         |   TAT->AAT    |          (1) JN563123_isolate_4059_P_23          |
|         82         |   AAC->AGC    |                    (1) Node15                    |
|         82         |   AAC->AAA    |(2) JN563120_isolate_4059_P_19, JN563123_isolat...|
|         84         |   AGG->AAG    |          (1) JN563134_isolate_4059_P_4           |
|         84         |   AGA->AGG    |                    (1) Node3                     |
|         85         |   TTG->TTA    |      (2) Node15, JN563146_isolate_4059_C_1       |
|         90         |   ACC->ACA    |          (1) JN563127_isolate_4059_P_29          |
|        105         |   ATT->ATC    |          (1) JN563142_isolate_4059_P_52          |
|        112         |   CCG->CCA    |(2) JN563128_isolate_4059_P_31, JN563127_isolat...|
|        118         |   CTA->TTA    |                 (2) Node8, Node3                 |
|        118         |   TTA->CTA    |                    (1) Node6                     |
|        120         |   TGT->TGC    |          (1) JN563128_isolate_4059_P_31          |
|        123         |   AAG->AGG    |          (1) JN563120_isolate_4059_P_19          |
|        123         |   AAG->AAA    |(1) JN562782_isolate_4059_C19_from_USA_envelope...|
|        124         |   AAG->AAA    |                    (1) Node33                    |
|        125         |   TTC->TTT    |          (1) JN563138_isolate_4059_P_47          |
|        127         |   GGA->GGT    |                    (1) Node25                    |
|        132         |   AAA->GGA    |                    (1) Node33                    |
|        133         |   AAT->AAC    |          (1) JN563137_isolate_4059_P_46          |
|        136         |   ACA->ACG    |(2) JN563137_isolate_4059_P_46, JN563121_isolat...|
|        142         |   GGA->GGG    |(1) JN562779_isolate_4059_C12_from_USA_envelope...|
|        149         |   ACC->ACT    |                    (1) Node25                    |
|        149         |   ACT->ACC    |                    (1) Node23                    |
|        151         |   CTG->TTG    |                    (1) Node33                    |
|        151         |   CTG->CTA    |          (1) JN563146_isolate_4059_C_1           |
|        151         |   TTG->CTG    |                    (1) Node6                     |
|        153         |   TTA->CTA    |          (1) JN563139_isolate_4059_P_48          |
|        155         |   GGC->GGT    |          (1) JN563134_isolate_4059_P_4           |
|        157         |   CTA->ATA    |                    (1) Node8                     |
|        161         |   GAG->GGG    |      (2) Node15, JN563138_isolate_4059_P_47      |
|        165         |   AGA->AAA    |                    (1) Node38                    |
|        166         |   TCT->TCC    |          (1) JN563138_isolate_4059_P_47          |
|        166         |   TCT->TCA    |          (1) JN563145_isolate_4059_P_9           |
|        167         |   GAT->GAC    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        170         |   TCT->TCA    |                    (1) Node25                    |
|        170         |   TCA->TCT    |                    (1) Node23                    |
|        171         |   AAC->AAT    |                (2) Node15, Node36                |
|        178         |   GTA->GTG    |                    (1) Node25                    |
|        179         |   CAG->CAA    |          (1) JN563132_isolate_4059_P_36          |
|        180         |   CTG->CTA    |          (1) JN563161_isolate_4059_C_4           |
|        180         |   CTG->TTG    |          (1) JN563141_isolate_4059_P_51          |
|        181         |   AAT->AAC    |          (1) JN563137_isolate_4059_P_46          |
|        183         |   ACT->ATT    |                    (1) Node15                    |
|        185         |   AAG->GAG    |          (1) JN563161_isolate_4059_C_4           |
|        185         |   AAG->AAA    |                    (1) Node33                    |
|        187         |   ATT->AAT    |                    (1) Node36                    |
|        187         |   AAT->ATT    |                    (1) Node20                    |
|        192         |   AAT->AAC    |          (1) JN563142_isolate_4059_P_52          |
|        197         |   AGA->ACA    |                    (1) Node15                    |
|        197         |   AGA->AAA    |                    (1) Node32                    |
|        197         |   AAA->ACA    |                    (1) Node36                    |
|        198         |   GGT->AGT    |          (1) JN563161_isolate_4059_C_4           |
|        198         |   AGT->GGT    |                    (1) Node6                     |
|        200         |   AAT->AGT    |          (1) JN563161_isolate_4059_C_4           |
|        200         |   AGT->AAT    |                    (1) Node6                     |
|        201         |   TTA->TTC    |          (1) JN563161_isolate_4059_C_4           |
|        201         |   TTC->TTA    |                    (1) Node6                     |
|        203         |   CAA->CCA    |          (1) JN563161_isolate_4059_C_4           |
|        203         |   CCA->CAA    |                    (1) Node6                     |
|        208         |   CAT->TAT    |   (3) Node6, Node79, JN563165_isolate_4059_C_9   |
|        212         |   GAC->GAT    |          (1) JN563134_isolate_4059_P_4           |
|        223         |   ACC->AAC    |(4) Node7, Node20, Node79, JN563165_isolate_405...|
|        223         |   AAC->AGC    |          (1) JN563137_isolate_4059_P_46          |
|        225         |   AAT->AGT    |(4) Node7, Node20, Node79, JN563165_isolate_405...|
|        226         |   AGA->ACA    |                    (1) Node7                     |
|        226         |   AGA->GGA    |                    (1) Node25                    |
|        226         |   AGA->AGG    |          (1) JN563119_isolate_4059_P_15          |
|        226         |   GGA->AGA    |                    (1) Node6                     |
|        226         |   GGA->ACA    |      (2) Node79, JN563165_isolate_4059_C_9       |
|        227         |   ACA->GCA    |          (1) JN563139_isolate_4059_P_48          |
|        228         |   GAA->CAA    |                    (1) Node25                    |
|        228         |   GAA->GAT    |(3) JN563144_isolate_4059_P_8, Node54, JN563123...|
|        228         |   GAT->GAA    |          (1) JN563143_isolate_4059_P_7           |
|        231         |   AAC->GAC    |                    (1) Node25                    |
|        232         |   ACT->ATT    |                    (1) Node24                    |
|        232         |   ATT->ACT    |   (3) Node6, Node79, JN563165_isolate_4059_C_9   |
|        234         |   AAA->AGA    |(4) Node7, Node25, Node79, JN563165_isolate_405...|
|        235         |   CTG->CTA    |      (2) JN563119_isolate_4059_P_15, Node90      |
|        237         |   GCT->GTT    |                    (1) Node51                    |
|        237         |   GTT->GCT    |                    (1) Node32                    |
|        238         |   AAG->AAC    |(3) JN563133_isolate_4059_P_39, Node51, JN56312...|
|        238         |   AAC->AAG    |                    (1) Node32                    |
|        238         |   AGC->AAC    |   (3) Node6, Node79, JN563165_isolate_4059_C_9   |
|        239         |   AAA->AAG    |          (1) JN563128_isolate_4059_P_31          |
|        241         |   CTA->AAA    |                    (1) Node51                    |
|        241         |   CAA->CTA    |                    (1) Node32                    |
|        242         |   GAA->GGA    |                    (1) Node25                    |
|        245         |   GGA->GGG    |      (2) JN563134_isolate_4059_P_4, Node51       |
|        245         |   GGA->GGT    |(4) JN563138_isolate_4059_P_47, JN563142_isolat...|
|        245         |   GGG->GGA    |                    (1) Node14                    |
|        246         |   ACT->GCT    |          (1) JN563134_isolate_4059_P_4           |
|        246         |   ACT->ATT    |                    (1) Node51                    |
|        251         |   AAA->AGA    |                    (1) Node23                    |
|        251         |   AAA->GTC    |                    (1) Node51                    |
|        251         |   AGA->AAA    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        253         |   GAG->GAA    |          (1) JN563121_isolate_4059_P_2           |
|        255         |   CCT->CCC    |                    (1) Node51                    |
|        255         |   CCT->CCA    |          (1) JN563130_isolate_4059_P_34          |
|        257         |   CAG->CCG    |          (1) JN563134_isolate_4059_P_4           |
|        257         |   CAG->GAG    |                    (1) Node24                    |
|        257         |   CAG->CAA    |          (1) JN563158_isolate_4059_C_21          |
|        257         |   GAG->GGG    |                    (1) Node25                    |
|        258         |   GGA->GGG    |          (1) JN563140_isolate_4059_P_50          |
|        261         |   CTA->CCA    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        262         |   GAA->AAA    |          (1) JN563146_isolate_4059_C_1           |
|        264         |   GTA->TTA    |      (2) Node7, JN563159_isolate_4059_C_22       |
|        264         |   GTA->GTG    |          (1) JN563137_isolate_4059_P_46          |
|        264         |   GTA->CTA    |                    (1) Node32                    |
|        264         |   TTA->CTA    |          (1) JN563121_isolate_4059_P_2           |
|        264         |   CTA->TTA    |                    (1) Node36                    |
|        269         |   AAT->AAC    |          (1) JN563141_isolate_4059_P_51          |
|        273         |   GAA->GAG    |          (1) JN563123_isolate_4059_P_23          |
|        274         |   TTT->TTC    |          (1) JN563123_isolate_4059_P_23          |
|        275         |   TTC->TTT    |          (1) JN563135_isolate_4059_P_41          |
|        276         |   TAC->TAT    |          (1) JN563161_isolate_4059_C_4           |
|        281         |   AAA->AAG    |          (1) JN563135_isolate_4059_P_41          |
|        281         |   AAA->CAA    |(1) JN562782_isolate_4059_C19_from_USA_envelope...|
|        282         |   CTG->ATG    |          (1) JN563139_isolate_4059_P_48          |
|        282         |   CTG->TTG    |          (1) JN563127_isolate_4059_P_29          |
|        286         |   ACT->ATT    |                (2) Node79, Node83                |
|        288         |   AAT->AAA    |  (3) Node15, JN563161_isolate_4059_C_4, Node32   |
|        288         |   AAT->ACT    |          (1) JN563134_isolate_4059_P_4           |
|        288         |   AAT->GAT    |                    (1) Node4                     |
|        288         |   GAT->AAT    |                (2) Node14, Node1                 |
|        289         |   AAT->AGT    |                (2) Node51, Node20                |
|        289         |   AGT->AAT    |                    (1) Node32                    |
|        291         |   ACT->AAT    |(2) JN563134_isolate_4059_P_4, JN563143_isolate...|
|        291         |   ACT->GCT    |      (2) JN563119_isolate_4059_P_15, Node1       |
|        291         |   GCT->ACT    |                    (1) Node4                     |
|        292         |   AGT->AAT    |          (1) JN563139_isolate_4059_P_48          |
|        293         |   ATT->ACT    |                    (1) Node51                    |
|        293         |   ATT->GTT    |          (1) JN563142_isolate_4059_P_52          |
|        293         |   ACT->ATT    |                    (1) Node40                    |
|        293         |   ACT->ACC    |                    (1) Node1                     |
|        293         |   ACC->ACT    |                    (1) Node4                     |
|        294         |   TGG->GGG    |  (3) Node51, JN563141_isolate_4059_P_51, Node4   |
|        294         |   GGG->TGG    |                (2) Node14, Node1                 |
|        295         |   GAT->GAA    |          (1) JN563133_isolate_4059_P_39          |
|        295         |   GAT->AAT    |(3) JN563128_isolate_4059_P_31, JN563139_isolat...|
|        295         |   GAT->GGT    |(2) JN563142_isolate_4059_P_52, JN563146_isolat...|
|        295         |   GAT->GAG    |(3) Node65, JN563140_isolate_4059_P_50, JN56312...|
|        295         |   GAG->GAA    |          (1) JN563129_isolate_4059_P_33          |
|        296         |   AAT->AAA    |(2) JN563145_isolate_4059_P_9, JN563130_isolate...|
|        296         |   AAT->GAT    |                    (1) Node36                    |
|        296         |   AAT->AGT    |          (1) JN563147_isolate_4059_C_10          |
|        296         |   GAT->AAT    |                    (1) Node40                    |
|        297         |   AAC->GAC    |          (1) JN563134_isolate_4059_P_4           |
|        297         |   AAC->AAT    |                    (1) Node25                    |
|        297         |   AAC->AGC    |      (2) Node33, JN563144_isolate_4059_P_8       |
|        297         |   GAC->AAC    |                (2) Node4, Node93                 |
|        298         |   ATG->ACT    |                    (1) Node49                    |
|        298         |   ATG->AGG    |          (1) JN563131_isolate_4059_P_35          |
|        298         |   ACT->ATG    |                    (1) Node40                    |
|        299         |   GAA->GCA    |                    (1) Node15                    |
|        300         |   GGA->GAA    |                (2) Node4, Node93                 |
|        301         |   GAT->AAT    |                    (1) Node18                    |
|        301         |   AAT->GAT    |                (2) Node4, Node93                 |
|        302         |   GGC->AGC    |                    (1) Node23                    |
|        303         |   ACA->ACG    |          (1) JN563152_isolate_4059_C_15          |
|        305         |   ACA->TCA    |          (1) JN563144_isolate_4059_P_8           |
|        306         |   CTC->CTT    |          (1) JN563137_isolate_4059_P_46          |
|        311         |   AGG->AGA    |(4) JN563134_isolate_4059_P_4, Node23, JN563144...|
|        311         |   AGA->AGG    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        313         |   TTT->ATT    |      (2) JN563134_isolate_4059_P_4, Node23       |
|        313         |   ATT->TTT    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        315         |   AAC->AAT    |                    (1) Node36                    |
|        315         |   AAT->AAC    |                    (1) Node20                    |
|        318         |   CAG->CAA    |(4) JN563146_isolate_4059_C_1, JN563132_isolate...|
|        318         |   CAA->CAG    |(4) Node46, JN563131_isolate_4059_P_35, Node6, ...|
|        319         |   AAG->AAA    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        329         |   ATC->ATT    |(2) JN563134_isolate_4059_P_4, JN563123_isolate...|
|        329         |   ATC->ATA    |          (1) JN563135_isolate_4059_P_41          |
|        330         |   AGC->GGA    |          (1) JN563134_isolate_4059_P_4           |
|        330         |   AGC->AGT    |(3) JN563143_isolate_4059_P_7, JN563129_isolate...|
|        330         |   AGT->AGC    |                    (1) Node36                    |
|        332         |   CGA->CGG    |      (2) Node7, JN563159_isolate_4059_C_22       |
|        334         |   AAC->AAT    |(6) JN563119_isolate_4059_P_15, Node51, Node59,...|
|        341         |   GGG->GGA    |          (1) JN563135_isolate_4059_P_41          |
|        341         |   GGA->GGG    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        342         |   CTG->TTG    |(5) Node41, Node51, Node57, JN563129_isolate_40...|
|        342         |   CTG->CTA    |          (1) JN563123_isolate_4059_P_23          |
|        343         |   CTA->TTA    |(2) JN563135_isolate_4059_P_41, JN563131_isolat...|
|        343         |   CTA->CTG    |          (1) JN563129_isolate_4059_P_33          |
|        345         |   ATG->ATA    |(4) Node7, JN563161_isolate_4059_C_4, Node22, J...|
|        350         |   AGT->GGT    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        351         |   GGT->GGA    |                    (1) Node33                    |
|        351         |   AAT->GGT    |                    (1) Node6                     |
|        352         |   GAC->GGT    |                    (1) Node22                    |
|        353         |   ACG->AGT    |      (2) Node7, JN563159_isolate_4059_C_22       |
|        353         |   ACG->GCG    |          (1) JN563134_isolate_4059_P_4           |
|        353         |   ACG->AAT    |                    (1) Node22                    |
|        353         |   AAT->GTT    |                    (1) Node36                    |
|        354         |   GAC->AAC    |                    (1) Node36                    |
|        354         |   AAC->GAC    |                    (1) Node22                    |
|        355         |   GAG->GAA    |          (1) JN563153_isolate_4059_C_16          |
|        355         |   GAG->GAT    |                    (1) Node15                    |
|        355         |   GAG->GGG    |(5) JN563134_isolate_4059_P_4, JN563133_isolate...|
|        355         |   GAG->GAC    |          (1) JN563130_isolate_4059_P_34          |
|        355         |   GAG->GGA    |          (1) JN563121_isolate_4059_P_2           |
|        355         |   GAG->ACG    |                    (1) Node22                    |
|        355         |   GGG->GAG    |                    (1) Node40                    |
|        355         |   ACG->GGG    |                    (1) Node36                    |
|        356         |   ACC->ACA    |(4) JN563145_isolate_4059_P_9, JN563129_isolate...|
|        356         |   CCC->ACC    |      (2) Node6, JN563159_isolate_4059_C_22       |
|        358         |   ATC->ACC    |(5) JN563145_isolate_4059_P_9, JN563129_isolate...|
|        358         |   ACC->ATC    |                    (1) Node36                    |
|        361         |   CCT->CCG    |          (1) JN563141_isolate_4059_P_51          |
|        367         |   AAG->AGG    |                (2) Node25, Node36                |
|        367         |   AGG->AAG    |(4) JN563119_isolate_4059_P_15, Node68, JN56314...|
|        368         |   GAC->GAA    |(5) JN563119_isolate_4059_P_15, JN563128_isolat...|
|        373         |   GAG->GAA    |      (2) JN563135_isolate_4059_P_41, Node36      |
|        373         |   GAA->GAG    |                    (1) Node22                    |
|        377         |   TAC->TAT    |          (1) JN563161_isolate_4059_C_4           |
|        380         |   ATC->ATT    |(2) JN563128_isolate_4059_P_31, JN563142_isolat...|
|        381         |   AAA->AGA    |(2) JN563157_isolate_4059_C_2, JN563159_isolate...|
|        387         |   ATA->ACA    |          (1) JN563134_isolate_4059_P_4           |
|        392         |   GCA->GCT    |          (1) JN563145_isolate_4059_P_9           |
|        393         |   AAG->AAA    |      (2) JN563151_isolate_4059_C_14, Node49      |
|        399         |   AGA->AAA    |          (1) JN563144_isolate_4059_P_8           |
|        399         |   AAA->AGA    |                    (1) Node14                    |
|        401         |   AAA->AAG    |                    (1) Node33                    |
|        406         |   TTA->TTT    |          (1) JN563164_isolate_4059_C_7           |
|        409         |   ATG->GTG    |(5) JN563151_isolate_4059_C_14, JN563134_isolat...|
|        414         |   TTG->CTG    |          (1) JN563131_isolate_4059_P_35          |
|        414         |   CTG->TTG    |                    (1) Node14                    |
|        426         |   CTG->CTA    |                    (1) Node20                    |
|        426         |   CTA->CTG    |                    (1) Node14                    |
|        435         |   TTA->TTG    |          (1) JN563133_isolate_4059_P_39          |
|        438         |   GGT->GGA    |                    (1) Node90                    |
|        443         |   CAG->CAA    |                (2) Node23, Node36                |
|        443         |   CAA->CAG    |      (2) Node48, JN563129_isolate_4059_P_33      |
|        444         |   AGC->AAC    |          (1) JN563134_isolate_4059_P_4           |
|        445         |   AAT->AAA    |      (2) JN563137_isolate_4059_P_46, Node36      |
|        445         |   AAA->AAG    |      (2) Node48, JN563129_isolate_4059_P_33      |
|        445         |   AAA->AAT    |                    (1) Node20                    |
|        445         |   AAG->AAA    |                    (1) Node14                    |
|        446         |   TTG->CTG    |      (2) Node48, JN563129_isolate_4059_P_33      |
|        446         |   CTG->TTG    |      (2) Node14, JN563152_isolate_4059_C_15      |
|        448         |   AGA->AGG    |                    (1) Node33                    |
|        454         |   CAA->CAG    |(4) Node57, JN563130_isolate_4059_P_34, JN56313...|
|        454         |   CAG->CAA    |                    (1) Node36                    |
|        457         |   TTG->CTG    |                    (1) Node33                    |
|        469         |   GCT->GCA    |                    (1) Node25                    |
|        469         |   GCA->GCT    |(4) Node23, JN563138_isolate_4059_P_47, JN56314...|
|        471         |   GTC->ATC    |      (2) Node15, JN563161_isolate_4059_C_4       |
|        474         |   GTG->GTA    |                    (1) Node7                     |


### Writing detailed analysis report to `/Users/sergei/Development/hyphy-analyses/AncestralSequences/4059.json'

## JSON output. 

### `tree`

Will store the phylogenetic tree with internal node labels in Newick format, which is useful for resolving auto-generated internal node names (like `Node8`)

### `labeled_tree`

The same phylogenetic tree, but with branch lengths used to denote the number of substitutions occuring on every branch

### `ancestral_sequences`

A dictionary that stores, for each internal node of the tree, the reconstructed ancestral sequence, e.g.:

```
"ancestral_sequences":{
   "Node1":"TGTGTCAATCTAAATTGCACTGATAAAGTGACAGTGAACAATTCGACCATAAAGAATACCACCAATGATGACGGGGTAGGAATGATGGATAAAGAAATGAAAAACTGCTCTTTCAATGTTACCACAAATGAAGGAAATAAGGTGCGGAAAGAATATGCACTTCTTTATAAACTTGATGTAGTATCAATAGATGGTAATAAAAATAATAAAGATAATAGTAATAGAAATAATAACACCTATAGTAACTATAGATTGATAAGTTGTAACACCTCAGTTATTACCCAGGCCTGTCCAAAAGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATGATAAGAAGTTCAATGGAAAAGGAGAGTGTAAAAATGTCAGCACAGTACAATGTACACATGGAATTAGACCAGTAGTATCAACTCACTTGTTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGATAATTTCTCAAACAATGCTAAAACTATAATAGTACAGCTGAATAAAACTGTAAAGATTAATTGCACAAGACCCAATAACAATACAAGAAGAAGTATAAGTTTCGGACCAGGGAGAGCATGGCATGCAACAACAGACATAGTAGGAGATATAAGACAAGCACATTGTACCATTAATGGAACAGAATGGAATAACATTTTAAAACTGGTAGTTAGCAAATTACAAGAACAATATGGGACTAATAAAACAATAAGATTTGAGCAACCTGTGCAGGGAGGGGACCTAGAAATTGTAATGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATACATCAAAACTGTTTAATAGTACTTGGAATAATACTGCTAGTACCTGGGATAATGACACTGAAGGAAATGGCACACTCACACTCCCATGTAAAATAAGACAAATTATAAATATGTGGCAAAAGGTAGGAAAAGCAATGTATGCCCCTCCCATCAGCGGACGAATTAACTGCTTATCAAATATTACTGGACTGCTATTAATGAGAGATGGTGGTAGTAATGACACGAACGAGCCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATACAAAGTAATCAAAATTGAACCATTAGGAATAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAAAGAAAAAAGAGCAGTGGGATTAGGAGCTATGTTCCTTGGGTTCCTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGTCACTAACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAACAGCAGAGCAAGCTGCTGAGAGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATAAAACAGCTCCAGGCAAGAGTCCTGGCTGTGGAAGCATACCTAAAGGATCAACAG",
   "Node11":"TGTGTACCTTTAAATTGCACTAAAGTGAAAAAAAAAAAAAACAAAACTAATCCTAATACCACCAATAATGGTGGGGTAGGAATGATGAATGAAGAAATAAAAAACTGCTCTTTCAATGTTACCACAAATGAAGGAAATAAGGTGAAAAAAGAATATGCACTTCTGTATAAACTTGATGTAGTATCAATAGATGGTAGTAGTAGTAAAAAAAAAAAAAAAAAAAAAAAGAATGAAACCTATAGTAACTATAGGTTGATAAGTTGTAACACCTCAGTTATTACCCAGGCCTGTCCAAAAGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATGATAAGAAGTTCAATGGAAAAGGAGAGTGTAAAAATGTCAGCACAGTACAATGTACACATGGAATTAGACCAGTAGTATCAACTCACCTGTTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGACAATTTCTCAAACAATGCTAAAACTATAATAGTACAGCTGAATAAAACTGTAAAGATTAATTGCACAAGACCCAATAACAATACAAGAAGAGGTATAAATTTAGGACAAGGGAGAGCATGGTATGCAACAACAGACATAGTAGGAGATATAAGACAAGCACATTGTAACATTAGTACAACAGAATGGAATAACACTTTAAGACTGGTAGTTAACAAATTACAAGAACAATATGGGACTAATAAAACAATAAAATTTGAGCAACCTGTGCAGGGAGGGGACCCAGAAATTTTAATGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATACATCAAAACTGTTTAATAGTACTTGGGATAATAAAACTAGTACTGGGAAAAATAACACTGAAGAAGATGGCACACTCACACTCCCATGTAAAATAAGGCAATTTATAAATATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGCGGACGGATTAACTGCTTATCAAATATTACTGGGCTGCTATTAATAAGAGATGGTGGTGGTAAAAAAAGTAACGAGACCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATACAAAGTAATCAAAATTGAACCATTAGGAATAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAAAGAAAAAAGAGCAGTGGGATTAGGAGCTATGTTCCTTGGGTTCCTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGTCACTAACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAACAGCAGAGCAAGCTGCTGAGAGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATAAAACAGCTCCAGGCAAGAGTCCTGGCTGTAGAAGCATACCTAAAGGATCAACAG",
   "Node14":"TGTGTACCTTTAAATTGCACTAAAGTGAAAAAAGTGAACAACAAAACTAATCCTAATACCACCAATAATGGTGGGGTAGGAATGATGAATGAAGAAATAAAAAACTGCTCTTTCAATGTTACCACAAATGAAGGAAATAAGGTGAAAAAAGAATATGCACTTCTGTATAAACTTGATGTAGTATCAATAGATGGTAGTAGTAGTAATAAAGATAATAGTAATAGAAAGAATGAAACCTATAGTAACTATAGGTTGATAAGTTGTAACACCTCAGTTATTACCCAGGCCTGTCCAAAAGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATGATAAGAAGTTCAATGGAAAAGGAGAGTGTAAAAATGTCAGCACAGTACAATGTACACATGGAATTAGACCAGTAGTATCAACTCACCTGTTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGACAATTTCTCAAACAATGCTAAAACTATAATAGTACAGCTGAATAAAACTGTAAAGATTAATTGCACAAGACCCAATAACAATACAAGAAGAGGTATAAATTTAGGACAAGGGAGAGCATGGTATGCAACAACAGACATAGTAGGAGATATAAGACAAGCACATTGTACCATTAATAGAACAGAATGGAATAACACTTTAAAACTGGTAGTTAACAAATTACAAGAACAATATGGAACTAATAAAACAATAAAATTTGAGCAACCTGTGCAGGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATACATCAAAACTGTTTAATAGTACTTGGAATAATACTACTAGTACTTGGGATAATAACACTGAAGAAGATGGCACACTCACACTCCCATGTAAAATAAGGCAATTTATAAATATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGCGGACGAATTAACTGCTTATCAAATATTACTGGGCTGCTATTAATGAGAGATGGTGGTGGTGGTGACACGAACGAGACCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATACAAAGTAATCAAAATTGAACCATTAGGAATAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGATTAGGAGCTATGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGTCACTGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAACAGCAGAGCAAATTGCTGAGAGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATAAAACAGCTCCAGGCAAGAGTCCTGGCTGTGGAAGCATACCTAAAGGATCAACAG",
...
```

### `substitution_counts`

An object with three keys

#### `Counts` 

An N (branches) x S (variable sites) matrix whose (i,j)-th entry stores the number of substitutions occuring along branch `i` at site `j`.

#### `Branches`

Row labels

#### `Sites`

Column labels.

For example:

```

"substitution_counts":{
   "Branches":{
     "0":"JN563162_isolate_4059_C_5",
     "1":"JN563149_isolate_4059_C_12",
     "2":"JN562779_isolate_4059_C12_from_USA_envelope_glycoprotein_env_gene_complete_cds",

     ...
     
   "Counts":    [
	 [0, 0, 0, 0, ...],
     [0, 0, 0, 0, ...],
     [1, 1, 0, 0, ...]
    ...
    ], 

	"Sites":{
     "0":1,
     "1":2,
     "2":3,
     ...
    }
```

Means that there is 1 substitution along branch `JN562779_isolate_4059_C12_from_USA_envelope_glycoprotein_env_gene_complete_cds` at sites 2 and 3 (the values in `Sites` are 0-based)

#### `substitution_map`

Stores all individual substitiutions, binned by site, type, and branch. For example, 

```
"substitution_map":{
   "1":{
     "GTA":{
       "GTC":{
         "0":"JN562779_isolate_4059_C12_from_USA_envelope_glycoprotein_env_gene_complete_cds",
         "1":"Node15",
         "2":"JN562781_isolate_4059_P28_from_USA_envelope_glycoprotein_env_gene_complete_cds",
         "3":"JN562780_isolate_4059_P21_from_USA_envelope_glycoprotein_env_gene_complete_cds",
         "4":"JN562784_isolate_4059_P26_from_USA_envelope_glycoprotein_env_gene_complete_cds",
         "5":"JN562783_isolate_4059_C6_from_USA_envelope_glycoprotein_env_gene_complete_cds"
        }
      },
     "GTC":{
       "GTA":{
         "0":"Node2"
        }
      }
    },
    ...
```

Means that at site 2 (again, 0-based), there are two types of substitutions, `GTA->GTC`, which occurs along **6** branches, listed there, and `GTA->GTC` occuring along a single branch: `Node2`.