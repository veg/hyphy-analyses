## Codon-aware MSA

This set of scripts facilitates creating in-frame multiple sequence alignments (MSA) from coding nucleotide sequences, only some of which are in frame. It **requires** the use of a general MSA program (e.g., `muscle` or `MAFFT`) that can handle aligning **protein** sequences. The workflow is as follows

1. Obtain coding sequences (nucleotides). **At least one of them must be in frame**
2. Run these sequences through `pre-msa.bf` in order to correct frame-shift mutations and translate the resulting sequences to proteins.
3. Take the output of step 2 and run in through the general MSA program to generate a **protein** MSA
4. Run the protein MSA and the frameshift corrected nucleotide sequences from step 2 through `post-msa.bf` to obtain a nucleotide msa. This step will also, optionally, compress all identical sequences, i.e. replace them with a single representative sequence.

For step 2, you can use `--E number` (where number is a small floating point value e.g. `0.05`) command line option to permit the alignment of sequences with low homology if you see warning messages about sequences that could not be aligned.

## Example

[1]. Coding sequences (unaligned) in `example.fas`. 

[2]. `hyphy pre-msa.bf --input example.fas` 

```
Analysis Description
--------------------
 Load nucleotide sequences (SOME of which **must be in frame**), perform
frameshift correction as needed, and translate to amino-acids for
subsequent alignment. 

- __Requirements__: Sequences

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.01

code: Universal
input: example.fas
[Data QC] Loaded 28 sequences on 2669 sites from **/Users/sergei/Development/hyphy-analyses/codon-msa/example.fas**
protein: /Users/sergei/Development/hyphy-analyses/codon-msa/example.fas_protein.fas
rna: /Users/sergei/Development/hyphy-analyses/codon-msa/example.fas_nuc.fas
[Data QC] Will write unaligned protein sequences for MSA to **/Users/sergei/Development/hyphy-analyses/codon-msa/example.fas_protein.fas**, and the corresponding nucleotide sequences to **/Users/sergei/Development/hyphy-analyses/codon-msa/example.fas_nuc.fas**
[Data QC] Found 1 unique sequences that were in frame
[Data QC] Correcting frame-shifts in the remaining reads
[Data QC] Checking for frame-preserving indels in other reads
[Next steps] Please run **/Users/sergei/Development/hyphy-analyses/codon-msa/example.fas_protein.fas** through an MSA program, and then run post-msa.bf on the output and **/Users/sergei/Development/hyphy-analyses/codon-msa/example.fas_nuc.fas** to recover the nucleotide MSA
```

[3]. `muscle -in example.fas_protein.fas -out example.fas_protein.msa`

[4].  `hyphy post-msa.bf --protein-msa example.fas_protein.msa --nucleotide-sequences example.fas_nuc.fas --output example.msa`

```
compress: Yes
code: Universal

Analysis Description
--------------------
 Map a protein MSA back onto nucleotide sequences 

- __Requirements__: A protein MSA and the corresponding nucleotide alignment

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.01

Load the protein MSA
protein-msa: example.fas_protein.msa
Load the unaligned in-frame sequences
nucleotide-sequences: example.fas_nuc.fas
output: example.msa
[UNIQUE SEQUENCES] Retained 27 unique  sequences
```


Or [to **not** compress duplicate sequences]

[4] `hyphy post-msa.bf --protein-msa example.fas_protein.msa --nucleotide-sequences example.fas_nuc.fas --output example-all.msa --compress No`

