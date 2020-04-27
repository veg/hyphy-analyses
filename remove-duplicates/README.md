## Remove duplicate sequences

This is a simple script to read in a multiple sequence alignment (and optionally a tree), remove all duplicate sequences (keep one copy), and trim the tree file accordingly. The resulting de-duplicated alignment (and trimmed tree, if provided) is saved to a new file.

### Example

```
hyphy remove-duplicates.bf --msa example.fas --tree example.nwk --output uniques.nxh
```

#### Output

Analysis Description
--------------------
 Read an alignment (and, optionally, a tree) remove duplicate sequences,
and prune the tree accordingly. 

- __Requirements__: An MSA and, optionally, a tree

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1

> Loaded an alignment with 90 sequences and 1260 sites from /Users/sergei/Development/hyphy-analyses/remove-duplicates/example.fas

There are **11** unique sequences in alignment 
An optional tree file to trim : tree: example.nwk


### Analysis arguments

```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

msa [required]
	The MSA to remove duplicate sequences from

tree
	An optional tree file to trim
	default value: None

output [required]
	Write de-duplicated MSA to
```

### Tips

To change the output format, add (in version 2.5.8 or later)

```
ENV="DATA_FILE_PRINT_FORMAT=N"
```

to the command line.

`N` is a number between 0 and 11 (4 is the default, see below). For example, to produce a FASTA output, use 

```
hyphy remove-duplicates.bf --msa example.fas --tree example.nwk --output uniques.fas ENV="DATA_FILE_PRINT_FORMAT=9"
```

```
kFormatMEGASequential             = 0,
kFormatMEGAInterleaved            = 1,
kFormatPHYLIPSequential           = 2,
kFormatPHYLIPInterleaved          = 3,
kFormatNEXUSLabelsSequential      = 4,
kFormatNEXUSLabelsInterleaved     = 5,
kFormatNEXUSSequential            = 6,
kFormatNEXUSInterleaved           = 7,
kFormatCharacterList              = 8,
kFormatFASTASequential            = 9,
kFormatFASTAInterleaved           = 10,
kFormatPAML                       = 11

```