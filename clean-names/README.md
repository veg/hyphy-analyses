## Convert sequence names

This is a simple script to read in a multiple sequence alignment and a tree, convert all the names to `HyPhy` valid identifiers and write out the alignment and the tree into a new file together. There's an option to match sequences and tips in the tree using a regular expression, so they don't have to be identical (just need to have unique substrings).

### Example

```
hyphy clean-names.bf --msa test/N1.HOG0007354_NT.fasta.fixed --tree test/N1.HOG0007354.treefile --regexp '([A-Z][A-Z0-9a-z_\\.\\-]+)$' --output test/N1.HOG0007354.nex
```

#### Output

``` 
Analysis Description
--------------------
 Read an alignment and a tree and rename all the sequences to conform
with HyPhy naming requirements. The result will be written as a combined
alignment, with the format specified by the DATA_FILE_PRINT_FORMAT HyPhy
environment variable 

- __Requirements__: An MSA and a tree.

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.0.1

Use the following regular expression to select a subset of leaves : regexp: ([A-Z][A-Z0-9a-z_\.\-]+)$

```

### Analysis arguments

```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

msa [required]
	The MSA to process

tree [required]
	The tree

regexp
	Use the following pattern to match sequences and tree labels
	default value: ^(.+)$

output [required]
	Write cleaned MSA+Tree to


```

### Tips

To change the output format, add (in version 2.5.8 or later)

```
ENV="DATA_FILE_PRINT_FORMAT=N"
```

to the command line.

`N` is a number between 0 and 11 (4 is the default, see below). For example, to produce a FASTA output, use 

```
hyphy remove-duplicates.bf --msa example.fas --tree example.nwk 
--output uniques.fas ENV="DATA_FILE_PRINT_FORMAT=9"
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

### Selecting how to match tips and sequences

By default, sequences in the alignments and tips in the tree must match 1-1. But if you have sequences named

```
>ID1_NNN
>ID2_NNN
>ID3_NNN
```

and a tree

```
((ID1|ABC,ID2|ABC), ID3|ABC)
```

supplying the `regexp` argument `^([a-zA-Z0-9]+)_+$` would match them on the prefix alphanumeric IDs for example 