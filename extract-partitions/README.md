## Extract partitons from a NEXUS file

This is a simple script to read in a multiple sequence alignment in NEXUS, and extract all `CHARSET` partitions into separate files.
### Example

```
hyphy extract-partitions.bf --msa example.nex --output output/parts 

```

#### Output

Analysis Description
--------------------
 Read in a multiple sequence alignment in NEXUS format, and extract all
CHARSET partitions into separate files. 

- __Requirements__: A NEXUS file with a CHARSET block

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1

> Loaded an alignment with 25 sequences and 1312 sites from /Users/sergei/Development/hyphy-analyses/extract-partitions/example.nex

There are **9** partitions in a file
The extension to use for split files : extension: nex

```
$ls output 
parts			parts_SPAN_2.nex	parts_SPAN_4.nex	parts_SPAN_6.nex	parts_SPAN_8.nex
parts_SPAN_1.nex	parts_SPAN_3.nex	parts_SPAN_5.nex	parts_SPAN_7.nex	parts_SPAN_9.nex

$more output/parts
{
 "partitions":  [
["SPAN_1", "0-112"],
  ["SPAN_2", "113-343"],
  ["SPAN_3", "344-525"],
  ["SPAN_4", "526-590"],
  ["SPAN_5", "591-723"],
  ["SPAN_6", "724-909"],
  ["SPAN_7", "910-1003"],
  ["SPAN_8", "1004-1195"],
  ["SPAN_9", "1196-1311"] 
  ],
 "trees":  [
 ...
 

```


### Analysis arguments

```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

msa [required]
	The NEXUS file to partition

extension
	The extension to use for split files
	default value: nex

output [required]
	Output parts here
```

### Tips

To change the output format, add (in version 2.5.8 or later)

```
ENV="DATA_FILE_PRINT_FORMAT=N"
```

to the command line.

`N` is a number between 0 and 11 (4 is the default, see below). For example, to produce a FASTA output, use 

```
hyphy extract-partitions.bf --msa example.nex --output output/parts  --extension .fas ENV="DATA_FILE_PRINT_FORMAT=9"
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