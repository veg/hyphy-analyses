RequireVersion ("2.3.12");


LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Map a protein MSA back onto nucleotide sequences
                            ",
                            terms.io.version :          "0.01",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A protein MSA and the corresponding nucleotide alignment"
                          };


KeywordArgument ("compress", "Only retain one copy of identical (at the nucleotide level sequences)", "Yes");

filter.compress = io.SelectAnOption ({
                                        "Yes" : "Replace identical sequences with a single copy; append copy number as :N", 
                                        "No"  : "Retain all sequences",
                                      }, 
                                      "Only retain one copy of identical (at the nucleotide level sequences)");


KeywordArgument ("code", "Genetic code", "Universal");
filter.code_info = alignments.LoadGeneticCode (null);

io.DisplayAnalysisBanner (filter.analysis_description);

KeywordArgument ("protein-msa", "Protein MSA");
console.log ("Load the protein MSA");

filter.protein = alignments.ReadProteinDataSet ("filter.protein_data", None);

KeywordArgument ("nucleotide-sequences", "In-frame nucleotide sequences");
//console.log ("Loading nucleotide unaligned sequences");
console.log ("Load the unaligned in-frame sequences");
filter.nuc = alignments.ReadNucleotideDataSet ("filter.nuc_data", None);

alignments.GetSequenceByName ("filter.nuc_data", None);

KeywordArgument ("output", "Write nucleotide MSA to");
filter.path = io.PromptUserForFilePath ("Save nucleotide in-frame MSA to ");
fprintf (filter.path, CLEAR_FILE, KEEP_OPEN);

filter.tags = {};

for (filter.id = 0; filter.id < filter.protein[terms.data.sequences]; filter.id += 1) {
    //filter.nuc_seq  = alignments.GetSequenceByName ("filter.nuc_data", _sequence_);
    filter.prot_seq = alignments.GetIthSequence ("filter.protein_data", filter.id);
    filter.nuc_seq = alignments.GetSequenceByName ("filter.nuc_data", filter.prot_seq[terms.id]);    
    filter.seq_tag = filter.prot_seq[terms.id];
    
       
    fprintf (filter.path, ">", filter.seq_tag, "\n", alignment.MapCodonsToAA(filter.nuc_seq, filter.prot_seq[terms.data.sequence] , 1, filter.code_info[terms.code.mapping]), "\n");
}

fprintf (filter.path, CLOSE_FILE);

if (filter.compress != "No") {
    utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT",9);
    filter.aligned = alignments.ReadNucleotideDataSet ("filter.dataset", filter.path);
    DataSetFilter filter.all = CreateFilter (filter.dataset, 1);
    io.ReportProgressMessage ("UNIQUE SEQUENCES", "Retained `alignments.CompressDuplicateSequences ('filter.all','filter.all.compressed', TRUE)` unique  sequences");
    fprintf (filter.path, CLEAR_FILE, filter.all.compressed);
}

