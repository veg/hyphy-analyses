RequireVersion ("2.4");


LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("lib/igscueal.bf");


utility.SetEnvVariable  ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable  ("ACCEPT_ROOTED_TREES", TRUE);


filter.analysis_description = {terms.io.info :
                            "
                                Load nucleotide sequences (SOME of which **must be in frame**),
                                perform frameshift correction as needed, and translate to amino-acids
                                for subsequent alignment.
                            ",
                            terms.io.version :          "0.01",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "Sequences"
                          };
                          
io.DisplayAnalysisBanner (filter.analysis_description);

KeywordArgument ("code", "Genetic code", "Universal");
filter.code_info = alignments.LoadGeneticCode (null);
KeywordArgument ("input", "A FASTA file with sequences (at least some must be full length ORF in one of three frames)");
filter.nuc_data = alignments.ReadNucleotideDataSet ("filter.raw_data", None);

io.ReportProgressMessage ("Data QC", "Loaded `filter.nuc_data[terms.data.sequences]` sequences on `filter.nuc_data[terms.data.sites]` sites from **`filter.nuc_data[terms.data.file]`**");

filter.all_sequences = alignments.GetSequenceNames ("filter.raw_data");


/* convert sequence names from RNA reads */

filter.RNA_reads     = {};
filter.lookup_cache  = {};
filter.seq_count     = 1;
filter.sequence_info = {};
filter.longest_seq   = "";
filter.longest_seq_L = 0;
filter.clean_seqs   = {};
filter.frameshifted = {};
filter.unique       = {};


KeywordArgument ("protein", "Translated sequences", filter.nuc_data[terms.data.file] + "_protein.fas");
filter.protein_path = io.PromptUserForFilePath ("Save translated RNA sequences file to");

KeywordArgument ("rna", "Corrected nucleotide sequences", filter.nuc_data[terms.data.file] + "_nuc.fas");
filter.nuc_path = io.PromptUserForFilePath ("Save reduced RNA sequences file to");

io.ReportProgressMessage ("Data QC", "Will write unaligned protein sequences for MSA to **`filter.protein_path`**, and the corresponding nucleotide sequences to **`filter.nuc_path`**");

fprintf (filter.protein_path, CLEAR_FILE, KEEP_OPEN);
fprintf (filter.nuc_path, CLEAR_FILE, KEEP_OPEN);
filter.sequences_with_copies = {};

alignments.GetSequenceByName ("filter.raw_data", null);

utility.ForEach (filter.all_sequences , "_seq_record_", 
'
    io.ReportProgressBar ("filter","Processing sequence " + filter.seq_count);
    filter.read_to_check = alignments.Strip (alignments.StripGaps (alignments.GetSequenceByName ("filter.raw_data", _seq_record_)));
    
    if (filter.unique / filter.read_to_check == FALSE) {    
        
        filter.RNA_reads[_seq_record_] = filter.read_to_check;
        
        filter.sequences_with_copies [filter.read_to_check] = {"0" : _seq_record_};
        filter.unique [filter.read_to_check] = _seq_record_;
    
        filter.sequence_info[_seq_record_] = alignments.TranslateCodonsToAminoAcidsWithAmbigsAllFrames (filter.RNA_reads[_seq_record_],
                               filter.code_info, filter.lookup_cache);
                               
                                                                                                 
        for (frame = 0; frame < 3; frame += 1) {

            if (((filter.sequence_info[_seq_record_])[frame])[terms.stop_codons] == 0) {
                if (((filter.sequence_info[_seq_record_])[frame])[terms.sense_codons] > filter.longest_seq_L) {
                    filter.longest_seq_L = ((filter.sequence_info[_seq_record_])[frame])[terms.sense_codons];
                    filter.longest_seq = ( filter.RNA_reads[_seq_record_]) [frame][ Abs(filter.RNA_reads[_seq_record_]) - 1];
                    filter.longest_seq_NL = Abs(filter.longest_seq);
                    if ( filter.longest_seq_NL % 3) {
                        filter.longest_seq_NL = filter.longest_seq_NL$3*3;
                        filter.longest_seq = filter.longest_seq[0][filter.longest_seq_NL-1];
                    }
                }
                filter.clean_seqs [_seq_record_] = ((filter.sequence_info[_seq_record_])[frame])[terms.data.sequence];
                break;
            }
        }   
                           
        if (frame == 3) {
            filter.frameshifted [_seq_record_] = 1;
        }  
    } else {
        filter.sequences_with_copies [filter.read_to_check] + _seq_record_;
    }
                           
    filter.seq_count += 1;
');

io.ClearProgressBar ();
io.ReportProgressMessage ("Data QC", "Found `Abs(filter.clean_seqs)` unique sequences that were in frame");
io.CheckAssertion ("filter.longest_seq_L>0", "There were no sequences that were in frame and had no stop codons");


filter.ref_seq = {"REF" : {'stripped' : filter.longest_seq}};
filter.options = IgSCUEAL.define_alignment_settings (filter.code_info);   
filter.options["code"] = filter.code_info;


function _write_to_file (key, value) {
    fprintf (filter.protein_path, ">", value, "\n",  filtered.aa_seq, "\n");
    fprintf (filter.nuc_path, ">", value, "\n", filtered.na_seq , "\n");
}

if (Abs(filter.frameshifted)) {
    io.ReportProgressMessage ("Data QC", "Correcting frame-shifts in the remaining reads");

   filter.seq_count = 1;
   utility.ForEachPair (filter.frameshifted, "_sequence_", "_value_",
    '
        io.ReportProgressBar ("filter","Processing sequence " + filter.seq_count);
        filter.cleaned = IgSCUEAL.align_sequence_to_reference_set (filter.RNA_reads[_sequence_], filter.ref_seq, filter.options);
        
        filtered.aa_seq = alignments.StripGaps(filter.cleaned["AA"]);
        filtered.na_seq = IgSCUEAL.strip_in_frame_indels(filter.cleaned["QRY"]);
                
        (filter.sequences_with_copies[filter.RNA_reads[_sequence_]])["_write_to_file"][""];
                
        filter.seq_count += 1;
    ');
    io.ClearProgressBar ();
    
} 



io.ReportProgressMessage ("Data QC", "Checking for frame-preserving indels in other reads");

filter.seq_count = 1;

utility.ForEachPair (filter.clean_seqs, "_sequence_", "_value_",
'
    io.ReportProgressBar ("filter","Processing sequence " + filter.seq_count);
    filter.cleaned = IgSCUEAL.align_sequence_to_reference_set (filter.RNA_reads[_sequence_], filter.ref_seq, filter.options);
    
    filtered.aa_seq = alignments.StripGaps(filter.cleaned["AA"]);
    filtered.na_seq = IgSCUEAL.strip_in_frame_indels(filter.cleaned["QRY"]);
    
    (filter.sequences_with_copies[filter.RNA_reads[_sequence_]])["_write_to_file"][""];
 
    filter.seq_count += 1;
    
');

io.ClearProgressBar ();

fprintf (filter.protein_path, CLOSE_FILE);
fprintf (filter.nuc_path, CLOSE_FILE);

io.ReportProgressMessage ("Next steps", "Please run **`filter.protein_path`** through an MSA program, and then run post-msa.bf on the output and **`filter.nuc_path`** to recover the nucleotide MSA");






