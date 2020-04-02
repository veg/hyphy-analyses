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

KeywordArgument ("duplicates", "Save JSON with compressed sequence IDs to", "/dev/null");
filter.duplicates = io.PromptUserForFilePath ("Save identities of sequences that have been collapsed to");
fprintf (filter.duplicates, CLEAR_FILE, KEEP_OPEN);

/*

{
 "SEQUENCE_MAP":  {
{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 4, 0, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 7, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 10, 0, 0, 0, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 0, 0, 0, 0, 0, 14, 0, 4, 15, 15, 16, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 0, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 19, 0, 0, 8, 0, 20, 20, 20, 0, 21, 4, 1, 0, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 22, 10, 10, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 23, 24, 10, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26, 26, 26, 26, 26, 26, 26, 27, 27, 10, 28, 28, 0, 0, 0, 0, 0, 0, 0, 29, 29, 29, 29, 29, 29, 29, 0, 30, 10, 0, 0, 0, 0, 0, 31, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 33, 10, 10, 10, 10, 10, 10, 34, 35, 36, 0, 0, 0, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38, 38, 0, 39, 39, 39, 39, 39, 0, 40, 0, 0, 0, 41, 42, 0, 43, 0, 0, 0, 0, 44, 0, 0, 0, 0, 45, 45, 45, 10, 10, 10, 10, 10, 10, 0, 46, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 10, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 0, 0, 0, 48, 10, 10, 49, 0, 0, 50, 10, 51, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 10, 10, 10, 10, 10, 10, 10, 52, 53, 54, 0, 0, 0, 0, 49, 0, 55, 56, 56, 56, 56, 56, 56, 56, 57, 58, 0, 59, 0, 60, 0, 0, 0, 61, 62, 62, 62, 63, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 64, 45, 65, 65, 0, 0, 0, 10, 0, 0, 0, 10, 66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 67, 68, 69, 10, 10, 70, 10, 0, 0, 0, 0, 0, 0, 10, 0, 57, 71, 72, 57, 10, 73, 73, 73, 74, 75, 75, 0, 0, 10, 76, 77, 78, 78, 78, 78, 10, 10, 10, 10, 10, 79, 79, 79, 79, 10, 0, 59, 61, 10, 10, 10, 10, 10, 10, 10, 80, 0, 81, 0, 82, 83, 10, 0, 0, 55, 84, 0, 39, 10, 0, 0, 0, 10, 85, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 86, 10, 10, 0, 0, 87, 88, 0, 10, 0, 0, 0, 0, 88, 88, 89, 10, 90, 90, 10, 10, 91, 92, 0, 0, 10, 93, 94, 95, 95, 95, 96, 0, 0, 97, 0, 0, 0, 92, 10, 98, 99, 4, 100, 0, 0, 0, 101, 102, 102, 103, 0, 5, 0, 5, 0, 0, 104, 105, 0, 0, 0, 106, 0, 0, 0, 0, 107, 108, 0, 0, 109, 0, 0, 0, 0, 109, 110, 111, 0, 101, 112, 113, 114, 115, 116, 117, 117, 0, 0, 0, 0, 118, 0, 119, 119, 10, 5, 10, 120, 120, 120, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 121, 121, 10, 0, 0, 0, 10, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 122, 54, 0, 123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 0, 0, 10, 0, 10, 124, 0, 10, 10, 0, 10, 10, 125, 125, 125, 125, 125, 125, 10, 126, 10, 0, 10, 0, 127, 0, 0, 128, 0, 129, 10, 10, 0, 0, 0, 10, 120, 10, 81, 47, 0, 10, 47, 10, 0, 0, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 130, 130, 130, 130, 130, 10, 10, 131, 132, 10, 133, 134, 0, 135, 49, 49, 10, 10, 136, 10, 10, 0, 137, 120, 0, 138, 0, 139, 140, 0, 0, 141, 0, 142, 0, 0, 0, 143, 0, 0, 0, 144, 144, 0, 10, 145, 146, 137, 10, 0, 0, 10, 0, 10, 147, 0, 0, 147, 147, 147, 10, 10, 148, 10, 149, 0, 150, 10, 0, 92, 10, 151, 152, 153, 154, 0, 0, 155, 155, 156, 157, 10, 10, 158, 131, 10, 0, 0, 0, 0, 159, 10, 0, 127, 10, 10, 0, 0, 159, 127, 127, 127, 127, 127, 10, 160, 120, 120, 120, 161, 161, 161, 0, 137, 120} 
  },
 "UNIQUE_COUNTS":  {
{500, 15, 2, 15, 5, 11, 13, 1, 9, 1, 189, 2, 25, 18, 1, 2, 1, 1, 9, 1, 3, 1, 1, 1, 1, 1, 7, 2, 2, 7, 1, 1, 1, 1, 1, 1, 1, 1, 2, 7, 1, 1, 1, 1, 1, 4, 1, 24, 1, 4, 1, 1, 1, 1, 2, 2, 7, 3, 1, 2, 1, 2, 3, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 4, 4, 1, 2, 1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 3, 1, 1, 3, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 9, 2, 1, 1, 1, 6, 1, 7, 1, 1, 5, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 2, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 3} 
  },
 "UNIQUE_INDICES":  {
{0, 9, 25, 28, 43, 86, 90, 107, 108, 129, 130, 134, 136, 161, 184, 187, 189, 190, 191, 230, 235, 239, 317, 330, 331, 333, 343, 350, 353, 362, 370, 377, 378, 433, 440, 441, 442, 446, 461, 464, 470, 474, 475, 477, 482, 487, 497, 517, 530, 533, 536, 538, 560, 561, 562, 569, 570, 577, 578, 580, 582, 586, 587, 590, 624, 626, 636, 659, 660, 661, 664, 675, 676, 679, 682, 683, 688, 689, 690, 699, 714, 716, 718, 719, 724, 732, 757, 762, 763, 772, 774, 778, 779, 783, 784, 785, 788, 791, 797, 798, 800, 804, 805, 807, 814, 815, 819, 824, 825, 828, 834, 835, 838, 839, 840, 841, 842, 843, 849, 851, 856, 870, 893, 896, 917, 924, 931, 936, 939, 941, 972, 979, 980, 982, 983, 985, 990, 994, 997, 999, 1000, 1003, 1005, 1009, 1013, 1017, 1018, 1026, 1034, 1036, 1038, 1043, 1044, 1045, 1046, 1049, 1051, 1052, 1055, 1062, 1077, 1081} 
  },
 "UNIQUE_SEQUENCES":162
}

*/

filter.tags = {};

for (filter.id = 0; filter.id < filter.protein[terms.data.sequences]; filter.id += 1) {
    //filter.nuc_seq  = alignments.GetSequenceByName ("filter.nuc_data", _sequence_);
    filter.prot_seq = alignments.GetIthSequence ("filter.protein_data", filter.id);
    filter.nuc_seq = alignments.GetSequenceByName ("filter.nuc_data", filter.prot_seq[terms.id]); 
    /*console.log (filter.prot_seq[terms.id]);
    console.log ( filter.nuc_seq);
    console.log ( filter.prot_seq[terms.data.sequence]); 
    console.log (alignment.MapCodonsToAA(filter.nuc_seq, filter.prot_seq[terms.data.sequence] , 1, filter.code_info[terms.code.mapping]));
    console.log ("\n\n");*/
    filter.seq_tag = filter.prot_seq[terms.id];
    
       
    fprintf (filter.path, ">", filter.seq_tag, "\n", alignment.MapCodonsToAA(filter.nuc_seq, filter.prot_seq[terms.data.sequence] , 1, filter.code_info[terms.code.mapping]), "\n");
}

fprintf (filter.path, CLOSE_FILE);

if (filter.compress != "No") {
    utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT",9);
    filter.aligned = alignments.ReadNucleotideDataSet ("filter.dataset", filter.path);
    DataSetFilter filter.all = CreateFilter (filter.dataset, 1);
    io.ReportProgressMessage ("UNIQUE SEQUENCES", "Retained `alignments.CompressDuplicateSequences ('filter.all','filter.all.compressed', TRUE)` unique  sequences");
    if (filter.duplicates != "/dev/null") {
        GetDataInfo (duplicate_info, filter.all, -2); 
        filter.names = alignments.GetSequenceNames ('filter.all');
        filter.map   = {};
        for (idx,id; in; duplicate_info["SEQUENCE_MAP"]) {
            filter.ref_name =  filter.names[(duplicate_info["UNIQUE_INDICES"])[id]];
            
            if ( filter.map / filter.ref_name == FALSE) {
                filter.map [filter.ref_name] = {};
            }
            filter.map [filter.ref_name] + filter.names[idx];
        }
        fprintf (filter.duplicates,  filter.map);
    }
    fprintf (filter.path, CLEAR_FILE, filter.all.compressed);
}

