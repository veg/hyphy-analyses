LoadFunctionLibrary ("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Read in a multiple sequence alignment in NEXUS format, and extract all CHARSET partitions into separate files.
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A NEXUS file with a CHARSET block"
                          };


io.DisplayAnalysisBanner (filter.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);
utility.SetEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", FALSE);

KeywordArgument                  ("msa", "The NEXUS file to partition");
SetDialogPrompt                  ("The NEXUS file to partition");

DataSet input_ds                 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter input              = CreateFilter (input_ds,1);

console.log ("> Loaded an alignment with `input.species` sequences and `input.sites` sites from `LAST_FILE_PATH`");

dfpm = utility.GetEnvVariable("DATA_FILE_PARTITION_MATRIX");

if (Type(dfpm) == "Matrix") {
    part_count = Columns(dfpm);
    if (part_count > 0) {
        console.log ("\nThere are **`part_count`** partitions in a file");

        KeywordArgument ("target","Which file to partition", "Self");
        part.self = io.SelectAnOption ({{"Self", "Partition this file"}, {"Other", "Partition a different file"}}, "Which file to partition");
        if (part.self == "Self") {
            trees = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");
            utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", FALSE);
            if (Type (trees) == "Matrix") {
                if (Rows (trees) != part_count) {
                    trees= null;
                } else {
                    utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", TRUE);
                 }
            } else {
                trees = null;
            }
        } else {
            KeywordArgument                  ("target-msa", "The file to partition");
            DataSet input_ds                 = ReadDataFile (PROMPT_FOR_FILE);
            DataSetFilter input              = CreateFilter (input_ds,1);
            trees = null;
        }

        KeywordArgument ("extension","The extension to use for split files", "nex");
        ext = io.PromptUserForString ("The extension to use for split files");
        
        KeywordArgument ("sites","Which sites to extract", "All");
        part.variable = io.SelectAnOption ({{"All", "Extract all sites"}, {"Variable", "Extract only variable sites"}}, "Which sites to extract");
        
        SetDialogPrompt ("The alignment file to partition");
        KeywordArgument ("output","Output parts here");    
        output = io.PromptUserForFilePath ("Output parts here");
        partitions = Transpose(dfpm);
        
        function variable_sites (s,freqs) {
            return (+freqs["_MATRIX_ELEMENT_VALUE_>0"]) > 1;
        }
        
        for (i = 0; i < part_count; i+=1) {
            json_out = {'partitions' : partitions, 'trees' : trees};
            io.SpoolJSON(json_out, output);
            if ( part.variable == "All") {
                DataSetFilter part_filter = CreateFilter (input_ds, 1, partitions[i][1]);
            } else {
                DataSetFilter part_filter1 = CreateFilter (input_ds, 1, partitions[i][1]); 
                Export (eds, part_filter1);
                DataSet export_ds = ReadFromString (eds);
                DataSetFilter part_filter = CreateFilter (export_ds, 1, "variable_sites");    
            }
            if (null != trees)  {
                utility.SetEnvVariable("DATAFILE_TREE", trees[i][1]);
            }
            fprintf (output + "_" + partitions[i][0] + "." + ext, CLEAR_FILE, part_filter);
        }
        return 0;

    }
}

fprintf (stdout, "No data partitions found.");
    
