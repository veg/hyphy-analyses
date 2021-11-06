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
        KeywordArgument ("extension","The extension to use for split files", "nex");
        ext = io.PromptUserForString ("The extension to use for split files");
        
        SetDialogPrompt ("The alignment file to partition");
        KeywordArgument ("output","Output parts here");    
        output = io.PromptUserForFilePath ("Output parts here");
        partitions = Transpose(dfpm);
        for (i = 0; i < part_count; i+=1) {
        json_out = {'partitions' : partitions, 'trees' : trees};
        io.SpoolJSON(json_out, output);
            DataSetFilter part_filter = CreateFilter (input_ds, 1, partitions[i][1]);
            if (null != trees)  {
                utility.SetEnvVariable("DATAFILE_TREE", trees[i][1]);
            }
            fprintf (output + "_" + partitions[i][0] + "." + ext, CLEAR_FILE, part_filter);
        }
        return 0;

    }
}

fprintf (stdout, "No data partitions found.");
    
