RequireVersion ("2.4.0");


LoadFunctionLibrary("libv3/all-terms.bf"); 
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


partitioned.analysis_description = {terms.io.info : "Fit several models of partition-level dN/dS variation to _partitioned_ data and perform inference on whether or not there is evidence of different selective pressures among partitions",
                               terms.io.version : "0.1",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment with at least two partitions and a phylogenetic tree"
                              };

io.DisplayAnalysisBanner (partitioned.analysis_description);

KeywordArgument ("code",        "Which genetic code should be used ", "Universal");  
KeywordArgument ("alignment",   "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("branches",   "The set of branches to test", "All");


partitioned.json = {};

namespace partitioned {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ("partitioned");
}

io.CheckAssertion ("partitioned.partition_count>=2","This analysis requires data with at least two partitions");

namespace partitioned {
    doGTR ("partitioned");
}

console.log (partitioned.gtr_results);

/*********************************************************************

Each partition gets its own copy of the **local** MG94xREV model
Models share nucleotide substitution rates, and base frequencies (CF3x4)

*********************************************************************/

partitioned.model_generator = "models.codon.MG94xREV";

namespace partitioned {
    models = {};
    utility.ForEach (filter_names, "_name_", "
        
    ");

}