RequireVersion ("2.5.44");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");


LoadFunctionLibrary     ("libv3/tasks/alignments.bf");

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */

LoadFunctionLibrary("TemplateModels/chooseGeneticCode.def");



outliers.analysis_description = {terms.io.info :
                            "
                            Read a SLAC result file, and report a heuristic scale with problematic sequence-level stretches.
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A SLAC JSON"
                          };


io.DisplayAnalysisBanner (outliers.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);

KeywordArgument ("slac", "A SLAC result file (JSON)");
outliers.slac.json = io.ParseJSON(io.PromptUserForFilePathRead ("A SLAC result file (JSON)"));
io.CheckAssertion("(outliers.slac.json[terms.json.input])[terms.json.partition_count] == 1", "Only single partition SLAC files are supported");

KeywordArgument ("window-width", "The width of a sliding window to use for filtering", 5);
outliers.min_window = io.PromptUser ("The width of a sliding window to use for filtering", 5, 3, 100, TRUE);

KeywordArgument ("window-fraction", "The fraction of multiple-hits in a window used to filtering", 0.4);
outliers.fraction = io.PromptUser ("The fraction of multiple-hits in a window used to filtering", 0.4, 0.05, 1.0, FALSE);

KeywordArgument ("hit-threshold", "Number of substitution on a branch to count towards 'filtering' (e.g. 2 for 2 or more changes)", 1.5);
outliers.slac.window_freq = io.PromptUser ("The fraction of multiple-hits in a window used to filtering", 1.5, 0, 3, FALSE);

KeywordArgument ("filter-partials", "Filter partially resolved codons", "No");

outliers.filter_partials = io.SelectAnOption (
    {"No" : "Keep partial (e.g. -GC) codons as they are", 
    "Yes" : "Replace partial codons with ---"}, 
    "Filter partial codons") == "Yes";

outliers.suspect_sites_by_seq = {};
outliers.treeS = ((outliers.slac.json[terms.json.input])[terms.json.trees])[0];

Topology T = outliers.treeS;
outliers.parents_by_seq = trees.ParentMap ("T");
outliers.internals = {};

outliers.descendants_by_node = {};

for (n, p; in; outliers.parents_by_seq) {
    outliers.internals[p] = 1;
}

outliers.seq_suspect_ranges = {};

outliers.tree_root = BranchName (T, BranchCount (T));
 
for (n; in; T) {
    if (n != outliers.tree_root) {
        if (outliers.internals[n]) {
            outliers.my_children = outliers.descendants_by_node[n];
        } else {
            outliers.my_children = {n : 1};
        }
        p = outliers.parents_by_seq [n];
        if (outliers.descendants_by_node / p == FALSE) {
            outliers.descendants_by_node [p] = {};
        }
        outliers.descendants_by_node [p] * outliers.my_children;
        outliers.seq_suspect_ranges [n] = {};
    }
}


for (branch, bdata; in; (outliers.slac.json[terms.json.branch_attributes])["0"]) {
    if (outliers.parents_by_seq[branch]) {
        outliers.suspect_sites_by_seq[branch] = {};
        for (index, count; in; bdata[terms.nonsynonymous_sub_count]) {
            outliers.total_count = count + (bdata[terms.synonymous_sub_count])[index];
            if (outliers.total_count >= outliers.slac.window_freq) {
                if (outliers.suspect_sites_by_seq / branch == 0) {
                    outliers.suspect_sites_by_seq[branch] = {};
                }
                (outliers.suspect_sites_by_seq[branch])[index] = 1;
            }
        }
    }
}


outliers.filter.sites = (outliers.slac.json[terms.json.input])[terms.json.sites];

//console.log (outliers.suspect_sites_by_seq["VS_BALACU1"]);

for (seq_idx, seq_mh; in; outliers.suspect_sites_by_seq) {
    counter      = seq_mh[0];
    window_start = 0;
    window_span  = 1;
    my_children = outliers.descendants_by_node[seq_idx];
    for (i = 1; i < outliers.filter.sites; i+=1) {
               
        if (window_span < outliers.min_window) {
            counter += seq_mh[i];
            window_span += 1;
        } else {
            if (counter / window_span >= outliers.fraction) {
               counter += seq_mh[i];
               window_span += 1;
            } else {
                if (window_span > outliers.min_window) {
                    for (j = window_start; j < window_start + window_span - 1; j+=1) {
                        (outliers.seq_suspect_ranges[seq_idx])[j] = 1;
                    }
                    
                    if (Abs (my_children)) {
                        for (c,ignore; in; my_children) {
                            for (j = window_start; j < window_start + window_span - 1; j+=1) {
                                (outliers.seq_suspect_ranges[c])[j] = 1;
                            }
                        }
                    }
                    
                    //outliers.seq_suspect_ranges[seq_idx] + {{window_start__, window_start__ + window_span__ - 2}};
                    window_start = i;
                    counter = seq_mh[i];
                    window_span = 1;
                } else {
                    counter += (seq_mh[i]-seq_mh[window_start]);
                    window_start += 1;
                }
            }
        }
    }
    
    window_span += -1;
    if (counter / window_span >= outliers.fraction) {
        for (j = window_start; j < window_start + window_span - 1; j+=1) {
            (outliers.seq_suspect_ranges[seq_idx])[j] = 1;
        }
    }
}



outliers.filtered = "";
outliers.filtered * 8192;

lfunction outliers.filter (code, filter) {    
    if (filter) { return "---" };
    return code;
}

lfunction outliers.filter_partials_callback (code) {    
    if (code[0] == '-' || code[1] == '-' || code [2] == '-') {
        return '---';
    }
    return code;
}

// Filter out internal nodes for now, but leave the possibility open to export a larger JSON
outliers.seq_suspect_ranges_tips = {};

for (outliers.seq_name, seq_ranges; in; outliers.seq_suspect_ranges) {
   if (outliers.internals [outliers.seq_name]) {
        continue;
   } 
   seq.info   = ((outliers.slac.json[terms.json.branch_attributes])["0"])[outliers.seq_name];
   seq.codons = seq.info[terms.codon];
   if (outliers.filter_partials) {
        seq.codons = seq.codons ["outliers.filter_partials_callback(_MATRIX_ELEMENT_VALUE_)"];
   }
   if (Abs (seq_ranges))  {
        outliers.seq_suspect_ranges_tips[outliers.seq_name] = utility.Keys(seq_ranges);
        seq.filter = {};      
        console.log ("Filtering " + (+seq_ranges) + " sites from `outliers.seq_name`");
        outliers.filtered * (">" + outliers.seq_name + "\n" + Join ("", seq.codons["outliers.filter(_MATRIX_ELEMENT_VALUE_,seq_ranges[_MATRIX_ELEMENT_COLUMN_])"]) + "\n");
     } else {
        outliers.filtered * (">" + outliers.seq_name + "\n" + Join ("", seq.codons) + "\n");
    }
}

outliers.filtered * 0;

KeywordArgument ("output", "Write filtering results to");
outliers.outpath = io.PromptUserForFilePath ("Write filtering results to");
fprintf (outliers.outpath, CLEAR_FILE, outliers.filtered);

DataSet outliers.check = ReadDataFile (outliers.outpath);

lfunction outliers.non_gap (site, freq) {
    return +freq ["_MATRIX_ELEMENT_VALUE_>0"];
}

utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT",9);
DataSetFilter outliers.check.filter = CreateFilter (outliers.check, 3, "outliers.non_gap");
fprintf (outliers.outpath, CLEAR_FILE, outliers.check.filter);

KeywordArgument ("outlier-coord-output", "Write outlier coordinates to");
outliers.coord_outpath = io.PromptUserForFilePath ("Write outlier coordinates to");
io.SpoolJSON(outliers.seq_suspect_ranges_tips, outliers.coord_outpath);
fprintf (outliers.coord_outpath, CLOSE_FILE);


