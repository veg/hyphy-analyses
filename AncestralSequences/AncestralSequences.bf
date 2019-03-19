RequireVersion ("2.4.0");


LoadFunctionLibrary("libv3/all-terms.bf"); 
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/ancestral.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

namespace terms {
    namespace ancestors {
        ancestral_sequences = "ancestral_sequences";
        substitution_counts = "substitution_counts";
        substitution_map    = "substitution_map";
        labeled_tree        = "labeled_tree";
    }
}

ancestors.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "left", terms.table_options.column_widths : {"0" : 20, "1" : 15, "2" : 50}};

ancestors.analysis_description = {terms.io.info : "Load a previously generated model fit file (NEXUS + HBL), reconstruct ancestral sequences, and map substitutions",
                               terms.io.version : "0.1",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "a fit file generated previously by a different HyPhy analysis"
                              };

io.DisplayAnalysisBanner (ancestors.analysis_description);

 
KeywordArgument ("fit",  "Load a previously saved HyPhy fit (NEXUS + HBL)" );  
ancestors.file = io.LoadFile     ("Load a previously saved HyPhy fit (NEXUS + HBL)");

ancestors.lf = utility.GetListOfLoadedLikelihoodFunctions(null);

io.CheckAssertion           ("utility.Array1D(ancestors.lf) == 1", "This analysis requires a fit file with exactly one likelihood function");
GetString                   (ancestors.lf_info, ^(ancestors.lf [0]), -1);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + '.ancestors.json')", ancestors.file + ".ancestors.json");
ancestors.json_file = io.PromptUserForFilePath ("Save the resulting JSON file to");

ancestors.trees           =  ancestors.lf_info["Trees"];
ancestors.partition_count =  utility.Array1D (ancestors.trees);
ancestors.json            = {};

io.ReportProgressMessageMD ("ancestral", "loading" , "Loaded a previously saved likelihood function with **`ancestors.partition_count`** partitions from \``ancestors.file`\`");

ancestors.cache           =  ancestral.build (ancestors.lf [0], 0, {});
io.ReportProgressMessageMD("ancestral", "reconstruction" , "Successfully reconstructed joint maximum likelihood ancestors");

ancestors.json [terms.json.tree_string]                 = Format (^(ancestors.trees[0]), 1, 1);
ancestors.json [terms.ancestors.ancestral_sequences]    = ancestral.Sequences (ancestors.cache);

ancestors.counts          =  ancestral.ComputeSubstitutionCounts (ancestors.cache, null, null, null);
ancestors.json [terms.ancestors.substitution_counts] = ancestors.counts;

ancestors.branch2id       =  utility.SwapKeysAndValues (ancestors.counts["Branches"]);

ancestors.branch_count   =  Rows (ancestors.counts["Counts"]);
ancestors.site_count     =  Columns (ancestors.counts["Counts"]);

io.ReportProgressMessageMD("ancestral", "reconstruction" , "**`ancestors.site_count`** sites have substitutions along at least one branch");


ancestors.json[terms.ancestors.labeled_tree] = tree.Annotate (ancestors.trees[0], {}, "[]", "ancestors.substitution_count");

io.ReportProgressMessageMD("ancestral", "branch" , "Phylogenetic tree annotated by the inferred number of substitutions");
io.ReportProgressMessageMD ("ancestral", "branch", "\`\`\`\n" + ancestors.json[terms.ancestors.labeled_tree] + "\n\`\`\`" );

ancestors.detailed_report_headers  = {{"Site", "Substitution", "Branches"}};

fprintf (stdout, "\n", io.FormatTableRow (ancestors.detailed_report_headers,ancestors.table_output_options));
ancestors.table_output_options[^"terms.table_options.header"] = FALSE;

ancestors.by_site = ancestors.by_site % 1;

ancestors.map = {}; // by site

utility.ForEach(ancestors.counts["Sites"], "site",
'
    ancestors.map [site] = (ancestral.ComputeSubstitutionBySite (ancestors.cache, site, null))["substitutions"];
    ancestors.report_site (site + 1, ancestors.map [site]);
    
');

ancestors.json[terms.ancestors.substitution_map] = ancestors.map;


io.ReportProgressMessageMD ("ancestor", "writing", "Writing detailed analysis report to \`" + ancestors.json_file  + "\'");
io.SpoolJSON (ancestors.json, ancestors.json_file );


//------------------------------------------------------------------------------------------------//

function ancestors.substitution_count (node) {
    return  "" + (+((ancestors.counts["Counts"])[+ancestors.branch2id[node["Name"]]][-1]));
}

//------------------------------------------------------------------------------------------------//

lfunction ancestors.report_site (site, counts) {
    row = {1,3};
    row [0] = "" + site;
    
    
    from = utility.Keys (counts);
    c1 = utility.Array1D (from);
     for (i = 0; i < c1; i += 1) {
        from_state = from [i];
        to = counts[from_state];
        to_states = utility.Keys (to);
        c2 = utility.Array1D (to_states);
        for (j = 0; j < c2; j += 1) {
            to_state = to_states[j];
            row[1] = from_state + "->" + to_state;
            row[2] = "(" + Abs (to[to_state]) + ") "+(Join (", ", to[to_state]));
            fprintf (stdout, io.FormatTableRow (row,^"ancestors.table_output_options"));
        }
    }
     
}



