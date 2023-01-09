RequireVersion ("2.5.46");


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
        joint               = "joint";
        marginal            = "marginal";
        support             = "marginal_support";
        sampled             = "sampled";
        count               = "samples";
        sampled_states      = "inferred";
        sampled_subs        = "sampled_substitution";
    }
}

ancestors.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "left", terms.table_options.column_widths : {"0" : 20, "1" : 20, "2" : 50}};

ancestors.analysis_description = {terms.io.info : "Load a previously generated model fit file (NEXUS + HBL), reconstruct ancestral sequences using several approaches, and map substitutions",
                               terms.io.version : "0.2",
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

KeywordArgument ("min-support", "Minimum marginal support for a character state to be reported", 1e-4);
ancestors.min_support = io.PromptUser ("Minimum marginal support for a character state to be reported", 1e-4, 0, 1, FALSE);

KeywordArgument ("samples", "Draw this many samples to estimate uncertainty (0 to skip)", 100);
ancestors.samples = io.PromptUser ("Draw this many samples to estimate uncertainty (0 to skip)", 100, 0, 10000, TRUE);



KeywordArgument ("format", "Output format", "JSON");
ancestors.format = io.SelectAnOption (
    {"JSON" : "Write detailed output suitable for visualization and exploration to a JSON file (default)",
    "CSV" : "A 'VCF-like' CSV file reporting position-by-position information"}, 
    "Output format");


if (ancestors.format == "JSON") {
    KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + '.ancestors.json')", ancestors.file + ".ancestors.json");
    ancestors.output_file = io.PromptUserForFilePath ("Save the resulting JSON file to");
} else {
   KeywordArgument ("output", "Write the resulting CSV to this file (default is to save to the same path as the alignment file + '.ancestors.csv')", ancestors.file + ".ancestors.csv");
   ancestors.output_file = io.PromptUserForFilePath ("Save the resulting CSV file to");

}

ancestors.trees           =  ancestors.lf_info["Trees"];
ancestors.partition_count =  utility.Array1D (ancestors.trees);
ancestors.json            = {};

io.ReportProgressMessageMD ("ancestral", "loading" , "Loaded a previously saved likelihood function with **`ancestors.partition_count`** partitions from \``ancestors.file`\`");

ancestors.cache           =  ancestral.build (ancestors.lf [0], 0, {});
io.ReportProgressMessageMD("ancestral", "reconstruction" , "Successfully reconstructed joint maximum likelihood ancestors");

ancestors.cache.marginal          =  ancestral.build (ancestors.lf [0], 0, {"marginal" : TRUE});
io.ReportProgressMessageMD("ancestral", "marginal" , "Successfully reconstructed marginal maximum likelihood ancestors");

if (ancestors.samples > 0) {
    io.ReportProgressBar("", "Performing ancestral state sampling");
    ancestors.sampled_states = {};
    ancestors.sampled           =  ancestral.build (ancestors.lf [0], 0, {"sample" : TRUE});
    ancestors.s_states =  ancestral.SequenceStates (ancestors.sampled);
    ancestors.sampled_substitutions = {};
    
    for (seq, states; in;  ancestors.s_states) {
        ancestors.sampled_states [seq] = {};
        for (state; in; states) {
            ancestors.sampled_states [seq] + {state : 1};
        }
    }
    
    ancestors.counts          =  ancestral.ComputeDetailedSubstitutionCounts  (  ancestors.sampled , null, null, null);

    
 
    if (ancestors.format == "JSON") {
        for (seq; in; ancestors.counts ["Branches"]) {
             ancestors.sampled_substitutions [seq] = {};
        }
        for (seq, site, state; in; ancestors.counts ["Substitutions"]) {
            if (state > 0) {
                seq_id = (ancestors.counts ["Branches"])[seq];
                utility.EnsureKey (ancestors.sampled_substitutions [seq_id], site);
                ((ancestors.sampled_substitutions [seq_id])[site])[(ancestors.counts ["Types"])[0][state] + ":" + (ancestors.counts ["Types"])[1][state]] = 1;
            }
        }
    }
    

    for (i = 1; i <  ancestors.samples; i+=1) {
        io.ReportProgressBar("", "Performing ancestral state sampling " + (i) + "/" + ancestors.samples);
        ancestors.sampled           =  ancestral.build (ancestors.lf [0], 0, {"sample" : TRUE});
        ancestors.s_states =  ancestral.SequenceStates (ancestors.sampled);
        for (seq, states; in;  ancestors.s_states) {
            for (s, state; in; states) {
                ((ancestors.sampled_states [seq])[s])[state] += 1;
            }
        }
        if (ancestors.format == "JSON") {
            ancestors.counts          =  ancestral.ComputeDetailedSubstitutionCounts  (  ancestors.sampled , null, null, null);
            for (seq, site, state; in; ancestors.counts ["Substitutions"]) {
                if (state > 0) {
                    seq_id = (ancestors.counts ["Branches"])[seq];
                    utility.EnsureKey (ancestors.sampled_substitutions [seq_id], site);
                    ((ancestors.sampled_substitutions [seq_id])[site])[(ancestors.counts ["Types"])[0][state] + ":" + (ancestors.counts ["Types"])[1][state]] += 1;
                }
            }
        }
        
    }

    io.ReportProgressMessageMD("ancestral", "reconstruction" , "Completed sampling");
    io.ReportProgressMessageMD("ancestral", "samples" , "Successfully generated `ancestors.samples` ancestor samples");

}

if (ancestors.format == "JSON") {
    ancestors.report_results (ancestors.cache, terms.ancestors.joint);
    ancestors.report_results (ancestors.cache.marginal, terms.ancestors.marginal);
    if (ancestors.samples > 0) {
        ancestors.json [terms.ancestors.sampled] = {};
        (ancestors.json [terms.ancestors.sampled])[terms.ancestors.count] = ancestors.samples;
        (ancestors.json [terms.ancestors.sampled])[terms.ancestors.sampled_states] = ancestors.sampled_states;
        (ancestors.json [terms.ancestors.sampled])[terms.ancestors.sampled_subs] = ancestors.sampled_substitutions;
    }
    
    io.ReportProgressMessageMD ("ancestor", "writing", "Writing detailed analysis report to \`" + ancestors.output_file  + "\'");
    io.SpoolJSON (ancestors.json, ancestors.output_file );

} else {
    ancestors.j_states =  ancestral.SequenceStates (ancestors.cache);
    ancestors.m_states =  ancestral.SequenceStates (ancestors.cache.marginal);
    
    io.ReportProgressMessageMD ("ancestor", "writing", "Writing CSV analysis report to \`" + ancestors.output_file  + "\'");
    if (ancestors.samples > 0) {
        fprintf ( ancestors.output_file, CLEAR_FILE, KEEP_OPEN, "Sequence,Site,JointML,MarginalML,MarginalSupport,MarginalAlternative,SampledML,SampledSupport,SampledAlternative\n");
        ancestors.row = {1,9};

    } else {
        fprintf ( ancestors.output_file, CLEAR_FILE, KEEP_OPEN, "Sequence,Site,JointML,MarginalML,MarginalSupport,MarginalAlternative\n");
        ancestors.row = {1,6};
    }
    

    anc.support = ancestral.Support (ancestors.cache.marginal, ancestors.min_support);
    
    
    for (seq, inferred_states; in; ancestors.j_states) {
         ancestors.marginal_support =  anc.support [seq];
         ancestors.row[0] = seq;
         for (site, state; in; inferred_states) {
            ancestors.row[1] = site+1;
            ancestors.row[2] = ""+state;
            ancestors.row[3] = (ancestors.m_states [seq])[site];
            ancestors.marginal_support =  (anc.support [seq])[site];
            ancestors.row[4] =  ancestors.marginal_support   [ancestors.row[3]];
            ancestors.marginal_alt = {};
            for (alt_state, alt_support; in;  ancestors.marginal_support) {
                if (alt_state !=  ancestors.row[3]) {
                    ancestors.marginal_alt + (alt_state + "(" + Format (alt_support, 0, 4) + ")");
                }
            }
            ancestors.row[5] = Join ("|", ancestors.marginal_alt );
            if (ancestors.samples > 0) {
                ancestor.s_max = Max(((ancestors.sampled_states [seq])[site]),1);
                ancestors.row[6] = ancestor.s_max["key"]; 
                ancestors.row[7] = ancestor.s_max["value"] / ancestors.samples; 
                ancestors.marginal_alt = {};
                for (alt_state, alt_support; in; (ancestors.sampled_states [seq])[site]) {
                    if (alt_state !=  ancestors.row[6]) {
                        ancestors.marginal_alt + (alt_state + "(" + Format (alt_support/ancestors.samples, 0, 4) + ")");
                    }
                }
                ancestors.row[8] = Join ("|", ancestors.marginal_alt );
             }
             

           
            
            fprintf ( ancestors.output_file, Join (",",ancestors.row),"\n");
        }
         
    }
    fprintf ( ancestors.output_file, CLOSE_FILE);
}



//------------------------------------------------------------------------------------------------//

function ancestors.substitution_count (node) {
    return  "" + (+(((ancestors.counts["Substitutions"])[+ancestors.branch2id[node["Name"]]][-1])["_MATRIX_ELEMENT_VALUE_>0"]));
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

function ancestors.report_results (anc, json_key) {
    ancestors.json [json_key] = {};
    (ancestors.json [json_key]) [terms.json.tree_string]                 = Format (^(ancestors.trees[0]), 1, 1);
    (ancestors.json [json_key])  [terms.ancestors.ancestral_sequences]    = ancestral.Sequences (anc);

    ancestors.counts          =  ancestral.ComputeDetailedSubstitutionCounts (anc, null, null, null);
    (ancestors.json [json_key]) [terms.ancestors.substitution_counts] = ancestors.counts;

    ancestors.branch2id       =  utility.SwapKeysAndValues (ancestors.counts["Branches"]);

    ancestors.branch_count   =  Rows (ancestors.counts["Substitutions"]);
    ancestors.site_count     =  Columns (ancestors.counts["Substitutions"]);

    io.ReportProgressMessageMD("ancestral", "reconstruction" , "**`ancestors.site_count`** sites have substitutions along at least one branch with `json_key` reconstruction");

    (ancestors.json [json_key]) [terms.ancestors.labeled_tree] = tree.Annotate (ancestors.trees[0], {}, "[]", "ancestors.substitution_count");

    io.ReportProgressMessageMD("ancestral", "branch" , "Phylogenetic tree annotated by the inferred (`json_key`) number of substitutions");
    io.ReportProgressMessageMD ("ancestral", "branch", "\`\`\`\n" +  (ancestors.json [json_key]) [terms.ancestors.labeled_tree] + "\n\`\`\`" );

    ancestors.detailed_report_headers  = {{"Site", "Substitution", "Branches"}};

    fprintf (stdout, "\n", io.FormatTableRow (ancestors.detailed_report_headers,ancestors.table_output_options));
    ancestors.table_output_options[^"terms.table_options.header"] = FALSE;

    ancestors.by_site = ancestors.by_site % 1;

    ancestors.map = {}; // by site
    
    for (site_id; in; Rows(ancestors.counts["Sites"])) {
        site = (ancestors.counts["Sites"])[site_id];
        ancestors.map [site] = (ancestral.ComputeSubstitutionBySite (anc, site, null))["substitutions"];
        ancestors.report_site ((+site) + 1, ancestors.map [site]);
    
    }
    
    anc.support = ancestral.Support (anc, ancestors.min_support);
    if (None != anc.support) {
         (ancestors.json [json_key])[terms.ancestors.support] = anc.support;
    }

    (ancestors.json [json_key])[terms.ancestors.substitution_map] = ancestors.map;
}


