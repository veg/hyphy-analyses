RequireVersion ("2.5.52");


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

fitter.analysis_description = {terms.io.info : "Fit an MG94xREV model with several selectable options frequency estimator and 
report the fit results including dN/dS ratios, and synonymous and non-synonymous branch lengths. v0.2 adds LRT test for dN/dS != 1. 
v0.3 adds LRT test support for dN/dS != 1 for local models. v0.4 adds the lineage option",
                               terms.io.version : "0.4",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree"
                              };

io.DisplayAnalysisBanner (fitter.analysis_description);

  
namespace fitter.terms {
    MG94 = "Standard MG94";
    LRT = "LRT";
}  

terms.fitter.ci = "Confidence Intervals";
terms.fitter.lrt = "LRT";
terms.fitter.partitioned = "partitioned";
terms.fitter.lineage = "lineage";
terms.fitter.MLE = "Lineage dN/dS";

 
KeywordArgument ("rooted", "Accept rooted trees", "No");
KeywordArgument ("type", "Model type: 
        - global (single dN/dS for all branches) 
        - local (separate dN/dS) 
        - partitioned (separate dN/dS for user-defined groups)
        - lineage (separate dN/dS for each tip's lineage)", terms.global, "Model Type");
KeywordArgument ("code",        "Which genetic code should be used", "Universal");  
KeywordArgument ("alignment",   "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",        "A phylogenetic tree", null, "Please select a tree file for the data:");
KeywordArgument ("frequencies", "Equilibrium frequency estimator", "CF3x4");

fitter.json    = { terms.json.analysis: fitter.analysis_description,
                   terms.json.input: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
               };
 
fitter.display_orders = {terms.original_name      :  -1,
                         terms.json.nucleotide_gtr: 0,
                         fitter.terms.MG94        : 1,
                         fitter.terms.dS          : 2,
                         fitter.terms.dN          : 3,
                         fitter.terms.LRT         : 2
                        };


selection.io.startTimer (fitter.json [terms.json.timers], "Overall", 0);

fitter.accept_rooted_trees = io.SelectAnOption ({"Yes" : "Accept rooted trees", 
                                            "No" : "Automatically reroot trees"}, "Accept rooted trees if present");

if (fitter.accept_rooted_trees == "Yes") {
    utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);
}

fitter.model_type = io.SelectAnOption ({
    terms.global : "Shared dN/dS for all branches", 
    terms.local : "Each branch has its own dN and dS", 
    terms.fitter.partitioned : "Each branch partition has its own dN/dS",
    terms.fitter.lineage : "Iterate over lineages (root-to-tip) for each tip, and estimate the dN/dS shared by lineage branches. Other branches have separate (nuisance) dN/dS values"}, 
"Model Type");


namespace fitter {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    if (model_type == ^"terms.fitter.partitioned") {
        load_file ({utility.getGlobalValue("terms.prefix"): "fitter", 
            utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranchSets"}});
    } else {
        load_file ({utility.getGlobalValue("terms.prefix"): "fitter", 
            utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
    }
}


fitter.frequency_type = io.SelectAnOption ({"CF3x4" : terms.frequencies.CF3x4, 
                                            "F3x4" : terms.frequencies.F3x4,
                                            "F1x4" : terms.frequencies.F1x4}, "Equilibrium frequency estimator");

KeywordArgument ("lrt",         "Perform LRT to test which for dN/dS == 1", "No");
fitter.compute_lrt = io.SelectAnOption ({"No"  : "Do not perform LRT",
                                        "Yes" : "Perform LRT to test omega == 1"}, "Perform LRT to test omega != 1") != "No";



KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')", fitter.codon_data_info [terms.json.json]);
fitter.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

namespace fitter {
    doGTR ("fitter");
}

io.ReportProgressMessageMD ("fitter", fitter.terms.MG94,  "Fitting `fitter.terms.MG94`");
selection.io.startTimer (fitter.json [terms.json.timers], fitter.terms.MG94 , fitter.display_order [fitter.terms.MG94 ]);

lfunction fitter.defineMG (type, code) {
    m = Call ("models.codon.MG_REV.ModelDescription", type, code);
    if (^"fitter.frequency_type" == "F3x4") {
        m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F3x4";
    } else {
        if (^"fitter.frequency_type" == "F1x4") {
            m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F1x4";
        } 
    }
    return m;
}

if (fitter.model_type == terms.fitter.partitioned) {
    fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, "fitter.defineMG", fitter.codon_data_info [utility.getGlobalValue("terms.code")],
        {
            terms.run_options.model_type            : terms.local,
            terms.run_options.retain_lf_object      : TRUE,
            terms.run_options.retain_model_object   : TRUE,
            terms.run_options.partitioned_omega     : fitter.selected_branches
        }, 
        fitter.gtr_results);
        
} else {
    
    fitter.mta = fitter.model_type;
    if (fitter.model_type == terms.fitter.lineage) {
        fitter.mta = terms.local;
    }
    
    fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, "fitter.defineMG", fitter.codon_data_info [utility.getGlobalValue("terms.code")],
        {
            terms.run_options.model_type: fitter.mta,
            terms.run_options.retain_lf_object: TRUE,
            terms.run_options.retain_model_object : TRUE
        }, 
        fitter.gtr_results);
}

io.ReportProgressMessageMD("fitter", fitter.terms.MG94 , "* " + selection.io.report_fit (fitter.results, 0, fitter.codon_data_info[terms.data.sample_size]));
fitter.global_dnds = selection.io.extract_global_MLE_re (fitter.results, terms.parameters.omega_ratio);

selection.io.json_store_lf (fitter.json,
                            fitter.terms.MG94 ,
                            fitter.results[terms.fit.log_likelihood],
                            fitter.results[terms.parameters],
                            fitter.sample_size,
                            utility.Map (fitter.results[terms.global], "_value_", '_value_ [terms.fit.MLE]'),
                            fitter.display_orders[fitter.terms.MG94 ]);


for (_value_; in; fitter.global_dnds) {
    fitter.omega_parameters = ((fitter.results[terms.global])[_value_[terms.description]])[terms.id];
    fitter.omega.CI = parameters.GetProfileCI(fitter.omega_parameters,fitter.results[terms.likelihood_function], 0.95);
    io.ReportProgressMessageMD ("fitter", fitter.terms.MG94, "* " + _value_[utility.getGlobalValue("terms.description")] + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4) + 
                    " (95% profile CI " + Format ((fitter.omega.CI )[terms.lower_bound],8,4) + "-" + Format ((fitter.omega.CI )[terms.upper_bound],8,4) + ")");
                    
    utility.EnsureKey(((fitter.json[terms.json.fits])[fitter.terms.MG94]), terms.fitter.ci);
    (((((fitter.json[terms.json.fits])[fitter.terms.MG94]))[ terms.fitter.ci]))[ _value_[utility.getGlobalValue("terms.description")]] = 
        {
            terms.lower_bound : (fitter.omega.CI )[terms.lower_bound],
            terms.upper_bound : (fitter.omega.CI )[terms.upper_bound]
        };
}

function fitter.correct_p_values () {
   fitter.corrected = (math.HolmBonferroniCorrection (
        utility.Map (fitter.lrt, "_value_", "_value_[terms.p_value]")
    ));

    fitter.fdr = (math.BenjaminiHochbergFDR (
        utility.Map (fitter.lrt, "_value_", "_value_[terms.p_value]")
    ));
    
    for (i,v; in; fitter.lrt) {
        v [terms.json.corrected_pvalue] = fitter.corrected[i];
        v ["FDR"] = fitter.fdr[i];
    }
        
    selection.io.json_store_branch_attribute(fitter.json, terms.fitter.lrt , terms.json.branch_attributes, fitter.display_orders[fitter.terms.MG94 ] + 4,
                                             "0",
                                             fitter.lrt);
}


if (fitter.model_type == terms.local) {

    fitter.table_screen_output  = {"0" : "Branch", "1" : "Length" , "2": "dN/dS", "3" : "Approximate dN/dS CI"};


    if (^"fitter.compute_lrt") {
        io.ReportProgressMessageMD("fitter", "LRT", "Running the likelihood ratio tests for dN/dS=1 and estimating confidence intervals for dN/dS along each branch");
        fitter.table_screen_output + "LRT p-value dN != dS";
	} else {
        io.ReportProgressMessageMD("fitter", "LRT", "Estimating confidence intervals for dN/dS along each branch");	
	}
	
	fitter.table_output_options = {terms.table_options.header : TRUE, 
                            terms.table_options.minimum_column_width: 16,
                            terms.table_options.column_widths: {
                                    "0" : 50,
                                    "1" : 10,
                                    "2" : 10,
                                    "3" : 20,
                                    "4" : 12}, 
                            terms.table_options.align : "center"};


    fprintf (stdout, "\n",
                    io.FormatTableRow (fitter.table_screen_output,fitter.table_output_options));
                    
    fitter.table_output_options [terms.table_options.header] = FALSE;            
    fitter.ci = {};
    fitter.lrt = {};
    fitter.save_lf = estimators.TakeLFStateSnapshot((^"fitter.results")[^"terms.likelihood_function"]);
    fitter.alphaP = "";
    fitter.betaP  = "";
       
    lfunction fitter.local_omega (set) {
        if (set) {
            parameters.SetConstraint (^"fitter.betaP", ^"fitter.alphaP", "");
        } else {
            parameters.ClearConstraint (^"fitter.betaP");
        }
        return 1;
    }
    
    lfunction profile_ci (tree_name, node_name, model_description, ignore) {
        omega = 1;
        omega :> 0;
        omega :< 10000;
        
        report_row = {};
        report_row + node_name;
        report_row + Format (((((^"fitter.results")[^"terms.branch_length"])["0"])[node_name])[^"terms.fit.MLE"], 0, 3);
          
        alphaName = tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")];
        betaName = tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")];
        saveAlpha = ^alphaName;
        saveBeta = ^betaName;
        
        ^"fitter.alphaP" = alphaName;
        ^"fitter.betaP" = betaName;
        
        omega = saveBeta/Max (saveAlpha, 1e-6);
        report_row + Format (parameters.NormalizeRatio (saveBeta, saveAlpha), 0, 3);
        
        ^betaName := omega * ^alphaName;
        ci_spec = {&omega : 1};
    
        (^"fitter.ci")[node_name] = parameters.GetProfileCI(ci_spec,(^"fitter.results")[^"terms.likelihood_function"], 0.95);
        
        report_row + (Format (((^"fitter.ci")[node_name])[^"terms.lower_bound"], 0, 3) + " - " + Format (((^"fitter.ci")[node_name])[^"terms.upper_bound"], 0, 3));
        
        ^alphaName = saveAlpha;
        ^betaName = saveBeta;  
        
        
        if (^"fitter.compute_lrt") {
            local.lrt = estimators.ConstrainAndRunLRT ((^"fitter.results")[^"terms.likelihood_function"], "fitter.local_omega");
            report_row + Format (local.lrt[^"terms.p_value"], 0, 4);
            (^"fitter.lrt")[node_name] = local.lrt;
        }  
        
        fprintf (stdout, io.FormatTableRow (report_row,^"fitter.table_output_options"));    
    }
    
    estimators.TraverseLocalParameters (fitter.results[terms.likelihood_function], fitter.results[utility.getGlobalValue("terms.model")], "profile_ci");
    
    
    
    selection.io.json_store_branch_attribute(fitter.json, terms.fitter.ci , terms.json.branch_attributes, fitter.display_orders[fitter.terms.MG94 ] + 3,
                                             "0",
                                             fitter.ci);
                        
    if (fitter.compute_lrt) {
        fitter.correct_p_values ();
    }
                                            
}

if (fitter.model_type == terms.fitter.lineage) {

    fitter.table_screen_output  = {"0" : "Lineage", "1" : "Root-to-tip" , "2": "dN/dS", "3" : "Approximate dN/dS CI"};

    if (^"fitter.compute_lrt") {
        io.ReportProgressMessageMD("fitter", "LRT", "Running the likelihood ratio tests for dN/dS=1 and estimating confidence intervals for dN/dS along each lineage");
        fitter.table_screen_output + "LRT p-value dN != dS";
	} else {
        io.ReportProgressMessageMD("fitter", "LRT", "Estimating confidence intervals for dN/dS along each lineage");	
	}
	
	fitter.table_output_options = {terms.table_options.header : TRUE, 
                            terms.table_options.minimum_column_width: 16,
                            terms.table_options.column_widths: {
                                    "0" : 50,
                                    "1" : 10,
                                    "2" : 10,
                                    "3" : 20,
                                    "4" : 12}, 
                            terms.table_options.align : "center"};
                            
     fprintf (stdout, "\n",
                    io.FormatTableRow (fitter.table_screen_output,fitter.table_output_options));
                    
    fitter.table_output_options [terms.table_options.header] = FALSE;            
    fitter.ci = {};
    fitter.lrt = {};
    fitter.MLE = {};
    fitter.lf_name = (^"fitter.results")[^"terms.likelihood_function"];
    fitter.save_lf = estimators.TakeLFStateSnapshot(fitter.lf_name);
    global fitter.lineage_omega = 1;
              
    
    GetString(fitter.lf_info, ^ fitter.lf_name, -1);
    fitter.tree_id = (fitter.lf_info[utility.getGlobalValue ("terms.fit.trees")])[0];
    fitter.parent_map = trees.ParentMap (fitter.tree_id);
    
    
    lfunction fitter.lineage_path (tree_name, node_name, model_description, ignore) {
          
        alphaName = tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")];
        betaName = tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")];
        if ((^"fitter.path_to_root")[node_name]) {
            saveAlpha = ^alphaName;
            saveBeta = ^betaName;
                        
            parameters.SetConstraint (betaName, "fitter.lineage_omega*" + alphaName, "");
            
        }
    }
    
    fitter.ci_spec = {"fitter.lineage_omega" : 1};
 
    
    function fitter.lrt_lineage_omega (set) {
        if (set) {
            fitter.SetToOne.stash = fitter.lineage_omega;
            parameters.SetConstraint ("fitter.lineage_omega", "1", "");

        } else {
            parameters.SetValue ("fitter.lineage_omega", fitter.SetToOne.stash);
        }
        return 1;
    }
    
    for (n; in; TipName(^fitter.tree_id, -1) ) {
        fitter.path_to_root = {};
        fitter.path_to_root [n] = 1;
        fitter.p = fitter.parent_map[n];
        fitter.ll = (((fitter.results[terms.branch_length])["0"])[n])[terms.fit.MLE];
        while (Abs (fitter.p)) {
            fitter.ll += (((fitter.results[terms.branch_length])["0"])[fitter.p])[terms.fit.MLE];
            fitter.path_to_root [fitter.p] = 1;
            fitter.p = fitter.parent_map[fitter.p];
        }
        fitter.path_mean = 1;
        estimators.TraverseLocalParameters (fitter.lf_name, fitter.results[utility.getGlobalValue("terms.model")], "fitter.lineage_path");
        Optimize (res, ^fitter.lf_name);
        
        

        report_row = {};
        report_row + n;
        report_row + Format (fitter.ll, 0, 3);
        report_row + Format (fitter.lineage_omega,0,3);
        fitter.MLE [n] = fitter.lineage_omega;
        
        fitter.ci[n] = parameters.GetProfileCI(fitter.ci_spec,fitter.lf_name, 0.95);

        report_row + (Format ((fitter.ci[n])[terms.lower_bound], 0, 3) + " - " + Format ((fitter.ci[n])[terms.upper_bound], 0, 3));
        
        if (fitter.compute_lrt) {
            local.lrt = estimators.ConstrainAndRunLRT (fitter.lf_name, "fitter.lrt_lineage_omega");
            report_row + Format (local.lrt[terms.p_value], 0, 4);
            fitter.lrt[n] = local.lrt;
        }         

        fprintf (stdout, io.FormatTableRow (report_row,^"fitter.table_output_options"));    

        estimators.RestoreLFStateFromSnapshot (fitter.lf_name, fitter.save_lf);
    }
    
    selection.io.json_store_branch_attribute(fitter.json, terms.fitter.ci , terms.json.branch_attributes, fitter.display_orders[fitter.terms.MG94 ] + 3,
                                         "0",
                                         fitter.ci);

    selection.io.json_store_branch_attribute(fitter.json, terms.fitter.MLE , terms.json.branch_attributes, fitter.display_orders[fitter.terms.MG94 ] + 4,
                                         "0",
                                         fitter.MLE);


    if (fitter.compute_lrt) {
        fitter.correct_p_values ();
    }

}


if (fitter.compute_lrt && fitter.model_type != terms.local && fitter.model_type != terms.fitter.lineage) {
    selection.io.startTimer (fitter.json [terms.json.timers], fitter.terms.LRT , fitter.display_orders [fitter.terms.LRT ]);
    io.ReportProgressMessageMD("fitter", "LRT", "Running the likelihood ratio tests for dN/dS=1");

    fitter.LRTs = {};
    
    function fitter.SetToOne (set) {
        if (set) {
            fitter.SetToOne.stash = Eval (fitter.parameter_id);
            parameters.SetConstraint (fitter.parameter_id, "1", "");
        } else {
            parameters.SetValue (fitter.parameter_id, fitter.SetToOne.stash);
        }
        return 1;
    }
    
    fitter.pvalues = {};

    for (parameter; in; fitter.global_dnds) {
        
        fitter.parameter_id = ((fitter.results[terms.global])[parameter[terms.description]])[terms.id];
        fitter.parameter_desc = parameter[terms.description];
        io.ReportProgressMessageMD("fitter", "LRT", "\n>Testing _`fitter.parameter_desc`_ == 1");
        fitter.LRTs     [fitter.parameter_desc] = (estimators.ConstrainAndRunLRT (fitter.results[terms.likelihood_function], "fitter.SetToOne"));
        io.ReportProgressMessageMD("fitter", "LRT", "\nLikelihood ratio test for _`fitter.parameter_desc` == 1_, **p = " + Format ((fitter.LRTs[fitter.parameter_desc])[terms.p_value], 8, 4) + "**.");
        fitter.pvalues  [fitter.parameter_desc] =  (fitter.LRTs     [fitter.parameter_desc])[terms.p_value];
    }
    
    
    fitter.json [terms.json.test_results] = fitter.LRTs;
    selection.io.stopTimer (fitter.json [terms.json.timers], fitter.terms.LRT);

}

utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, fitter.terms.MG94 , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ],
                                             _key_,
                                             selection.io.extract_branch_info((fitter.results[terms.branch_length])[_key_], "selection.io.branch.length"));');
                                             
fitter.ESEN_trees = estimators.FitMGREVExtractComponentBranchLengths (fitter.codon_data_info , fitter.results );

fitter.stree_info  = trees.ExtractTreeInfo ((fitter.ESEN_trees [terms.fit.synonymous_trees])[0]);
fitter.nstree_info = trees.ExtractTreeInfo ((fitter.ESEN_trees [terms.fit.nonsynonymous_trees])[0]);

fitter.codon_counts = genetic_code.ComputePairwiseDifferencesAndExpectedSites (fitter.codon_data_info[terms.code],{});

fitter.efv = (fitter.results[terms.efv_estimate])[utility.Keys(fitter.results[terms.efv_estimate])[0]];



fitter.S = +(fitter.efv $ fitter.codon_counts[terms.genetic_code.SS]); 
fitter.NS = +(fitter.efv $ fitter.codon_counts[terms.genetic_code.NS]); 


utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, terms.genetic_code.synonymous , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ] + 1,
                                             _key_,
                                             fitter.stree_info[terms.branch_length])');
                                   

utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, terms.genetic_code.nonsynonymous , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ] +2,
                                             _key_,
                                             fitter.nstree_info[terms.branch_length])');

fitter.dS =  fitter.stree_info;                                     
for (fitter.n, fitter.b; in; fitter.stree_info[terms.branch_length]) {
    (fitter.dS [terms.branch_length])[fitter.n] = fitter.b * (fitter.S+fitter.NS)/fitter.S;
}

fitter.dN =  fitter.nstree_info;                                     
for (fitter.n, fitter.b; in; fitter.nstree_info[terms.branch_length]) {
    (fitter.dN [terms.branch_length])[fitter.n] = fitter.b * (fitter.S+fitter.NS)/fitter.NS;
}


utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, terms.json.dS , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ] + 1,
                                             _key_,
                                             fitter.dS[terms.branch_length])');
                                             
utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, terms.json.dN, terms.branch_length, fitter.display_orders[fitter.terms.MG94 ] + 1,
                                             _key_,
                                             fitter.dN[terms.branch_length])');


io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 + terms.genetic_code.synonymous, "**Synonymous tree** \n" + (fitter.ESEN_trees [terms.fit.synonymous_trees])[0]);
io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 + terms.genetic_code.nonsynonymous, "**Non-synonymous tree** \n" + (fitter.ESEN_trees [terms.fit.nonsynonymous_trees])[0]);
io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 + terms.genetic_code.nonsynonymous, "**Combined tree** \n" + (fitter.ESEN_trees [terms.json.trees])[0]);

KeywordArgument ("save-fit", "Save MG94 model fit to this file (default is not to save)", "/dev/null");
io.SpoolLFToPath(fitter.results[terms.likelihood_function], io.PromptUserForFilePath ("Save MG94 model fit to this file ['/dev/null' to skip]"));

selection.io.stopTimer (fitter.json [terms.json.timers], fitter.terms.MG94);

selection.io.stopTimer (fitter.json [terms.json.timers], "Overall");
io.ReportProgressMessageMD ("fitter", "writing", "Writing detailed analysis report to \`" + fitter.codon_data_info [terms.json.json] + "\'");
io.SpoolJSON (fitter.json, fitter.codon_data_info [terms.json.json]);

return fitter.results;

