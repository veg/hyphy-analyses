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


fitter.analysis_description = {terms.io.info : "Fit an MG94xREV model with several selectable options frequency estimator and report the fit results including dN/dS ratios, and synonymous and non-synonymous branch lengths",
                               terms.io.version : "0.1",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree"
                              };

io.DisplayAnalysisBanner (fitter.analysis_description);

  
namespace fitter.terms {
    MG94 = "Standard MG94";
}  
 
KeywordArgument ("code",        "Which genetic code should be used", "Universal");  
KeywordArgument ("alignment",   "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",        "A phylogenetic tree", null, "Please select a tree file for the data:");
KeywordArgument ("type",        "Model type: global (single dN/dS for all branches) or local (separate dN/dS)", terms.global, "Model Type");
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
                         fitter.terms.dN          : 3
                        };


selection.io.startTimer (fitter.json [terms.json.timers], "Overall", 0);


namespace fitter {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fitter", 
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}


fitter.model_type = io.SelectAnOption ({terms.global : "Shared dN/dS for all branches", terms.local : "Each branch has its own dN and dS"}, "Model Type");
fitter.frequency_type = io.SelectAnOption ({"CF3x4" : terms.frequencies.CF3x4, 
                                            "F3x4" : terms.frequencies.F3x4,
                                            "F1x4" : terms.frequencies.F1x4}, "Equilibrium frequency estimator");


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

fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, "fitter.defineMG", fitter.codon_data_info [utility.getGlobalValue("terms.code")],
    {
        terms.run_options.model_type: fitter.model_type,
        terms.run_options.retain_lf_object: TRUE,
        terms.run_options.retain_model_object : TRUE
    }, 
    fitter.gtr_results);


io.ReportProgressMessageMD("fitter", fitter.terms.MG94 , "* " + selection.io.report_fit (fitter.results, 0, fitter.codon_data_info[terms.data.sample_size]));
fitter.global_dnds = selection.io.extract_global_MLE_re (fitter.results, terms.parameters.omega_ratio);

utility.ForEach (fitter.global_dnds, "_value_", 'io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 , "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

selection.io.json_store_lf (fitter.json,
                            fitter.terms.MG94 ,
                            fitter.results[terms.fit.log_likelihood],
                            fitter.results[terms.parameters],
                            fitter.sample_size,
                            utility.Map (fitter.results[terms.global], "_value_", '_value_ [terms.fit.MLE]'),
                            fitter.display_orders[fitter.terms.MG94 ]);


utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, fitter.terms.MG94 , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ],
                                             _key_,
                                             selection.io.extract_branch_info((fitter.results[terms.branch_length])[_key_], "selection.io.branch.length"));');
                                             
fitter.ESEN_trees = estimators.FitMGREVExtractComponentBranchLengths (fitter.codon_data_info , fitter.results );

fitter.stree_info  = trees.ExtractTreeInfo ((fitter.ESEN_trees [terms.fit.synonymous_trees])[0]);
fitter.nstree_info = trees.ExtractTreeInfo ((fitter.ESEN_trees [terms.fit.nonsynonymous_trees])[0]);


utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, terms.genetic_code.synonymous , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ] + 1,
                                             _key_,
                                             fitter.stree_info[terms.branch_length])');

utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fitter.json, terms.genetic_code.nonsynonymous , terms.branch_length, fitter.display_orders[fitter.terms.MG94 ] +2,
                                             _key_,
                                             fitter.nstree_info[terms.branch_length])');


io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 + terms.genetic_code.synonymous, "**Synonymous tree** \n" + (fitter.ESEN_trees [terms.fit.synonymous_trees])[0]);
io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 + terms.genetic_code.nonsynonymous, "**Non-synonymous tree** \n" + (fitter.ESEN_trees [terms.fit.nonsynonymous_trees])[0]);

KeywordArgument ("save-fit", "Save MG94 model fit to this file (default is not to save)", "/dev/null");
io.SpoolLFToPath(fitter.results[terms.likelihood_function], io.PromptUserForFilePath ("Save MG94 model fit to this file ['/dev/null' to skip]"));

selection.io.stopTimer (fitter.json [terms.json.timers], fitter.terms.MG94);

selection.io.stopTimer (fitter.json [terms.json.timers], "Overall");
io.ReportProgressMessageMD ("fitter", "writing", "Writing detailed analysis report to \`" + fitter.codon_data_info [terms.json.json] + "\'");
io.SpoolJSON (fitter.json, fitter.codon_data_info [terms.json.json]);

return fitter.results;

