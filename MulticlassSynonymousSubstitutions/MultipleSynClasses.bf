RequireVersion ("2.5.28");


LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("lib/mss.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/models/codon/MSS.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


fitter.analysis_description = {terms.io.info : "Fit an MG94xREV model where synonymous substitutions are partitioned into several classes and within- and between-class rates are estimated. There are several selectable options for the frequency 
estimators and report the fit results including dN/dS ratios, and synonymous and non-synonymous branch lengths. v0.2 adds the ability to compute confidence intervals. v0.3 adds the ability to perform LRTs.",
                               terms.io.version : "0.3",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree, and a TSV file with codon class partitioning"
                              };

io.DisplayAnalysisBanner (fitter.analysis_description);


namespace fitter.terms {
    MG94 = "MG94 with multiple classes of synonymous substitutions";
    LRT = "LRT testing;"
}

KeywordArgument ("code",        "Which genetic code should be used", "Universal");
KeywordArgument ("alignment",   "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",        "A phylogenetic tree", null, "Please select a tree file for the data:");
KeywordArgument ("classes",     "A TSV file with three columns (AA, Codon, Class) which is used to partition synonymous substitutions into groups");
KeywordArgument ("neutral",     "Neutral reference class");
KeywordArgument ("type",        "Model type: global (single dN/dS for all branches) or local (separate dN/dS)", terms.global, "Model Type");
KeywordArgument ("frequencies", "Equilibrium frequency estimator", "CF3x4");
KeywordArgument ("ci",          "Compute profile confidence intervals", "No");
KeywordArgument ("lrt",         "Perform LRT to test which rates are different from the neutral rate", "No");
KeywordArgument ("F81",         "Use the F81 nucleotide component to look for the effect of nucleotide biases on rate estimates", "No");

fitter.json    = {
                    terms.json.analysis: fitter.analysis_description,
                    terms.json.input: {},
                    terms.json.fits : {},
                    terms.json.timers : {},
                  };

fitter.display_orders = {terms.original_name      :  -1,
                         terms.json.nucleotide_gtr: 0,
                         fitter.terms.MG94        : 1,
                         fitter.terms.dS          : 2,
                         fitter.terms.dN          : 3,
                         fitter.terms.LRT: 2,
                        };
                        
terms.fitter.ci = "Confidence Intervals";


selection.io.startTimer (fitter.json [terms.json.timers], "Overall", 0);


namespace fitter {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fitter",
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}

fitter.codons_by_class = models.codon.MSS.LoadClasses (null);
fitter.model_type = io.SelectAnOption ({terms.global : "rates shared by all branches", terms.local : "separate rates for each branch"}, "Model Type");
fitter.frequency_type = io.SelectAnOption ({"CF3x4" : terms.frequencies.CF3x4,
                                            "F3x4" : terms.frequencies.F3x4,
                                            "F1x4" : terms.frequencies.F1x4}, "Equilibrium frequency estimator");

fitter.compute_ci = io.SelectAnOption ({"No"  : "Do not compute profile confidence intervals for substitution rates",
                                        "Yes" : "Compute profile confidence intervals for substitution rates"}, "Compute profile confidence intervals") != "No";

fitter.compute_lrt = io.SelectAnOption ({"No"  : "Do not perform LRT",
                                        "Yes" : "Perform LRT to compare synonymous rates and omega == 1"}, "Perform LRT to test rate equality and omega != 1") != "No";

fitter.useF81 = io.SelectAnOption ({"No"  : "Use GTR (all nucleotide rates are different)",
                                        "Yes" : "Use F81 (all nucleotide rates are the same)"}, "Use the F81 nucleotide component to look for the effect of nucleotide biases on rate estimates") != "No";



KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')", fitter.codon_data_info [terms.json.json]);
fitter.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

namespace fitter {
    doGTR ("fitter");
}

if (fitter.useF81) {
    estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results[terms.global]);
}

io.ReportProgressMessageMD ("fitter", fitter.terms.MG94,  "Fitting `fitter.terms.MG94`");
selection.io.startTimer (fitter.json [terms.json.timers], fitter.terms.MG94 , fitter.display_orders [fitter.terms.MG94 ]);

function fitter.defineMG (type,code) {
    return models.codon.MSS.ModelDescription (type,code,fitter.codons_by_class);
}

fitter.opt.options = {
        terms.run_options.model_type: fitter.model_type,
        terms.run_options.retain_lf_object: TRUE,
        terms.run_options.retain_model_object : TRUE
        
    };

if (fitter.useF81) {
    fitter.opt.options [terms.run_options.apply_user_constraints] = "fitter.fix_to_F81";
}

fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, "fitter.defineMG", fitter.codon_data_info [utility.getGlobalValue("terms.code")],
    fitter.opt.options,
    fitter.gtr_results);


io.ReportProgressMessageMD("fitter", fitter.terms.MG94 , "* " + selection.io.report_fit (fitter.results, 0, fitter.codon_data_info[terms.data.sample_size]));
fitter.global_dnds = selection.io.extract_global_MLE_re (fitter.results, terms.parameters.omega_ratio + "|" + terms.parameters.synonymous_rate);


selection.io.json_store_lf (fitter.json,
                            fitter.terms.MG94 ,
                            fitter.results[terms.fit.log_likelihood],
                            fitter.results[terms.parameters],
                            fitter.sample_size,
                            utility.Map (fitter.results[terms.global], "_value_", '_value_ [terms.fit.MLE]'),
                            fitter.display_orders[fitter.terms.MG94 ]);


if (fitter.compute_ci) {
    utility.ForEach (fitter.global_dnds, "_value_", 
    '
        
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
    ');

} else {
    utility.ForEach (fitter.global_dnds, "_value_", 'io.ReportProgressMessageMD ("fitter", fitter.terms.MG94 , "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');
}




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

if (fitter.compute_lrt) {
    selection.io.startTimer (fitter.json [terms.json.timers], fitter.terms.LRT , fitter.display_orders [fitter.terms.LRT ]);

    io.ReportProgressMessageMD("fitter", "LRT", "Running likelihood ratio tests to compare all rates to the neutral rate (=1)");

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
        io.ReportProgressMessageMD("fitter", "LRT", "\nLikelihood ratio test for _`fitter.parameter_desc` == 1_, uncorrected **p = " + Format ((fitter.LRTs[fitter.parameter_desc])[terms.p_value], 8, 4) + "**.");
        fitter.pvalues  [fitter.parameter_desc] =  (fitter.LRTs     [fitter.parameter_desc])[terms.p_value];
    }
    
    fitter.corrected = math.HolmBonferroniCorrection (fitter.pvalues);
    
    io.ReportProgressMessageMD("fitter", "LRT", "#### Holm-Bonferroni Corrected p-values");
    
    for (parameter, info; in; fitter.LRTs) {
        (fitter.LRTs[parameter])[terms.json.corrected_pvalue] = fitter.corrected[parameter];
        (fitter.LRTs[parameter])[terms.json.uncorrected_pvalue] = (fitter.LRTs[parameter])[terms.p_value];
        (fitter.LRTs[parameter]) - terms.p_value;
        io.ReportProgressMessageMD("fitter", "LRT", "- Likelihood ratio test for _`parameter`_ == 1, corrected **p = " + Format ((fitter.LRTs[parameter])[terms.json.corrected_pvalue], 8, 4) + "**.");
    }
    
    fitter.json [terms.json.test_results] = fitter.LRTs;
    selection.io.stopTimer (fitter.json [terms.json.timers], fitter.terms.LRT);

}

selection.io.stopTimer (fitter.json [terms.json.timers], "Overall");
io.ReportProgressMessageMD ("fitter", "writing", "Writing detailed analysis report to \`" + fitter.codon_data_info [terms.json.json] + "\'");
io.SpoolJSON (fitter.json, fitter.codon_data_info [terms.json.json]);

return fitter.results;

lfunction fitter.fix_to_F81 (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    c = 0;
    re = terms.nucleotideRate("[ACGT]","[ACGT]");
    for (m; in; model_objects) {
        for (d, p; in; (m[^"terms.parameters"])[^"terms.global"]) {
            if (regexp.Find (d, re)) {
                parameters.SetConstraint (p, "1", "global");
                c+=1;
            }
        }
    }
    return c;
}


