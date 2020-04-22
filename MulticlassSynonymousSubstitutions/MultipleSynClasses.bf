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

SetDialogPrompt ("A TSV file with three columns (AA, Codon, Class) which is used to partition synonymous substitutions into groups");
fitter.classes = io.ReadDelimitedFile (null, "\t", TRUE);

io.CheckAssertion("utility.Array1D(fitter.classes[terms.io.header])==3", "Expected a TSV file with 3 columns");

fitter.codons_by_class = {};
utility.ForEach (fitter.classes[terms.io.rows], "_record_",
'
    fitter.codons_by_class[_record_[1]] = _record_[2];
');

fitter.classes = utility.Values(fitter.codons_by_class);
fitter.class_count = utility.Array1D(fitter.classes);
io.CheckAssertion("fitter.class_count>=2", "Expected at least 2 codon classes");

fitter.choices = {fitter.class_count,2};
for (i = 0; i < fitter.class_count; i += 1) {
    fitter.choices[i][0] = fitter.classes[i];
    fitter.choices[i][1] = "Codon class " + fitter.classes[i];
}

fitter.neutral_reference = io.SelectAnOption  (fitter.choices, "Select the codon class which will serve as the neutral rate reference (relative rate = 1)");

fitter.model_type = io.SelectAnOption ({terms.global : "rates shared by all branches", terms.local : "separate rates for each branch"}, "Model Type");
fitter.frequency_type = io.SelectAnOption ({"CF3x4" : terms.frequencies.CF3x4,
                                            "F3x4" : terms.frequencies.F3x4,
                                            "F1x4" : terms.frequencies.F1x4}, "Equilibrium frequency estimator");

fitter.compute_ci = io.SelectAnOption ({"No"  : "Do not compute profile confidence intervals for substitution rates",
                                        "Yes" : "Compute profile confidence intervals for substitution rates"}, "Compute profile confidence intervals") != "No";

fitter.compute_lrt = io.SelectAnOption ({"No"  : "Do not perform LRT",
                                        "Yes" : "Perform LRT to compare synonymous rates and omega == 1"}, "Perform LRT to test rate equality and omega != 1") != "No";


KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')", fitter.codon_data_info [terms.json.json]);
fitter.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

namespace fitter {
    doGTR ("fitter");
}

io.ReportProgressMessageMD ("fitter", fitter.terms.MG94,  "Fitting `fitter.terms.MG94`");
selection.io.startTimer (fitter.json [terms.json.timers], fitter.terms.MG94 , fitter.display_orders [fitter.terms.MG94 ]);

fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, "fitter.defineMG", fitter.codon_data_info [utility.getGlobalValue("terms.code")],
    {
        terms.run_options.model_type: fitter.model_type,
        terms.run_options.retain_lf_object: TRUE,
        terms.run_options.retain_model_object : TRUE
    },
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

//----------------------------------------------------------------------------------------------------------------

lfunction fitter.defineMG (type, code) {
    m = Call ("models.codon.MG_REV.ModelDescription", type, code);
    if (^"fitter.frequency_type" == "F3x4") {
        m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F3x4";
    } else {
        if (^"fitter.frequency_type" == "F1x4") {
            m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F1x4";
        }
    }
    m[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which allows multiple classes of synonymous substitution rates";
    m[utility.getGlobalValue("terms.model.q_ij")] = "fitter.codon.MG_REV_MC._GenerateRate";
    return m;
}

lfunction models.codon.MG_REV_MH.ModelDescription(type, code) {

    // piggy-back on the standard MG_REV model for most of the code

    mg_base = models.codon.MG_REV.ModelDescription (type, code);
    mg_base[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which allows for two-hit substitutions";
    mg_base[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.MG_REV_MH._GenerateRate";

    return mg_base;
}


lfunction fitter.codon.MG_REV_MC._GenerateRate (fromChar, toChar, namespace, model_type, model) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    omega_term = utility.getGlobalValue ("terms.parameters.omega_ratio");
    alpha_term = utility.getGlobalValue ("terms.parameters.synonymous_rate");
    beta_term  = utility.getGlobalValue ("terms.parameters.nonsynonymous_rate");
    omega      = "omega";
    alpha      = "alpha";
    beta       = "beta";

    _tt = model[utility.getGlobalValue("terms.translation_table")];

    if (diff_count == 1) {

        _GenerateRate.p[model_type] = {};
        _GenerateRate.p[utility.getGlobalValue("terms.global")] = {};

        nuc_rate = "";

        for (i = 0; i < diff_count; i += 1) {
            if ((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] > (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")]) {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")] + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")];
            } else {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] +(_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")];
            }
            nuc_p = parameters.ApplyNameSpace(nuc_p, namespace);
            (_GenerateRate.p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")], (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")])] = nuc_p;

            nuc_rate = parameters.AppendMultiplicativeTerm (nuc_rate, nuc_p);
       }


        rate_entry = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {

            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (_GenerateRate.p[model_type])[omega_term] = aa_rate;
            } else {
                aa_rate = beta;
                (_GenerateRate.p[model_type])[beta_term] = aa_rate;
            }
            rate_entry += "*" + aa_rate;
        } else {

            class_from = (^"fitter.codons_by_class")[fromChar];
            class_to   = (^"fitter.codons_by_class")[toChar];

            if (class_from == class_to) {
                if (class_from == ^"fitter.neutral_reference") {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from;
                        (_GenerateRate.p[model_type])[alpha_term + " within codon class " + class_from] = codon_rate;
                        rate_entry += "*" + codon_rate;
                    } else {
                        rate_entry = nuc_rate;
                    }
                } else {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from;
                    } else {
                        codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from, namespace);
                    }
                    (_GenerateRate.p[model_type])[alpha_term + " within codon class " + class_from] = codon_rate;
                    rate_entry += "*" + codon_rate;
                }
            } else {
                if (class_from > class_to) {
                    class_from = (^"fitter.codons_by_class")[toChar];
                    class_to   = (^"fitter.codons_by_class")[fromChar];
                }
                if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from + "_" + class_to;
                } else {
                    codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from + "_" + class_to, namespace);
                }
                (_GenerateRate.p[model_type])[alpha_term + " between codon classes " + class_from + " and "  + class_to] = codon_rate;
                rate_entry += "*" + codon_rate;
            }
        }

        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
    }

    return _GenerateRate.p;
}

