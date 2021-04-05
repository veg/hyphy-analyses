
RequireVersion("2.5.0");

LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4

// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/estimators.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");



utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);
utility.SetEnvVariable ("USE_MEMORY_SAVING_DATA_STRUCTURES", 1e8);
utility.SetEnvVariable ("LF_SMOOTHING_SCALER",0.05);


/*------------------------------------------------------------------------------*/

relax.analysis_description = {
                               terms.io.info : "RELAX (a random effects test of selection relaxation) uses a random effects branch-site model framework to test whether a set of 'Test' branches evolves under relaxed selection relative to a set of 'Reference' branches (R), as measured by the relaxation parameter (K).
                                                This version is customized to fit **only** the general exploratory model and then examine whether or not **k** is different from 1 for every individual branch",
                               terms.io.version : "3.1",
                               terms.io.reference : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                               terms.io.authors : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"
                              };

relax.json    = { terms.json.analysis: relax.analysis_description,
                  terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : {}
                  };

relax.relaxation_parameter        = "relax.K";
relax.rate_classes                = 3;

relax.initial_ranges              = {};
relax.initial_grid.N              = 500;

relax.MG94_name = terms.json.mg94xrev_sep_rates;
relax.general_descriptive_name = "General descriptive";


terms.relax.k          = "relaxation or intensification parameter";
terms.relax.k_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "50"
    };

terms.relax.k_range1    = {
        terms.lower_bound: "1",
        terms.upper_bound: "50"
    };

terms.relax.t_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "1e26"
    };


relax.p_threshold = 0.05;


relax.display_orders = {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        relax.MG94_name: 1,
                        relax.general_descriptive_name: 2
                       };

relax.table_screen_output  = {{"Branch name", "k MLE", "LRT", "Uncorrected p-value"}};
relax.table_screen_output2  = {{"Branch name", "k MLE", "LRT", "Corrected p-value"}};
relax.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "center",
                              terms.table_options.column_widths: {
            "0": 40,
            "1": 12,
            "2": 12,
            "3": 16
        }};
        
relax.report.header_done   = False;


/*------------------------------------------------------------------------------*/


KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
    /**
        keyword, description (for inline documentation and help messages), no default value,
        meaning that it will be required
    */

KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
    /** the use of null as the default argument means that the default expectation is for the
        argument to be missing, i.e. the tree is expected to be in the file
        the fourth, optional argument, can match this keyword with the dialog prompt / choice list title,
        meaning that it can only be consumed when this dialog prompt / choice list is invoked
        This allows handling some branching logic conditionals
    */


KeywordArgument ("rates", "The number omega rate classes to include in the model [2-10, default 3]", relax.rate_classes);

KeywordArgument ("fast-test", "Use a faster test (could have some false positives)", "Yes");


io.DisplayAnalysisBanner ( relax.analysis_description );

selection.io.startTimer (relax.json [terms.json.timers], "Overall", 0);

namespace relax {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "relax", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
    LoadFunctionLibrary ("SelectionAnalyses/modules/grid_compute.ibf");
}



relax.rate_classes = io.PromptUser ("The number omega rate classes to include in the model", relax.rate_classes, 2, 10, TRUE);

relax.fast_test = io.SelectAnOption (
    {"Yes" : "When scanning individual branches, fix global model parameters at their MLEs from the alternative; could create some false positives",
    "No" : "When scanning individual branches, re-estimate global model parameters at their MLEs; slower but more accurate"}, 
    "Use a faster test for branch scanning") == "Yes";

selection.io.startTimer (relax.json [terms.json.timers], "Preliminary model fitting", 1);

namespace relax {
    doGTR ("relax");
}

estimators.fixSubsetOfEstimates(relax.gtr_results, relax.gtr_results[terms.global]);

namespace relax {
    scaler_prefix = "relax.scaler";
    doPartitionedMG ("relax", FALSE);
}


io.ReportProgressMessageMD ("RELAX", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");


relax.final_partitioned_mg_results = estimators.FitMGREV (relax.filter_names, relax.trees, relax.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: relax.selected_branches,
}, relax.partitioned_mg_results);

io.ReportProgressMessageMD("RELAX", "codon-refit", "* " + selection.io.report_fit (relax.final_partitioned_mg_results, 0, relax.codon_data_info[terms.data.sample_size]));

relax.global_dnds  = selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
relax.report_dnds = {};

utility.ForEach (relax.global_dnds, "_value_", '
    io.ReportProgressMessageMD ("RELAX", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));
    relax.report_dnds [(regexp.FindSubexpressions (_value_[terms.description], "^" + terms.parameters.omega_ratio + ".+\\*(.+)\\*$"))[1]] = {"0" : {terms.json.omega_ratio : _value_[terms.fit.MLE], terms.json.proportion : 1}};
');




//Store MG94 to JSON
selection.io.json_store_lf_withEFV (relax.json,
                            relax.MG94_name,
                            relax.final_partitioned_mg_results[terms.fit.log_likelihood],
                            relax.final_partitioned_mg_results[terms.parameters],
                            relax.sample_size,
                            utility.ArrayToDict (utility.Map (relax.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (relax.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            relax.display_orders[relax.MG94_name]);
//single partition only for relax, but can't hurt .
utility.ForEachPair (relax.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(relax.json,relax.MG94_name, terms.branch_length, relax.display_orders[relax.MG94_name],
                                             _key_,
                                             selection.io.extract_branch_info((relax.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');




selection.io.stopTimer (relax.json [terms.json.timers], "Preliminary model fitting");
parameters.DeclareGlobalWithRanges (relax.relaxation_parameter, 1, 0, 50);

relax.ge_guess = None;

while (1) {

    relax.ge.bsrel_model =  model.generic.DefineMixtureModel("relax.BS_REL.ModelDescription",
            "relax.ge", {
                "0": parameters.Quote(terms.local),
                "1": relax.codon_data_info[terms.code],
                "2": parameters.Quote (relax.rate_classes) // the number of rate classes
            },
            relax.filter_names,
            None);

    for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
        parameters.SetRange (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range_almost_01);
    }
    parameters.SetRange (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);


    relax.model_object_map = { "relax.ge" :       relax.ge.bsrel_model };

    io.ReportProgressMessageMD ("RELAX", "gd", "Fitting the general descriptive (separate k per branch) model");
    selection.io.startTimer (relax.json [terms.json.timers], "General descriptive model fitting", 2);

    relax.distribution = models.codon.BS_REL.ExtractMixtureDistribution(relax.ge.bsrel_model);
    PARAMETER_GROUPING = {};
    PARAMETER_GROUPING + relax.distribution["rates"];
    PARAMETER_GROUPING + relax.distribution["weights"];

    if (Type (relax.ge_guess) != "Matrix") {
        // first time in
        relax.initial.test_mean    =

        math.Mean (
            utility.Map (selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+"), "_v_", "_v_[terms.fit.MLE]"));

        //console.log (relax.initial.test_mean);
        relax.init_grid_setup        (relax.distribution);
        relax.initial_grid         = estimators.LHC (relax.initial_ranges,relax.initial_grid.N);
        relax.initial_grid = utility.Map (relax.initial_grid, "_v_",
            'relax._renormalize (_v_, "relax.distribution", relax.initial.test_mean)'
        );
        relax.nm.precision = -0.00025*relax.final_partitioned_mg_results[terms.fit.log_likelihood];

        parameters.DeclareGlobalWithRanges ("relax.bl.scaler", 1, 0, 1000);

        //VERBOSITY_LEVEL = 10;

        relax.grid_search.results =  estimators.FitLF (relax.filter_names, relax.trees,{ "0" : {"DEFAULT" : "relax.ge"}},
                                    relax.final_partitioned_mg_results,
                                    relax.model_object_map,
                                    {
                                        "retain-lf-object": TRUE,
                                        terms.run_options.apply_user_constraints: "relax.init.k",
                                        terms.run_options.proportional_branch_length_scaler :
                                                                                {"0" : "relax.bl.scaler"},

                                        terms.run_options.optimization_settings :
                                            {
                                                "OPTIMIZATION_METHOD" : "nedler-mead",
                                                "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
                                                "OPTIMIZATION_PRECISION" : relax.nm.precision
                                            } ,

                                        terms.search_grid : relax.initial_grid
                                    }
        );


        relax.general_descriptive.fit =  estimators.FitLF (relax.filter_names,
                                    relax.trees,
                                    { "0" : {"DEFAULT" : "relax.ge"}},
                                    relax.grid_search.results,
                                    relax.model_object_map,
                                    {
                                        terms.run_options.apply_user_constraints: "relax.init.k",
                                        terms.run_options.retain_lf_object : TRUE

                                    });


  } else {
        parameters.SetStickBreakingDistribution (relax.distribution, relax.ge_guess);
        relax.general_descriptive.fit =  estimators.FitLF (relax.filter_names,
                                        relax.trees,
                                        { "0" : {"DEFAULT" : "relax.ge"}},
                                        relax.final_partitioned_mg_results,
                                        relax.model_object_map,
                                        {
                                            terms.run_options.apply_user_constraints: "relax.init.k",
                                            terms.run_options.retain_lf_object : TRUE

                                        });
   }






    estimators.TraverseLocalParameters (relax.general_descriptive.fit [terms.likelihood_function], relax.model_object_map, "relax.set.k2");

    relax.ge.likelihood_function = relax.general_descriptive.fit [terms.likelihood_function];
    relax.general_descriptive.fit = estimators.FitExistingLF (relax.ge.likelihood_function, relax.model_object_map);

    selection.io.stopTimer (relax.json [terms.json.timers], "General descriptive model fitting");

    io.ReportProgressMessageMD("RELAX", "ge", "* " + selection.io.report_fit (relax.general_descriptive.fit, 9, relax.codon_data_info[terms.data.sample_size]));
    io.ReportProgressMessageMD("RELAX", "ge", "* The following baseline rate distribution for branch-site combinations was inferred");
    relax.inferred_ge_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistributionFromFit (relax.ge.bsrel_model, relax.general_descriptive.fit)) % 0;
    selection.io.report_dnds (relax.inferred_ge_distribution);

 

    if (relax.rate_classes > 2) {
        relax.same_rate = -1;
        for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
            if (Abs(relax.inferred_ge_distribution[relax.i][0] - relax.inferred_ge_distribution[relax.i-1][0]) < 1e-5) {
                relax.same_rate = relax.i;
                break;
            }
        }
       
        relax.cutoff = 0.1 / relax.codon_data.sites;
        if (relax.same_rate >= 0 || Min (relax.inferred_ge_distribution[-1][1],0) < relax.cutoff) {
            io.ReportProgressMessageMD("RELAX", "ge", "\n##### Because some of the rate classes were collapsed to 0, the model is likely overparameterized. RELAX will reduce the number of site rate classes by one and repeat the fit now.\n----\n");
            relax.rate_classes = relax.rate_classes - 1;
            relax.ge_guess = {relax.rate_classes, 2};
            relax.shift    = 0;
            for (relax.i = 0; relax.i < relax.rate_classes; relax.i += 1) {
                if ((relax.i == relax.same_rate || relax.inferred_ge_distribution[relax.i][1] < relax.cutoff) && relax.shift == 0) {
                    relax.shift += 1;
                    continue;
                }
                relax.ge_guess[relax.i][0] = relax.inferred_ge_distribution[relax.i + relax.shift][0];
                relax.ge_guess[relax.i][1] = relax.inferred_ge_distribution[relax.i + relax.shift][1];
            }
            //console.log (relax.ge_guess);
            continue;
        }
    }


    relax.distribution_for_json = {'Shared' : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : relax.inferred_ge_distribution [_index_][0],
                                                           terms.json.proportion  : relax.inferred_ge_distribution [_index_][1]}")
                                   };
    selection.io.json_store_lf (relax.json,
                                relax.general_descriptive_name,
                                relax.general_descriptive.fit[terms.fit.log_likelihood],
                                relax.general_descriptive.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                                relax.codon_data_info[terms.data.sample_size],
                                relax.distribution_for_json,
                                relax.display_orders[relax.general_descriptive_name]
                            );

    selection.io.json_store_branch_attribute(relax.json, relax.general_descriptive_name, terms.branch_length, relax.display_orders[relax.general_descriptive_name],
                                                 0,
                                                 selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "selection.io.branch.length"));

    relax.k_estimates = selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "relax.extract.k");
    relax.k_variables = selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "relax.extract.k.name");

    relax.k_stats = math.GatherDescriptiveStats (utility.Map (utility.UniqueValues (relax.k_estimates), "_value_", "0+_value_"));

    io.ReportProgressMessageMD("RELAX", "ge", "* Branch-level `terms.relax.k` distribution has mean " + Format (relax.k_stats[terms.math.mean], 5,2) + ", median " +
                                                 Format (relax.k_stats[terms.math.median], 5,2) + ", and 95% of the weight in " + Format (relax.k_stats[terms.math._2.5], 5,2) + " - " + Format (relax.k_stats[terms.math._97.5], 5,2));


    selection.io.json_store_branch_attribute(relax.json, "k (general descriptive)", terms.json.branch_label, relax.display_orders[relax.general_descriptive_name],
                                                 0,
                                                 relax.k_estimates);

    break;
}

selection.io.startTimer (relax.json [terms.json.timers], "Branch scan", 3);

relax.branch_tests = {};
relax.branch.to.constrain = "";
relax.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{relax.ge.likelihood_function}},
                               terms.mpi.Models : {{"relax.ge.bsrel_model"}},
                               terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
                               terms.mpi.Variables : {{"terms.relax.k","relax.k_estimates","relax.model_object_map","relax.general_descriptive.fit","relax.branch.to.constrain","relax.fast_test"}},
                               terms.mpi.Functions : {{"relax.set.k_is_one"}}
                             });

utility.ForEachPair (relax.k_estimates, "_branch_", "_estimate_",
        '
            mpi.QueueJob (relax.queue , "relax.handle_a_branch", {
                                                                "0" : relax.ge.likelihood_function,
                                                                "1" : _branch_,
                                                                "2" : _estimate_
                                                             },
                                      "relax.store_results");

        '
    );


mpi.QueueComplete (relax.queue);

relax.corrected_p = math.HolmBonferroniCorrection(utility.Map(relax.branch_tests,"_v_","_v_['p-value']"));

utility.ForEachPair (relax.corrected_p, "_key_", "_value_", "(relax.branch_tests[_key_])['Corrected p-value']=_value_");

(relax.table_output_options)[utility.getGlobalValue("terms.table_options.header")] = TRUE;

console.log ("");

io.ReportProgressMessageMD("RELAX", "Results", "Tests of individual branches for relaxation/intensification of selection");
fprintf (stdout, "\n",
    io.FormatTableRow (relax.table_screen_output2,relax.table_output_options));
        
(relax.table_output_options)[utility.getGlobalValue("terms.table_options.header")] = FALSE;

relax.sig_count = 0;

utility.ForEachPair (relax.branch_tests, "_key_", "_record_", 
"
    relax.padder = '';
    if (_record_['Corrected p-value'] <= relax.p_threshold) {
        relax.sig_count += 1;
        relax.padder = ' (*)';
    }
    result = {{_key_, Format (_record_[terms.fit.MLE],8,2), Format (_record_[terms.LRT], 8,2), Format (_record_['Corrected p-value'], 8,2) + relax.padder}};
    fprintf (stdout,io.FormatTableRow (result, relax.table_output_options));
");

selection.io.stopTimer (relax.json [terms.json.timers], "Branch scan");

console.log ("\n\n## RELAX scan result\n**" + relax.sig_count + "** branches had significant relaxation/intensification of selection at significance level of " + relax.p_threshold);

selection.io.stopTimer (relax.json [terms.json.timers], "Overall");

relax.json [terms.json.test_results] = relax.branch_tests;

io.SpoolJSON (relax.json, relax.codon_data_info [terms.json.json]);

return relax.json;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

lfunction relax.extract.k(branch_info) {
    return (branch_info[utility.getGlobalValue("terms.relax.k")])[utility.getGlobalValue("terms.fit.MLE")];
}

lfunction relax.extract.k.name (branch_info) {
    return (branch_info[utility.getGlobalValue("terms.relax.k")])[utility.getGlobalValue("terms.id")];
}

//------------------------------------------------------------------------------

lfunction relax.set.k (tree_name, node_name, model_description, ignore) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        t = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")];
        parameters.SetConstraint (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("relax.relaxation_parameter"), "");
        parameters.SetRange (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("terms.relax.k_range"));
        parameters.SetRange (tree_name + "." + node_name + "." + t, utility.getGlobalValue ("terms.relax.t_range"));
    }
    return tree_name + "." + node_name + "." + k;
}

//------------------------------------------------------------------------------

lfunction relax.set.k2 (tree_name, node_name, model_description, ignore) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        parameters.RemoveConstraint (tree_name + "." + node_name + "." + k);
    }
    return tree_name + "." + node_name + "." + k;
}

//------------------------------------------------------------------------------

lfunction relax.set.k_is_one (tree_name, node_name, model_description, ignore) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        if (^'relax.branch.to.constrain' == node_name) {
            parameters.SetConstraint (tree_name + "." + node_name + "." + k, "1", "");
        }
    }
    return tree_name + "." + node_name + "." + k;
}

//----------------------------------------------------------------------------------------
lfunction relax.handle_a_branch (lf, branch, estimate) {
    ^'relax.branch.to.constrain' = branch;
    save = estimators.TakeLFStateSnapshot(lf);
    if (^"relax.fast_test") { 
        parameters.FixParameterSet (utility.Map ((^'relax.general_descriptive.fit')[^'terms.global'], "_value_", "_value_[terms.id]"));
    }
    estimators.TraverseLocalParameters (lf, ^'relax.model_object_map', "relax.set.k_is_one");
    relax.general_descriptive.null = estimators.FitExistingLF (lf, ^'relax.model_object_map');
    estimators.RestoreLFStateFromSnapshot (lf, save);
    return relax.general_descriptive.null;
}

//----------------------------------------------------------------------------------------

lfunction relax.store_results (node, result, arguments) {

    if ( ^'relax.report.header_done' == FALSE) {
        io.ReportProgressMessageMD("RELAX", "Branch", "Testing individual branches for relaxation/intensification of selection");
        fprintf (stdout, "\n",
            io.FormatTableRow (^'relax.table_screen_output',^'relax.table_output_options'));
        ^'relax.report.header_done' = TRUE;
        (^'relax.table_output_options')[utility.getGlobalValue("terms.table_options.header")] = FALSE;
    }

    branch = arguments [1];
    estimate = arguments[2];
    (^'relax.branch_tests')[branch] =
        math.DoLRT (result[^'terms.fit.log_likelihood'], (^'relax.general_descriptive.fit')[^'terms.fit.log_likelihood'],1);
    ((^'relax.branch_tests')[branch])[^"terms.fit.MLE"] = estimate;
    result = {{branch, Format (estimate,8,2), Format (((^'relax.branch_tests')[branch])[^'terms.LRT'], 8,2), Format (((^'relax.branch_tests')[branch])[^'terms.p_value'], 8,4)}};
    fprintf (stdout,
                io.FormatTableRow (result,  ^'relax.table_output_options'));
}

//------------------------------------------------------------------------------

lfunction relax.init.k (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {

    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, "relax.set.k");
    rc = utility.getGlobalValue ("relax.rate_classes");
    /*if (rc > 2) {
        parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("relax.ge.bsrel_model") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),rc-1)), utility.getGlobalValue("terms.parameters.one"), utility.getGlobalValue("terms.global"));
    }*/
    /*parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("relax.ge.bsrel_model") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),utility.getGlobalValue ("relax.rate_classes"))),
                             "1/(" +
                                Join ("*", utility.Map (
                                    utility.Range (utility.getGlobalValue ("relax.rate_classes") - 1, 1, 1),
                                    "_value_",
                                    'model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,_value_))'
                                    ))
                             + ")",
                            "global");*/

    return 0;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL.ModelDescription (type, code, components) {
    model = models.codon.BS_REL.ModelDescription(utility.getGlobalValue ('terms.global'), code, components);
    model [utility.getGlobalValue("terms.model.defineQ")] = "relax.BS_REL._DefineQ";
    return model;
}

//------------------------------------------------------------------------------


lfunction relax.DistributionGuess (mean) {
    rc = utility.getGlobalValue ("relax.rate_classes");

    guess = {rc,2};

    guess[rc-1][0] = 5;
    guess[rc-1][1] = 0.1;

    for (k = 0; k < rc - 1; k += 1) {
        guess[k][0] = 0.1 ^ (1 / (1 + k));
        guess[k][1] = (0.9) / (rc-1) ;
    }

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {


    p = {};
    diff = models.codon.diff(fromChar, toChar);

    if (None != diff) {
        p[model_type] = {};
        p[utility.getGlobalValue("terms.global")] = {};

        if (diff[utility.getGlobalValue("terms.diff.from")] > diff[utility.getGlobalValue("terms.diff.to")]) {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.to")] + diff[utility.getGlobalValue("terms.diff.from")];
        } else {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.from")] + diff[utility.getGlobalValue("terms.diff.to")];
        }
        nuc_rate = parameters.ApplyNameSpace(nuc_rate, namespace);
        (p[utility.getGlobalValue("terms.global")])[terms.nucleotideRate(diff[utility.getGlobalValue("terms.diff.from")], diff[utility.getGlobalValue("terms.diff.to")])] = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (p[model_type])[omega_term] = aa_rate;
                utility.EnsureKey (p, utility.getGlobalValue("terms.local"));
                 (p[utility.getGlobalValue("terms.local")])[utility.getGlobalValue ("terms.relax.k")] = "k";
                 aa_rate += "^k";
            } else {
                aa_rate = beta;
                (p[model_type])[beta_term] = aa_rate;
            }
            p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + aa_rate;
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                (p[model_type])[alpha_term] = alpha;
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + alpha;
            } else {
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }
    }


    return p;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._DefineQ (bs_rel, namespace) {
    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};

    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
           return relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       } else {

       }
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}



//------------------------------------------------------------------------------

lfunction relax.grid.MatrixToDict (grid) {
    return utility.Map (utility.MatrixToListOfRows (grid), "_value_",
                                                                '{  terms.relax.k : {
                                                                            terms.id : relax.relaxation_parameter,
                                                                            terms.fit.MLE : _value_[1]
                                                                        }

                                                                 }');
}

//------------------------------------------------------------------------------

function relax.init_grid_setup (omega_distro) {
    utility.ForEachPair (omega_distro[terms.parameters.rates], "_index_", "_name_",
        '
            if (_index_[0] < relax.rate_classes - 1) { // not the last rate
                  relax.initial_ranges [_name_] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
            }  else {
                relax.initial_ranges [_name_] = {
                    terms.lower_bound : 1,
                    terms.upper_bound : 10
                };
            }
        '
    );


    utility.ForEachPair (omega_distro[terms.parameters.weights], "_index_", "_name_",
        '
             relax.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : 1
            };
        '
    );

}

//------------------------------------------------------------------------------

lfunction relax._renormalize (v, distro, mean) {

    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    d = Rows (m);
    m = +(m[-1][0] $ m[-1][1]); // current mean
    for (i = 0; i < d; i+=1) {
        (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] = (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] / m * mean;
    }
    return v;

}
