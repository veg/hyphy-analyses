RequireVersion ("2.4.0");


LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");

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

partitioned.rate_classes = 2;


KeywordArgument ("code",        "Which genetic code should be used ", "Universal");
KeywordArgument ("alignment",   "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("branches",   "The set of branches to test", "All");


partitioned.json = {};

namespace partitioned {
    MG94 = "MG94";
    BSREL = "BSREL";
};

partitioned.json    = { terms.json.analysis: partitioned.analysis_description,
                       terms.json.input: {},
                       terms.json.fits : {},
                       terms.json.timers : {},
                  };

selection.io.startTimer (partitioned.json [terms.json.timers], "Overall", 0);

partitioned.display_orders = {terms.original_name: -1,
                         terms.json.nucleotide_gtr: 0,
                         partitioned.MG94: 1,
                         partitioned.BSREL: 2
                        };

namespace partitioned {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ("partitioned");
}

io.CheckAssertion ("partitioned.partition_count>=2","This analysis requires data with at least two partitions");

KeywordArgument ("rates", "The number omega rate classes to include in the model [1-10, default 2]", partitioned.rate_classes);
partitioned.rate_classes = io.PromptUser ("The number omega rate classes to include in the model", partitioned.rate_classes, 1, 10, TRUE);

KeywordArgument ("branch-lengths", "Proportional branch lengths across partitions", "Proportional");
partitioned.proportional_lengths = io.SelectAnOption ({"Proportional" : "Branch lengths are proportional across partitions", 
                                    "Unlinked"  : "Branch lengths are unlinked across partitions"},
                                    "Branch lengths across partitions"
                                    ) == "Proportional";
                                    
                                    
KeywordArgument ("bs-rel", "Run BS-REL analyses", "Yes");
partitioned.run_bs_rel = io.SelectAnOption ({"Yes" : "Run BS-REL analyses", 
                                    "No"  : "Skip BS-REL analyses"},
                                    "Run BS-REL (rate distribution)  analyses"
                                    ) == "Yes";

KeywordArgument ("mg", "Run MG94 analyses", "Yes");
partitioned.run_mg_rev = io.SelectAnOption ({"Yes" : "Run MG94 analyses", 
                                    "No"  : "Skip MG94 analyses"},
                                    "Run MG (mean dN/dS) analyses"
                                    ) == "Yes";
                                    

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'partitioned.json')", partitioned.codon_data_info [terms.json.json]);
partitioned.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");


selection.io.startTimer (partitioned.json [terms.json.timers], "Preliminary model fitting", 1);

namespace partitioned {
    doGTR ("partitioned");
}

selection.io.stopTimer (partitioned.json [terms.json.timers], "Preliminary model fitting");
partitioned.p_values = {};
// stick raw p-values here for correcting later


if (partitioned.run_mg_rev) {
    /*********************************************************************

    Build a model where each partition and each branch set receive their
    own dN/dS values. Start with the local model, and then define the constraints

    Branch lengths are proportional between the partitions

    *********************************************************************/

    partitioned.model_generator = "models.codon.MG_REV.ModelDescription";
    partitioned.MG_models = {};
    partitioned.model_map = {};

    partitioned.omega_ratios = {};
    // maps {"partition", "branch set"} pair to the name of the omega parameter

    partitioned.model_id = "partitioned.MG94";
    partitioned.MG_models[partitioned.model_id] =  model.generic.DefineModel(partitioned.model_generator, partitioned.model_id , {
                "0": "terms.local",
                "1" : partitioned.codon_data_info [utility.getGlobalValue("terms.code")]
            }, partitioned.filter_names, None);


    utility.ForEachPair (partitioned.filter_names, "_index_", "_value_",
    '
         partitioned.model_map [_index_] = {"DEFAULT" : partitioned.model_id};
    ');

    io.ReportProgressMessageMD ("partitioned", "MG94", "Fitting MG94xREV with separate dN/dS ratios for each partition");

    partitioned.mg_fit = estimators.FitLF(
                                            partitioned.filter_names,
                                            partitioned.trees,
                                            partitioned.model_map,
                                            partitioned.gtr_results,
                                            partitioned.MG_models,
                                            {
                                                terms.run_options.apply_user_constraints : "partitioned.define_omega_constraints",
                                                terms.run_options.retain_lf_object : TRUE
                                            }
                                          );
                                          
    partitioned.bs_rel_init = partitioned.mg_fit;

                                      
    /*fprintf ("/Users/sergei/Desktop/fit.dump", CLEAR_FILE, partitioned.mg_fit);
    Export (lfe, ^(partitioned.mg_fit[terms.likelihood_function]));
    fprintf ("/Users/sergei/Desktop/lf.dump", CLEAR_FILE, lfe);

    fscanf ("/Users/sergei/Desktop/fit.dump", "Raw", partitioned.mg_fit); partitioned.mg_fit = Eval (partitioned.mg_fit);
    ExecuteAFile ("/Users/sergei/Desktop/lf.dump");*/

    io.ReportProgressMessageMD("partitioned", "MG94", "* " + selection.io.report_fit (partitioned.mg_fit, 0, partitioned.codon_data_info[terms.data.sample_size]));

    partitioned.global_dnds = selection.io.extract_global_MLE_re (partitioned.mg_fit, "^" + utility.getGlobalValue("terms.parameters.omega_ratio"));

    partitioned.omega.CI = {};
    partitioned.omega_parameters = {};

    utility.ForEach (partitioned.global_dnds, "_value_", '
        if (null != regexp.Find (_value_[utility.getGlobalValue("terms.description")],"\|" + terms.tree_attributes.test)) {
           partitioned.omega_parameters + ((partitioned.mg_fit[terms.global])[_value_[terms.description]])[terms.id];
           partitioned.omega.CI [_value_[terms.description]] = parameters.GetProfileCI(partitioned.omega_parameters[Abs(partitioned.omega_parameters)-1],partitioned.mg_fit[terms.likelihood_function], 0.95);
            io.ReportProgressMessageMD ("partitioned", "MG94", "* " + _value_[utility.getGlobalValue("terms.description")] + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4) + 
                    " (95% profile CI " + Format ((partitioned.omega.CI [_value_[terms.description]])[terms.lower_bound],8,4) + "-" + Format ((partitioned.omega.CI [_value_[terms.description]])[terms.upper_bound],8,4) + ")");
       
        } else {
            io.ReportProgressMessageMD ("partitioned", "MG94", "* " + _value_[utility.getGlobalValue("terms.description")] + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4));
        }
    ');

    io.ReportProgressMessageMD ("partitioned", "MG94-null", "Checking for difference in partition-level dN/dS using a likelihood ratio test");
    parameters.ConstrainParameterSet (partitioned.omega_parameters,null);
    partitioned.mg_fit.null = estimators.FitExistingLF (partitioned.mg_fit[terms.likelihood_function], partitioned.MG_models);

    io.ReportProgressMessageMD("partitioned", "MG94-null", "* " + selection.io.report_fit (partitioned.mg_fit.null  , 0, partitioned.codon_data_info[terms.data.sample_size]));
    io.ReportProgressMessageMD("partitioned", "MG94-null", "* Test for dN/dS equality across partitions:  " + selection.io.report_lrt (partitioned.ComputeLRT (partitioned.mg_fit, partitioned.mg_fit.null, "Test for dN/dS equality across partitions")));
} else {
    partitioned.bs_rel_init = partitioned.gtr_results;
}


if (partitioned.run_bs_rel) {
    partitioned.has_background = FALSE;

    utility.ForEachPair (partitioned.selected_branches, "_partition_", "_selection_",
        "_selection_ = utility.Filter (_selection_, '_value_', '_value_ != terms.tree_attributes.test');
         if (utility.Array1D (_selection_)) { partitioned.has_background = TRUE;} ");


    

    partitioned.BSREL.model_generator = "models.codon.BS_REL.ModelDescription";
    partitioned.BSREL_models = {};
    partitioned.BSREL_model_map = {};
    partitioned.BSREL_distributions = {};
    partitioned.BSREL_bg_distributions = {};


    for (partitioned.k = 0; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
        partitioned.model_id = partitioned.bsrel_model_id (partitioned.k);
        partitioned.model_id_bg = partitioned.bsrel_model_id_bg (partitioned.k);
        partitioned.BSREL_models[partitioned.model_id] =  model.generic.DefineMixtureModel(partitioned.BSREL.model_generator,
            partitioned.model_id, {
                "0": parameters.Quote(terms.global),
                "1": partitioned.codon_data_info[terms.code],
                "2": parameters.Quote (partitioned.rate_classes) // the number of rate classes
            },
            partitioned.filter_names,
            None);
        
        if (partitioned.has_background) {
            partitioned.BSREL_models[partitioned.model_id_bg] =  model.generic.DefineMixtureModel(partitioned.BSREL.model_generator,
                partitioned.model_id_bg, {
                    "0": parameters.Quote(terms.global),
                    "1": partitioned.codon_data_info[terms.code],
                    "2": parameters.Quote (partitioned.rate_classes) // the number of rate classes
                },
                partitioned.filter_names,
                None);
            
            partitioned.BSREL_bg_distributions[partitioned.k] = models.codon.BS_REL.ExtractMixtureDistribution(partitioned.BSREL_models[partitioned.model_id_bg]);
        }
        
         partitioned.BSREL_model_map [partitioned.k] = { partitioned.model_id : utility.Filter (partitioned.selected_branches[partitioned.k], '_value_', '_value_ == terms.tree_attributes.test'),
                                                         partitioned.model_id_bg : utility.Filter (partitioned.selected_branches[partitioned.k], '_value_', '_value_ != terms.tree_attributes.test')};

         partitioned.BSREL_distributions [partitioned.k] = models.codon.BS_REL.ExtractMixtureDistribution(partitioned.BSREL_models[partitioned.model_id]);
     
         for (partitioned.i = 1; partitioned.i < partitioned.rate_classes; partitioned.i += 1) {
             parameters.SetRange (model.generic.GetGlobalParameter (partitioned.BSREL_models[partitioned.model_id] , terms.AddCategory (terms.parameters.omega_ratio,partitioned.i)), terms.range01);
         }
         parameters.SetRange (model.generic.GetGlobalParameter (partitioned.BSREL_models[partitioned.model_id] , 
                                                              terms.AddCategory (terms.parameters.omega_ratio, partitioned.rate_classes)), 
                                                              terms.range_gte1);
    }



    partitioned.model_keys = utility.Keys (partitioned.BSREL_models);

    for (partitioned.k = 1; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
        models.BindGlobalParameters ({"0" : partitioned.BSREL_models[partitioned.model_keys[0]], "1" : partitioned.BSREL_models[partitioned.model_keys[partitioned.k]]}, terms.nucleotideRate("[ACGT]","[ACGT]"));
    }



    io.ReportProgressMessageMD ("partitioned", "BSREL", "Fitting partition-specific BS-REL models with " + partitioned.rate_classes + " rate classes");


    //utility.ToggleEnvVariable("VERBOSITY_LEVEL", 10);


    partitioned.bsrel_fit = estimators.FitLF(
                                            partitioned.filter_names,
                                            partitioned.trees,
                                            partitioned.BSREL_model_map,
                                            partitioned.bs_rel_init ,
                                            partitioned.BSREL_models,
                                            {
                                                terms.run_options.apply_user_constraints : "partitioned.define_bsrel_constraints",
                                                terms.run_options.retain_lf_object : TRUE
                                            }
                                          );
    KeywordArgument ("save-fit", "Save BUSTED model fit to this file (default is not to save)", "/dev/null");
    io.SpoolLFToPath(partitioned.bsrel_fit[terms.likelihood_function], io.PromptUserForFilePath ("Save BS-REL model fit to this file ['/dev/null' to skip]"));

    /*fprintf ("/Users/sergei/Desktop/fit.dump", CLEAR_FILE, partitioned.bsrel_fit);
    Export (lfe, ^(partitioned.bsrel_fit[terms.likelihood_function]));
    fprintf ("/Users/sergei/Desktop/lf.dump", CLEAR_FILE, lfe);

    return 0;*/

    //fscanf ("/Users/sergei/Desktop/fit.dump", "Raw", partitioned.bsrel_fit); partitioned.bsrel_fit = Eval (partitioned.bsrel_fit);
    //ExecuteAFile ("/Users/sergei/Desktop/lf.dump");

    io.ReportProgressMessageMD("partitioned", "BSREL", "* " + selection.io.report_fit (partitioned.bsrel_fit, 0, partitioned.codon_data_info[terms.data.sample_size]));

    partitioned.has_selection = FALSE;

    for (partitioned.k = 0; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
        io.ReportProgressMessageMD("partitioned", "BSREL", "* Inferred rate distribution for \`" +  (partitioned.partitions_and_trees [partitioned.k])[terms.data.name] + "\`");
        partitioned.inferred_test_distribution = parameters.GetStickBreakingDistribution (partitioned.BSREL_distributions [partitioned.k]) % 0;
        selection.io.report_dnds (partitioned.inferred_test_distribution);
        partitioned.has_selection = partitioned.has_selection || partitioned.inferred_test_distribution[partitioned.rate_classes-1][0] > 1 &&  partitioned.inferred_test_distribution[partitioned.rate_classes-1][1] > 0;
        if (partitioned.has_background) {
            io.ReportProgressMessageMD("partitioned", "BSREL", "* _background_ rate distribution for \`" +  (partitioned.partitions_and_trees [partitioned.k])[terms.data.name] + "\`");
            partitioned.inferred_test_distribution = parameters.GetStickBreakingDistribution (partitioned.BSREL_bg_distributions[partitioned.k]) % 0;
            selection.io.report_dnds (partitioned.inferred_test_distribution);
   
        }
    } 

    if (partitioned.has_selection) {
        io.ReportProgressMessageMD ("partitioned", "BSREL-ns", "Testing for positive selection on any of the partitions");
        for (partitioned.k = 0; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
         parameters.SetConstraint (model.generic.GetGlobalParameter (partitioned.BSREL_models[partitioned.bsrel_model_id (partitioned.k)] , 
                                                              terms.AddCategory (terms.parameters.omega_ratio, partitioned.rate_classes)), 
                                                              "1","");
        }

        partitioned.bsrel_fit.no_selection = estimators.FitExistingLF (partitioned.bsrel_fit[terms.likelihood_function], partitioned.BSREL_models);
        io.ReportProgressMessageMD("partitioned", "BSREL-ns", "* " + selection.io.report_fit (partitioned.bsrel_fit.no_selection, 0, partitioned.codon_data_info[terms.data.sample_size]));
    
        for (partitioned.k = 0; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
            io.ReportProgressMessageMD("partitioned", "BSREL-ns", "* Inferred rate distribution for \`" +  (partitioned.partitions_and_trees [partitioned.k])[terms.data.name] + "\`");
            partitioned.inferred_test_distribution = parameters.GetStickBreakingDistribution (partitioned.BSREL_distributions [partitioned.k]) % 0;
            selection.io.report_dnds (partitioned.inferred_test_distribution);
            if (partitioned.has_background) {
                io.ReportProgressMessageMD("partitioned", "BSREL", "* _background_ rate distribution for \`" +  (partitioned.partitions_and_trees [partitioned.k])[terms.data.name] + "\`");
                partitioned.inferred_test_distribution = parameters.GetStickBreakingDistribution (partitioned.BSREL_bg_distributions[partitioned.k]) % 0;
                selection.io.report_dnds (partitioned.inferred_test_distribution);
   
            }
        }

        io.ReportProgressMessageMD("partitioned", "BSREL-ns", "* Test for positive selection on any of the partitions:  " + selection.io.report_lrt (partitioned.ComputeLRT (partitioned.bsrel_fit, partitioned.bsrel_fit.no_selection, "Test for positive selection on any of the partitions")));
    
        for (partitioned.k = 0; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
            parameters.RemoveConstraint (model.generic.GetGlobalParameter (partitioned.BSREL_models[partitioned.bsrel_model_id (partitioned.k)] , 
                                                              terms.AddCategory (terms.parameters.omega_ratio, partitioned.rate_classes)));
        }
    
        estimators.ApplyExistingEstimates(partitioned.bsrel_fit[terms.likelihood_function], partitioned.BSREL_models, partitioned.bsrel_fit, None);

    }


    io.ReportProgressMessageMD ("partitioned", "BSREL-null", "Testing for equality of rate distributions across partitions");

    for (partitioned.k = 1; partitioned.k < partitioned.partition_count; partitioned.k+=1) {
        models.BindGlobalParameters ({"0" : partitioned.BSREL_models[partitioned.model_keys[0]], "1" : partitioned.BSREL_models[partitioned.model_keys[partitioned.k]]}, terms.parameters.omega_ratio);
        models.BindGlobalParameters ({"0" : partitioned.BSREL_models[partitioned.model_keys[0]], "1" : partitioned.BSREL_models[partitioned.model_keys[partitioned.k]]}, terms.mixture.mixture_aux_weight);
    }

    partitioned.bsrel_fit.null = estimators.FitExistingLF (partitioned.bsrel_fit[terms.likelihood_function], partitioned.BSREL_models);

    io.ReportProgressMessageMD("partitioned", "BSREL-null", "* " + selection.io.report_fit (partitioned.bsrel_fit.null, 0, partitioned.codon_data_info[terms.data.sample_size]));
    io.ReportProgressMessageMD("partitioned", "BSREL-null", "* Test for selective differences between partitions:  " + selection.io.report_lrt (partitioned.ComputeLRT (partitioned.bsrel_fit, partitioned.bsrel_fit.null, "Test for selective differences between partitions")));

    io.ReportProgressMessageMD("partitioned", "BSREL-null", "* Inferred joint rate distribution");
    partitioned.inferred_test_distribution = parameters.GetStickBreakingDistribution (partitioned.BSREL_distributions [0]) % 0;
    selection.io.report_dnds (partitioned.inferred_test_distribution);
}


console.log ("----\n## Partitioned analyses test summary (corrected p-values)");
partitioned.p_values.corrected = math.HolmBonferroniCorrection (partitioned.p_values);
utility.ForEachPair (partitioned.p_values.corrected, "_name_", "_p_", 
'
console.log ( "* Likelihood ratio test for " + _name_ + ", **p = " + Format (_p_, 8, 4) + "**.");
'
);



selection.io.stopTimer (partitioned.json [terms.json.timers], "Overall");
io.SpoolJSON (partitioned.json, partitioned.codon_data_info [terms.json.json]);
return partitioned.json;


//---------

function partitioned.make_key (filter_index, branch_set) {
    return (partitioned.partitions_and_trees [filter_index])[terms.data.name] + "|" + branch_set;
}


function partitioned.define_omega_constraints (lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects) {
    // obtain the list of branch partitions
    partitioned.define_omega_constraints.branch_sets = {};
    utility.ForEach (partitioned.selected_branches, "_list_", '
         utility.ForEachPair (utility.UniqueValues (_list_), "_ignore_", "_set_",
            "partitioned.define_omega_constraints.branch_sets[_set_] = TRUE");
    ');
    partitioned.define_omega_constraints.branch_sets = utility.Keys (partitioned.define_omega_constraints.branch_sets);
    partitioned.define_omega_constraints.model_key = (utility.Keys (model_objects))[0];
    partitioned.define_omega_constraints.i = 0;
    partitioned.define_omega_constraints.scalers = {};
    
    utility.ForEachPair (data_filter, "_index_", "_value_",
    '
        for (partitioned.define_omega_constraints.k = 0; partitioned.define_omega_constraints.k < utility.Array1D (partitioned.define_omega_constraints.branch_sets); partitioned.define_omega_constraints.k+=1) {
            partitioned.define_omega_constraints.omega = "partitioned.omega._" + partitioned.define_omega_constraints.i;
            partitioned.define_omega_constraints.key = partitioned.make_key (_index_,partitioned.define_omega_constraints.branch_sets[partitioned.define_omega_constraints.k]);

            partitioned.omega_ratios [partitioned.define_omega_constraints.key] = partitioned.define_omega_constraints.omega;

            model.generic.AddGlobal (model_objects[partitioned.define_omega_constraints.model_key],
                                     partitioned.define_omega_constraints.omega ,
                                     (utility.getGlobalValue("terms.parameters.omega_ratio")) + " for *" + partitioned.define_omega_constraints.key + "*");

            parameters.DeclareGlobal    (partitioned.define_omega_constraints.omega, None);

            
            
            partitioned.define_omega_constraints.i += 1;
        }
    ');
    for (partitioned.define_omega_constraints.k = 0; partitioned.define_omega_constraints.k < partitioned.partition_count; partitioned.define_omega_constraints.k += 1) {
        // iterate over branches in the tree
        partitioned.define_omega_constraints.model_id = (model_map[partitioned.define_omega_constraints.k])["DEFAULT"];
        
        partitioned.define_omega_constraints.alpha = model.generic.GetLocalParameter(partitioned.MG_models[partitioned.define_omega_constraints.model_id], utility.getGlobalValue("terms.parameters.synonymous_rate"));
        partitioned.define_omega_constraints.beta  = model.generic.GetLocalParameter(partitioned.MG_models[partitioned.define_omega_constraints.model_id], utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
        partitioned.define_omega_constraints.treeID = (lf_components[2*partitioned.define_omega_constraints.k+1]);
        assert (Abs( partitioned.define_omega_constraints.alpha) > 0 && Abs ( partitioned.define_omega_constraints.beta), "Local parameters are not defined for `partitioned.define_omega_constraints.model_id`");

        utility.ForEach (^partitioned.define_omega_constraints.treeID, "_branch_name_", '
             parameters.SetConstraint("`partitioned.define_omega_constraints.treeID`."+_branch_name_+".`partitioned.define_omega_constraints.beta`", 
                                      "`partitioned.define_omega_constraints.treeID`."+_branch_name_+".`partitioned.define_omega_constraints.alpha`*" + 
                                        partitioned.omega_ratios [partitioned.make_key  (partitioned.define_omega_constraints.k, (partitioned.selected_branches[partitioned.define_omega_constraints.k])[_branch_name_])], 
                                      "");
            
        ');
 
         if (partitioned.proportional_lengths && partitioned.define_omega_constraints.k) {
            partitioned.define_omega_constraints.scaler = "partitioned.tree_scaler._" + partitioned.define_omega_constraints.k;
            parameters.DeclareGlobalWithRanges (partitioned.define_omega_constraints.scaler  , 1, 0, 1000);
          
            model.generic.AddGlobal (model_objects[partitioned.define_omega_constraints.model_key],
                                     partitioned.define_omega_constraints.scaler ,
                                     (utility.getGlobalValue("terms.model.branch_length_scaler")) + " for partition *" + partitioned.make_key (partitioned.define_omega_constraints.k, "*"));
            ReplicateConstraint ("this1.?.`partitioned.define_omega_constraints.alpha`:=`partitioned.define_omega_constraints.scaler`*this2.?.`partitioned.define_omega_constraints.alpha`", ^partitioned.define_omega_constraints.treeID, ^partitioned.define_omega_constraints.treeID.ref)
 
         } else {
            partitioned.define_omega_constraints.treeID.ref = (lf_components[2*partitioned.define_omega_constraints.k+1]);
         }
   }
    
    // apply proportional branch constraints
    
    
    
    return 0;
}

//------------------------------------------------------------------------------

function partitioned.bsrel_model_id (component) {
    return "partitioned.BSREL.model" + component;
}

//------------------------------------------------------------------------------

function partitioned.bsrel_model_id_bg (component) {
    return "partitioned.BSREL.background_model" + component;
}

//------------------------------------------------------------------------------

function partitioned.define_bsrel_constraints (lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects) {

    
    
    for (partitioned.define_bsrel_constraints.k = 0; partitioned.define_bsrel_constraints.k < partitioned.partition_count; partitioned.define_bsrel_constraints.k += 1) {
          if (partitioned.proportional_lengths && partitioned.define_bsrel_constraints.k) {
            partitioned.define_bsrel_constraints.scaler = "partitioned.bsrel.tree_scaler._" + partitioned.define_bsrel_constraints.k;
            parameters.DeclareGlobalWithRanges (partitioned.define_bsrel_constraints.scaler , 1, 0, 1000);
            
            partitioned.define_bsrel_constraints.treeID = (lf_components[2*partitioned.define_bsrel_constraints.k+1]);
            model.generic.AddGlobal (model_objects[partitioned.bsrel_model_id (partitioned.define_bsrel_constraints.k)],
                                     partitioned.define_bsrel_constraints.scaler ,
                                     (utility.getGlobalValue("terms.model.branch_length_scaler")) + " for partition *" + partitioned.make_key (partitioned.define_bsrel_constraints.k, "*"));
            ReplicateConstraint ("this1.?.`partitioned.define_bsrel_constraints.time`:=`partitioned.define_bsrel_constraints.scaler`*this2.?.`partitioned.define_bsrel_constraints.time`", ^partitioned.define_bsrel_constraints.treeID, ^partitioned.define_bsrel_constraints.treeID.ref)
 
         } else {
            partitioned.define_bsrel_constraints.treeID.ref = (lf_components[2*partitioned.define_bsrel_constraints.k+1]);
            partitioned.define_bsrel_constraints.time = Call ((model_objects[partitioned.bsrel_model_id (partitioned.define_bsrel_constraints.k)])[terms.model.time], 0);
         }
    }
    
    
    return 0;
}

//------------------------------------------------------------------------------

function partitioned.ComputeLRT (ha, h0, tag) {
    
    partitioned.p_values [tag] = 1-CChi2 (2*(ha[terms.fit.log_likelihood]-h0[terms.fit.log_likelihood]),ha[terms.parameters]-h0[terms.parameters]);
    
    return {terms.LRT : 2*(ha[terms.fit.log_likelihood]-h0[terms.fit.log_likelihood]),
            terms.p_value : partitioned.p_values [tag]};
}