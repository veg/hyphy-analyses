RequireVersion ("2.4.0");


LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/models/model_functions.bf");

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

partitioned.mg_fit = estimators.FitLF(
                                        partitioned.filter_names,
                                        partitioned.trees,
                                        partitioned.model_map,
                                        partitioned.gtr_results,
                                        partitioned.MG_models,
                                        {
                                            terms.run_options.apply_user_constraints : "partitioned.define_omega_constraints"
                                        }
                                      );

console.log (partitioned.mg_fit);

function partitioned.make_key (filter_index, branch_set) {
    return (partitioned.partitions_and_trees [filter_index])[terms.data.name] + "|" + branch_set;
}

/**
new_globals = {};
utility.ForEachPair(partition_omega, "_key_", "_value_",
    '`&new_globals` [_key_] = (`&name_space` + ".omega_" + Abs (`&new_globals`)); model.generic.AddGlobal (`&mg_rev`, `&new_globals` [_key_] , (utility.getGlobalValue("terms.parameters.omega_ratio")) + " for *" + _key_ + "*")');
parameters.DeclareGlobal(new_globals, None);
**/

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
    utility.ForEachPair (data_filter, "_index_", "_value_",
    '
        for (partitioned.define_omega_constraints.k = 0; partitioned.define_omega_constraints.k < utility.Array1D (partitioned.define_omega_constraints.branch_sets); partitioned.define_omega_constraints.k+=1) {
            partitioned.define_omega_constraints.omega = "partitioned.omega._" + partitioned.define_omega_constraints.i;
            partitioned.define_omega_constraints.key = partitioned.make_key (_index_,partitioned.define_omega_constraints.branch_sets[partitioned.define_omega_constraints.k]);

            partitioned.omega_ratios [partitioned.define_omega_constraints.key] = partitioned.define_omega_constraints.omega;

            model.generic.AddGlobal (model_objects[partitioned.define_omega_constraints.model_key],
                                     partitioned.define_omega_constraints.omega ,
                                     (utility.getGlobalValue("terms.parameters.omega_ratio")) + " for partition|branchset *" + partitioned.define_omega_constraints.key + "*");

            parameters.DeclareGlobal    (partitioned.define_omega_constraints.omega, None);
            partitioned.define_omega_constraints.i += 1;
        }
    ');
    console.log (partitioned.omega_ratios);
    assert (0);
    return 0;
}
