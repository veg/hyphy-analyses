
function pogofit.fitGTR_fixalpha (current_results) {

    filter_info    = {};
    trees = {};
    initial_values = {terms.global : {}, terms.branch_length : {}};
    index_to_file_name   = {};
    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
        file_path = pogofit.file_list [file_index];
        dataset_name = "pogofit.msa.part" + file_index;
        data_info = alignments.ReadNucleotideDataSet (dataset_name, file_path);
        data_info = alignments.EnsureMapping (dataset_name, data_info);

        partition_specification = { "0" : {terms.data.name : "all", terms.data.filter_string : "", term.data.tree : ((current_results[file_index])[terms.fit.trees])[0]}};


        this_bl = (current_results[terms.branch_length])[file_index];
        this_tree = (current_results[terms.fit.trees])[file_index];

        filter_info [file_index] = (alignments.DefineFiltersForPartitions (partition_specification,
                                                                            dataset_name ,
                                                                            dataset_name,
                                                                            data_info))[0];
        trees [file_index] = {terms.trees.newick :  this_tree};
        (initial_values[terms.branch_length])[file_index] = this_bl;
    }
    filter_names = utility.Map (filter_info, "_value_", "_value_[terms.data.name]");

    utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
    
    // run options including fix branch lengths
    fit_options = {
        terms.run_options.proportional_branch_length_scaler: {
        },
        terms.run_options.optimization_settings : 
        {
            "OPTIMIZATION_METHOD" : "coordinate-wise",
            "OPTIMIZATION_PRECISION" : pogofit.precision
        },
        terms.run_options.retain_lf_object : TRUE
    };
    //(fit_options[terms.run_options.proportional_branch_length_scaler])[0] = terms.model.branch_length_constrain;

    
    // set intial values to user chosen matrix (same as baseline)
    for (l1 = 0; l1 < 20; l1 += 1) {
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
            rate_term = terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]);
            (initial_values[terms.global]) [rate_term] = {terms.fit.MLE : (pogofit.initial_rates[models.protein.alphabet[l1]])[models.protein.alphabet[l2]]}; 
        }
    }
    // TODO: SOMETHING CHANGED WITH THE CODE REFACTOR??
    alpha_term = "Gamma distribution shape parameter";
    alpha = ((current_results[terms.global])[alpha_term])[terms.fit.MLE];
    (initial_values[terms.global]) [alpha_term] = {terms.fit.MLE : alpha , terms.fix : TRUE};

    
    pogofit.rev_mle = estimators.FitSingleModel_Ext (
                                        filter_names,
                                        trees,
                                        pogofit.rev_model,
                                        initial_values,
                                        fit_options
                                   );                         

                  

    /*   
    // Uncomment these lines if you'd like to save the NEXUS LF.                        
    lf_id = pogofit.rev.mle[terms.likelihood_function];
    Export(pogofit.finalphase_LF, ^lf_id);
    fprintf(pogofit.final_likelihood_function, pogofit.finalphase_LF);
    */
    
    // Save the rev.mle into the analysis_results, and cache it.
    (^"pogofit.analysis_results")[pogofit.final_phase] = pogofit.rev_mle;

    console.log (""); // clear past the optimization progress line
    utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);


    // Trees as dictionary for compatibility with rest of the output.
    pogofit.rev_mle[terms.fit.trees] = utility.SwapKeysAndValues(utility.MatrixToDict(pogofit.rev_mle[terms.fit.trees]));
    
    return pogofit.rev_mle;

}



function pogofit.fitBaselineTogether () {


    filter_info    = {};
    trees = {};
    index_to_file_name   = {};
    
    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
        file_path = pogofit.file_list [file_index];
        dataset_name = "pogofit.msa.part" + file_index;
        
        file_info = alignments.ReadNucleotideDataSet (dataset_name, file_path);
        file_info = alignments.EnsureMapping(dataset_name, file_info);

        utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
        ExecuteCommands ('partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (file_info[terms.data.partitions], pogofit.name_mapping)',
                         {"0" : "Y"});
        utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);

        tree = utility.Map (partitions_and_trees, "_value_", '_value_[terms.data.tree]');

        tree_with_lengths = (tree["0"])[terms.trees.newick_with_lengths];
        
        partition_specification = { "0" : {
                                            terms.data.name : "all", 
                                            terms.data.filter_string : "", 
                                            terms.data.tree : tree_with_lengths}};
             
        filter_info [file_index] = (alignments.DefineFiltersForPartitions (partition_specification,
                                                                            dataset_name ,
                                                                            dataset_name,
                                                                            file_info))[0];                                                        
        trees [file_index] = {terms.trees.newick : tree_with_lengths};        
    }

    utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", 0);

    pogofit.baseline_mle = estimators.FitSingleModel_Ext (
                                            utility.Map (filter_info, "_value_", "_value_[terms.data.name]"),
                                            trees,
                                            pogofit.baseline_model_desc_gamma,
                                            None,
                                            None
                                           );
    
    // Save the rev.mle into the analysis_results, and cache it.
    (^"pogofit.analysis_results")[pogofit.baseline_phase] = pogofit.baseline_mle;

    console.log (""); // clear past the optimization progress line
    utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);


    // Trees as dictionary for compatibility with rest of the output.
    pogofit.baseline_mle[terms.fit.trees] = utility.SwapKeysAndValues(utility.MatrixToDict(pogofit.baseline_mle[terms.fit.trees]));
    
    return pogofit.baseline_mle;

}











/********************************************************************************************************************/
/********************************************* MODEL DEFINITIONS ****************************************************/
/********************************************************************************************************************/

/**
 * @name models.protein.Baseline.ModelDescription.withGamma
 * @description Define baseline (standard matrix) model w/ +F and *no* four-category gamma rate variation
 */
function pogofit.Baseline.ModelDescription(type){
    def = Call( models.protein.empirical.plusF_generators[pogofit.baseline_model], type);
    return def;
}
/**
 * @name models.protein.Baseline.ModelDescription.withGamma
 * @description Define baseline (standard matrix) model w/ +F and *yes* four-category gamma rate variation
 */
function pogofit.Baseline.ModelDescription.withGamma(type){
    def = pogofit.Baseline.ModelDescription(type);
	def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.Gamma.factory ({utility.getGlobalValue("terms.rate_variation.bins") : 4});
    return def;
}


/**
 * @name pogofit.REV.ModelDescription.freqs
 * @description Define REV model frequencies as empirical
 */
function pogofit.REV.ModelDescription.freqs (model, namespace, datafilter) {
    model[terms.efv_estimate] = pogofit.shared_EFV;
    model[terms.model.efv_estimate_name] = terms.frequencies.predefined;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}

/**
 * @name pogofit.REV.ModelDescription
 * @description Define a REV model with constant site rates
 */
function pogofit.REV.ModelDescription (type) {
    def = models.protein.REV.ModelDescription(type);
    if (Type (pogofit.shared_EFV) == "Matrix") {
        def [terms.model.frequency_estimator] = "pogofit.REV.ModelDescription.freqs";
    }
    return def;
}

/**
 * @name pogofit.REV.ModelDescription.withGamma
 * @description Define a REV model with Gamma rate variation
 */
function pogofit.REV.ModelDescription.withGamma (type) {
    def = pogofit.REV.ModelDescription(type);
    def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.Gamma.factory ({utility.getGlobalValue("terms.rate_variation.bins"): 4});
    return def;
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/



/********************************************************************************************************************/
/********************************************** FITTING FUNCTIONS ***************************************************/
/********************************************************************************************************************/

/**
 * @name pogofit.fitBaselineToFile
 * @description Fits an empirical amino acid model to dataset for branch length optimization
 * @param {String} filename - The name of the file containing the dataset to which the amino acid model will be fitted
 * @return the fitted MLE
 */
function pogofit.fitBaselineToFile (filename) {
    
    
   // utility.EnsureKey(pogofit.analysis_results, pogofit.index_to_filename[filename]);  
    
    pogofit.file_info = alignments.ReadNucleotideDataSet ("pogofit.msa",
                                                              filename);
    pogofit.name_mapping = pogofit.file_info[utility.getGlobalValue("terms.data.name_mapping")];
    if (None == pogofit.name_mapping) { /** create a 1-1 mapping if nothing was done */
        pogofit.name_mapping = {};
        utility.ForEach (alignments.GetSequenceNames ("pogofit.msa"), "_value_", "`&pogofit.name_mapping`[_value_] = _value_");
    }
    utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
    ExecuteCommands ('pogofit.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (pogofit.file_info[terms.data.partitions], pogofit.name_mapping)',
                     {"0" : "Y"});
    utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);



    pogofit.partition_count      = Abs (pogofit.partitions_and_trees);
    io.CheckAssertion ("pogofit.partition_count==1", "This analysis can only handle a single partition");



    pogofit.filter_specification = alignments.DefineFiltersForPartitions (pogofit.partitions_and_trees,
                                                                            "pogofit.msa" ,
                                                                            "pogofit.filter.",
                                                                            pogofit.file_info);


    pogofit.full_trees = utility.Map (pogofit.partitions_and_trees, "_value_", '_value_[terms.data.tree]');
    pogofit.full_data_filter = utility.Map (pogofit.filter_specification, "_value_", "_value_[terms.data.name]");


    /********** Store dataset information *************/
    /* CURRENTLY DOESN'T WORK IN MPI FOR REASONS TBD */
//     pogofit.output_data_info = { utility.getGlobalValue("terms.json.sequences"): pogofit.file_info[utility.getGlobalValue("terms.data.sequences")],
//                                      utility.getGlobalValue("terms.json.sites"): pogofit.file_info[utility.getGlobalValue("terms.data.sites")],
//                                      utility.getGlobalValue("terms.json.trees"): (pogofit.full_trees["0"])[utility.getGlobalValue("terms.trees.newick_with_lengths")],
//                                      utility.getGlobalValue("terms.json.file"): filename
//                                    };
// 
// 
    // In case there were no branch lengths
//     if (  Abs( (pogofit.full_trees["0"])[utility.getGlobalValue("terms.branch_length")] ) == 0 ){
//         pogofit.output_data_info[ utility.getGlobalValue("terms.json.trees") ] = (pogofit.full_trees["0"])[utility.getGlobalValue("terms.trees.newick")];
//     }
//     utility.ForEach (utility.Keys (pogofit.name_mapping), "branch_name",
//                              "utility.EnsureKey (pogofit.output_data_info[terms.original_name], branch_name)");
// 
//     utility.ForEach (utility.Keys (pogofit.name_mapping), "branch_name",
//                              "(pogofit.output_data_info[terms.original_name])[branch_name] = pogofit.name_mapping[branch_name]");
// 
// 
//     (pogofit.analysis_results[pogofit.index_to_filename[filename]])[utility.getGlobalValue("terms.json.input")] = pogofit.output_data_info;
// 
//     
    /****************************************************/

    
    pogofit.baseline_mle = estimators.FitSingleModel_Ext(pogofit.full_data_filter,
                                                        pogofit.full_trees,
                                                        pogofit.baseline_model_desc,
                                                        None,
                                                        None);
                                                        //
                                                        
    pogofit.baseline_mle - terms.global; // delete empty key
    return pogofit.baseline_mle;
}



/**
 * @name pogofit.handle_baseline_callback
 * @param node - node name which processed the given data
 * @param {Dict} result - Dictionary of fitted information for given data
 * @param {Dict} arguments - Dictionary with single key:value :: 0:datafile name
 * @description Handle MPI callback after fitting a baseline amino acid model (for initial branch length optimization)
 */
function pogofit.handle_baseline_callback (node, result, arguments) {

    savekey = pogofit.index_to_filename[arguments[0]];
    
    utility.EnsureKey(pogofit.analysis_results, savekey);
    utility.EnsureKey(pogofit.analysis_results[savekey], pogofit.baseline_phase);
    (pogofit.analysis_results[savekey])[pogofit.baseline_phase] = result;
     
    io.ReportProgressMessageMD ("Protein GTR Fitter", "Initial branch length fit",
                                "Received file '" + arguments[0] + "' from node " + node + ". LogL = " + result[terms.fit.log_likelihood]);
}





function pogofit.fitGTR_gamma (current_bl, current_gtr, phase) {

    //file_list = utility.Keys (current_results); ---> pogofit.file_list
    //file_count = utility.Array1D (file_list);   ---> pogofit.file_list_count
    // NOTE: pogofit.index_to_filename is {filename:0, filename:1}

    partition_info = {};
    filter_info    = {};
    trees = {};
    initial_values = {terms.global : {}, terms.branch_length : {}};
    index_to_file_name   = {};



    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
        file_path = pogofit.file_list [file_index];
        dataset_name = "pogofit.msa.part" + file_index;
        partition_info [file_index] = alignments.ReadNucleotideDataSet (dataset_name, file_path);
        partition_specification = { "0" : {terms.data.name : "all", terms.data.filter_string : "", terms.data.tree : ((current_bl[file_index])[terms.fit.trees])[0]}};

        filter_info [file_index] = (alignments.DefineFiltersForPartitions (partition_specification,
                                                                            dataset_name ,
                                                                            dataset_name,
                                                                            partition_info [file_index]))[0];
        trees [file_index] = {terms.trees.newick :  ((current_bl[file_index])[terms.fit.trees])[0]};
        (initial_values[terms.branch_length])[file_index] = ((current_bl[file_index])[terms.branch_length])[0];
    }
    initial_values[terms.global] = current_gtr[terms.global];



    utility.SetEnvVariable    ("VERBOSITY_LEVEL", 1);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
 

   pogofit.rev_mle = estimators.FitSingleModel_Ext (
                                       utility.Map (filter_info, "_value_", "_value_[terms.data.name]"),
                                        trees,
                                        pogofit.rev_model_gamma,
                                        initial_values,
                                       {terms.run_options.retain_lf_object : TRUE,
                                       terms.run_options.optimization_settings : 
                                            {
                                                "OPTIMIZATION_METHOD" : "nedler-mead",
                                                "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
                                                "OPTIMIZATION_PRECISION" : 1
                                            }
                                         }
                                       );
    /*   
    // Uncomment these lines if you'd like to save the NEXUS LF.                        
    lf_id = pogofit.rev.mle[terms.likelihood_function];
    Export(pogofit.finalphase_LF, ^lf_id);
    fprintf(pogofit.final_likelihood_function, pogofit.finalphase_LF);
    */
    pogofit.rev_mle - terms.likelihood_function;
    
    // Save the rev.mle into the analysis_results, and cache it.
    (^"pogofit.analysis_results")[this_phase] = pogofit.rev_mle;

    console.log (""); // clear past the optimization progress line
    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);


    // Trees as dictionary for compatibility with rest of the output.
    pogofit.rev_mle[terms.fit.trees] = utility.SwapKeysAndValues(utility.MatrixToDict(pogofit.rev_mle[terms.fit.trees]));
    
    return pogofit.rev_mle;

}





function pogofit.fitGTR_twophase(current_results, phase, isfinalphase) {

    //file_list = utility.Keys (current_results); ---> pogofit.file_list
    //file_count = utility.Array1D (file_list);   ---> pogofit.file_list_count
    // NOTE: pogofit.index_to_filename is {filename:0, filename:1}

    partition_info = {};
    filter_info    = {};
    trees = {};
    initial_values = {terms.global : {}, terms.branch_length : {}};
    proportional_scalers = {};
    index_to_file_name   = {};

    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
        file_path = pogofit.file_list [file_index];
        dataset_name = "pogofit.msa.part" + file_index;
        partition_info [file_index] = alignments.ReadNucleotideDataSet (dataset_name, file_path);
        partition_specification = { "0" : {terms.data.name : "all", terms.data.filter_string : "", terms.data.tree : ((current_results[file_index])[terms.fit.trees])[0]}};


        filter_info [file_index] = (alignments.DefineFiltersForPartitions (partition_specification,
                                                                            dataset_name ,
                                                                            dataset_name,
                                                                            partition_info [file_index]))[0];
        trees [file_index] = {terms.trees.newick :  ((current_results[file_index])[terms.fit.trees])[0]};
        (initial_values[terms.branch_length])[file_index] = ((current_results[file_index])[terms.branch_length])[0];
        if (phase == pogofit.prefinal_phase) {
            scaler = "pogofit.gtr_scaler_" + file_index;
            parameters.DeclareGlobalWithRanges (scaler, 1, 0, 1000);
            proportional_scalers[file_index] = scaler;
        }

    
    }

    utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", 0);

    if (phase == pogofit.prefinal_phase) {
        // Set initial values 
        for (l1 = 0; l1 < 20; l1 += 1) {
            for (l2 = l1 + 1; l2 < 20; l2 += 1) {
                (initial_values[terms.global]) [terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2])] = {terms.fit.MLE : 0.1}; // set all to 1
            }
        }
        // fit the model
        pogofit.mle = estimators.FitSingleModel_Ext (
                                            utility.Map (filter_info, "_value_", "_value_[terms.data.name]"),
                                            trees,
                                            pogofit.rev_model,
                                            initial_values,
                                            {terms.run_options.proportional_branch_length_scaler : proportional_scalers}
                                       );

    } else
    {
        // FINAL TUNING
        pogofit.mle = estimators.FitSingleModel_Ext (
                                            utility.Map (filter_info, "_value_", "_value_[terms.data.name]"),
                                            trees,
                                            pogofit.rev_model,
                                            initial_values,
                                             {terms.run_options.retain_lf_object : TRUE}
                                       );
    }

                                      
    // Save the rev.mle into the analysis_results, and cache it.
    (^"pogofit.analysis_results")[phase] = pogofit.mle;

    console.log (""); // clear past the optimization progress line
    utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);


    // Trees as dictionary for compatibility with rest of the output.
    pogofit.mle[terms.fit.trees] = utility.SwapKeysAndValues(utility.MatrixToDict(pogofit.mle[terms.fit.trees]));
    
    return pogofit.mle;

}




/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/



/********************************************************************************************************************/
/********************************************** UTILITY FUNCTIONS ***************************************************/
/********************************************************************************************************************/

function pogofit.save_options() {
    pogofit.analysis_results[utility.getGlobalValue("terms.json.options")] = {utility.getGlobalValue("pogofit.options.frequency_type"): pogofit.frequency_type,
                                                                              utility.getGlobalValue("pogofit.options.baseline_model"): pogofit.baseline_model,
                                                                              utility.getGlobalValue("pogofit.options.imputation"): pogofit.imputation};

    pogofit.analysis_results[utility.getGlobalValue("terms.json.input")] = {utility.getGlobalValue("terms.json.file"): pogofit.input_file,
                                                                            pogofit.options.number_of_datasets: pogofit.file_list_count,
                                                                            pogofit.options.dataset_information: {}};


    /* Temporarily, we save input file information here in a highly redundant fashion, but doesn't seem possible to do in MPI...? */
    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
    
        filename = pogofit.file_list[file_index];
        utility.EnsureKey(pogofit.analysis_results, file_index);  //TODO POSSIBLE REMOVE
    
        pogofit.file_info = alignments.ReadNucleotideDataSet ("pogofit.msa",
                                                                      filename);
        pogofit.name_mapping = pogofit.file_info[utility.getGlobalValue("terms.data.name_mapping")];
        if (None == pogofit.name_mapping) { /** create a 1-1 mapping if nothing was done */
            pogofit.name_mapping = {};
            utility.ForEach (alignments.GetSequenceNames ("pogofit.msa"), "_value_", "`&pogofit.name_mapping`[_value_] = _value_");
        }
        utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
        ExecuteCommands ('pogofit.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (pogofit.file_info[terms.data.partitions], pogofit.name_mapping)',
                          {"0" : "Y"});
        utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);


        pogofit.filter_specification = alignments.DefineFiltersForPartitions (pogofit.partitions_and_trees,
                                                                                    "pogofit.msa" ,
                                                                                    "pogofit.filter.",
                                                                                    pogofit.file_info);
        pogofit.tree = utility.Map (pogofit.partitions_and_trees, "_value_", '_value_[terms.data.tree]');


        pogofit.output_data_info = { utility.getGlobalValue("terms.json.sequences"): pogofit.file_info[utility.getGlobalValue("terms.data.sequences")],
                                         utility.getGlobalValue("terms.json.sites"): pogofit.file_info[utility.getGlobalValue("terms.data.sites")],
                                         utility.getGlobalValue("terms.json.trees"): (pogofit.tree["0"])[utility.getGlobalValue("terms.trees.newick_with_lengths")],
                                         utility.getGlobalValue("terms.json.file"): filename,
                                         utility.getGlobalValue("terms.original_name"): {}
                                       };
        
        
        //In case there were no branch lengths
        if (  Abs( (pogofit.tree["0"])[utility.getGlobalValue("terms.branch_length")] ) == 0 ){
            pogofit.output_data_info[ utility.getGlobalValue("terms.json.trees") ] = (pogofit.full_trees["0"])[utility.getGlobalValue("terms.trees.newick")];
        }
        utility.ForEach (utility.Keys (pogofit.name_mapping), "branch_name",
                                 "utility.EnsureKey (pogofit.output_data_info[terms.original_name], branch_name)");

        utility.ForEach (utility.Keys (pogofit.name_mapping), "branch_name",
                                 "(pogofit.output_data_info[terms.original_name])[branch_name] = pogofit.name_mapping[branch_name]");


        ((pogofit.analysis_results[terms.json.input])[pogofit.options.dataset_information])[file_index] = pogofit.output_data_info;
    }
}


lfunction pogofit.startTimer(timers, key) {
    timers[key] = {
        utility.getGlobalValue("terms.timers.timer"): Time(1),
    };

}
lfunction pogofit.stopTimer(timers, key) {
    (timers[key])[utility.getGlobalValue("terms.timers.timer")] = Time(1) - (timers[key])[utility.getGlobalValue("terms.timers.timer")];
}







// From the fitted model results, create EFV array which can be used for custom model 
function pogofit.extract_efv () {
    if (pogofit.frequency_type == pogofit.ml_freq){
        fitted_efv = {20, 1};
        for (i = 0; i < 20; i+=1)
        {
            efv_search = terms.characterFrequency(models.protein.alphabet[i]);  
            fitted_efv[i] = ((pogofit.gtr_fit[terms.global])[efv_search])[terms.fit.MLE];
    
        }
        norm =  +fitted_efv;
        fitted_efv = fitted_efv * (1/norm);
    }
    if (pogofit.frequency_type == pogofit.emp_freq) {
        fitted_efv = (pogofit.gtr_fit[terms.efv_estimate])["VALUEINDEXORDER"][0];
    }
    return fitted_efv
}




// From the fitted model results and create rate dictionary which can be used for .fitted_model
function pogofit.extract_rates() {

    rij = {};
    for (l1 = 0; l1 < 20 - 1; l1 += 1) {
        rij[models.protein.alphabet[l1]] = {};
        for (l2 = l1 + 1; l2 < 20; l2 += 1) { 
            rate_search = terms.aminoacidRate(models.protein.alphabet[l1], models.protein.alphabet[l2]);       
            (rij[models.protein.alphabet[l1]])[models.protein.alphabet[l2]] = ((pogofit.gtr_fit[terms.global])[rate_search])[terms.fit.MLE];
        }
    }
    return rij
}


// From the fitted model results and create rate dictionary which can be used for .fitted_model
function pogofit.extract_rates_imputation() {


    pogofit.thekeys =  utility.Keys(pogofit.site_counts);
    summation = 0.0;
    for (x = 0; x < utility.Array1D(pogofit.thekeys); x+=1)
    {
        this_key = pogofit.thekeys[x];
        summation += ((1. + pogofit.site_counts[this_key]) * pogofit.tree_lengths[this_key] );            
    }
    
    rij = {};
    for (l1 = 0; l1 < 20 - 1; l1 += 1) {
        rij[models.protein.alphabet[l1]] = {};
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
            rate_search = terms.aminoacidRate(models.protein.alphabet[l1], models.protein.alphabet[l2]);
            this_rate = ((pogofit.gtr_fit[terms.global])[rate_search])[terms.fit.MLE];
            if (this_rate == 0.0) 
            {
                efv_sum = pogofit.final_efv[l1] + pogofit.final_efv[l2];
                this_rate = 1./ ((pogofit.final_efv[l1] + pogofit.final_efv[l2]) * summation);
            }
            (rij[models.protein.alphabet[l1]])[models.protein.alphabet[l2]] = this_rate;
        }
    }
    return rij
}




function pogofit.write_model_to_file() {

    pogofit.external_order   = {{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}};
    pogofit.hyphy_order_dict = utility.MatrixToDict(models.protein.alphabet); 

    if (pogofit.output_format == pogofit.output_hyphy)
    {
        pogofit.save_hyphy_model();
    }

    if (pogofit.output_format == pogofit.output_paml)
    {
        pogofit.save_paml_model();
    }
    if (pogofit.output_format == pogofit.output_raxml)
    {
        pogofit.save_raxml_model();
    }
    if (pogofit.output_format == pogofit.output_all)
    {
        pogofit.save_raxml_model();
        pogofit.save_paml_model();
        pogofit.save_hyphy_model();
    }
}


function pogofit.save_hyphy_model(){
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, CLEAR_FILE, "Rij = " + pogofit.final_rij + ";");
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, "\n\n\n");
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, "EFV = " + pogofit.final_efv + ";");
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, "\n\n\nCI = " + pogofit.final_ci + ";");
}

function pogofit.save_raxml_model(){

    pogofit.raxml_output = "";
    pogofit.raxml_file = pogofit.output_model_prefix + pogofit.raxml_model_ext;
    
    for (i = 0; i < 20; i +=1)
    {
        for (j = 0; j < 20; j += 1)
        {
        
            if (i == j)
            {
                pogofit.raxml_output += "0.0\n";
            }
            else
            {
                pogofit.aa1 = pogofit.external_order[i];
                pogofit.aa2 = pogofit.external_order[j];
        
                if (pogofit.aa1 < pogofit.aa2)
                {
                    pogofit.rate = (pogofit.final_rij[pogofit.aa1])[pogofit.aa2];
                }
                else
                {
                    pogofit.rate = (pogofit.final_rij[pogofit.aa2])[pogofit.aa1];
                }
                pogofit.raxml_output += pogofit.rate;
                pogofit.raxml_output += "\n"; 
            }
        }    
    }
    for (i = 0; i < 20; i += 1)
    {
        pogofit.raxml_output += pogofit.final_efv[ pogofit.hyphy_order_dict[pogofit.external_order[i]] ];
        // strip hack
        if (i <= 18){
            pogofit.raxml_output += "\n";
        }
    }

    fprintf(pogofit.raxml_file, CLEAR_FILE, pogofit.raxml_output);
}




function pogofit.save_paml_model(){

    pogofit.paml_output = "";
    pogofit.paml_file = pogofit.output_model_prefix + pogofit.paml_model_ext;

    for (i = 1; i < 20; i +=1)
    {
        pogofit.row = "";
        for (j = 0; j < i; j += 1)
        {
        
            pogofit.aa1 = pogofit.external_order[i];
            pogofit.aa2 = pogofit.external_order[j];
        
            if (pogofit.aa1 < pogofit.aa2){
                pogofit.rate = (pogofit.final_rij[pogofit.aa1])[pogofit.aa2];
            }
            else{
                pogofit.rate = (pogofit.final_rij[pogofit.aa2])[pogofit.aa1];
            }
            
            pogofit.row += pogofit.rate;
            // strip hack
            if (j != (i-1)){
                pogofit.row += " ";
            }
        }
    
        pogofit.paml_output += pogofit.row + "\n";
    }

    pogofit.paml_output += "\n";
    for (i = 0; i < 20; i += 1)
    {
        pogofit.paml_output += pogofit.final_efv[ pogofit.hyphy_order_dict[pogofit.external_order[i]] ];
        // strip hack
        if (i <= 18){
            pogofit.paml_output += " ";
        }
    }
    fprintf(pogofit.paml_file, CLEAR_FILE, pogofit.paml_output);
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/





/**
 * @name pogofit.run_gtr_iteration_branch_lengths
 * @description Optimizes branch lengths for all datasets using the REV model fitted in the current iteration
 * @return Dictionary containing summed LogL values from branch length optimizations and the phase index for this iteration
 */
function pogofit.run_gtr_iteration_branch_lengths () {

    pogofit.queue = mpi.CreateQueue ({  "Headers"   : utility.GetListOfLoadedModules ("libv3/") ,
                                            "Functions" :
                                            {
                                                {"pogofit.REV.ModelDescription",
                                                 "pogofit.REV.ModelDescription.withGamma",
                                                 "pogofit.REV.ModelDescription.freqs",
                                                 "models.protein.REV.ModelDescription.withGamma"
                                                }
                                            },
                                            "Variables" : {{
                                                "pogofit.shared_EFV",
                                                "pogofit.rev_model_gamma",
                                                "pogofit.phase3",
                                                "pogofit.file_list_count"
                                            }}
                                         });


//console.log(pogofit.current_gtr_fit);
// global: sub rates
// EFV
// branch length
//// 0
////// [all the branches]
//// 1
//////  branch lengthssssss
// Trees (dict now!)
// logl
// parameters 


    io.ReportProgressMessageMD ("Protein GTR Fitter", pogofit.phase3, "Retuning branch lengths (" + pogofit.phase3 + ")");

    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
            
        io.ReportProgressMessageMD ("Protein GTR Fitter", pogofit.phase_key,
                                    "Dispatching file '" + pogofit.file_list[file_index] + "' " + (file_index+1) + "/" + pogofit.file_list_count);
        
        mpi.QueueJob (pogofit.queue, "pogofit.UpdateBLWithREV", {"0" : pogofit.file_list[file_index],
                                                                      "1" : pogofit.gtr_fit[terms.global],
                                                                      "2" : (pogofit.gtr_fit[terms.branch_length])[file_index]},
                                                                    "pogofit.handle_branch_length_callback");
    }
    mpi.QueueComplete (pogofit.queue);

    pogofit.run_gtr_iteration_branch_lengths.logL = math.Sum (utility.Map (utility.Filter (pogofit.analysis_results, "_value_", "_value_/pogofit.phase_key"), "_value_", "(_value_[pogofit.phase_key])[terms.fit.log_likelihood]"));

    io.ReportProgressMessageMD ("Protein GTR Fitter", pogofit.phase_key,
                            "Overall Log(L) = " + pogofit.run_gtr_iteration_branch_lengths.logL);
}

/**
 * @name pogofit.handle_gtr_callback
 * @description Handle MPI callback after fitting a REV model 
 */
function pogofit.handle_branch_length_callback (node, result, arguments) {


    savekey = pogofit.index_to_filename[arguments[0]];


    utility.EnsureKey(pogofit.analysis_results, pogofit.phase3);
    utility.EnsureKey(pogofit.analysis_results[savekey], pogofit.phase3);
    (pogofit.analysis_results[pogofit.phase3])[savekey] = result;

    io.ReportProgressMessageMD ("Protein GTR Fitter", "* " + ^"pogofit.phase3",
                                "Received file '" + arguments[0] + "' from node " + node + ". LogL = " + result[terms.fit.log_likelihood]);


}
/**
 * @name pogofit.UpdateBLWithREV
 * @description Use a previously-fitted average REV amino acid model to a file, specifically for branch length optimization under this model
 * @param {String} filename - the filename of the dataset to be fitted
 * @param {Dict} rates - the rates for the GTR model used in fitting
 * @param {Dict} branch_lengths - the current branch length values for this dataset
 * @return the fitted MLE
 */
function pogofit.UpdateBLWithREV (filename, rates, branch_lengths) {


    pogofit.file_info = alignments.ReadNucleotideDataSet ("pogofit.msa",
                                                              filename);
    pogofit.name_mapping = pogofit.file_info[utility.getGlobalValue("terms.data.name_mapping")];
    if (None == pogofit.name_mapping) { /** create a 1-1 mapping if nothing was done */
        pogofit.name_mapping = {};
        utility.ForEach (alignments.GetSequenceNames ("pogofit.msa"), "_value_", "`&pogofit.name_mapping`[_value_] = _value_");
    }

    utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");

    ExecuteCommands ('pogofit.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (pogofit.file_info [utility.getGlobalValue("terms.data.partitions")], pogofit.name_mapping)',
                     {"0" : "Y"});

    utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);
    pogofit.filter_specification = alignments.DefineFiltersForPartitions (pogofit.partitions_and_trees,
                                                                            "pogofit.msa" ,
                                                                            "pogofit.filter.",
                                                                            pogofit.file_info);

    pogofit.rev_file_mle = {terms.global : {}};

    for (l1 = 0; l1 < 20; l1 += 1) {
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
            rate_term = terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]);
            (pogofit.rev_file_mle[terms.global]) [rate_term] =
                {terms.fit.MLE : (rates[rate_term])[terms.fit.MLE] , terms.fix : TRUE};
        }
    }

    pogofit.rev_file_mle [terms.branch_length] = { "0" : branch_lengths };


    utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);
    pogofit.rev_file_mle = estimators.FitSingleModel_Ext (
                                        utility.Map (pogofit.filter_specification, "_value_", "_value_[terms.data.name]"), // value => value['name']
                                        utility.Map (pogofit.partitions_and_trees, "_value_", "_value_[terms.data.tree]"), // value => value['tree']
                                        pogofit.rev_model_gamma, 
                                        pogofit.rev_file_mle,
                                        None
                                   );
    utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    console.log (""); // clear past the optimization progress line

    pogofit.rev_file_mle - terms.global; // delete redundant keys

    return pogofit.rev_file_mle;

}


