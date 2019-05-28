RequireVersion("2.3.9");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("libv3/models/rate_variation.bf");
LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein/REV.bf");
LoadFunctionLibrary("libv3/models/protein.bf");
LoadFunctionLibrary("pogofit_helper_fixgamma.bf"); // Functions, model definitions used for this batchfile.


/*------------------------------------------------------------------------------*/

//utility.ToggleEnvVariable ("OPTIMIZATION_TIME_HARD_LIMIT", 15);

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

pogofit.analysis_banner = {
    terms.io.info: "PogoFit, *P*rotein *G*TR *Fit*ter: Fit a general time reversible (GTR) model to a collection of training protein sequence alignments.",
    terms.io.version: "0.01",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "spond@temple.edu; spielman@rowan.edu",
    terms.io.requirements: "All alignments must be in HyPhy-format: Each file must contain a protein multiple sequence alignment and newick phylogeny. NEXUS input is not accepted."
};
io.DisplayAnalysisBanner(pogofit.analysis_banner);

pogofit.baseline_phase   = "Baseline Fit";
pogofit.final_phase      = "GTR Fit";

pogofit.options.imputation = "Impute zero rates";
pogofit.options.dataset_information = "Dataset information";
pogofit.options.number_of_datasets = "Number of datasets";
pogofit.options.frequency_type = "frequency estimation";
pogofit.options.baseline_model   = "baseline model";


pogofit.ml_freq                = "ML";
pogofit.emp_freq               = "Emp";

pogofit.single   = "Single";
pogofit.multiple = "Multiple";

pogofit.impute = "Yes";
pogofit.no_impute = "No";


pogofit.output_hyphy    = "HyPhy";
pogofit.output_paml     = "PAML";
pogofit.output_raxml    = "RAxML";
pogofit.output_all      = "All";
pogofit.hyphy_model_ext = ".POGOFIT.fitted_model";
pogofit.paml_model_ext  = ".POGOFIT.paml";
pogofit.raxml_model_ext = ".POGOFIT.raxml";

pogofit.analysis_results = {terms.json.analysis: pogofit.analysis_banner,
                                terms.json.input: {},
                                terms.json.timers: {}};
pogofit.timers = {};

/********************************************** MENU PROMPTS ********************************************************/
/********************************************************************************************************************/

KeywordArgument ("mode",        "Single or Multiple files ", pogofit.multiple);  


// Prompt for number of files to analyze, and read file list accordingly //
pogofit.one_or_many  = io.SelectAnOption ({{pogofit.multiple, "Infer a protein model from multiple training datasets (this is more common)."}, 
                                               {pogofit.single, "Infer a protein model from a single training datasets."}}, 
                                                "How many datasets will be used to fit the protein model?");                                    
if (pogofit.one_or_many == pogofit.single) {
    KeywordArgument ("alignment",        "The alignment to analyze");  
    pogofit.input_file = io.PromptUserForString ("Provide the filename of the alignment to analyze");
    pogofit.file_list           = {{pogofit.input_file}};
} else {
    if (pogofit.one_or_many == pogofit.multiple) {
        KeywordArgument ("list",        "The list of alignments (one per line) to analyze");  
        SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
        fscanf (PROMPT_FOR_FILE, "Lines", pogofit.file_list);
        pogofit.input_file  = utility.getGlobalValue("LAST_FILE_PATH");
        // convert paths to globals relative to the input file
        pogofit.input_file_dir = (io.splitFilePath (   pogofit.input_file)) [0];   
        pogofit.file_list = utility.Map (pogofit.file_list, "_file_name_", 'pogofit.input_file_dir + _file_name_');
    }
}


pogofit.output_model_prefix = pogofit.input_file;
pogofit.json_file           = pogofit.input_file  + ".POGOFIT.json";
pogofit.file_list           = io.validate_a_list_of_files (pogofit.file_list);
pogofit.file_list_count     = Abs (pogofit.file_list);
pogofit.index_to_filename   = utility.SwapKeysAndValues(pogofit.file_list);


KeywordArgument ("baseline-model",        "The empirical protein model to use for optimizing branch lengths", "WAG");  

pogofit.baseline_model  = io.SelectAnOption (models.protein.empirical_models,
                                                "Select an empirical protein model to use for optimizing the provided branch lengths:");

// Prompt for F inference //

KeywordArgument ("frequencies",      "Equilibrium frequency estimator", pogofit.emp_freq);  

pogofit.frequency_type  = io.SelectAnOption ({{pogofit.emp_freq, "+F (Empirical)"}, {pogofit.ml_freq, "Maximum likelihood"}},
                                                "Select an frequency specification:");
 
KeywordArgument ("output-format",      "Output format for the fitted model", pogofit.output_all);  
                     
// Prompt for output format //
pogofit.output_format  = io.SelectAnOption ({
                                                  {pogofit.output_hyphy, "HyPhy-formatted model (extension `.fitted_model`)"},
                                                  {pogofit.output_paml, "PAML-formatted model (extension `.paml`)"},
                                                  {pogofit.output_raxml, "RAXML-formatted model (extension `.raxml`)"},
                                                  {pogofit.output_all, "Output all file formats"}},
                                                 "Select an output format for the fitted model:");

// Prompt for zero-rate imputation //

KeywordArgument ("zero-rates",      "Should zero rates be imputed or left at 0", pogofit.impute);  

pogofit.imputation  = io.SelectAnOption ({{pogofit.impute, "Impute zero rates as in Nickle et al. 2007 (Recommended)"},
                                          {pogofit.no_impute, "Leave zero rates at zero"}},
                                           "Impute zero rates for final model files (*excluding* JSON)?:");

pogofit.use_rate_variation = "Yes"; 

KeywordArgument ("precision",        "Optimization precision", 0.001);  

pogofit.precision = io.PromptUser ("Optimization precision", 0.001, 1e-5, 1, FALSE);

pogofit.save_options();

/*
pogofit.baseline_model_name = pogofit.baseline_model + "+F, with 4 category Gamma rates";
pogofit.baseline_model_desc = "pogofit.Baseline.ModelDescription.withGamma";
pogofit.initial_rates       = Eval("models.protein." + pogofit.baseline_model + ".Rij");

if (pogofit.frequency_type == pogofit.emp_freq){
    pogofit.rev_model = "models.protein.REV.ModelDescription";
}
if (pogofit.frequency_type == pogofit.ml_freq){
    pogofit.rev_model = "models.protein.REVML.ModelDescription";
}
*/

pogofit.baseline_model_name = pogofit.baseline_model + "+F, with 4 category Gamma rates";
pogofit.baseline_model_desc_nogamma = "pogofit.Baseline.ModelDescription";
pogofit.baseline_model_desc_gamma = "pogofit.Baseline.ModelDescription.withGamma";

pogofit.initial_rates       = Eval("models.protein." + pogofit.baseline_model + ".Rij");
pogofit.rev_model           = "models.protein.REV.ModelDescription.withGamma";


/********************************************************************************************************************/
/********************************************* ANALYSIS BEGINS HERE *************************************************/
/********************************************************************************************************************/


pogofit.startTimer (pogofit.timers, "Total time");


pogofit.queue = mpi.CreateQueue ({  utility.getGlobalValue("terms.mpi.Headers")   : utility.GetListOfLoadedModules ("libv3/") ,
                                        utility.getGlobalValue("terms.mpi.Functions") :
                                        {
                                            {"models.protein.REV.ModelDescription.withGamma",
                                             "models.protein.REV.ModelDescription.withGDD4",
                                             "pogofit.REV.ModelDescription",
                                             "pogofit.REV.ModelDescription.withGamma",
                                             "pogofit.REV.ModelDescription.freqs",
                                             "pogofit.Baseline.ModelDescription.withGamma",
                                             "pogofit.Baseline.ModelDescription",
                                             "pogofit.fitBaselineToFile"
                                            }
                                        },
                                        utility.getGlobalValue("terms.mpi.Variables") : {{
                                            "pogofit.baseline_model_desc_gamma",
                                            "pogofit.baseline_model_desc_nogamma",
                                            "pogofit.rev_model",
                                            "pogofit.baseline_model",
                                            "pogofit.index_to_filename",
                                            "pogofit.analysis_results",
                                            "pogofit.baseline_phase",
                                            "pogofit.shared_EFV"
                                        }}
                                     });
  


/******************************************* STEP ONE *******************************************************
        Perform an initial fit of the Baseline model+F+4G for full data. Estimates shared alpha as well.
*************************************************************************************************************/
console.log("\n\n[PHASE 1] Performing initial branch length optimization using " + pogofit.baseline_model);

pogofit.startTimer (pogofit.timers, pogofit.baseline_phase);

pogofit.baseline_fit = pogofit.fitBaselineTogether();

pogofit.stopTimer (pogofit.timers, pogofit.baseline_phase);
/*************************************************************************************************************/
/*************************************************************************************************************/



/******************************************* STEP TWO *******************************************************
        Fit a full GTR model to all dataset(s) jointly, using the baseline model as initial rates
*************************************************************************************************************/
console.log("\n\n[PHASE 2] Optimizing protein model");

pogofit.startTimer (pogofit.timers, pogofit.final_phase);

//pogofit.baseline_fit = utility.Map (utility.Filter (pogofit.analysis_results, "_value_", "_value_/'" + pogofit.baseline_phase + "'"), "_value_", "_value_['" + pogofit.baseline_phase + "']");

pogofit.gtr_fit = pogofit.fitGTR_fixalpha(pogofit.baseline_fit);
                                                                                                                              
pogofit.stopTimer (pogofit.timers, pogofit.final_phase);
/*************************************************************************************************************/
/*************************************************************************************************************/




/*********************** Save custom model to file(s) as specified **************************/
pogofit.final_efv = pogofit.extract_efv(); // fitted frequencies

if (pogofit.imputation == pogofit.impute) {
    // Tree length for each alignment
    pogofit.tree_lengths = {};
    utility.ForEachPair (pogofit.gtr_fit[terms.branch_length], "_part_", "_value_",
    '
        pogofit.tree_lengths[_part_] = math.Sum(utility.Map (_value_, "_data_",
        "
            _data_ [terms.fit.MLE]
        " 
        ))
    '
    );


    // Site counts for each alignment
    pogofit.site_counts = {};
    utility.ForEachPair( (pogofit.analysis_results[terms.json.input])[pogofit.options.dataset_information], "_key_", "_value_",
    '
        pogofit.site_counts[_key_] = _value_[terms.json.sites];
    '
    );
    pogofit.final_rij = pogofit.extract_rates_imputation();
}
else  {
    pogofit.final_rij = pogofit.extract_rates();
}

io.ReportProgressMessageMD ("Protein GTR Fitter", "Confidence intervals", "Computing confidence intervals for individual rate parameters");

namespace pogofit {
    final_ci    = {};
    lf_id = gtr_fit[^"terms.likelihood_function"];
    for (l1 = 0; l1 < 20; l1 += 1) {
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
            rate_term = terms.aminoacidRate ((^"models.protein.alphabet")[l1],(^"models.protein.alphabet")[l2]);
            parameter_name = ((gtr_fit[^"terms.global"])[rate_term])[^"terms.id"];
            if (utility.Has ((gtr_fit[^"terms.global"])[rate_term], ^"terms.constraint", "String")) {
               console.log ((^"models.protein.alphabet")[l1] + "--" + (^"models.protein.alphabet")[l2] + ": Constrained");
            } else {
                ci = parameters.GetProfileCI (parameter_name, lf_id, 0.95);
                
                console.log ((^"models.protein.alphabet")[l1] + "--" + (^"models.protein.alphabet")[l2] + ": " + Format (ci [^("terms.fit.MLE")], 8, 4) + " [" + Format (ci [^("terms.lower_bound")],8,4) + ", " + Format (ci [^("terms.upper_bound")],8,4) + "]");
                
                ci["From"] = (^"models.protein.alphabet")[l1];
                ci["To"] = (^"models.protein.alphabet")[l2];
                final_ci [rate_term] = ci;
            }
        }
    }
}


console.log("\n\n Saving results");

pogofit.write_model_to_file();

/************************************* Save analysis JSON ***********************************/
pogofit.stopTimer (pogofit.timers, "Total time");
pogofit.analysis_results[terms.json.timers] = pogofit.timers;
io.SpoolJSON(pogofit.analysis_results, pogofit.json_file);

console.log("\n\nAnalysis complete!");

