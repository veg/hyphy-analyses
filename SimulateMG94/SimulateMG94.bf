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


simulator.analysis_description = {terms.io.info         : "Simulate codon data using the MG94 model of sequence evolution",
                               terms.io.version      : "0.1",
                               terms.io.authors      : "Sergei L Kosakovsky Pond",
                               terms.io.contact      : "spond@temple.edu",
                               terms.io.requirements : "a phylogenetic tree with branch lengths and other annotations"
                              };

io.DisplayAnalysisBanner (simulator.analysis_description);


KeywordArgument ("seed",                 "Random seed (0 to use default initialization)", "0");
KeywordArgument ("code",                 "Which genetic code should be used", "Universal");
KeywordArgument ("tree",                 "A phylogenetic tree with branch lengths or annotations");
KeywordArgument ("sites",                "How many codon sites to simulate", 500);
KeywordArgument ("replicates",           "How many replicates", 1);
KeywordArgument ("base-frequencies",     "Base frequencies to use. 'equal' or 9 comma-separated values [A in first codon position, C-1, G-1, A-2, C-2, G-2...] or 12 comma-separated values [A in first codon position, C-1, G-1, T-1, A-2, C-2, G-2, T-2...] to specify positional nucleotide frequencies]", "equal");
KeywordArgument ("frequency-estimator",  "Equilibrium frequency estimator", "CF3x4");
KeywordArgument ("AC",                   "The AC substitution rate relative to the AG rate (=1)", "0.5");
KeywordArgument ("AT",                   "The AT substitution rate relative to the AG rate (=1)", "0.5");
KeywordArgument ("CG",                   "The CG substitution rate relative to the AG rate (=1)", "0.5");
KeywordArgument ("CT",                   "The CT substitution rate relative to the AG rate (=1)", "1.0");
KeywordArgument ("GT",                   "The GT substitution rate relative to the AG rate (=1)", "0.5");


simulator.seed = +io.PromptUserForString ("Random generator seed (0 to use default initialization)");
if (simulator.seed != 0) {
    SetParameter (RANDOM_SEED, simulator.seed, 0);
}

simulator.code = alignments.LoadGeneticCode  (null);
simulator.tree = trees.LoadAnnotatedTopology (TRUE);

simulator.sites      = io.PromptUser ("The number of codons per alignment", 300, 1, 1e7, TRUE);
simulator.replicates = io.PromptUser ("The number of replicate alignments to generate", 1, 1, 1e7, TRUE);
simulator.efv        = io.PromptUserForString ("Base frequencies specification");

if (simulator.efv  == "equal") {
    simulator.efv = {4,3} ["0.25"];
} else {
     if (simulator.efv == "HIV") {
        simulator.efv = {
            {0.41, 0.34, 0.41} 
            {0.16, 0.20, 0.12} 
            {0.25, 0.17, 0.14} 
            {0.18, 0.28, 0.33} 
        };
     } else {
         simulator.efv = Eval ("{{" +  simulator.efv + "}}");
         if (utility.Array1D (simulator.efv) == 9) {
                //simulator.t   = {3,3}["simulator.efv[_MATRIX_ELEMENT_COLUMN_*4+_MATRIX_ELEMENT_ROW_]"];
                simulator.efv4 = {4,3};
                for (simulator.c = 0; simulator.c < 3; simulator.c += 1) {
                    for (simulator.r = 0; simulator.r < 3; simulator.r += 1) {
                        simulator.efv4[simulator.r][simulator.c] =  simulator.efv [simulator.c * 3 + simulator.r];
                    }
                    simulator.efv4[3][simulator.c] = 1 - (+simulator.efv4[-1][simulator.c]);
                }
            
                simulator.efv = simulator.efv4 $ Eval({{1/(+simulator.efv4[-1][0]),
                                                        1/(+simulator.efv4[-1][1]),
                                                        1/(+simulator.efv4[-1][2])}});
        
         } else {
             if (utility.Array1D (simulator.efv) == 12) {
                simulator.efv = {4,3}["simulator.efv[_MATRIX_ELEMENT_COLUMN_*4+_MATRIX_ELEMENT_ROW_]"];
                simulator.efv = simulator.efv $ Eval({{1/(+simulator.efv[-1][0]),
                                                  1/(+simulator.efv[-1][1]),
                                                  1/(+simulator.efv[-1][2])}});
             } else {
                io.ReportAnExecutionError ("Incorrect dimensions for the base frequency argument (9 or 12 comma separated terms)");
             }
         }
    }
}

// TODO check that simulator.efv; no negative or 0 entries (for CF3x4)

console.log (">Nucleotide frequencies used for simulator\n" + simulator.efv);

simulator.frequency_type = io.SelectAnOption ({"CF3x4" : terms.frequencies.CF3x4, 
                                               "F3x4" : terms.frequencies.F3x4,
                                               "F1x4" : terms.frequencies.F1x4}, "Equilibrium frequency estimator");


simulator.model = model.generic.DefineModel ("simulator.defineMG.frequencies" , "simulator.MG94", {"0" : parameters.Quote(terms.local), "1" : simulator.code[terms.code]}, null, null);

parameters.SetValue (((simulator.model [terms.parameters])[terms.global])[terms.nucleotideRateReversible("A","C")],io.PromptUser ("Relative AC rate", 0.5, 0, 1000, FALSE));
parameters.SetValue (((simulator.model [terms.parameters])[terms.global])[terms.nucleotideRateReversible("A","T")],io.PromptUser ("Relative AT rate", 0.5, 0, 1000, FALSE));
parameters.SetValue (((simulator.model [terms.parameters])[terms.global])[terms.nucleotideRateReversible("C","G")],io.PromptUser ("Relative CG rate", 0.5, 0, 1000, FALSE));
parameters.SetValue (((simulator.model [terms.parameters])[terms.global])[terms.nucleotideRateReversible("C","T")],io.PromptUser ("Relative CT rate", 1.0, 0, 1000, FALSE));
parameters.SetValue (((simulator.model [terms.parameters])[terms.global])[terms.nucleotideRateReversible("G","T")],io.PromptUser ("Relative GT rate", 0.5, 0, 1000, FALSE));

KeywordArgument ("branch-variation",     "The model for describing branch-to-branch variation in omega ratios","constant");

simulator.module.branch = io.PromptUserForString ('Branch variation module');

ExecuteAFile (PATH_TO_CURRENT_BF + "modules/branch-variation/" + simulator.module.branch);

model.ApplyModelToTree ("simulator.T", simulator.tree, {"0" : simulator.model}, null); 

/** validate the tree **/

simulator.validation_error                     = simulator.validate_tree (simulator.tree);
io.CheckAssertion ("simulator.validation_error == ''", simulator.validation_error);


utility.ForEachPair (simulator.tree [terms.trees.partitioned], "_name_", "_value_", '
    simulator.set_branch_rates (simulator.model, "simulator.T", _name_, 
        {
            terms.trees.model_map : (simulator.tree[terms.trees.model_map])[_name_],
            terms.trees.meta : (simulator.tree[terms.trees.meta])[_name_],
            terms.trees.partitioned : _value_,
            terms.branch_length : (simulator.tree[terms.branch_length ])[_name_]
        });
');

KeywordArgument ("site-variation",       "The model for describing site-to-site variation in relative alpha (dS) and beta (dN) rates", "constant");
simulator.module.site = io.PromptUserForString ('Site variation module');

ExecuteAFile (PATH_TO_CURRENT_BF + "modules/site-variation/"   + simulator.module.site);


simulator.site_profile = simulator.prepare_site_distribution (simulator.model, simulator.sites, "simulator.T", simulator.tree);

console.log (simulator.site_profile);

simulator.sites_by_profile = {
};

for (simulator.i = 0; simulator.i < simulator.sites; simulator.i += 1) {
    simulator.site_profile_value = "" + simulator.set_site_omega (simulator.model, simulator.i, null);
    if (utility.Has (simulator.sites_by_profile, simulator.site_profile_value, "AssociativeList") == FALSE) {
        simulator.sites_by_profile [simulator.site_profile_value] = {};
    }
    simulator.sites_by_profile [simulator.site_profile_value] + simulator.i;
}

simulator.matrix = {2,4};
simulator.matrix [0][0] = "A"; simulator.matrix [0][1] = "C"; simulator.matrix [0][2] = "G"; simulator.matrix [0][3] = "T"; 
simulator.matrix [1][0] = "3"; simulator.matrix [1][1] = simulator.code[terms.code.stops];

simulator.root_freqs = simulator.model[terms.efv_estimate];

KeywordArgument ("output",       "Write simulated alignments (as FASTA) to the following prefix path, using the syntax ${path}.replicate.index");
simulator.path = io.PromptUserForFilePath ("Save simulator settings to this path, and replicates to ${path}.replicate.index");


fprintf (simulator.path, CLEAR_FILE, {
    terms.model : simulator.model,
    terms.data.tree  : simulator.tree,
    "simulator.site.profile" : simulator.site_profile
});


utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT",9);
utility.SetEnvVariable ("DATAFILE_TREE",simulator.tree[terms.trees.newick_annotated]);
utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA",TRUE);

simulator.rate_type   = 0;
simulator.inverse_map = {};

utility.ForEachPair (simulator.sites_by_profile, "_rate_distribution_", "_site_counts_", '
    simulator.apply_site_distribution (simulator.model, _rate_distribution_,  "simulator.T");
    
    utility.ForEach (_site_counts_, "_site_id_", "
        simulator.inverse_map + (\\"\\" + 3*(_site_id_) + \\"-\\" + (3*_site_id_+2));
    ");
    
    for (simulator.i = 0; simulator.i < simulator.replicates; simulator.i += 1) {
        DataSet simulated_data = Simulate (simulator.T, simulator.root_freqs, simulator.matrix, utility.Array1D (_site_counts_));
        if (simulator.rate_type == 0) {
            DataSetFilter all      = CreateFilter (simulated_data, 1);
        } else {
            DataSet existing_data  = ReadDataFile (simulator.path + ".replicate." + (1+simulator.i));
            DataSet combined_data  = Concatenate (existing_data, simulated_data);
            DataSetFilter all      = CreateFilter (combined_data, 1);            
        }
        fprintf (simulator.path + ".replicate." + (1+simulator.i) , CLEAR_FILE, all);
    }
    simulator.rate_type += 1;
');

simulator.inverse_map = Join ("," ,simulator.inverse_map);

for (simulator.i = 0; simulator.i < simulator.replicates; simulator.i += 1) {
    DataSet existing_data  = ReadDataFile (simulator.path + ".replicate." + (1+simulator.i));
    DataSetFilter all      = CreateFilter (existing_data, 1, simulator.inverse_map);   
    fprintf (simulator.path + ".replicate." + (1+simulator.i) , CLEAR_FILE, all);
}
//----------------------------------------------------------------------------------------
  
lfunction simulator.defineMG.frequencies (type,code) {
    model = Call ("models.codon.MG_REV.ModelDescription", type, code);
    if (^"simulator.frequency_type" == "F3x4") {
        model[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F3x4";
    } else {
        if (^"simulator.frequency_type" == "F1x4") {
            model[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F1x4";
        } 
    }
    model [^"terms.efv_estimate"] = ^"simulator.efv";
    return model;
}


