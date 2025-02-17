RequireVersion ("2.5.65");


LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);


simulator.analysis_description = {terms.io.info      : "Simulate protein data using available evolutionary models",
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
KeywordArgument ("root-seq",             "Use a specific root sequence to simulate from (overrides --sites)", "None");
KeywordArgument ("extension",            "Use the following file extension for replicates", ".nex");
KeywordArgument ("frequencies",          "The substitution model to use", "empirical");
KeywordArgument ("model",                "The substitution model to use", "empirical");


simulator.seed = +io.PromptUserForString ("Random generator seed (0 to use default initialization)");
if (simulator.seed != 0) {
    SetParameter (RANDOM_SEED, simulator.seed, 0);
}

simulator.code = alignments.LoadGeneticCode  (null);
simulator.tree = trees.LoadAnnotatedTopology (FALSE);

simulator.sites      = io.PromptUser ("The number of residues per alignment", 300, 1, 1e7, TRUE);
simulator.replicates = io.PromptUser ("The number of replicate alignments to generate", 1, 1, 1e7, TRUE);

simulator.root_seq     = io.PromptUserForString ("Use a specific root sequence to simulate from (overrides --sites)");
simulator.extension     = io.PromptUserForString ("Use the following file extension");

if (simulator.root_seq != "None") {
    simulator.sites = Abs (simulator.root_seq);
} else {
    simulator.root_seq = null;
}

simulator.efv        = io.PromptUserForString ("Base frequencies specification");

if (simulator.efv == "empirical") {} else {
    if (simulator.efv  == "equal") {
        simulator.efv = {20,1} ["0.05"];
    } else {
        if (io.FileExists (simulator.efv)) {
            DataSet ds = ReadDataFile (simulator.efv);
            HarvestFrequencies (simulator.efv, ds, 1, 1, 1);
        } else {
             simulator.efv = Eval ("{{" +  simulator.efv + "}}");
             if (utility.Array1D (simulator.efv) == 20) {
                    simulator.efv = {20,1};
                    for (simulator.c = 0; simulator.c < 20; simulator.c += 1) {
                        simulator.efv [simulator.c] =  simulator.efv [simulator.c];
                    }
    
                    simulator.efv = simulator.efv $ Eval(1/(simulator.efv));
    
             } else {
                io.ReportAnExecutionError ("Incorrect dimensions for the base frequency argument (20 comma separated terms)");
             }
        }
    }
}

// TODO check that simulator.efv; no negative or 0 entries (for CF3x4)


console.log (">Protein frequencies used for simulator\n" + simulator.efv);

simulator.module.model = io.PromptUserForString ('Substitution model module');
ExecuteAFile (PATH_TO_CURRENT_BF + "modules/model/" + simulator.module.model);


simulator.model = simulator.define_model ();

KeywordArgument ("branch-variation",     "The model for describing branch-to-branch variation in substitution rates","constant");

simulator.module.branch = io.PromptUserForString ('Branch variation module');

ExecuteAFile (PATH_TO_CURRENT_BF + "modules/branch-variation/" + simulator.module.branch);

model.ApplyModelToTree ("simulator.T", simulator.tree, {"0" : simulator.model}, null);

simulator.has_categories = utility.Array1D ((simulator.model[terms.parameters])[terms.category]) > 0;

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

KeywordArgument ("site-variation",       "The model for describing site-to-site variation in rates", "constant");
simulator.module.site = io.PromptUserForString ('Site variation module');

ExecuteAFile (PATH_TO_CURRENT_BF + "modules/site-variation/"   + simulator.module.site);

simulator.site_profile = simulator.prepare_site_distribution (simulator.model, simulator.sites, "simulator.T", simulator.tree);

simulator.sites_by_profile = {
};

for (simulator.i = 0; simulator.i < simulator.sites; simulator.i += 1) {
    simulator.site_profile_value = "" + simulator.set_site_rate (simulator.model, simulator.i, null);
    if (utility.Has (simulator.sites_by_profile, simulator.site_profile_value, "AssociativeList") == FALSE) {
        simulator.sites_by_profile [simulator.site_profile_value] = {};
    }
    simulator.sites_by_profile [simulator.site_profile_value] + simulator.i;
}

simulator.matrix = {2,20};

simulator.aa = "ACDEFGHIKLMNPQRSTVWY";

for (i = 0; i < 20; i += 1) {
    simulator.matrix[0][i] = simulator.aa[i];
}

simulator.matrix [1][0] = "1";

simulator.root_freqs = simulator.model[terms.efv_estimate];

KeywordArgument ("output",       "Write simulated alignments (as FASTA) to the following prefix path, using the syntax ${path}.replicate.index");
simulator.path = io.PromptUserForFilePath ("Save simulator settings to this path, and replicates to ${path}.replicate.index.extension");


utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT",9);
utility.SetEnvVariable ("DATAFILE_TREE",simulator.tree[terms.trees.newick_annotated]);
utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA",TRUE);

simulator.rate_type   = 0;
simulator.inverse_map = {};

simulator.string_buffer = {};
for (simulator.i = 0; simulator.i < simulator.replicates; simulator.i += 1) {
    simulator.string_buffer[simulator.i] = {};
}

simulator.mode = 1;
simulator.counter = 0;
simulator.BL = {};

for (_rate_distribution_, _site_counts_; in; simulator.sites_by_profile) {
    simulator.apply_site_distribution (simulator.model, _rate_distribution_,  "simulator.T");
    
    simulator.bl_code = simulator.model[terms.model.get_branch_length];
    if (Abs (simulator.bl_code)) {
        simulator.tree_length = 0;
        simulator.tree_root = BranchName (simulator.T, BranchCount (simulator.T));
        for (b; in; simulator.T) {
            if (b != simulator.tree_root) {
                simulator.tree_length += Call (simulator.bl_code ,simulator.model,"simulator.T", b);
            }
        }
    } else {
        simulator.tree_length = +BranchLength (simulator.T, -1);
    }
    
    io.ReportProgressBar ("SIMULATING", "Rate regime " + simulator.mode + " of " + utility.Array1D (simulator.sites_by_profile) + " (branch length = " + Format (simulator.tree_length, 8, 3) + ")");
    simulator.mode   += 1;
    
    

    for (_site_id_; in; utility.DictToArray (_site_counts_)) {
        simulator.inverse_map [_site_id_] = ("" + (simulator.counter));
        simulator.BL[simulator.counter] = simulator.tree_length;
        simulator.counter += 1;
    
    }
    
        
    simulator.site_block = utility.Array1D (_site_counts_);
    
    if (None != simulator.root_seq) {
        simulator.template = {utility.Array1D (_site_counts_), 1};
        simulator.template[0] = "";
        for (k, vl; in; _site_counts_) {
             simulator.template[+k] = simulator.root_seq[vl];
        }
        simulator.start_from_seq_seed = Join ("", simulator.template);

        simulator.start_from_seq = "";  simulator.start_from_seq * (Abs (simulator.start_from_seq_seed) * simulator.replicates);
        for (i = 0; i < simulator.replicates; i+=1) {
            simulator.start_from_seq * simulator.start_from_seq_seed;
        }

        simulator.start_from_seq * 0;
        //console.log (simulator.start_from_seq);
        DataSet simulated_data = Simulate (simulator.T, simulator.root_freqs, simulator.matrix, simulator.start_from_seq);

    } else {
        //console.log (simulator.matrix);
        DataSet simulated_data = Simulate (simulator.T, simulator.root_freqs, simulator.matrix, simulator.site_block*simulator.replicates);
    }
 
    // simulate ALL sites from one scenario here
    if (simulator.rate_type == 0) {
         GetString (simulator.sim_names, simulated_data, -1);
         for (simulator.i = 0; simulator.i < simulator.replicates; simulator.i += 1) {
            simulator.j = 0;
            for (_value_; in; simulator.sim_names) {
                (simulator.string_buffer[simulator.i])+"";
                (simulator.string_buffer[simulator.i])[simulator.j] * (simulator.sites);
                simulator.j += 1;
            }
         }
    }

    for (simulator.i = 0; simulator.i < simulator.replicates; simulator.i += 1) {
    
        DataSetFilter all      = CreateFilter (simulated_data, 1, siteIndex%simulator.replicates == simulator.i);
        for (simulator.j = 0; simulator.j < all.species; simulator.j+=1) {
            GetDataInfo (sim.string, all, simulator.j);
            (simulator.string_buffer[simulator.i])[simulator.j] * sim.string;
        }
    }

    simulator.rate_type += 1;
}

if (Type (simulator.report) == "AssociativeList") {
    fprintf (simulator.path, CLEAR_FILE, {
        terms.model : simulator.model,
        terms.data.tree  : simulator.tree,
        terms.json.branch_lengths : simulator.BL,
        "simulator.site.profile" : simulator.site_profile,
        "simulator.additional_settings" : simulator.report
    });
} else {
    fprintf (simulator.path, CLEAR_FILE, {
        terms.model : simulator.model,
        terms.data.tree  : simulator.tree,
        terms.json.branch_lengths : simulator.BL,
        "simulator.site.profile" : simulator.site_profile
    });
}

io.ClearProgressBar ();


io.ReportStatsMD ("Branch length statistics (per site)", math.GatherDescriptiveStats(Transpose(utility.DictToArray (simulator.BL))));

simulator.inverse_map = Join ("," ,simulator.inverse_map);

for (simulator.i = 0; simulator.i < simulator.replicates; simulator.i += 1) {
    simulator.j = 0;
    utility.ForEach (simulator.sim_names,"_value_", "
                (simulator.string_buffer[simulator.i])[simulator.j] * 0;
                (simulator.string_buffer[simulator.i])[simulator.j] = '>' + _value_ + '\n' + (simulator.string_buffer[simulator.i])[simulator.j];
                simulator.j += 1;
            ");


    simulator.replicate_file = simulator.path + ".replicate." + (1+simulator.i) + simulator.extension;
    fprintf (simulator.replicate_file, CLEAR_FILE, Join ("\n",(simulator.string_buffer[simulator.i])));
    DataSet existing_data  = ReadDataFile (simulator.replicate_file);
    utility.SetEnvVariable ("DATAFILE_TREE",simulator.tree[terms.trees.newick_annotated]);
    utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA",TRUE);
    DataSetFilter all      = CreateFilter (existing_data, 1, simulator.inverse_map);
    fprintf (simulator.replicate_file, CLEAR_FILE, all);
}
