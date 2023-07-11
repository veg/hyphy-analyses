LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/models/parameters.bf");

KeywordArgument  ("fit", "Load the RELAX fit file");
sim.filepath = io.LoadFile("Load the RELAX fit file");

KeywordArgument  ("skip-tree", "Do not include the generating tree in replicate files","False");
sim.skip_tree = io.SelectAnOption ({
                                        {"False", "Include the tree"}
                                        {"True", "Skip the tree"}
                                  }, "Do not include the generating tree in replicate files");


console.log (">Loaded fit file from `sim.filepath`");

sim.lf = (utility.GetListOfLoadedLikelihoodFunctions (None));
// list of likelihood function IDs

assert (utility.Array1D(sim.lf) == 1, "The fit file must contain a single likelihood function");

sim.lf  = sim.lf [0];  
// the ID of the single loaded likelihood function

GetString (sim.parameters, ^sim.lf,-1);


sim.globals = sim.parameters [terms.parameters.global_independent];
sim.locals = sim.parameters [terms.parameters.local_independent];
console.log (">The likelihood function `sim.lf` has `utility.Array1D(sim.globals)` global parameters and `utility.Array1D(sim.locals)` local parameters");



//console.log ("\nGlobal parameters and their MLEs (ranges)");

for (sim.param; in; sim.globals) {
    sim.param.MLE = Eval (sim.param);
    sim.param.range = parameters.GetRange(sim.param);
    //console.log ("\t" + sim.param + " = " + sim.param.MLE + " [" + sim.param.range [terms.lower_bound] + ", " +sim.param.range [terms.upper_bound] + "]");
    
    ExecuteCommands ('
        KeywordArgument  ("`sim.param`", "Value for `sim_param` to use during simulations", `sim.param.MLE`);
        sim.use_this_value = io.PromptUser ("Value for `sim.param` to use during simulations",
                                            `sim.param.MLE`,
                                            `sim.param.range [terms.lower_bound]`,
                                            `sim.param.range [terms.upper_bound]`,
                                            0);
    ');
    
    ^sim.param = sim.use_this_value;
}

if (utility.Array1D(sim.locals) > 0) {
     ExecuteCommands ('
        KeywordArgument  ("local-scaler", "Scale local parameters", 1.);
        sim.use_this_value = io.PromptUser ("Scale local parameters (multiplicative)",
                                            1,
                                            0,
                                            10000,
                                            0);
    ');
    
    for (sim.param; in; sim.locals) {
        ^sim.param = ^sim.param * sim.use_this_value;
    }
}

KeywordArgument  ("replicates", "How many replicates to generate", 10);
sim.replicates = io.PromptUser ("How many replicates to generate", 10, 1, 100000, TRUE);

console.log ("Will generate `sim.replicates` replicates");

KeywordArgument ("output", "Write the function used for data simulation and replicates to this path");
sim.path = io.PromptUserForFilePath ( "Write the function used for data simulation and replicates to this path");

Export (sim.lf_export, ^sim.lf);

fprintf (sim.path, CLEAR_FILE, sim.lf_export);

function not_internal (name, sequence) {
    return None == regexp.Find (name, "^(Node|NODE)");
}

function is_internal (name, sequence) {
    return None != regexp.Find (name, "^(Node|NODE)");
}

if (sim.skip_tree == "True") {
    utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA", FALSE);
} else {
    utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA", TRUE);
}


for (i = 1; i <= sim.replicates; i+=1) {
    io.ReportProgressBar("", "Generating replicate " + i + "/" +  sim.replicates);
    DataSet sim.simulated = SimulateDataSet (^sim.lf, "", sim.cat_states, sim.cat_names, sim.path + ".replicate." + i + ".ancestral");
    DataSetFilter sim.simulated.filter = CreateFilter (sim.simulated,1,"","not_internal");
    DataSetFilter sim.simulated.filter.ancestor = CreateFilter (sim.simulated,1,"","is_internal");
    
    if (sim.skip_tree != "True") {
       utility.SetEnvVariable ("DATAFILE_TREE",  Format (^((sim.parameters["Trees"])[0]), 1,1));
    } 
    
    fprintf (sim.path + "." + i + ".nex", CLEAR_FILE, sim.simulated.filter);
    fprintf (sim.path + "." + i + ".internal-nex", CLEAR_FILE, sim.simulated.filter.ancestor);
}
io.ClearProgressBar();

