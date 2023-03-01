RequireVersion ("2.5.47");

LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");
LoadFunctionLibrary     ("libv3/convenience/random.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/models/parameters.bf");
LoadFunctionLibrary     ("libv3/models/rate_variation.bf");
LoadFunctionLibrary     ("libv3/tasks/estimators.bf");
LoadFunctionLibrary     ("SelectionAnalyses/modules/io_functions.ibf");

utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);



namespace terms.relax_joint {
    mapping = "mapping";
    files = "files";
};


relax_joint.analysisDescription = {terms.io.info : "Load a collection of RELAX models fitted to individual alignments, and perform a joint test of relaxation/intensification, by using a single 'K' parameter for all genes, or a mixture of several K parameters (a random effects model)",
                           terms.io.version : "0.1.0",
                           terms.io.reference : "TBA",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "A collection of RELAX alternative model fit files"
                          };


io.DisplayAnalysisBanner (relax_joint.analysisDescription);

            
/* 1b. User Input and data load
------------------------------------------------------------------------------*/

KeywordArgument ("filelist","List of files to include in this analysis");
relax_joint.file_list = io.get_a_list_of_files(io.PromptUserForFilePathRead ("List of files to include in this analysis"));

KeywordArgument ("output", "Write comparison JSON to");
relax_joint.output_path = io.PromptUserForFilePath ("Save the resulting JSON file to");

KeywordArgument ("components", "How many random effect components (1-4, '1' for shared K)? ", 1);
relax_joint.components = io.PromptUser ("How many random effect components ('1' for shared K)", 1, 1, 4, TRUE);

relax_joint.file_count = utility.Array1D (relax_joint.file_list);
io.CheckAssertion("relax_joint.file_count >= 1", "A non-empty file list is required");

io.ReportProgressMessageMD("relax_joint", "data" , "* Loaded a list with **" + relax_joint.file_count  + "** files");

relax_joint.path_ordering = {};
relax_joint.trees = {};
relax_joint.file_prefix = {};
relax_joint.order2path = {};
relax_joint.lf_info = {};
relax_joint.distributions = {};
relax_joint.distribution_estimates_ind = {};
relax_joint.K = {};
relax_joint.K_estimates = {};

relax_joint.likelihoodFunctionComponents = {};


for (relax_joint.counter = 0; relax_joint.counter < relax_joint.file_count; relax_joint.counter += 1  ) {
     relax_joint.path = relax_joint.file_list[relax_joint.counter ];
     io.ReportProgressMessageMD("relax_joint", "data" , "* Loading file \`" + relax_joint.path + "\`");
     relax_joint.namespace = "relax_joint_" + relax_joint.counter;
     ExecuteCommands ('
        namespace `relax_joint.namespace` {
            ExecuteAFile (^"relax_joint.path");
            lfSpec  = _extractLFInfo (^"relax_joint.counter");
            (^"relax_joint.lf_info")[^"relax_joint.path"] = lfSpec;
            for (i, n; in; lfSpec["Datafilters"] ) {
                (^"relax_joint.likelihoodFunctionComponents") + n;
                (^"relax_joint.likelihoodFunctionComponents") + (lfSpec["Trees"])[i];
            }
            (^"relax_joint.distributions")[^"relax_joint.path"] =  _extractRateDistribution (lfSpec);
            (^"relax_joint.K")[^"relax_joint.path"] =  _extractKParameter (lfSpec);
            (^"relax_joint.K_estimates")[^"relax_joint.counter"] =  Eval (((^"relax_joint.K")[^"relax_joint.path"]));
            (^"relax_joint.distribution_estimates_ind")[^"relax_joint.counter"] = parameters.GetStickBreakingDistribution ( (^"relax_joint.distributions")[^"relax_joint.path"] ) % 0;
            
            io.ReportProgressMessageMD("relax_joint", "data" , "**K** = " + Format ((^"relax_joint.K_estimates")[^"relax_joint.counter"], 6, 2) );
            selection.io.report_dnds ((^"relax_joint.distribution_estimates_ind")[^"relax_joint.counter"] );

            
        }
     ');
     relax_joint.file_prefix  [relax_joint.path] = relax_joint.namespace;
     relax_joint.order2path [Abs (relax_joint.path_ordering)] = relax_joint.path;
     relax_joint.path_ordering [relax_joint.path] = Abs (relax_joint.path_ordering);
}   

ExecuteCommands ('
    LikelihoodFunction relax_joint.LF = (' + Join (",",relax_joint.likelihoodFunctionComponents) + ');
');

LFCompute (relax_joint.LF, LF_START_COMPUTE);
LFCompute (relax_joint.LF, independentLL);
LFCompute (relax_joint.LF, LF_DONE_COMPUTE);

relax_joint.json = {
    terms.json.analysis: relax_joint.analysisDescription,
    'files' : relax_joint.file_list,
    'distributions' : {
        "independent" : relax_joint.distribution_estimates_ind
    },
    'K' : {
        "independent" : relax_joint.K_estimates,
    },
    'fits' : {
        'unconstrained' : independentLL
    }
};

// apply the constraint 


relax_joint.free = (Columns (relax_joint.K))[0];
relax_joint.geo_mean = Eval (relax_joint.free);


terms.relax.k_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "50"
    };
    
for (relax_joint.counter, relax_joint.k_parameter; in; relax_joint.K) {
    if (relax_joint.k_parameter != relax_joint.free) {      
        relax_joint.geo_mean  = relax_joint.geo_mean  * Eval (relax_joint.k_parameter);
        parameters.SetConstraint (relax_joint.k_parameter, relax_joint.free,"");
    } 
    parameters.SetRange(relax_joint.k_parameter, terms.relax.k_range);
}

relax_joint.geo_mean = relax_joint.geo_mean ^ (1/relax_joint.file_count );
^relax_joint.free  = relax_joint.geo_mean;
 

io.ReportProgressMessageMD("relax_joint", "optimize", "Fitting the joint RELAX model");
io.ReportProgressMessageMD("relax_joint", "optimize", "\n- Independent model likelihood " + Format (independentLL, 10, 3));
utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

Optimize (relax_joint.MLE, relax_joint.LF);
io.ReportProgressMessageMD("relax_joint", "optimize", "\n- Joint model likelihood " + Format (relax_joint.MLE[1][0], 10, 3));

relax_joint.ci = parameters.GetProfileCI(relax_joint.free, "relax_joint.LF", 0.95);

io.ReportProgressMessageMD ("relax_joint", "optimize", "\n- Shared K = " + Format (Eval (relax_joint.free), 10, 3) + 
                    " (95% profile CI " + Format ((relax_joint.ci )[terms.lower_bound],8,4) + "-" + Format ((relax_joint.ci )[terms.upper_bound],8,4) + ")");
                    
relax_joint.lrt_joint = math.DoLRT (relax_joint.MLE[1][0], independentLL, relax_joint.file_count-1);
io.ReportProgressMessageMD("relax_joint", "results", "- p-value for file-level K (vs a single K) " + relax_joint.lrt_joint[terms.p_value]);

KeywordArgument ("save-fit", "Save RELAX alternative model fit to this file (default is not to save)", "/dev/null");
relax.save_fit_path = io.PromptUserForFilePath ("Save the joint RELAX model fit to this file ['/dev/null' to skip]");
io.SpoolLFToPath("relax_joint.LF", relax.save_fit_path);

                
relax_joint.distributions_alt = {};

for (relax_joint.counter, relax_joint.path; in; relax_joint.file_list) {
    io.ReportProgressMessageMD("relax_joint", "optimize" , "* File \`" + relax_joint.path + "\`");
    relax_joint.distributions_alt [relax_joint.counter] = parameters.GetStickBreakingDistribution ( relax_joint.distributions[relax_joint.path] ) % 0;
    selection.io.report_dnds (relax_joint.distributions_alt [relax_joint.counter]);
}   


(relax_joint.json["distributions"])["RELAX"] = relax_joint.distributions_alt;
(relax_joint.json["fits"])["RELAX"] = relax_joint.MLE[1][0];
(relax_joint.json["K"])["RELAX"] = Eval (relax_joint.free);
(relax_joint.json["K"])["CI"] = 
        {
            terms.lower_bound : (relax_joint.ci )[terms.lower_bound],
            terms.upper_bound : (relax_joint.ci )[terms.upper_bound]
        };

^relax_joint.free  := 1; 
io.ReportProgressMessageMD("relax_joint", "optimize-null", "Fitting the joint RELAX NULL model");
utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

Optimize (relax_joint.MLE_null, relax_joint.LF);
io.ReportProgressMessageMD("relax_joint", "optimize", "\n- Joint NULL model likelihood " + Format (relax_joint.MLE_null[1][0], 10, 3));

relax_joint.distributions_alt = {};

for (relax_joint.counter, relax_joint.path; in; relax_joint.file_list) {
    io.ReportProgressMessageMD("relax_joint", "optimize" , "* File \`" + relax_joint.path + "\`");
    relax_joint.distributions_alt [relax_joint.counter] = parameters.GetStickBreakingDistribution ( relax_joint.distributions[relax_joint.path] ) % 0;
    selection.io.report_dnds (relax_joint.distributions_alt [relax_joint.counter]);
}   


(relax_joint.json["distributions"])["RELAX-null"] = relax_joint.distributions_alt;
(relax_joint.json["fits"])["RELAX-null"] = relax_joint.MLE_null[1][0];


relax_joint.lrt_relax = math.DoLRT (relax_joint.MLE_null[1][0], relax_joint.MLE[1][0], 1);

relax_joint.json["test"] = relax_joint.lrt_relax ;
relax_joint.json["test-joint"] = relax_joint.lrt_joint ;

io.ReportProgressMessageMD("relax_joint", "results", "- p-value for K!=1 " + relax_joint.lrt_relax[terms.p_value]);

io.SpoolJSON (relax_joint.json, relax_joint.output_path);


// ----------------------------------------------------------------
// HELPER FUNCTIONS
// ----------------------------------------------------------------


lfunction _extractLFInfo (id) {
    GetString (lf, LikelihoodFunction, id);
    GetString (lfInfo, ^lf, -1);
    return lfInfo;
}

lfunction _extractRateDistribution (info) {

   patterns = {
        "rates" : "reference.omega[0-9]+$",
        "weights" : "reference.bsrel_mixture_aux_[0-9]+$",
        
   };
   matches = regexp.PartitionByRegularExpressions (info["Global Independent"], patterns);
   N = Abs (matches[patterns["rates"]]);
   M = Abs (matches[patterns["weights"]]);
   assert (N > 1 && N == M+1, "Missing distribution parameters");
   
   
   
   return {
    "rates"   : utility.sortStrings (Columns (matches[patterns["rates"]])),
    "weights" : utility.sortStrings (Columns (matches[patterns["weights"]]))   
   };
}

lfunction _extractKParameter (info) {

   patterns = {
        "K" : "\\.K$",
        
   };
   matches = regexp.PartitionByRegularExpressions (info["Global Independent"], patterns);
   N = Abs (matches[patterns["K"]]);
   
   assert (N == 1, "Missing K parameter");
   
   return (utility.sortStrings (Columns (matches[patterns["K"]])))[0];
}

//------------------------------------------------------------------------------

function _init_grid_setup (omega_distro,omega_distro2) {

    _initial_ranges = {};
    _N = utility.Array1D (omega_distro[terms.parameters.rates]);
    
    for (_index_, _name_; in; omega_distro[terms.parameters.rates]) {
        _r1 = Eval (_name_);
        _r2 = Eval((omega_distro2[terms.parameters.rates])[_index_]);
        
        if (_index_ < _N - 1) {
            _initial_ranges[_name_] = {
                terms.lower_bound : Min(_r1,_r2)/2,
                terms.upper_bound : Min(1,Max(_r1,_r2)*2)
            };
        } else {
             _initial_ranges[_name_] = {
                terms.lower_bound : Max(1,Min(_r1,_r2)/2),
                terms.upper_bound : Max(_r1,_r2)
            };       
        }
    }

    for (_index_, _name_; in; omega_distro[terms.parameters.weights]) {
        _p1 = Eval (_name_);
        _p2 = Eval((omega_distro2[terms.parameters.weights])[_index_]);
        _initial_ranges[_name_] = {
            terms.lower_bound : Min(_p1,_p2)/2,
            terms.upper_bound : Min(1,Max(_p1,_p2)*2)
        };
    }
    return _initial_ranges;
}
