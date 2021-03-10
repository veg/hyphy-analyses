RequireVersion ("2.5.28");

LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");
LoadFunctionLibrary     ("libv3/convenience/random.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/models/parameters.bf");
LoadFunctionLibrary     ("libv3/tasks/estimators.bf");
LoadFunctionLibrary     ("SelectionAnalyses/modules/io_functions.ibf");



filter.analysis_description = {terms.io.info :
                            "
                            Read two fit files from BUSTED analyses and conduct an LRT test to determine if the inferred dN/dS rate distributions are different between them
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "Two BUSTED fit files"
                          };


io.DisplayAnalysisBanner (filter.analysis_description);


KeywordArgument     ("fit1", "First fit file");
namespace first {
    SetDialogPrompt ("Load second fit File");
    ExecuteAFile (PROMPT_FOR_FILE);
    path = ^"LAST_FILE_PATH";

}

KeywordArgument     ("fit2", "Second fit file");

SetDialogPrompt ("Load second fit File");
ExecuteAFile (PROMPT_FOR_FILE);
second.path = ^"LAST_FILE_PATH";

KeywordArgument ("output", "Write comparison JSON to");
json.path = io.PromptUserForFilePath ("Save the resulting JSON file to");

KeywordArgument ("grid-size", "The number of points in the initial distributional guess for likelihood fitting", 250);
initial_grid.N = io.PromptUser ("The number of points in the initial distributional guess for likelihood fitting", 250, 1, 10000, TRUE);

lfSpec1 = _extractLFInfo (0);
lfSpec2 = _extractLFInfo (1);
distribution1 = _extractRateDistribution (lfSpec1);
distribution2 = _extractRateDistribution (lfSpec2);

assert (Abs (distribution1["rates"])==Abs (distribution2["rates"]), "Unequal number of dN/dS rate classes in the models");

io.ReportProgressMessageMD("BUSTED", "load", "Loaded distributions report");

io.ReportProgressMessageMD("BUSTED", "load", "\n>dN/dS distribution for `first.path`");
selection.io.report_dnds (parameters.GetStickBreakingDistribution (distribution1) % 0);

io.ReportProgressMessageMD("BUSTED", "load", "\n>dN/dS distribution for `second.path`");
selection.io.report_dnds (parameters.GetStickBreakingDistribution (distribution2) % 0);


likelihoodFunctionComponents = {};

for (i, n; in; lfSpec1["Datafilters"] ) {
    likelihoodFunctionComponents + n;
    likelihoodFunctionComponents + (lfSpec1["Trees"])[i];
}

for (i, n; in; lfSpec2["Datafilters"] ) {
    likelihoodFunctionComponents + n;
    likelihoodFunctionComponents + (lfSpec2["Trees"])[i];
}

ExecuteCommands ('
LikelihoodFunction composite = (' + Join (",",likelihoodFunctionComponents) + ');
');


LFCompute (composite, LF_START_COMPUTE);
LFCompute (composite, independentLL);
LFCompute (composite, LF_DONE_COMPUTE);

json = {
    'files' : {{first.path, second.path}},
    'distributions' : {
        "0" : parameters.GetStickBreakingDistribution (distribution1) % 0, 
        "1" : parameters.GetStickBreakingDistribution (distribution2) % 0
    },
    'fits' : {
        'unconstrained' : independentLL
    }
};

io.ReportProgressMessageMD("BUSTED", "optimize", "Fitting the constrained model");
io.ReportProgressMessageMD("BUSTED", "optimize", "\n- Independent model likelihood " + Format (independentLL, 10, 3));

utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
df = 0;

initial_grid         = estimators.LHC (_init_grid_setup (distribution1,distribution2),initial_grid.N$2);



for (p = 0; p < initial_grid.N$2; p+=1) {
    if (Random (0,1) < 0.5) {
        _d1 = distribution1;
    } else {
        _d1 = distribution2;
    }
    point = {};
    for (i, d; in; distribution1["rates"]) {
        _r1 = Eval ((_d1["rates"])[i]);        
        if (Random (0,1) < 0.5) {
            _var = _r1 * 0.1;
        } else {
            _var = 0;
        }
        point[d] = {terms.id: d, terms.fit.MLE : _r1 + random.normal.standard () * _var};
    }
    for (i, d; in; distribution1["weights"]) {
        _r1 = Eval ((_d1["weights"])[i]);        
        if (Random (0,1) < 0.5) {
            _var = _r1 * 0.1;
        } else {
            _var = 0;
        }
        point[d] = {terms.id: d, terms.fit.MLE : _r1 + random.normal.standard () * _var};
    }
    initial_grid + point;
    
}

extra_grid_points = {
    "1" : {},
    "2" : {},
    "3" : {}
};

for (i, d; in; distribution1["rates"]) {
    (extra_grid_points [1])[d] = {terms.id: (distribution1["rates"])[i], terms.fit.MLE : Eval ((distribution1["rates"])[i])};
    (extra_grid_points [2])[d] = {terms.id: (distribution1["rates"])[i], terms.fit.MLE : Eval ((distribution2["rates"])[i])};
    (extra_grid_points [3])[d]= {terms.id: (distribution1["rates"])[i], terms.fit.MLE : 0.5 *(Eval ((distribution1["rates"])[i]) + Eval ((distribution2["rates"])[i]))};
    parameters.SetConstraint ((distribution2["rates"])[i], d, "");
    df+=1;
}

for (i, d; in; distribution1["weights"]) {
    (extra_grid_points [1])[d] = {terms.id: (distribution1["weights"])[i], terms.fit.MLE : Eval ((distribution1["weights"])[i])};
    (extra_grid_points [2])[d] = {terms.id: (distribution1["weights"])[i], terms.fit.MLE : Eval ((distribution2["weights"])[i])};
    (extra_grid_points [3])[d]= {terms.id: (distribution1["weights"])[i], terms.fit.MLE : 0.5 *(Eval ((distribution1["weights"])[i]) + Eval ((distribution2["weights"])[i]))};
    parameters.SetConstraint ((distribution2["weights"])[i], d, "");
    df+=1;
}

initial_grid + extra_grid_points[1];
initial_grid + extra_grid_points[2];
initial_grid + extra_grid_points[3];


io.ReportProgressMessageMD("BUSTED", "optimize", "- Degrees of freedom = " + df);
io.ReportProgressMessageMD("BUSTED", "optimize", "- Computing the likelihood function on `Abs(initial_grid)` initial points");

grid_results = mpi.ComputeOnGrid (&composite, initial_grid, "mpi.ComputeOnGrid.SimpleEvaluator", "mpi.ComputeOnGrid.ResultHandler");

grid_min = Min (grid_results,1);
grid_max = Max (grid_results,1);

io.ReportProgressMessageMD("BUSTED", "optimize", "- log(L) ranges from " + Format (grid_min[terms.data.value], 8, 3) + " to " + grid_max[terms.data.value]);

parameters.SetValues (initial_grid[grid_max["key"]]);

io.ReportProgressMessageMD("BUSTED", "optimize", "\n>dN/dS distribution used to initialize the joint optimization");
selection.io.report_dnds (parameters.GetStickBreakingDistribution (distribution2) % 0);


Optimize (null_MLE, composite);

io.ReportProgressMessageMD("BUSTED", "results", "Model fit results");

lrt = math.DoLRT (null_MLE[1][0], independentLL, df);

(json['fits'])['constrained'] = null_MLE[1][0];
(json['distributions'])['shared'] =parameters.GetStickBreakingDistribution (distribution2) % 0;

json['test']  = lrt;

io.ReportProgressMessageMD("BUSTED", "results", "\n- Constrained model likelihood " + Format (null_MLE[1][0], 10, 3));
io.ReportProgressMessageMD("BUSTED", "results", "- p-value " + lrt[terms.p_value]);

io.ReportProgressMessageMD("BUSTED", "results", "\nInferred shared dN/dS distribution");
selection.io.report_dnds ((json['distributions'])['shared']);

io.SpoolJSON (json, json.path);


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
        "rates" : "test.omega[0-9]+$",
        "weights" : "test.bsrel_mixture_aux_[0-9]+$",
        
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
