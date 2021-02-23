RequireVersion ("2.5.28");

LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/models/parameters.bf");
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

for (i, d; in; distribution1["rates"]) {
    parameters.SetConstraint ((distribution2["rates"])[i], d, "");
    df+=1;
}

for (i, d; in; distribution1["weights"]) {
    parameters.SetConstraint ((distribution2["weights"])[i], d, "");
    df+=1;
}

io.ReportProgressMessageMD("BUSTED", "optimize", "- Degrees of freedom = " + df);


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

