/*
 This batch implements a model test of a 12 rate NREV model versus a standard GTR and strandGTR model. 
 The strand GTR model constrains rates so that AG := TC, not AG := GA as in the GTR model.

 Code by Wayne Delport, inspired by Darren P. Martin

  wdelport@mac.com
  19 February 2009
  
 Refreshed for HyPhy 2.5x by Sergei L Kosakovsky Pond (Nov 9th, 2021)
 
*/

RequireVersion ("2.5.35");

LoadFunctionLibrary ("libv3/all-terms.bf");
LoadFunctionLibrary ("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
LoadFunctionLibrary ("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary ("libv3/convenience/random.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);


/*---------------------Display analysis information-------------------------------*/

nrm.analysis_description = {
    terms.io.info: "
    Perform a fit of GTR, strand non-reversible, and fully non-reversible model (with Gamma rate variation) to a nucleotide alignment.
    Report estimated rate matrices, and perform nested model fits.
    ",
    terms.io.version: "0.1",
    terms.io.reference: "TBD",
    terms.io.authors: "Wayne Delport and Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "a nucleotide alignment and a phylogenetic tree (rooted)"
};
io.DisplayAnalysisBanner(nrm.analysis_description);

nrm.json = {
    terms.json.analysis: nrm.analysis_description,
    terms.json.input: {},
    terms.json.fits: {},
    terms.json.timers: {},
};


/*---------------------Get data and tree-------------------------------*/

KeywordArgument ("alignment",   "Sequence alignment to screen for recombination");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");

namespace nrm {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_nuc_file ("nrm");
}


if (((nrm.partitions_and_trees[0])[terms.data.tree])[terms.trees.rooted] == FALSE) {
    io.ReportWarning ("The provided tree was NOT rooted; this may return incorrect non-reversible model results");
};

/* ---------------------------- Model stuff -------------------------------------------- */

/* user defined functions */
function PopulateNucleotideModelMatrix ( ModelMatrixName, NREVBiasTerms ) {
  
    utility.ExecuteInGlobalNamespace (ModelMatrixName + "=None");
	^ModelMatrixName = { 4, 4 };

	hv = 0;
	for ( h = 0; h < 4; h += 1 ) {
		for ( v = 0; v < 4; v += 1 ) {
			if ( h != v ) {
				modelString = (ModelMatrixName + "["+h+"]["+v+"] := " + NREVBiasTerms [ hv ] + "t;\n");
				//fprintf ( stdout, modelString ); 
				ExecuteCommands ( modelString );
				hv += 1;
			}
		}		
	}
}

lfunction GetEFV ( model ) {
    GetInformation (mi, ^model);
    for (v; in; mi) {
        ^v = 1;
    }
    
    GetString(mi, ^model, -2);
    rm  = mi["RATE_MATRIX"];
    mbf = mi["MULT_BY_FREQ"];
    fm  = mi["EQ_FREQS"];
    
    
    mm = ^rm;
        
    d = Rows (mm);
    for (i = 0; i < d; i+=1) {
        de = 0;
        for (j = 0; j < d; j+=1) {
            if (i != j) {
               if (mbf) {
                 mm[i][j] = mm[i][j] * (^fm)[j];
               }
               de += mm[i][j];
            }
        }
        mm[i][i] = -de;
    }
         
    for (i = 0; i < d; i+=1) {
        mm[i][d-1] = 1;
    }
    
    return (Inverse (mm))[d-1][-1];
    
}


lfunction ReportMatrixAndFreqs ( rates, frequencies ) {
  nucs = "ACGT";
  
  console.log ("\n#### Rate and equilibrium frequency estimates.\n");
  
  fprintf (stdout, "| From/To  |");
  for (i = 0; i < 4; i+=1) {
        fprintf (stdout, "     " ,nucs[i], "    |");
  }
  fprintf (stdout, " Frequency| \n");
  for (i = 0; i < 6; i+=1) {
     fprintf (stdout, "|:--------:");
  }
  fprintf (stdout, "|\n");
  
  res = {4,4};
  
  rc = 0;
  for (i = 0; i < 4; i+=1) {
     fprintf (stdout, "|     " ,nucs[i], "    |");
     for (j = 0; j < 4; j+=1) {
        if (i == j) {
            fprintf (stdout, "     *    |");
        } else {
            if (Abs(rates[rc])) {
                rv = Eval (rates[rc]);
            } else {
                rv = 1;
            }
            res[i][j] = rv;
            rc += 1;
            fprintf (stdout, Format (rv, 9,4), " |");
        }
     }
     fprintf (stdout, Format (frequencies [i], 9, 4), " |\n");
  }
  
  return res;
  
}

/* user defined functions */
lfunction nrm.extractLFInfo ( res, lf_id, xp, samplesize ) {
    
}

/* end of user defined functions */

/* ---------------------------- Model stuff -------------------------------------------- */
/* define a nucleotide bias correction matrix AG := 1 */

/* built like this incase we want to scale to codon models with nuc bias terms */

global alpha = .35; alpha:>0.001;alpha:<100;
category c = (4, EQUAL, MEAN, 
				GammaDist(_x_,alpha,alpha), 
				CGammaDist(_x_,alpha,alpha), 
				0 , 
		  	    1e25,
		  	    CGammaDist(_x_,alpha+1,alpha)
		  	 );
	
NREVBiasTerms = { { "c*AC*", "c*", "c*AT*", "c*CA*", "c*CG*", "c*CT*", "c*GA*", "c*GC*", "c*GT*", "c*TA*", "c*TC*", "c*TG*" } };

ratesArray = { { "AC", "", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG" } };

global AC = 1;
global AT = 1;
global CA = 1;
global CG = 1;
global CT = 1;
global GA = 1;
global GC = 1;
global GT = 1;
global TA = 1;
global TC = 1;
global TG = 1;

nrm.filter_name = (nrm.filter_specification [0])[terms.data.name];
HarvestFrequencies (nrm.dataFrequencies,^nrm.filter_name,1,1,1);
io.ReportProgressMessageMD ("nrm", "gtr", "Fitting the GTR + G model with empirical base frequencies");



/* standard GTR */
CA := AC;
GA := 1;
GC := CG;
TA := AT;
TC := CT;
TG := GT;

PopulateNucleotideModelMatrix ( "GTRMatrix", NREVBiasTerms);
Model GTRModel = ( GTRMatrix, nrm.dataFrequencies, 1 );


utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", FALSE);
ExecuteCommands ("Tree T = " + (nrm.trees[0])[terms.trees.newick_with_lengths]);
LikelihoodFunction lf_gtr = ( ^nrm.filter_name, T );

Optimize ( res_gtr, lf_gtr );
nrm.gtr_fit = estimators.ExtractMLEFromObject ("lf_gtr");

nrm.EFV = GetEFV ("GTRModel");

io.ReportProgressMessageMD ("nrm", "gtr", "\n>" + selection.io.report_fit  (nrm.gtr_fit, 3, nrm.nuc_data_info[terms.data.sample_size]) + "\n" + selection.io.report_fit_secondary_stats (nrm.gtr_fit));
io.ReportProgressMessageMD ("nrm", "gtr", "_Gamma shape parameter_ = " + Format (alpha, 8, 4));

nrm.gtr_rates = ReportMatrixAndFreqs (ratesArray, nrm.EFV);


selection.io.json_store_lf_withEFV(nrm.json, "GTR" ,nrm.gtr_fit [terms.fit.log_likelihood],
                            nrm.gtr_fit [terms.parameters],
                            nrm.nuc_data_info[terms.data.sample_size], 
                            nrm.gtr_rates, 
                            nrm.EFVs,
                            3);
                            

io.ReportProgressMessageMD ("nrm", "stGTR", "Fitting the NREV6 + G model with empirical base frequencies");

parameters.RemoveConstraint (ratesArray);
CA := GT;
GA := CT;
GC := CG;
TA := AT;
TC := 1;
TG := AC;

ExecuteCommands ("Tree T = " + (nrm.trees[0])[terms.trees.newick_with_lengths]);
LikelihoodFunction lf_gtr = ( ^nrm.filter_name, T );
Optimize ( res_gtr, lf_gtr );


nrm.EFV = GetEFV ("GTRModel");
nrm.stgtr_fit = estimators.ExtractMLEFromObject ("lf_gtr");

io.ReportProgressMessageMD ("nrm", "NREV6", "\n>" + selection.io.report_fit  (nrm.stgtr_fit, 3, nrm.nuc_data_info[terms.data.sample_size]) + "\n" + selection.io.report_fit_secondary_stats (nrm.stgtr_fit));
io.ReportProgressMessageMD ("nrm", "NREV6", "_Gamma shape parameter_ = " + Format (alpha, 8, 4));

nrm.stgtr_rates = ReportMatrixAndFreqs (ratesArray, nrm.EFV);

selection.io.json_store_lf_withEFV(nrm.json, "NREV6" ,nrm.stgtr_fit [terms.fit.log_likelihood],
                            nrm.stgtr_fit [terms.parameters],
                            nrm.nuc_data_info[terms.data.sample_size], 
                            nrm.stgtr_rates, 
                            nrm.EFV,
                            3);

io.ReportProgressMessageMD ("nrm", "NREV12", "Fitting the NREV12 + G model with empirical root frequencies");

parameters.RemoveConstraint (ratesArray);
Model GTRModel = ( GTRMatrix, nrm.dataFrequencies, 0 );
ExecuteCommands ("Tree T = " + (nrm.trees[0])[terms.trees.newick_with_lengths]);
LikelihoodFunction lf_gtr = ( ^nrm.filter_name, T );
Optimize ( res_gtr, lf_gtr );


nrm.nrm_fit = estimators.ExtractMLEFromObject ("lf_gtr");

io.ReportProgressMessageMD ("nrm", "NREV12", "\n>" + selection.io.report_fit  (nrm.nrm_fit, 3, nrm.nuc_data_info[terms.data.sample_size]) + "\n" + selection.io.report_fit_secondary_stats (nrm.nrm_fit));
io.ReportProgressMessageMD ("nrm", "NREV12", "_Gamma shape parameter_ = " + Format (alpha, 8, 4));

nrm.EFV = GetEFV ("GTRModel");
nrm.nrm_rates = ReportMatrixAndFreqs (ratesArray, nrm.EFV);

selection.io.json_store_lf_withEFV(nrm.json, "NREV12" ,nrm.nrm_fit [terms.fit.log_likelihood],
                            nrm.nrm_fit [terms.parameters],
                            nrm.nuc_data_info[terms.data.sample_size], 
                            nrm.nrm_rates, 
                            nrm.EFV,
                            3);

io.ReportProgressMessageMD ("nrm", "NREV12F", "Fitting the NREV12F + G model with estimated root frequencies");

global fA = 1/3; fA :< 1;
global fC = 1/3; fC :< 1;
global fG = 1/3; fG :< 1;

nrm.rootFreqs = {{fA}
                 {(1-fA)*fC}, 
                 {(1-fA)*(1-fC)*fG}, 
                 {(1-fA)*(1-fC)*(1-fG)}
                };


ExecuteCommands ("Tree T = " + (nrm.trees[0])[terms.trees.newick_with_lengths]);
LikelihoodFunction3 lf_gtr = ( ^nrm.filter_name, T, nrm.rootFreqs);


Optimize ( res_gtr, lf_gtr );

nrm.nrmf_fit = estimators.ExtractMLEFromObject ("lf_gtr");

io.ReportProgressMessageMD ("nrm", "NREV12F", "\n>" + selection.io.report_fit  (nrm.nrmf_fit, 3, nrm.nuc_data_info[terms.data.sample_size]) + "\n" + selection.io.report_fit_secondary_stats (nrm.nrmf_fit));
io.ReportProgressMessageMD ("nrm", "NREV12F", "_Gamma shape parameter_ = " + Format (alpha, 8, 4));


io.ReportProgressMessageMD ("nrm", "NREV12F", "\n#### Estimates for root frequencies");
io.ReportProgressMessageMD ("nrm", "NREV12F", "- A : " + nrm.rootFreqs[0]);
io.ReportProgressMessageMD ("nrm", "NREV12F", "- C : " + nrm.rootFreqs[1]);
io.ReportProgressMessageMD ("nrm", "NREV12F", "- G : " + nrm.rootFreqs[2]);
io.ReportProgressMessageMD ("nrm", "NREV12F", "- T : " + nrm.rootFreqs[3]);

nrm.EFV = GetEFV ("GTRModel");
nrm.nrmf_rates = ReportMatrixAndFreqs (ratesArray, nrm.EFV);

selection.io.json_store_lf_withEFV(nrm.json, "NREV12F" ,nrm.nrmf_fit [terms.fit.log_likelihood],
                            nrm.nrmf_fit [terms.parameters],
                            nrm.nuc_data_info[terms.data.sample_size], 
                            nrm.nrmf_rates, 
                            nrm.EFV,
                            0);


KeywordArgument ("save-fit", "Save NRM+F model fit to this file (default is not to save)", "/dev/null");
io.SpoolLFToPath("lf_gtr", io.PromptUserForFilePath ("Save NRM+F model fit to this file ['/dev/null' to skip]"));


KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'NRM.json')", nrm.nuc_data_info [terms.json.json]);
nrm.nuc_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");


nrm.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): TRUE,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 10,
        utility.getGlobalValue("terms.table_options.align"): "center",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 12,
            "1": 12,
            "2": 15,
            "3": 15,
            "4": 15,
            "5": 20
        }
    };

io.ReportProgressMessageMD ("nrm", "tests", "Model comparison results");

nrm.report = {{"Null","Alt.","LRT","Deg. freedom","p","Delta c-AIC"}};

nrm.tests = {};
nrm.pairs = {
    "0" : {"Null" : "GTR", "Alt" : "NREV6", "nested" : FALSE},
    "1" : {"Null" : "GTR", "Alt" : "NREV12", "nested" : TRUE},
    "2" : {"Null" : "NREV6", "Alt" : "NREV12", "nested" : TRUE},
    "3" : {"Null" : "NREV12", "Alt" : "NREV12F", "nested" : TRUE}
};

fprintf(stdout, "\n", io.FormatTableRow(nrm.report , nrm.table_output_options));
nrm.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

for (p; in; nrm.pairs ) {  

   nrm.row = {6,1};
   nrm.row[0] = p["Null"];
   nrm.row[1] = p["Alt"];
    
   if (p["nested"]) {
        nrm.df = ((nrm.json[terms.json.fits])[p["Alt"]])[terms.json.parameters] - ((nrm.json[terms.json.fits])[p["Null"]])[terms.json.parameters];
        nrm.row[3] = nrm.df;
        nrm.lrt = math.DoLRT (((nrm.json[terms.json.fits])[p["Null"]])[terms.json.log_likelihood],((nrm.json[terms.json.fits])[p["Alt"]])[terms.json.log_likelihood], nrm.df);
        nrm.row[2] = Format (nrm.lrt[terms.LRT], 10, 4);
        nrm.row[4] = Format (nrm.lrt[terms.p_value], 10, 4);
        nrm.tests [p["Null"] + " vs "+ p["Alt"]] = nrm.lrt;
    } else {
        nrm.row[2] = null;
        nrm.row[3] = null;
        nrm.row[4] = null;
   }
   nrm.row[5] = Format (((nrm.json[terms.json.fits])[p["Null"]])[terms.json.AICc]-((nrm.json[terms.json.fits])[p["Alt"]])[terms.json.AICc], 10, 4);
   fprintf(stdout, io.FormatTableRow(nrm.row, nrm.table_output_options));
}


nrm.json[terms.json.test_results] =nrm.tests;

io.SpoolJSON (nrm.json, nrm.nuc_data_info [terms.json.json]);




