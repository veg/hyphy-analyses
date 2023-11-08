LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/mpi.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");



/* 1. include a file to define the genetic code
   Note the use of base directory and path forming variables to make this analysis 
   independent of directory placement
 */

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

paml.branch_site.json = {};

namespace paml.branch_site {
    LoadFunctionLibrary (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles/SelectionAnalyses/modules/shared-load-file.bf");
    load_file (
                {
                    utility.getGlobalValue("terms.prefix"): "paml.branch_site", 
                    utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}
                }
            );
}

paml.filter_name = paml.branch_site.filter_names["0"];
tree_annotation = (paml.branch_site.partitions_and_trees["0"])[terms.data.tree];

if (utility.Array1D (tree_annotation[terms.trees.model_list]) < 2) {
    tree_annotation = trees.extract_paml_annotation ((tree_annotation[terms.trees.newick_annotated])[0][-1]);
}

treeString = tree_annotation[terms.trees.newick];

paml.ds_name = (paml.branch_site.codon_data_info)[terms.data.dataset];

/* 4. Compute nucleotide counts by position for the F3x4 estimator */

COUNT_GAPS_IN_FREQUENCIES = 0;
HarvestFrequencies (baseFreqs,^paml.filter_name,3,1,1);
COUNT_GAPS_IN_FREQUENCIES = 1;

_Genetic_Code = paml.branch_site.codon_data_info[terms.code];

fprintf (stdout, "\nBase composition:\n\tA: ", Format (baseFreqs[0][0],10,5),",",Format (baseFreqs[0][1],10,5),",",Format (baseFreqs[0][2],10,5),
								    "\n\tC: ", Format (baseFreqs[1][0],10,5),",",Format (baseFreqs[1][1],10,5),",",Format (baseFreqs[1][2],10,5), 
									"\n\tG: ", Format (baseFreqs[2][0],10,5),",",Format (baseFreqs[2][1],10,5),",",Format (baseFreqs[2][2],10,5), 
									"\n\tT: ", Format (baseFreqs[3][0],10,5),",",Format (baseFreqs[3][1],10,5),",",Format (baseFreqs[3][2],10,5), "\n");
										  


/* 6. define the 'site_kind' variable as a discrete category variable; 
	  it decides which class a site belongs to, but does not
	  determine omega ratios directly (see below for this) */


global P_0     = 0.5;
P_0 :< 1;

global P_1     = 0.5;
P_1 :< 1;


rateClasses = 4;

categFreqMatrix = {{P_0,
                    (1-P_0)*P_1,
                    (1-P_0)*(1-P_1)/(P_0+(1-P_0)*P_1)*P_0,
                    (1-P_0)*(1-P_1)/(P_0+(1-P_0)*P_1)*(1-P_0)*P_1}} ;


categRateMatrix = {{1,2,3,4}};
		
	
category site_kind = (rateClasses, categFreqMatrix , MEAN, ,categRateMatrix, 1, 4);

/* 7. define the GY94 rate matrix; for now each branch will have it's own
   dS and dN, we will constrain them later */

global kappa_inv = 1;

ModelMatrixDimension = 64;
for (h = 0; h<64; h += 1)  {
	if (_Genetic_Code[h]==10) /* stop codon */
	{
		ModelMatrixDimension = ModelMatrixDimension-1;
	}
}

GY_Matrix = {ModelMatrixDimension,ModelMatrixDimension};
hshift = 0;
for (h=0; h<64; h=h+1) {
	if (_Genetic_Code[h]==10) {
		hshift += 1;
	}
	else {
		vshift = hshift;
		for (v = h+1; v<64; v=v+1) {
			diff = v-h;
			if (_Genetic_Code[v]==10) {
				vshift = vshift+1;
			}
			else {
			  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* one step */
			  	{
			  		if (h$4==v$4)
			  		{
			  			transition = v%4;
			  			transition2= h%4;
			  		}
			  		else
			  		{
			  			if(diff%16==0)
			  			{
			  				transition = v$16;
			  				transition2= h$16;
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := kappa_inv*synRate;
			  				GY_Matrix[v-vshift][h-hshift] := kappa_inv*synRate;
			  			}
			  			else
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := synRate;
			  				GY_Matrix[v-vshift][h-hshift] := synRate;			  			
			  			}
				  	}
			  		else
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := kappa_inv*nonSynRate;
			  				GY_Matrix[v-vshift][h-hshift] := kappa_inv*nonSynRate;
			  			}
			  			else
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := nonSynRate;
			  				GY_Matrix[v-vshift][h-hshift] := nonSynRate;			  			
			  			}
		  			}
			  	}
			 }
		 }
	}	
}

/*8. build codon frequencies (use the F3x4 estimator) */

PIStop = 1.0;
codonFreqs = {ModelMatrixDimension,1};
hshift = 0;

for (h=0; h<64; h=h+1) {
	first  = h$16;
	second = h%16$4;
	third  = h%4;
	if (_Genetic_Code[h]==10) {
		hshift = hshift+1;
		PIStop = PIStop-baseFreqs[first][0]*baseFreqs[second][1]*baseFreqs[third][2];
		continue; 
	}
	codonFreqs[h-hshift]=baseFreqs[first][0]*baseFreqs[second][1]*baseFreqs[third][2];
}

codonFreqs = codonFreqs*(1.0/PIStop);

/*9. define the codon model */

Model GY_Model = (GY_Matrix,codonFreqs,1);

/*10. Define the tree and pick the foreground branch, displaying a tree window to facilitate selection; 
the latter step is executed for 2 of 3 model choices */

Tree 	  givenTree = treeString;

USE_LAST_RESULTS    = 0;
OPTIMIZATION_METHOD = 4;

/* Approximate kappa and branch lengths using an HKY85 fit */

HKY85_Matrix = {{*,t*kappa_inv,t,t*kappa_inv}
				{t*kappa_inv,*,kappa_inv*t,t}
				{t,t*kappa_inv,*,kappa_inv*t}
				{t*kappa_inv,t,kappa_inv*t,*}};
			
HarvestFrequencies (nucFreqs,^paml.ds_name,1,1,1);
Model HKY85_Model = (HKY85_Matrix,nucFreqs,1);

Tree		  nucTree = treeString;
DataSetFilter nucData = CreateFilter (^paml.ds_name,1);

fprintf (stdout, "Obtaining nucleotide branch lengths and kappa to be used as starting values...\n");
LikelihoodFunction	nuc_lf = (nucData,nucTree);
Optimize(nuc_mle,nuc_lf);
fprintf (stdout, "\n", Format (nucTree,1,1), "\nkappa=", Format (1/kappa_inv,8,3), "\n");

USE_LAST_RESULTS = 1;

choices = {};
for (m;in;tree_annotation[terms.trees.model_list]) {
    if (Abs (m)) {
        choices[m] = "Set foreground branches to group `m`";
    }
}

assert (Abs (choices) >= 1, "MUST have at least once labeled set of branches");

if (Abs (choices) > 1) {
    choice_made = io.SelectAnOption (choices, "Use these branches as foreground");
} else {
    choice_made = utility.Keys (choices)[0];
}

selected_branches = {};
for (b,m; in; tree_annotation[terms.trees.model_map]) {
    if (m == choice_made) {
        selected_branches + b;
    }
}

fprintf (stdout, "\n\n", Abs (selected_branches)," foreground branch(es) set to: ", Join(", ", selected_branches), "\n");
	
/* 15. Constrain dS and dN in the tree to based upon different models */

global omega_0 = 0.25;
omega_0 :< 1;
ClearConstraints (givenTree);


global omega_2 = 2.0;
omega_2:>1;
global omega_FG := ((site_kind==1)*omega_0+(site_kind==2)+(site_kind>2)*omega_2); 			/* foreground model */
global omega_BG := (((site_kind==1)+(site_kind==3))*omega_0+(site_kind==2)+(site_kind==4)); /* background model */

for (b; in; selected_branches) {
	ExecuteCommands ("givenTree.`b`.nonSynRate:=omega_FG*givenTree.`b`.synRate;");
}

	/* constrain other branches next */
ReplicateConstraint ("this1.?.nonSynRate:=omega_BG*this2.?.synRate",givenTree,givenTree);

/* 16. define and optimize the likelihood function */

bNames = BranchName   (givenTree,-1);
nucBL  = BranchLength (nucTree,-1);

for (bc=0; bc<Columns(bNames)-1; bc += 1) {
	ExecuteCommands ("givenTree."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t;");
}

codBL  = BranchLength (givenTree,-1);

for (bc=0; bc<Columns(bNames)-1; bc += 1) {
	if (nucBL[bc]>0) {
		ExecuteCommands ("givenTree."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t*"+nucBL[bc]/codBL[bc]+";");
	}
}

LikelihoodFunction lf = (^paml.filter_name, givenTree);
start.grid  =  {    
    "0" : {
        "P_0"        : 0.8,
        "P_1"        : 0.5,
        "omega_0"    : 0.1,
        "omega_2"    : 1.0,
    },
    "1" : {
        "P_0"        : 0.8,
        "P_1"        : 0.5,
        "omega_0"    : 0.1,
        "omega_2"    : 2.0
    },
    "2" : {
        "P_0"        : 0.8,
        "P_1"        : 0.5,
        "omega_0"    : 0.1,
        "omega_2"    : 5.0
    },
    "3" : {
        "P_0"        : 0.5,
        "P_1"        : 0.5,
        "omega_0"    : 0.1,
        "omega_2"    : 2.0
    },
    "4" : {
        "P_0"        : 0.5,
        "P_1"        : 0.5,
        "omega_0"    : 0.1,
        "omega_2"    : 2
    }
};
                        

while (1) {
	Optimize 		   (mles,lf,  {
            "OPTIMIZATION_START_GRID" : start.grid             
    });
    
	fprintf (stdout, "\nLog(L) = ", mles[1][0]);
	mle_saved = estimators.TakeLFStateSnapshot("lf");
	
	samples = 500;
	fprintf (stdout, "\nChecking for convergence by Latin Hypercube Sampling (this may take a bit of time...)\n");
		
    sampling_ranges = {
                        "P_0" : {
                            terms.lower_bound : P_0/2,
                            terms.upper_bound : Min (P_0*2,1)
                        },
                        "P_1" : {
                            terms.lower_bound : P_1/2,
                            terms.upper_bound : Min (P_1*2,1).
                        },  
                        "omega_0" : {
                            terms.lower_bound : omega_0 / 3,
                            terms.upper_bound : omega_0 * 3
                        },  
                        "omega_2" : {
                            terms.lower_bound : 1,
                            terms.upper_bound : 5
                        }
                    };              
                        
 
    samples = estimators.LHC (sampling_ranges, samples);
    grid_results = mpi.ComputeOnGrid (&lf,samples, "mpi.ComputeOnGrid.SimpleEvaluator", "mpi.ComputeOnGrid.ResultHandler");
	grid_max = Max (grid_results,0);
	
	
	estimators.RestoreLFStateFromSnapshot("lf", mle_saved);
	 
	if (grid_max["value"] > mles[1][0]) {
	    
	    parameters.SetValues(samples[grid_max["key"]]);
	    start.grid  = null;
		fprintf (stdout, "\nConvergence checks FAILED : a better score found. Restarting the optimization procedure");
		continue;
	} 

	
	fprintf (stdout, "\nThe estimation procedure appears to have converged.\n");
	
	cmc_ll = mles[1][0];
	
	break;
}

cmc.json = {
    terms.json.fits : {
        
    },
    terms.json.test_results : {
    
    }
};


/* 17. Report inferred rate distribition to screen */

function print_fit (model, ll) {
    
    (cmc.json[terms.json.fits])[model] = {
        terms.json.log_likelihood : ll,
        "Class 1 omega" : omega_0,
        "Class 1 weight" : categFreqMatrix[0],
        "Class 2 omega" : 1,
        "Class 2 weight" : categFreqMatrix[1],
        "Class 3 omega FG" : omega_2,
        "Class 3 omega BG" : omega_0,
        "Class 3 weight" : categFreqMatrix[2],
        "Class 4 omega FG" : omega_2,
        "Class 4 omega BG" : 1,
        "Class 4 weight" : categFreqMatrix[3]
    };
    
    fprintf (stdout, "** MODEL `model` **\nLog (L) = `ll`\n",
            "\nInferred rate distribution:",
					 "\n\tClass 0.  omega_0 = ", Format (omega_0, 5,3), " weight = ", Format (categFreqMatrix[0],8,5),
					 "\n\tClass 1.  omega  := ", Format (1, 5,3), " weight = ", Format (categFreqMatrix[1],8,5),
					 "\n\tClass 2a. Background omega_0 = ", Format (omega_0, 5,3), " foreground omega_2 = ", Format (omega_2, 5,3), " weight = ", Format (categFreqMatrix[2],8,5), 
					 "\n\tClass 2b. Background omega  := ", Format (1, 5,3), " foreground omega_2 = ", Format (omega_2, 5,3), " weight = ", Format (categFreqMatrix[3],8,5), "\n");
}

print_fit ("Branch-site model A ", cmc_ll);

fprintf (stdout, 
"
===================================
** Fitting Model A omega2:=1     **
===================================
");

omega_2 := 1;
Optimize (mle_MA2, lf);

print_fit ("Model A omega2:=1", mle_MA2[1][0]);
(cmc.json [terms.json.test_results]) ["Model A|Model A neutral"] =  math.DoLRT (mle_MA2[1][0], cmc_ll,1);

fprintf (stdout, "\n p-value, Model A vs Model A neutral= ", Format (((cmc.json [terms.json.test_results]) ["Model A|Model A neutral"])[terms.p_value], 5, 10));

fprintf (stdout, 
"
===============================
** Fitting null model M1a    **
===============================
");


P_1 := 1;
Optimize (mle_M1a, lf);
print_fit ("M1a model", mle_M1a[1][0]);
(cmc.json [terms.json.test_results]) ["Model A|M1a"] =  math.DoLRT (mle_M1a[1][0], cmc_ll , 2);

fprintf (stdout, "\n p-value, Model A vs M1a = ", Format (((cmc.json [terms.json.test_results]) ["Model A|M1a"])[terms.p_value], 5, 10));


io.SpoolJSON (cmc.json, io.PromptUserForFilePath("Write model fit result to"));

