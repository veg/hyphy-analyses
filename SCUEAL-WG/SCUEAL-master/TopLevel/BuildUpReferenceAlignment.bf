_maximumPureRecombLevel = 0.05;


ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"GrabBag.bf");
ExecuteAFile ("../HBF/utils.ibf");

SetDialogPrompt 					("Existing reference file:");
DataSet reference_seq 				= ReadDataFile (PROMPT_FOR_FILE);
referenceFile    					= LAST_FILE_PATH;
referenceAlignmentPath				= referenceFile;
fprintf (stdout, 					"Read ",reference_seq.species, " baseline reference sequences\n");

SetDialogPrompt 					("Sequence file to add to the reference from:");
DataSet ds_in 						= ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, 					"Read ",ds_in.species, " potential new reference sequences\n");
queryFileLog					    = LAST_FILE_PATH + ".log";
fprintf (queryFileLog, CLEAR_FILE, KEEP_OPEN, "screen_results = {};");

/* nucleotide distance threshold */



DISTANCE_PROMPTS				  = 1;
DataSetFilter		filteredData  = CreateFilter (reference_seq,1);
ExecuteAFile 					  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"chooseDistanceFormula.def");
nuc_threshold					  = prompt_for_a_value ("Include only sequences at least this distant from the nearest reference sequence",_maximumPureRecombLevel, 0, 5, 0);

runInMPIMode  = 1;
DataSetFilter ds_fil = CreateFilter (ds_in,1);

DO_NOT_RELOAD_ESTIMATES = 1;
KEEP_ALL_GAPS_IN		= 1;
HAVE_PRESET_VALUES		= 1;


ChoiceList (runType,"Pure/CRF",1,SKIP_NONE,
			"Pure","Pure subtypes only (terminate as soon as recombination is found)",
			"CRF","Accurate recombination mapping"
		    );

if (runType < 0)
{
	return -1;
}

if (runType == 0)
{
	noMoreBPsThan     = 2;
	initPopSize       = 64;
	stoppingCriterion = 50;
	_BICMinLength     = 300;
}
else
{
	noMoreBPsThan 	  = 50;
	initPopSize   	  = 96;
	stoppingCriterion = 50;
	BICMinLength      = 100;
}

alignmentType = 2;

for (seqLoop = 0; seqLoop < ds_in.species; seqLoop = seqLoop + 1)
{
	inOptions2 = {};
	inOptions2["1"] = referenceFile;
	inOptions2["0"] = "Universal";
	GetString (sName, ds_in, seqLoop);
	inOptions2["2"] = sName;
	GetDataInfo (sData, ds_fil, seqLoop);
	inOptions2["3"] = sData;
	
	fprintf (stdout, "Processing ", sName, "\n");
	
	ExecuteAFile (referenceFile + ".labels");
	ExecuteAFile ("../HBF/I_am_the_aligner.bf", inOptions2);

	if (Abs(outputAlignment))
	{
		fprintf (stdout, "\tFinished aligning. Checking minimum distance threshold\n");
		DataSet	 			distanceCheck = ReadFromString (outputAlignment);
		DataSetFilter		filteredData  = CreateFilter (distanceCheck,1);
		InitializeDistances				  (0);
		GetDataInfo 		(alStr,filteredData,distanceCheck.species-1);
		fprintf 			(stdout, "\nAligned string\n", alStr, "\n");
		
		nonGap = 0;
		for (k = 0; k < Abs(alStr); k = k+1)
		{
			if (alStr[k] != "-")
			{
				nonGap = nonGap+1;
			}
		}
		
		if (nonGap/Abs(alStr) < 0.85)
		{
			fprintf (stdout, "[WARNING: low sequence coverage, ", nonGap, "/", Abs(alStr), "]\n");
		}
		
		for (seqCheck = 0; seqCheck < distanceCheck.species-1; seqCheck = seqCheck + 1)
		{
			d = ComputeDistanceFormula (seqCheck, distanceCheck.species-1);
			GetString (closeRelative, distanceCheck, seqCheck);
			fprintf (stdout, "Distance to ", closeRelative, " = ", d, "\n");
			if (d < nuc_threshold)
			{
				break;
			}
		}
		if (seqCheck < distanceCheck.species-1)
		{
			fprintf (stdout, "\tDistance check failed (vs ",closeRelative,"). Skipping this sequence...\n");
			continue;
		}
		
		fprintf (stdout, "\tRunning subtyping on the new sequence\n");
		speciesIndex = 0;
		populationSize  		= initPopSize;
		/*stoppingCriterion		= 50;*/
		rvChoice				= 0;
		siteRateClasses			= 3;
		
		BICMinLength			= _BICMinLength;
		ExecuteAFile ("../HBF/SingleSequenceScan2.bf");
		
		/* how many breakpoints */

		bestLineageAssignments = ConvertToPart(overallBestFound);
		fprintf (stdout, "\nPredicted subtype     : ", bestAssignment, 
						 "\nModel score			  : ", returnAVL["SCORE"],
				 	 	 "\nModel averaged support: ", Format(matchingSum/totalSum*100,8,4), "%\n",
				 	 	 ConvertToPartString(overallBestFound),"\n");
				 	 	 
		if (runType == 0 && returnAVL["RECOMB"] > _maximumPureRecombLevel)
		{
			fprintf (stdout, "[NON-PURE SUBTYPE SKIPPED]\n");
			continue;
		}
		
		fittedModels 	= Abs(MasterList);
		fittedModelsD	= Rows(MasterList);
		
		models			= {fittedModels,1};
		scores			= {fittedModels,1};
		
		for (k = 0; k < fittedModels; k = k+1)
		{
			thisModel = fittedModelsD[k];
			models[k] = thisModel^{{"\\\n",";"}};
			scores[k] = MasterList[thisModel];
		}
				
		fprintf (queryFileLog, returnAVL, "_hyphyAssociativeArray[\"Tree\"] = \"",Format(referenceTopology,1,1),"\";",
							   "_hyphyAssociativeArray[\"Models\"] = ",models,";",
							   "_hyphyAssociativeArray[\"Scores\"] = ",scores,";",
					       	   "\nscreen_results[\"", sName, "\"] = _hyphyAssociativeArray;");
		
		fittedModelsD = 0; models = 0; scores = 0;
		
		bp_count = Rows(bestPC);
		fprintf (stdout, "Add ", sName, " to the reference alignment (y/n)?");
		fscanf  (stdin,"String", doAdd);
		if (doAdd == "y")
		/*add this sequence to the reference alignment*/
		{
			GetDataInfo (alignedQuery,filteredDataBig,filteredDataBig.species-1);
			
			alreadyDoneTop		   = {};
			
			refTopAVL				= Format (referenceTopology,1,0);
			INTERNAL_NODE_PREFIX	= "Node";
			
			Tree				   refTreeBase = refTopAVL;
			refTopAVL			   = refTreeBase^0;
			currentTreeString	   = SCUEALTreeAVL2String(refTopAVL);	
			sequencesToAdd		   = {};
			newString = ""; newString * 128; 

			sequenceLabels = _subtypeAssignmentByNode;
			fprintf 		(stdout, "\nLabel for tip ", sName, " : ");
			fscanf			(stdin, "String", nodeLabel);
			
			for (_bp = 0; _bp <= bp_count; _bp = _bp+1)
			{
				thisLineage			  = bestLineageAssignments[_bp][0];
				if (alreadyDoneTop[thisLineage] == 0)
				{
					thisLineageName		  = (refTopAVL[thisLineage])["Name"];
					alreadyDoneTop	[thisLineage] = 1;
					Topology cT			  = currentTreeString;
					newTopTopAVL		  = cT^0;
					if (bp_count)
					{
						taxonName = sName + "_CRF_" + Abs(alreadyDoneTop);
					}
					else
					{
						taxonName = sName;
					}
					sequenceLabels [taxonName&&1] = nodeLabel;
					for (_ni = 1; _ni < Abs(newTopTopAVL); _ni = _ni + 1)
					{
						if (((newTopTopAVL[_ni])["Name"]&&1) == thisLineageName)
						{
							break;
						}
					}
					SCUEALInsertANode 			 ("newTopTopAVL",_ni,taxonName,"NODE"+Abs(newTopTopAVL));
					currentTreeString	   = SCUEALTreeAVL2String(newTopTopAVL);	
					/* generate a new reference string */
					
					bpF 	 = -1;
					bpF2	 = -1;
					
					newString * ("\n>"+taxonName+"\n");
					
					for (h=1; h<=bp_count; h=h+1)
					{
						bpF2 = bppMap[sortedBP[h][1]];
						curSegLength = bpF2-bpF;
						if (thisLineage != bestLineageAssignments[h-1][0])
						{
							for (z=0; z < curSegLength; z = z + 1)
							{
								newString * "-";
							}
						}
						else
						{
							newString * (alignedQuery[bpF+1][bpF2]);
						}
						bpF = bpF2;
					}
					
					if (bpF2<ds.sites)
					{
						curSegLength = ds.sites-bpF2;
						if (thisLineage != bestLineageAssignments[bp_count][0])
						{
							for (z=0; z < curSegLength-1; z = z + 1)
							{
								newString * "-";
							}
						}
						else
						{
							newString * (alignedQuery[bpF+1][Abs(alignedQuery)-1]);
						}
					}
				}
			}
			newString * 0; 
			DeleteObject			(exportT);
			Tree	exportT			= currentTreeString;
			currentTreeString		= Format(exportT,0,0);

			DataSet addOn 			= ReadFromString (newString);
			DataSet reference_seq 	= Combine (ds,addOn);
			DataSetFilter exporter	= CreateFilter (reference_seq,1,"",speciesIndex!=ds.species-1);
			IS_TREE_PRESENT_IN_DATA	= 1;
			DATA_FILE_PRINT_FORMAT  = 6;
			DATAFILE_TREE		    = currentTreeString;
			fprintf 				(referenceFile,CLEAR_FILE,exporter);
				
			fprintf (stdout, 		"\nAuto-generating internal node labels\n"); 
			
			_doLabelGeneration	();
			
			outLabels = "_subtypeAssignmentByNode  = " + sequenceLabels + ";\n\n_crfEquiv = " + _crfEquiv + ";";
			
			/*
			outLabels = ("" + sequenceLabels)^{{"_hyphyAssociativeArray","_subtypeAssignmentByNode"}};
			*/
			
			fprintf (stdout, "\nUpdating the reference alignment...\n");
			fprintf (stdout, "\nFitting a codon model to estimate the ancestor...");
			
			codonOptions = {};
			codonOptions ["0"] = "Universal";
			codonOptions ["1"] = referenceFile;
			codonOptions ["2"] = "MG94CUSTOM";
			codonOptions ["3"] = "Global";
			codonOptions ["4"] = "012345";
			codonOptions ["6"] = "y";
			codonOptions ["7"] = "Estimate";
			
			INTERNAL_NODE_PREFIX	= "NODE";

			/*
			USE_DISTANCES			= 0;
			*/
			
			VERBOSITY_LEVEL        = 10;
			OPTIMIZATION_PRECISION = 1;
			
			ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
											   + "AnalyzeCodonData.bf", codonOptions);
			
			OPTIMIZATION_PRECISION = 0.001;
			VERBOSITY_LEVEL 	   = 0;
			fprintf (stdout, "\nReconstructing ancestors...");		
						
			DataSet			ancestralSequences = ReconstructAncestors (lf);
			DataSet			jointDS			   = Combine (ds,ancestralSequences);
			DataSetFilter	referenceFilter	   = CreateFilter (jointDS,1,"",speciesIndex <= filteredData.species);
			
			IS_TREE_PRESENT_IN_DATA 		   = 1;
			DATA_FILE_PRINT_FORMAT			   = 6;
			DATAFILE_TREE					   = currentTreeString;
			
			fprintf (referenceFile,CLEAR_FILE,referenceFilter);
			outLabelFile = referenceFile + ".labels";
			fprintf (outLabelFile, CLEAR_FILE, outLabels);

		}
	}
	else
	{
		fprintf (stdout, "\tAlignment failed. Skipping sequence...\n");
	}
}
fprintf (queryFileLog, CLOSE_FILE);

