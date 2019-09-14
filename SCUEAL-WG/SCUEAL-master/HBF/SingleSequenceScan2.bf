RequireVersion			("2.5");
LoadFunctionLibrary		("PS_Plotters.bf");
LoadFunctionLibrary		("GrabBag.bf");
LoadFunctionLibrary		("ReadDelimitedFiles.bf");

verboseFlag				   = 0;
timeLimit				   = 1000000;


if (runInMPIMode && MPI_NODE_COUNT>1)
{
	verboseFlag = 0;
	timeLimit   = 3600;
}

partCount				   = 2;
VERBOSITY_LEVEL			   = -1;
DO_NOT_RELOAD_ESTIMATES    = 1;
USE_ADAPTIVE_VARIABLE_STEP = 1;

stepByStepLogging		= 1;

settingsFilePath 		=  PATH_TO_CURRENT_BF || ("\\"+DIRECTORY_SEPARATOR);
settingsFilePath		=  PATH_TO_CURRENT_BF[0][settingsFilePath[Rows(settingsFilePath)-3]] + 
						   "Configs"+DIRECTORY_SEPARATOR+"settings.ibf";

if (runInMPIMode == 0 && !settingsFilePath)
{
	if (verboseFlag)
	{
		fprintf (stdout, "Loading settings...\n");
	}
	ExecuteAFile ("../Configs/settings.ibf");
	
	referenceAlignmentPath 		=  PATH_TO_CURRENT_BF || ("\\"+DIRECTORY_SEPARATOR);
	referenceAlignmentPath		=  PATH_TO_CURRENT_BF[0][referenceAlignmentPath[Rows(referenceAlignmentPath)-3]] + 
							   		"data"+DIRECTORY_SEPARATOR+referenceAlignmentFileName;

}
else
{
	if (HAVE_PRESET_VALUES == 0)
	{
		populationSize  		= 64;
		stoppingCriterion		= 50;
		rvChoice				= 0;
		siteRateClasses			= 3;
		noMoreBPsThan			= 10;
		BICMinLength			= 100;
	}
}


startBPC				= 1;

convergenceBump			= 5;
populationSizeBump		= 8;
_contigGapThresh		= Min(50,BICMinLength);
_randomizedCRFGroups	= 5000;

//fprintf (stdout, "\nBICMinLength:", BICMinLength, "\n_contigGapThresh:", _contigGapThresh, "\n");


produceOffspring		= 3*populationSize$4;
incestDistance  		= 0;
generationCount		  	= 5000;
maxSampleTries			= populationSize*10;
mutationThreshhold		= 0.001;
mutationProbStart		= 0.20;
mutationProbDecrease    = 0.95;
mutationProbMin			= 0.05;
annealingPhase			= 3*stoppingCriterion$2;
SHORT_MPI_RETURN		= 1;
totalSampleCounter		= 0;
localMutationRate		= 0.05;
localMutationMultiplier = 0.45;
localMutationMultiplierDecrease = 0.03;
localMutationMultiplierMin		= 0.10;
localMutationInterval	= Min(stoppingCriterion-1,stoppingCriterion * localMutationMultiplier $ 1);
shortestAllowedSegment  = 0;


prematureTermination	= 0;
scuealTimer				= Time(1);

sampleCount				= 0;
familyControlSize		= produceOffspring$8;

rateClassesCount		= 2;
MESSAGE_LOGGING			= 0;
cAICPenaltyPerSite		= 50;
adjustAICScores			= 1;
matingChoice			= 0;
strictViableCheck		= 0;
excludeSameSequence		= 0;
mildMatingProcess		= 0;
useGrayCodes			= 1;
swapBeforeMateProb		= 0.10;
mutateAfterMate			= 0.01;

runSubpopulations		= 3;
tryBreakpointFusion		= 1;

initialPops				= {};

if (verboseFlag)
{
	fprintf (stdout, "\nRun options: ",
					 "\n\tPopulation size: ", populationSize,
					 "\n\tCovergence criterion: ", stoppingCriterion,
					 "\n\tRate variation option: ", rvChoice,
					 "\n\tMaximum number of breakpoints: ", noMoreBPsThan,
					 "\n\tMinimum fragment length: ", BICMinLength,
					 "\n\tLocal Mut. Interval: ", localMutationInterval,
					 "\n");
}



/* ________________________________________________________________________________________________*/

global AC 				= 1;
global AT 				= 1;
global CG 				= 1;
global CT 				= 1;
global GT 				= 1;

INTERNAL_NODE_PREFIX	= "INODE";

/* ________________________________________________________________________________________________*/

MasterList						= {};
REPLACE_TREE_STRUCTURE  		= 1;
bppMap							= {};
bppMapLength					= 0;
SHORT_MPI_RETURN				= 1;
totalBitSize					= 0;
branchBits						= 0;
bitsPerPart						= 0;
FILE_SEPARATOR			   		= "__FILE_SEPARATOR__";
treeStringCache					= {};
treeConstraintCache				= {};

analysisPrefix					= "";

/*************************************************************************************/

function reduceSubtype (in)
{
	if (Abs(_subtypeReductions[in]))
	{
		return _subtypeReductions[in];
	}
	_locST = _subtypeAssignmentByNode[in];
	if (Abs(_locST)==0)
	{
		_locST = _subtypeAssignmentByNode[in[1][Abs(in)-1]];
		if (Abs(_locST)==0)
		{
			return "Unlabeled";
		}
	}
	return _locST;
}

/*************************************************************************************/

function	AssembleSubtypeAssignment (theDatum,doConversion)
{
	if (doConversion)
	{
		theDatumString = ConvertToPartString (theDatum);
	}
	else
	{
		theDatumString = theDatum;
	}
	
	_s = 0;
	while (theDatumString[_s] != "\n")
	{
		_s = _s+1;
	}
	theDatumString = splitOnRegExp(
					 theDatumString[_s+1][Abs(theDatumString)-1],",");
					 
	subtypeCounter  = {};				 
	crfEquivCounter = {};
	
	_hasCRFs   = Abs(_crfEquiv);
	_foundCRFs = {}; 
	_howManyParts =  Abs (theDatumString);
	crfEquivCounter = {_howManyParts,1};
	_crfReducedCounter = {};
	
	for (_s=0; _s < _howManyParts; _s=_s+1)
	{
		thisNodeSubtype  = reduceSubtype(theDatumString[_s]);
		if (_s > 0)
		{
			subtypeString    = subtypeString + "," + thisNodeSubtype;
		}
		else
		{
			subtypeString  = thisNodeSubtype;
		}
		subtypeCounter  [thisNodeSubtype] = 1;
		if (_hasCRFs)
		{
			_thisCRF = _crfEquiv[theDatumString[_s]];
			if (Abs(_thisCRF))
			{
				_foundCRFs[thisNodeSubtype] = 1;
				crfEquivCounter[_s] = _thisCRF;
				_crfReducedCounter[_thisCRF] = 1;
			}
			else
			{
				crfEquivCounter[_s] = thisNodeSubtype;
				_crfReducedCounter[thisNodeSubtype] = 1;
			}
		}
		subtypeCounter  [thisNodeSubtype] = 1;
	}
	
	classType = 0;
	if (Abs (theDatumString) > 1 )
	{
		if (Abs (subtypeCounter) > 1)
		{
			subtypeString = subtypeString + " inter-subtype recombinant";
			classType = 1;
		}
		else
		{
			subtypeString = reduceSubtype(theDatumString[0]) + " intra-subtype recombinant (" + 
							 (Abs (theDatumString)-1) + " breakpoints)";
			classType = 2;
		}
	}
	
	return subtypeString;
}

/*************************************************************************************/

function isEqualOrSuperset (profBits, referenceBits)
{
	size1 = Rows (profBits);
	size2 = Rows (referenceBits);
	if (size1 >= size2)
	{
		currentRefBit  = 0;
		currentProfBit = 0;
		while (currentProfBit < size1 && currentRefBit < size2 )
		{
			if (profBits[currentProfBit] == referenceBits[currentRefBit][0])
			{
				currentProfBit  = currentProfBit + 1;
				currentRefBit   = currentRefBit + 1;
			}
			else
			{
				if (profBits[currentProfBit] == referenceBits[currentRefBit][0])
				{
					currentProfBit = currentProfBit + 1;	
				}
				else
				{
					return 0;
				}
			}
		}
		if (currentProfBit == size1 && currentRefBit == size2)
		{
			/*fprintf (stdout, profile, ">=", reference, "\n");*/
			return 1;
		}
	}
	return 0;
}

/*************************************************************************************/

function	AssembleSubtypeAssignmentSimple (theDatum,doConversion)
{
	fullSubtypeString = AssembleSubtypeAssignment(theDatum, doConversion);
	if (classType > 0)
	{
		if (classType == 2)
		{
			return (Rows(subtypeCounter))[0];
		}
		else
		{
			_fcrfC = Abs(_foundCRFs);
			if (_fcrfC)
			{
				/* check for CRF-like structure */
				tryTheseCRFs = Rows (_foundCRFs);
				_crfLike	 = {};
				
				for (crfID = 0; crfID < _fcrfC; crfID = crfID+1)
				{
					_crfToTry = tryTheseCRFs[crfID];
					crfStructureCounts = Abs(CRFMosaics[_crfToTry]);
					for (crfID2 = 0; crfID2 < crfStructureCounts; crfID2 = crfID2+1)
					{
						if (verboseFlag)
						{
							fprintf (stdout, crfEquivCounter, "\n", ((CRFMosaics[_crfToTry])[crfID2]), "\n");
						}
						if (isEqualOrSuperset (crfEquivCounter, (CRFMosaics[_crfToTry])[crfID2]))
						{
							return _crfToTry + "-like";
						}
					}
				}
				if (Abs(_crfReducedCounter) == 1)
				{
					return (Rows(_crfReducedCounter))[0];
				}
				subtypeCounter = _crfReducedCounter;
			}	
		}
		
		if (Abs(subtypeCounter) > 2)
		{
			return "Complex";
		}
		if ((Rows(subtypeCounter))[0] < (Rows(subtypeCounter))[1])
		{
			return (Rows(subtypeCounter))[0] + "," + (Rows(subtypeCounter))[1] + " recombinant";
		}
		return (Rows(subtypeCounter))[1] + "," + (Rows(subtypeCounter))[0] + " recombinant";
	}
	return fullSubtypeString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StringToMatrix (zz)
{
	return zz;
}

_psTreePlots    = {};
treeImageHeight = 0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function DefineAPart (filterString,daTree,index,bls,mrca_name)
{
	ExecuteCommands ("DataSetFilter filter_" + index + " = CreateFilter (filteredDataBig,1,\""+filterString+"\");");
	ExecuteCommands ("Tree 			tree_ " + index + " = " + daTree + ";treeAVL_export=tree_" + index + "^0;");
	_interiorBranchCount = 0;
	for (_k = 1; _k < Abs (treeAVL_export); _k = _k+1)
	{
		_k2 = (treeAVL_export[_k])["Name"];
		if (_k2 == mrca_name)
		{
			plotRootLoc = _interiorBranchCount;
		}
		else
		{
			if (Abs((treeAVL_export[_k])["Children"]))
			{
				_interiorBranchCount = _interiorBranchCount + 1;
			}
		}
		(TREE_OUTPUT_OPTIONS [_k2]) = {};
		if (Abs(_subtypeAssignmentByNode[_k2]))
		{
			(TREE_OUTPUT_OPTIONS [_k2])["TREE_OUTPUT_BRANCH_TLABEL"] = reduceSubtype(_k2);
		}
		else
		{
			(TREE_OUTPUT_OPTIONS [_k2])["TREE_OUTPUT_BRANCH_TLABEL"] = "";		
		}
		(TREE_OUTPUT_OPTIONS [_k2])["TREE_OUTPUT_BRANCH_COLOR"]  = {{0.8,0.8,0.8}};
		if (_k2 == "_SPLICE_NODE_")
		{
			ExecuteCommands("tree_" + index + "." + _k2 + ".t = " + bls[1]);
			ExecuteCommands("tree_" + index + ".QUERY.t = " + bls[2]);
			splicedAt	 = (treeAVL_export[((treeAVL_export[_k])["Children"])[0]])["Name"];
			(TREE_OUTPUT_OPTIONS [splicedAt])["TREE_OUTPUT_BRANCH_DASH"] = {{2,2,0}};
			ExecuteCommands("tree_" + index + "."+splicedAt+".t = " + bls[0]);
			
		}
		else
		{
			if (nodeNameToAVL[_k2])
			{
				ExecuteCommands ("tree_" + index + "." + _k2 + ".t = baselineTree." + _k2 + ".t");
			}
			_k3 = branchColorMap[_k2]-1;
			if (_k3>=0)
			{
				(TREE_OUTPUT_OPTIONS [_k2])["TREE_OUTPUT_BRANCH_COLOR"]  = plotColors[_k3][-1];
			}
		}
		
	}
	(TREE_OUTPUT_OPTIONS ["QUERY"])["TREE_OUTPUT_BRANCH_DASH"] = {{2,2,0}};
	(TREE_OUTPUT_OPTIONS ["_SPLICE_NODE_"])["TREE_OUTPUT_BRANCH_DASH"] = {{2,2,0}};
	(TREE_OUTPUT_OPTIONS ["QUERY"])["TREE_OUTPUT_BRANCH_LABEL"] = "0 0 3 0 360 arc fill";
	ACCEPT_ROOTED_TREES = 1;
	ExecuteCommands ("plotTree=tree_"+index+"[plotRootLoc][\"EXPECTED_NUMBER_OF_SUBSTITUTIONS\"];");
	
	
	if (Type (plotTree) != "Unknown")
	{
		Tree toPlot = plotTree;
		treeImageHeight = Max(200,10*TipCount(toPlot));
		_psTreePlots[index-1]=PSTreeString (toPlot,"STRING_SUPPLIED_LENGTHS",{{treePlotWidth,treeImageHeight}});
	}
	else
	{
		treeImageHeight = Max(200,10*TipCount(referenceTopology));
		ExecuteCommands ("_psTreePlots[index-1]=PSTreeString (tree_"+index+",\"EXPECTED_NUMBER_OF_SUBSTITUTIONS\",{{treePlotWidth,treeImageHeight}});");
	}
	
	if (index>1)
	{
		return ",filter_" + index + ",tree_" + index;
	}
	return "filter_" + index + ",tree_" + index;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (rateMatrix,bestModelBL,plotRoot)
{
	TREE_OUTPUT_OPTIONS = {};
	TREE_OUTPUT_OPTIONS ["TREE_OUTPUT_LAYOUT"] = 0;
	TREE_OUTPUT_OPTIONS ["TREE_OUTPUT_EMBED"]  = 1;
	
	SetParameter(ds,ds.species-1,"QUERY");
	UseModel	(nucModel);
	sortedBP 	= ConvertToPart (rateMatrix);
	
	theAVL		= {};
	theAVL 		["BP"]    = {};
	theAVL 		["Trees"] = {};
	partSpecs		 	  = {};
	
	bpF 			= -1;
	bpF2			= -1;
	
	v 		  = Rows (sortedBP);
	treePlotWidth = totalPlotWidth/v;
	spliceLoc = sortedBP[0][0];
	LF_DEF	  = "LikelihoodFunction lf_export = (";

	for (h=1; h<v; h=h+1)
	{
		bpF2 					= bppMap[sortedBP[h][1]];
		LF_DEF					= LF_DEF + DefineAPart	(""+(bpF+1)+"-"+bpF2, treeStringCache[sortedBP[h-1][0]], h, bestModelBL[h-1][-1],plotRoot);
	}
	
	if (bpF2<ds.sites)
	{
		LF_DEF					= LF_DEF + DefineAPart	( ""+(bpF2+1)+"-"+(ds.sites-1), treeStringCache[sortedBP[h-1][0]], h, bestModelBL[h-1][-1],plotRoot);
	}
	return LF_DEF+")";
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function CleanUpMPI (dummy)
{
	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
				{
					fromNode = ReceiveJobs (0,0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}			
	}
	return 0;
}
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function adjustAICScore (theLL,bpMatrix)
{
	myDF	     = baseParams + tbc2 + (Rows (bpMatrix)-1); 
	if (icOption == 0)
	{
		daAICScore   = 2*(myDF*(baseSites/(baseSites-myDF-1)) - theLL) ;
		allDF	     = myDF + 3*Rows(bpMatrix);
		baseAICScore = 2*(allDF*(baseSites/(baseSites-allDF-1)) - theLL);
	}
	else
	{
		daAICScore   = myDF*Log(baseSites) - 2*theLL ;
		baseAICScore = (myDF + 3*Rows(bpMatrix))*Log(baseSites) - 2*theLL;	
	}
	lastBpos     = 0;
	
	allTheSame 	 = 1;
	
	for (_pid = 1; _pid < Rows(bpMatrix); _pid = _pid+1)
	{
		thisSpan = bppMap[bpMatrix[_pid][1]] - lastBpos+1;
		lastBpos = bppMap[bpMatrix[_pid][1]];
		allTheSame = allTheSame && (bpMatrix[_pid][0] == bpMatrix[_pid-1][0]);
		if (icOption == 0)
		{
			if (thisSpan > tbc2)
			{
				daAICScore = daAICScore + 6*(thisSpan/(thisSpan-tbc2-1));
			}
			else
			{
				daAICScore = daAICScore + 6*cAICPenaltyPerSite;
			}
		}
		else
		{
			if (thisSpan < BICMinLength)
			{
				daAICScore = daAICScore + 1000000;
			}
			else
			{
				if (excludeSameSequence && bpMatrix[_pid][0] == bpMatrix[_pid-1][0])
				{
					daAICScore = daAICScore + 1000000;
				}
				else
				{
					daAICScore = daAICScore + 3*Log(thisSpan);
				}
			}
		}
	}
	savebpm = bpMatrix;
	thisSpan = baseSites-lastBpos;
	if (icOption == 0)
	{
		if (thisSpan > tbc2)
		{
			daAICScore = daAICScore + 6*(thisSpan/(thisSpan-tbc2-1));
		}
		else
		{
			daAICScore = daAICScore + 6*cAICPenaltyPerSite;
		}
		return -Max(daAICScore,baseAICScore);
	}
	if (thisSpan < BICMinLength)
	{
		daAICScore = daAICScore + 1000000;
	}
	else
	{
		daAICScore = daAICScore + 3*Log(thisSpan) /*- allTheSame * Log(baseSites)*/;
	}
	
	/*
	if (currentBPC == 2)
	{
		fprintf (stdout, sortedScores, "\n");	
	}
	*/
	return -daAICScore;
}



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji)
{
	myLL  = lf_MLES[0];
	myDF  = baseParams+lf_MLES[1]+tbc-1;
	
	if (icOption == 1)
	{
		myAIC = 2*myLL-myDF*Log(baseSites);	
	}
	else
	{
		myAIC = 2*(myLL-myDF*(baseSites/(baseSites-myDF-1)));
	}
	sortedBP = {{-1}};
	
	if (resultProcessingContext==0)
	{
		sortedScores[ji][0] = myAIC;
		if (ji>=0)
		{
			jobPrint = ConvertToPartString (currentPopulation[ji]);
			if (adjustAICScores)
			{
				myAIC	 = adjustAICScore (myLL, sortedBP);
			}
			v 		 = Rows (sortedBP);
			sortedScores[ji][0] = myAIC;
		}
		if (verboseFlag>5)
		{
			fprintf (stdout, "Individual ",ji," AIC-c = ",-myAIC,"\n");
			//fprintf (stdout, ConvertToPartString (currentPopulation[ji]), "\n");
		}
	}
	else
	{
		intermediateProbs[ji][0] = myAIC;	
		if (ji>=0)
		{
			jobPrint = ConvertToPartString (children[ji-populationSize]);
			if (adjustAICScores)
			{
				myAIC	 = adjustAICScore (myLL, sortedBP);
			}
			v = Rows (sortedBP);
			intermediateProbs[ji][0] = myAIC;	
		}
		if (verboseFlag>5)
		{
			fprintf (stdout, "Offspring ",ji," AIC-c = ",-myAIC,"\n");
		}
	}
		
	if (Columns (sortedBP)>1)
	{
		MasterList [_ModelKeyString] = myAIC;
	}
	return fromNode-1;
}



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

compressedString = {{1,1}};

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MakeStringCanonical (someString, dummy)
{
	localPartCount	 = 		(Columns (someString) * Rows (someString) - branchBits) $ bitsPerPart;
	if (localPartCount > 1)
	{
		bpOrdering = decodeIndividual (someString);
		bpOrdering = (bpOrdering)["_MATRIX_ELEMENT_VALUE_*_MATRIX_ELEMENT_COLUMN_+(_MATRIX_ELEMENT_COLUMN_==0)*(_MATRIX_ELEMENT_ROW_-1)"];
		bpOrdering = bpOrdering%1;
		
		reordered  = someString;
		h 		   = branchBits;
		for (mpiNode=1; mpiNode<=localPartCount; mpiNode=mpiNode+1)
		{
			readFrom 		= bpOrdering[mpiNode][0] * bitsPerPart + branchBits;
			if (readFrom < 0)
			{
				fprintf (stdout, someString, "\n", branchBits, "\n", bitsPerPart, "\n", bpOrdering, "\n");
				return someString;
			}
			readTo	 		= readFrom+bitsPerPart;
			for (v=readFrom; v<readTo; v=v+1)
			{
				reordered[h]	= someString[v];
				h=h+1;
			}
		}
		return reordered;
		
	}
	return someString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertToPartStringInt (sortedBP)
{
	bpF 	 = -1;
	bpF2	 = -1;
	
	minPartLength 	  = 1e100;
	
	_ConstraintString = "";
	_ConstraintString * 256;
	_SpliceString	  = "";
	_SpliceString	  * 256;
	
	v 			= Rows (sortedBP); /* number of partitions */
	spliceName 	= (refTopAVL[sortedBP[0][0]])["Name"];

	for (h=1; h<v; h=h+1)
	{
		bpF2 = bppMap[sortedBP[h][1]];
	
		if (h>1)
		{
			_ConstraintString * ",";
			_SpliceString 	  * ",";
		}
		
		_SpliceString * spliceName;
		_ConstraintString * (""+(bpF+1)+"-"+bpF2);		
		curSegLength = bpF2-bpF;

		bpF = bpF2;
		
		if (curSegLength < minPartLength && curSegLength>0)
		{
			minPartLength = curSegLength;
		}
		spliceName 	= (refTopAVL[sortedBP[h][0]])["Name"];
	}
	
	if (bpF2<ds.sites)
	{
		if (v > 1)
		{
			_ConstraintString * (","+(bpF2+1)+"-"+(ds.sites-1));		
			_SpliceString * ("," + spliceName);
		}
		else
		{
			_ConstraintString * (""+(bpF2+1)+"-"+(ds.sites-1));		
			_SpliceString * (spliceName);
		
		}
		curSegLength = ds.sites-bpF2;

		if (curSegLength < minPartLength && curSegLength>0)
		{
			minPartLength = curSegLength;
		}
	}

	_ConstraintString * 0;
	_SpliceString	  * 0;
	_ModelKeyString = _ConstraintString + "\n" + _SpliceString;
	return _ModelKeyString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertToPartString (pString)
{
	return  ConvertToPartStringInt(ConvertToPart (pString));
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertFromPartString (modelString)
{
	s = 0;
	while (modelString[s] != "\n")
	{
		s=s+1;
	}
	bps_s 				= splitOnRegExp(modelString[0][s],"\\,");
	node_s				= splitOnRegExp(modelString[s+1][Abs(modelString)-1],"\\,");
	part_c				= Abs(bps_s);
	matrix_form			= {part_c,2};
	matrix_form			[0][1] = -1;
	matrix_form			[0][0] = nodeNameToAVL[node_s[0]]-1;
	for (s = 1; s < part_c; s=s+1)
	{
		matrix_form[s][0]	= nodeNameToAVL[node_s[s]]-1;
		matrix_form[s][1]	= inverseBppMap[(-1)+(splitOnRegExp(bps_s[s],"\\-"))[0]];
	}
	return encodeIndividual(matrix_form);
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

/* assumes that inVec is a _column_ */
function binaryToDecimal (inVec, inFrom, inStep)
{
	if (useGrayCodes)
	{
		auxVec = {inStep,1};
		auxVec [0] = inVec[inFrom];
		for (_k=1; _k<inStep; _k=_k+1)
		{
			if (auxVec[_k-1])
			{
				auxVec[_k] = 1-inVec[inFrom+_k];
			}
			else
			{
				auxVec[_k] = inVec[inFrom+_k];			
			}
		}
	}
	else
	{
		auxVec = inVec[{{inFrom,0}}][{{inFrom+inStep-1,0}}];
	}
	binaryVec = {1,inStep};
	return    (binaryVec["2^(inStep-1-_MATRIX_ELEMENT_COLUMN_)"]*auxVec)[0];
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

/* assumes that inVec is a _column_ */
function decimalToBinary (inVec&, inFrom, inStep, inValue)
{
	for (_k=inFrom+inStep-1; _k>=inFrom; _k=_k-1)
	{
		inVec[_k] = inValue%2;
		inValue = inValue$2;
	}
	if (useGrayCodes)
	{
		auxVec = inVec;
		for (_k=inFrom+1; _k<=inFrom+inStep-1; _k=_k+1)
		{
			toggle = inVec[_k-1];
			if (toggle == 1)
			{
				auxVec[_k] = 1-inVec[_k];
			}
		}	
		inVec = auxVec;
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function  swapBreakpoints (indString)
{
	ind_matrix = decodeIndividual (indString);
	_hh = Rows(ind_matrix);
	if (_hh > 1)
	{
		swapPoint   			   = Random (1,_hh)$1;
		t						   = ind_matrix[swapPoint][0];
		ind_matrix[swapPoint][0]   = ind_matrix[swapPoint-1][0];
		ind_matrix[swapPoint-1][0] = t;
		
		if (swapPoint > 1)
		{
			from_seg				= ind_matrix[swapPoint-1][1];
		}
		else
		{
			from_seg				= 0;
		}
		if (swapPoint < _hh-2)
		{
			to_seg					= ind_matrix[swapPoint+1][1];
		}
		else
		{
			to_seg					=  Abs(bppMap);
		}
		ind_matrix[swapPoint][1]	=  from_seg + (to_seg-ind_matrix[swapPoint][1]);
		for (_i2=0; _i2<_hh; _i2=_i2+1)
		{
			ind_matrix [_i2][0] = ind_matrix [_i2][0]-1;
			if (ind_matrix[_i2][1] < 0 && _i2)
			{
				return indString;
			}
		}
		return encodeIndividual (ind_matrix);
	}
	return indString;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function decodeIndividual (indString)
{
	_hh 			 = 		0; 
	localPartCount	 = 		(Columns (indString) * Rows (indString) - branchBits) $ bitsPerPart;
	spliceLocations  = 		{localPartCount+1,2};
	
	/* map splice branches to indices */
		
	aBP = binaryToDecimal(indString,0,branchBits)%tbc;
	
	spliceLocations [0][0] = aBP+1;
	spliceLocations [0][1] = -1;
	
	_hh = branchBits;
	for (_i2=0; _i2<localPartCount; _i2=_i2+1)
	{
		spliceLocations[_i2+1][1] = binaryToDecimal(indString,_hh,bppSize)%bppMapSize;
		spliceLocations[_i2+1][0] = binaryToDecimal(indString,_hh+bppSize,branchBits)%tbc+1;
		_hh 				      = _hh + bitsPerPart;
	}

	return spliceLocations;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function encodeIndividual (indMatrix)
{
	_vectDim		 = 		branchBits + (Rows(indMatrix)-1) * bitsPerPart;
	
	_hh 			 = 		0; 
	newIndiv		 = 		{_vectDim,1};
	decimalToBinary	 ("newIndiv",_hh,branchBits, indMatrix[0][0]);
	
	/* map splice branches to indices */
		
	_hh = branchBits;
	for (_i2=1; _i2<Rows(indMatrix); _i2=_i2+1)
	{
		decimalToBinary	("newIndiv",_hh,bppSize, indMatrix[_i2][1]);
		decimalToBinary	("newIndiv",_hh+bppSize,branchBits, indMatrix[_i2][0]);
		_hh 				      = _hh + bitsPerPart;
	}

	return newIndiv;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertToPart (pString)
{
	sortedBP  	= decodeIndividual (pString)%1;
	return 	   sortedBP;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PrepareSampleForARun (sortedBP, is_banned&)
{
	paramSpec = sortedBP;
	is_banned = 0;
	my_size   = Rows(sortedBP);
	
	for (mpiNode=0; mpiNode < my_size; mpiNode = mpiNode+1)
	{
		seq = nodeIDMap[sortedBP[mpiNode][0]];
		if (mpiNode < my_size-1)
		{
			loc = bppMap[sortedBP[mpiNode+1][1]];
		}
		else
		{
			loc = filteredData.sites;
		}
		
		if (hasBannedBP[seq])
		{
			spanEnd = loc;
			if (mpiNode)
			{
				spanStart   = bppMap[sortedBP[mpiNode][1]+1];
			}
			else
			{
				spanStart = 0;
			}
			
			bannedSpan = +(bpStencil ["_MATRIX_ELEMENT_COLUMN_>=spanStart__&&_MATRIX_ELEMENT_COLUMN_<=spanEnd__"]$(bannedBreakpointLocations[seq]));
			if (bannedSpan > 0)
			{
				is_banned = 1;
				if (verboseFlag > 5)
				{
					fprintf (stdout, sortedBP, " is banned!\n");
					fprintf (stdout,  spanStart, "-", spanEnd, "\n");
					fprintf (stdout, bannedBreakpointLocations[seq], "\n");
				}
				break;
			}
		}

		paramSpec[mpiNode][0] = loc;
		paramSpec[mpiNode][1] = seq;
		
		if (debug)
		{
			visitStats [paramSpec[mpiNode][1]][sortedBP[mpiNode+1][1]] = visitStats [paramSpec[mpiNode][1]][sortedBP[mpiNode+1][1]] + 1;
		}
	}
	
	if (paramSpec[0][0] == 0)
	{
		paramSpec[0][0] = 1;
	}
	if (is_banned == 0)
	{
		mpiNode = my_size-1;
		seq = nodeIDMap[sortedBP[mpiNode][0]];
		if (hasBannedBP[seq])
		{
			if ((bannedBreakpointLocations[seq])[filteredData.sites-1])
			{
				is_banned = 1;
			}
			else
			{
				if (mpiNode)
				{
					loc = bppMap[sortedBP[mpiNode][1]+1];
				}
				else
				{
					loc = 0;
				}				
				if ((bannedBreakpointLocations[seq])[loc])
				{
					is_banned = 1;
				}
			}
		}
		paramSpec[mpiNode][1] = seq;
		paramSpec[mpiNode][0] = filteredData.sites;
	}
	
	if (Abs(_validBreakpointStructure))
	{
		if (Abs (_validBreakpointStructure[my_size-1]))
		{
			ss2 = AssembleSubtypeAssignment(ConvertToPartStringInt (sortedBP),0);
			
			if (matchStringToSetOfPatterns(ss2,_validBreakpointStructure[my_size-1]) < 0)
			{
				is_banned = 1;
			}
		}
	}
	
	if (verboseFlag > 5)
	{
		if (is_banned == 0)
		{
			fprintf (stdout, sortedBP, " is OK!\n");
		}
	}
	
	return paramSpec;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample (dummy,jobIndex)
{	
	myAIC	 = MasterList[ConvertToPartString (cString)];

	if (myAIC<0)
	{		
		if (resultProcessingContext==0)
		{
			sortedScores[jobIndex][0] = myAIC;
			if (verboseFlag>5)
			{
				fprintf (stdout, "Individual ",jobIndex," AIC-c = ",-myAIC, "\n");
			}
		}
		else
		{
			intermediateProbs[jobIndex][0] = myAIC;	
			if (verboseFlag>5)
			{
				fprintf (stdout, "Offspring ",jobIndex," AIC-c = ",-myAIC,"\n");
			}
		}	
		return 0;
	}

	paramSpec = PrepareSampleForARun(sortedBP, "is_banned");
	
	if (is_banned)
	{
		lf_MLES = {{-1e10, 3, 0}{0,0,0}};
	}
	else
	{
		lf_MLES = runAModel (paramSpec,branchOptionValue);
	}
	mpiNode = ReceiveJobs (0,jobIndex);
	
	return 0;	
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {totalBitSize,1};
	for (h=0; h<totalBitSize; h=h+1)
	{
		rModel[h] = Random(0,2)$1;
	}
	return MakeStringCanonical(rModel,0);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function IsChildViable (putativeChild)
{
	sampleString = 	ConvertToPartString (putativeChild);
	if (strictViableCheck)
	{
		myAIC		 = MasterList[sampleString];
	}
	else
	{
		myAIC 		 = 0;
	}
	testChild 	 = putativeChild;
	
	myAIC = -populationStrings [sampleString];
	
	
	if (myAIC<(-0.1))
	{
		ordering = {1,stateVectorDimension};
		ordering = Random(ordering["_MATRIX_ELEMENT_COLUMN_"],0);
		
		for (_lc2 = 0; _lc2 < stateVectorDimension && myAIC<(-0.1); _lc2 = _lc2+1)
		{
			testChild [ordering[_lc2]]	 = 1-testChild [ordering[_lc2]];
			testChild					 = MakeStringCanonical(testChild,1);
			sampleString 				 = ConvertToPartString (testChild);
			if (strictViableCheck)
			{
				myAIC 						 = MasterList		   [sampleString];			
			}
			else
			{
				myAIC 						 = 0;
			}	
			if (myAIC > (-0.1))
			{
				myAIC = -populationStrings [sampleString];
			}
		}
	}

	return testChild;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


dataType = 0;
debug	 = 0;

/*

fscanf	 	(stdin, "String", baseFilePath);
fscanf	 	(stdin, "String", referenceTopologyString);
fscanf	 	(stdin, "String", modelDesc);

fscanf		(stdin, "Number", querySequenceID);
fscanf		(stdin, "Number", branchOptionValue);
fscanf		(stdin, "Number", icOption);
fscanf		(stdin, "String", correctModelFile);

*/


if (EXISTING_ALIGNMENT_RUN == 1)
{
	querySequenceID 			= whichSeq;
}
else
{
	DataSet 	  ds 		    = ReadFromString (outputAlignment);
	querySequenceID 			= ds.species-1;
}

referenceTopologyString 	= DATAFILE_TREE;
modelDesc					= "012345";
icOption					= 1; /* BIC */
branchOptionValue			= 0; /* 0-3,1-2 branches */
/*
correctModelFile			= "";
*/

if (Abs(correctModelFile))
{
	ExecuteAFile(correctModelFile);
}



/* END DEBUG BITS */

Tree			referenceTopology = referenceTopologyString;
refTopAVL						  = referenceTopology^0;
_subtypeReductions				  = computeSubtypeReductions ();

if (EXISTING_ALIGNMENT_RUN == 1)
{
	_subtypeAssignmentByNode = {};
	
	tc = TipCount (referenceTopology);
	isALeaf = {};
	for (k=0; k < tc; k = k+1)
	{
		nodeName 	  = TipName (referenceTopology,k);
		_subtypeAssignmentByNode [nodeName] = nodeName;
		isALeaf [nodeName] = 1;
	}
	
	tc = BranchCount (referenceTopology);
	for (k=0; k < tc; k = k+1)
	{
		nodeName 	  = BranchName (referenceTopology,k);
		subtreeAVL	  = referenceTopology[nodeName];
		haveLabels	  = {};
		nodeNames	  = Rows(subtreeAVL);
		for (k2 = 0; k2 < Abs (subtreeAVL); k2 = k2+1)
		{
			if (isALeaf[nodeNames[k2]])
			{
				haveLabels [_subtypeAssignmentByNode[nodeNames[k2]]] = 1;
			}
		}
		haveSubtypes = Rows (haveLabels);
		nodeLabel = haveSubtypes[0];
		for (k2 = 1; k2 < Columns (haveSubtypes); k2 = k2 + 1)
		{
			nodeLabel = nodeLabel + "/" + haveSubtypes[k2];
		}
		_subtypeAssignmentByNode [nodeName&&1] = nodeLabel;
	}
}

/* make a list of valid sequence names */
tbc								  =   TipCount (referenceTopology);
validTaxonNames					  =   {};
CRFGroups						  =   {};
CRFTags							  =   {};


if (Abs(crf_for_grep) == 0)
{
	crf_for_grep					  =   "(.+)_CRF_([0-9]+)(.*)$";
}

for (k=0; k<tbc; k=k+1)
{
	taxonName 				   = TipName(referenceTopology,k)&&1;
	is_crf 					   = taxonName$crf_for_grep;
	if (is_crf[0] >= 0) /* this is a CRF */
	{
		crf_tag = taxonName[is_crf[2]][is_crf[3]];
		if (CRFTags[crf_tag] == 0)
		{
			CRFTags[crf_tag] = 1+Abs(CRFTags);
		}
		crf_tag = 0+CRFTags[crf_tag];
		crf_bit = 0+taxonName[is_crf[4]][is_crf[5]];
		CRFGroups[taxonName] = {{crf_tag__,crf_bit__}};
	}
	validTaxonNames[taxonName] = 1;
}



if (verboseFlag)
{
	fprintf (stdout, "Fitting a baseline nucleotide model\n");
}

dataNames = {ds.species,1};
GetString (h, ds,querySequenceID);
validTaxonNames [h&&1] = 1;

for (k=0; k<ds.species; k=k+1)
{
	GetString (h, ds,k);
	dataNames [k] = h&&1;
}


DataSetFilter filteredData  = CreateFilter (ds,1,"",validTaxonNames[dataNames[speciesIndex]]);

//PATHtosave = ""+MPI_NODE_ID + ".seq";
//fprintf (PATHtosave, CLEAR_FILE, filteredData);

DataSetFilter referenceData = CreateFilter (ds,1,"",validTaxonNames[dataNames[speciesIndex]] && (speciesIndex != querySequenceID));

GetString	  (qsName,ds,querySequenceID);

/* find "informative sites" */

inverseBppMap = {filteredData.sites,1};

for (h=0; h<filteredData.sites; h=h+1)
{
	filterString 			 = Format(h,20,0);
	DataSetFilter siteFilter = CreateFilter (filteredData,1,filterString);
	HarvestFrequencies 		   (f1, siteFilter, 1, 1, 0);
	if (((Transpose(f1))["1"]*f1["_MATRIX_ELEMENT_VALUE_>0"])[0]>1)
	{
		inverseBppMap[h] = Abs(bppMap);
		bppMap[Abs(bppMap)] = h;
	}
	else
	{
		inverseBppMap[h] = -1;
	}
}


bppMapSize		= Abs(bppMap);
bppSize 		= (Log(Abs(bppMapSize))/Log(2)+1)$1;
branchBits		= (Log(Abs(refTopAVL)-2)/Log(2)+1)$1;

bitsPerPart		= bppSize + branchBits;

if (verboseFlag)
{
	fprintf (stdout, "There are ",Abs(bppMap)," potential breakpoints and ", Abs(refTopAVL)-2," candidate branches. ",
						   "\nBit size of the sample (per partition) is ", bitsPerPart,"\n");
}

partCount = 2;
h 		  = Abs(bppMap);

if (h <= partCount)
{
	fprintf (stdout,   "ERROR: \nThere are too few potential break points to support ", partCount-1, " recombination events.\n");
	return 1;
}

maxBPP 	   = Min(h-1,noMoreBPsThan);
maxBPP	   = Min(maxBPP, filteredData.sites$BICMinLength);

ModelTitle = ""+modelDesc[0];
			
rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
paramCount	  = 0;

modelConstraintString = "";

for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
		{
			ModelTitle  = ModelTitle+modelDesc[customLoopCounter2];	
			if (rateBiasTerms[customLoopCounter2] == "1")
			{
				modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter]+":="+rateBiasTerms[customLoopCounter2]+";";
			}
			else
			{
				modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter2]+":="+rateBiasTerms[customLoopCounter]+";";			
			}
			break;
		}
	}
	if (customLoopCounter==customLoopCounter2)
	{
		ModelTitle = ModelTitle+modelDesc[customLoopCounter2];	
	}
}	

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
	analysisPrefix * modelConstraintString;
}

HarvestFrequencies 		  (nucEFV, filteredData, 1, 1, 1);


if (rvChoice)
{
	siteRateClasses = Min(8,Max(2,siteRateClasses $ 1));
	
	if (rvChoice == 1)
	{
		gdDefString = "";
		gdDefString * 1024;
		for (mi=1; mi<siteRateClasses; mi=mi+1)
		{
			gdDefString*("global PS_"+mi+" = 1/"+((siteRateClasses+1)-mi)+";\nPS_"+mi+":<1;\n");
		}
		
		gdDefString*("\n\nglobal RS_1 = .3;\nRS_1:<1;RS_1:>0.000000001;\n");

		for (mi=3; mi<=siteRateClasses; mi=mi+1)
		{
			gdDefString*("global RS_"+mi+" = 1.5;"+"\nRS_"+mi+":>1;RS_"+mi+":<100000;\n");
		} 

		rateStrMx    = {siteRateClasses,1};
		rateStrMx[0] = "RS_1";
		rateStrMx[1] = "1";

		for (mi=3; mi<=siteRateClasses; mi=mi+1)
		{
			rateStrMx[mi-1] = rateStrMx[mi-2]+"*RS_"+mi;
		} 	

		freqStrMx    = {siteRateClasses,1};
		freqStrMx[0] = "PS_1";

		for (mi=1; mi<siteRateClasses-1; mi=mi+1)
		{
			freqStrMx[mi] = "";
			for (mi2=1;mi2<=mi;mi2=mi2+1)
			{
				freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
			}
			freqStrMx[mi] = freqStrMx[mi]+"PS_"+(mi+1);	
		}	

		freqStrMx[mi] = "";
		for (mi2=1;mi2<mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi+")";	


		gdDefString*("\n\nglobal c_scale:="+rateStrMx[0]+"*"+freqStrMx[0]);

		for (mi=1; mi<siteRateClasses; mi=mi+1)
		{
			gdDefString*("+"+rateStrMx[mi]+"*"+freqStrMx[mi]);
		}

		gdDefString*(";\ncategFreqMatrix={{"+freqStrMx[0]);

		for (mi=1; mi<siteRateClasses; mi=mi+1)
		{
			gdDefString*(","+freqStrMx[mi]);
		}

		gdDefString*("}};\ncategRateMatrix={{"+rateStrMx[0]+"/c_scale");

		for (mi=1; mi<siteRateClasses; mi=mi+1)
		{
			gdDefString*(","+rateStrMx[mi]+"/c_scale");
		}

		gdDefString*("}};\n\ncategory c  = ("+siteRateClasses+", categFreqMatrix , MEAN, ,categRateMatrix, 0, 1e25);\n\n");
		gdDefString*0;
		ExecuteCommands (gdDefString);	
	}
	else
	{
		global betaP = 1;
		global betaQ = 1;
		betaP:>0.05;betaP:<85;
		betaQ:>0.05;betaQ:<85;
		category pc = (siteRateClasses-1, EQUAL, MEAN, 
						_x_^(betaP-1)*(1-_x_)^(betaQ-1)/Beta(betaP,betaQ), /* density */
						IBeta(_x_,betaP,betaQ), /*CDF*/
						0, 				   /*left bound*/
						1, 			   /*right bound*/
						IBeta(_x_,betaP+1,betaQ)*betaP/(betaP+betaQ)
					   );
		
		global alpha = .5;
		alpha:>0.01;alpha:<100;
		category c = (siteRateClasses, pc, MEAN, 
						GammaDist(_x_,alpha,alpha), 
						CGammaDist(_x_,alpha,alpha), 
						0 , 
				  	    1e25,
				  	    CGammaDist(_x_,alpha+1,alpha)
				  	 );
			

	}
	
	NucleotideMatrix	 = {{*,c*AC*t,c*t,c*AT*t}
							{c*AC*t,*,c*CG*t,c*CT*t}
							{c*t,c*CG*t,*,c*GT*t}
							{c*AT*t,c*CT*t,c*GT*t,*}};
													
}
else
{
	NucleotideMatrix	 = {{*,AC*t,t,AT*t}{AC*t,*,CG*t,CT*t}{t,CG*t,*,GT*t}{AT*t,CT*t,GT*t,*}};
}

Model nucModel   		= (NucleotideMatrix, nucEFV, 1);

/* check parameter counts */

if (verboseFlag)
{
	fprintf (stdout, "FITTING THE REFERENCE TREE TO OBTAIN INITIAL PARAMETER ESTIMATES\n");
}

Tree	baselineTree	= referenceTopologyString;
LikelihoodFunction		baselineLF = (referenceData, baselineTree);

GetString 				(varList, baselineLF, -1);
baseParams 		   		= Columns(varList["Global Independent"])+3;
perPartParameterCount	= Columns(varList["Local Independent"]);
baseSites		   		= filteredData.sites;

if (icOption == 0)
{
	if (baseParams + 2 * (perPartParameterCount+2) >= baseSites - 1)
	{
		fprintf (stdout,   "ERROR: \nToo few sites for reliable c-AIC inference.\n");
		return 0;
	}
}

if (EXISTING_ALIGNMENT_RUN == 0)
{
	if (DO_NOT_RELOAD_ESTIMATES == 1)
	{
		checkForSavedOptions = 1;
	}
	else
	{
		checkForSavedOptions = referenceAlignmentPath + ".params";
	}
	if (!checkForSavedOptions)
	{
		fscanf (checkForSavedOptions, "NMatrix",baseRes);
		ExecuteAFile (checkForSavedOptions);
		for (k=0; k<Columns(baseRes); k=k+1)
		{
			SetParameter (baselineLF, k, baseRes[0][k]);
		}
		LFCompute (baselineLF,LF_START_COMPUTE);
		LFCompute (baselineLF,res);
		LFCompute (baselineLF,LF_DONE_COMPUTE);	
		baseRes[1][0] = res;
	}
	else
	{
		if (verboseFlag)
		{
			VERBOSITY_LEVEL = 1;
		}
		OPTIMIZATION_PRECISION = 0.01;
		Optimize (baseRes, baselineLF);
		OPTIMIZATION_PRECISION = 0.001;
		if (verboseFlag)
		{
			VERBOSITY_LEVEL = 0;
		}
		if (DO_NOT_RELOAD_ESTIMATES == 0)
		{
			GetString (_lfInfo,baselineLF,-1);
			saveStr = ""; saveStr * 128;
			varList = _lfInfo["Global Independent"];
			for (_gb_idx = 0; _gb_idx < Columns (varList); _gb_idx = _gb_idx + 1)
			{
				ExecuteCommands ("pv="+varList[_gb_idx]);
				saveStr * (varList[_gb_idx] + "=" + pv + ";");
			} 	
			varList = _lfInfo["Local Independent"];
			for (_gb_idx = 0; _gb_idx < Columns (varList); _gb_idx = _gb_idx + 1)
			{
				ExecuteCommands ("pv="+varList[_gb_idx]);
				saveStr * (varList[_gb_idx] + "=" + pv + ";\n");
			} 	
			saveStr * 0;
			fprintf  (checkForSavedOptions, CLEAR_FILE, saveStr, "\nbaseRes=", baseRes, ";\n");
			return 0;
		}
	}
}
else
{
	Optimize (baseRes, baselineLF);
}

if (verboseFlag)
{
	fprintf  (stdout, "Baseline log-L: ", baseRes[1][0], "\n");
}

baseLineBranchNames = BranchName (baselineTree,-1);

if (baseParams>3)
{
	ConstraintString = "";
	ConstraintString *256;
	for (h=0; h<baseParams-3; h=h+1)
	{
		GetString (v,baselineLF,h);
	    ExecuteCommands ("vlue = " + v + ";");
		if (verboseFlag)
		{
	    	fprintf  (stdout, v, " = ", vlue, "\n");
	    }
		ConstraintString * (v+":="+vlue+";\n");
	}
	ConstraintString*0;
	ExecuteCommands (ConstraintString);
}



if (verboseFlag)
{
	fprintf (stdout, "\nDETERMINING THE BEST QUERY SEQUENCE PLACEMENT FOR A SINGLE PARTITION\n");
}

bestAIC 		= 1e100;
bestSBP 		= 0;
MPINodeState 	= {MPI_NODE_COUNT-1,3};

tbc 	= BranchCount (referenceTopology) + TipCount (referenceTopology);
tbc2 	= tbc + 2;
treeS   = tbc+1;

if (debug)
{
	visitStats = {tbc,bppMapSize};
}

DataSetFilter filteredDataBig  = CreateFilter (ds,1,"",validTaxonNames[dataNames[speciesIndex]]);
DataSetFilter filteredDataSM   = CreateFilter (ds,1,"",validTaxonNames[dataNames[speciesIndex]] && speciesIndex != querySequenceID);

GetDataInfo 								  (map1, filteredDataSM);
GetDataInfo 								  (map2, filteredDataBig);

mapBigToSmall 				  				= {Columns(filteredDataBig.site_freqs),1};

seen1 = {};
seen2 = {};

u1 = 1;
u2 = 1;

for (k=0; k<filteredDataBig.sites; k=k+1)
{
	t1 = map1[k];
	t2 = map2[k];
	
	if (seen2[t2] == 0)
	{
		if (seen1[t1] > 0) 
		{	
			mapBigToSmall[u2-1] = seen1[t1]-1;
			seen2		 [t2]   = u2;
			u2 					= u2+1;
		}
		else
		{
			mapBigToSmall[u2-1] = u1-1;
			seen2[t2] = u2;
			seen1[t1] = u1;
			u2 = u2+1;
			u1 = u1+1;
		}
	}
}

rawToUniqueMap 		= {filteredDataBig.sites,1};
rawToUniqueMap		= rawToUniqueMap["mapBigToSmall[map2[_MATRIX_ELEMENT_ROW_]]"];	

TRY_NUMERIC_SEQUENCE_MATCH	 = 1;


	ConstructCategoryMatrix 		(probMatrix, baselineTree);
	/*saveRawValues = probMatrix["Values"];

	minV = 1e100;
	for (k = 0; k < Rows(saveRawValues); k=k+1)
	{
		myMin = saveRawValues[k][-1]*({4,1}["1"]);
		minV = Min(myMin[0], minV);
	}*/
	
	if (rvChoice)
	{
		categoryShifter = Rows (probMatrix["Values"])/siteRateClasses;
	}
	nodeToID = {};
	for (k=0; k<treeS; k=k+1)
	{
		nodeToID[(probMatrix["Nodes"])[k]] = k;
	}

	prepNumericArray = {};
	
	prepNumericArray 		["FILTER_NAMES"]      = {3,1};
	(prepNumericArray 		["FILTER_NAMES"])[0] = "Clade1";
	(prepNumericArray 		["FILTER_NAMES"])[1] = "Clade2";
	(prepNumericArray 		["FILTER_NAMES"])[2] = "Query";
	prepNumericArray 		["FILTER_ARRAYS"]	   = {};


	if (rvChoice)
	{
		NucleotideMatrixF	 = {{*,c*AC__*t,c*t,c*AT__*t}
								{c*AC__*t,*,c*CG__*t,c*CT__*t}
								{c*t,c*CG__*t,*,c*GT__*t}
								{c*AT__*t,c*CT__*t,c*GT__*t,*}};
								
							
	
	}
	else
	{
		NucleotideMatrixF	 = {{*,AC__*t,t,AT__*t}{AC__*t,*,CG__*t,CT__*t}{t,CG__*t,*,GT__*t}{AT__*t,CT__*t,GT__*t,*}};
	}

	Model	nucFixed			 = (NucleotideMatrixF,nucEFV,1);
	Tree tryTree 				 = (Clade1,Clade2,Query);


nodeIDMap 		 = {};
baseBranchValues = {tbc,1};
nodeNameToAVL	 = {};

for (k=1; k<=tbc; k=k+1)
{
	nodeIDMap[k] 					= nodeToID[(refTopAVL[k])["Name"]];
	nodeNameToAVL[(refTopAVL[k])["Name"]] = k;
	ExecuteCommands 				("val = baselineTree." +(refTopAVL[k])["Name"] + ".t;");
	baseBranchValues [nodeIDMap[k]] = val;
}

/* find "banned" ranges for CRFs and pre-populate initial solutions based on CRFs */

GetString 	(allSeqs, filteredData, -1);

bannedBreakpointLocations = {};
hasBannedBP				  = {};
CRFStructure			  = {};

gapStructureRegex		  = "(\\-){"+_contigGapThresh+"}\\-+";
stencil     			  = {1,filteredData.sites};

for (h = 0; h < filteredData.species; h = h+1)
{
	crf_group = CRFGroups[allSeqs[h]];
	
	if (Rows(crf_group))
	{
		GetDataInfo (qS, filteredData, h);
		if (Abs(CRFStructure[crf_group[0]]) == 0)
		{
			CRFStructure[crf_group[0]] = {};
		}
		
		haveGaps    = qS||gapStructureRegex;
		seq = (-1) + nodeNameToAVL[allSeqs[h]];
		
		if (haveGaps[0]>=0)
		{
			bannedSites = stencil;
			for (k = Rows (haveGaps)-1; k > 0 ; k=k-2)
			{
				bannedSites = bannedSites + stencil["_MATRIX_ELEMENT_COLUMN_>=haveGaps[k__-1]&&_MATRIX_ELEMENT_COLUMN_<=haveGaps[k__]"]; 
				if (haveGaps[k] != filteredData.sites - 1) // NOT the last site
				{
					loc = 0 + inverseBppMap[haveGaps[k]];
					(CRFStructure[crf_group[0]])[Abs(CRFStructure[crf_group[0]])] = {{seq__,loc__+1}};
				}
				else
				{
					if (_bannedEndAdjust)
					{
						bannedSites = bannedSites - stencil["_MATRIX_ELEMENT_COLUMN_>=haveGaps[k__-1]&&_MATRIX_ELEMENT_COLUMN_<haveGaps[k__-1]+_contigGapThresh-1"];
					}
				}
				
				if (haveGaps[k-1] == 0)
				{
					if (_bannedStartAdjust)
					{
						bannedSites = bannedSites - stencil["_MATRIX_ELEMENT_COLUMN_<=haveGaps[k__]&&_MATRIX_ELEMENT_COLUMN_>=haveGaps[k__]-_contigGapThresh+1"];
					}
				
				}
			}
			if (haveGaps[0] > 0)
			{
				(CRFStructure[crf_group[0]])[Abs(CRFStructure[crf_group[0]])] = {{seq__,1}};
			}
			if (verboseFlag)
			{
				fprintf (stdout, allSeqs[h], " has ", +bannedSites, " banned sites\n");
			}
			seq 							= nodeIDMap[seq+1];
			bannedBreakpointLocations [seq] = bannedSites;
			hasBannedBP				  [seq] = 1;
		}
		else
		{
			(CRFStructure[crf_group[0]])[0] = {{seq__,1}};
		}
	}
}


CRFMosaics				  = {};
if (Abs(_crfEquiv))
{
	for (h = 1; h <= Abs (CRFStructure); h=h+1)
	{
		thisStructure = CRFStructure[h];
		thisCRF		  = _subtypeAssignmentByNode[(refTopAVL[(thisStructure[0])[0]+1])["Name"]];
		k = Abs(thisStructure);
		
		toSort	   = {k,2};
		toSort	   = toSort["_MATRIX_ELEMENT_ROW_"];

		for (ci = 0; ci < k; ci = ci+1)
		{
			toSort[ci][0] = inverseBppMap[(thisStructure[ci])[1]-1];
		}
		toSort = toSort%0;
		
		thisMosaic = {k,2};
		for (ci2 = 0; ci2 < k; ci2 = ci2+1)
		{
			ci = toSort[ci2][1];
			thisMosaic[ci2][0] = _crfEquiv[(refTopAVL[(thisStructure[ci])[0]+1])["Name"]];
			if (thisMosaic[ci2][0] == 0)
			{
				break;
			}
			thisMosaic[ci2][1] = ""+inverseBppMap[(thisStructure[ci])[1]-1];
		}
		if (ci2 == k)
		{
			if (Abs(CRFMosaics[thisCRF]) == 0)
			{
				(CRFMosaics[thisCRF]) = {};
			}
			(CRFMosaics[thisCRF])[Abs(CRFMosaics[thisCRF])] = thisMosaic;
		}
	}
}




/* now add banned breakpoints for internal nodes */
bpStencil     = {1,filteredData.sites}["1"];

for (k=1; k<=tbc; k=k+1)
{
	cc = Abs((refTopAVL[k])["Children"]);
	if (cc>0)
	{
		bannedSites = bpStencil;
		hasBanned	= 1;
		for (ci = 0; ci < cc; ci = ci+1)
		{
			cci = nodeIDMap[((refTopAVL[k])["Children"])[ci]];
			if (hasBannedBP[cci])
			{
				bannedSites	  = bannedSites$bannedBreakpointLocations[cci];
			}
			else
			{
				hasBanned = 0;
				break;
			}	
		}
		
		if (hasBanned)
		{
			bannedSiteCount = (bpStencil * Transpose(bannedSites))[0];
			if (bannedSiteCount)
			{
				if (verboseFlag)
				{
					fprintf (stdout, (refTopAVL[k])["Name"], " has ", bannedSiteCount, " banned sites\n");
				}
				cci = nodeIDMap[k];
				hasBannedBP[cci] = 1;
				bannedBreakpointLocations[cci] = bannedSites;
			}
		}
	}
}

crf_count = Abs (CRFStructure);

for (k = 1; k <= crf_count; k=k+1)
{
	thisCRF 	   = CRFStructure [k];
	h 			   = Abs (thisCRF)-1;
	if (h>0)
	{
		individual	   = {h+1,2};
		for (v = 0; v <= h; v=v+1)
		{
			individual[v][0] = (thisCRF[v])[0];
			individual[v][1] = (thisCRF[v])[1]-1;
		}
		if (Abs(initialPops[h])==0)
		{
			initialPops[h] = {};
		}
		(initialPops[h]) + encodeIndividual(individual%1);
		
	}
}


totalAddedCRFGroups = 0;

if (Abs(_additionalStartingBreakpoints))
{
	upto = Columns (allSeqs);
	
	for (k = 0; k < Abs(_additionalStartingBreakpoints); k+=1)
	{
		mxDef 		= _additionalStartingBreakpoints[k];
		howManyABP	= Columns (mxDef);
		
		bySegment = {};
		
		for (whichBP = 0; whichBP < howManyABP; whichBP += 1)
		{
			bySegment + {};
			for (s = 0; s < upto; s+=1)
			{
				if ((allSeqs[s]$mxDef[whichBP])[0]>=0)
				{
					(bySegment      [whichBP]) + allSeqs[s];
				}
			}
		}
		
		totalCombinations = 1;
		for (whichBP = 0; whichBP < howManyABP; whichBP += 1)
		{
			totalCombinations = totalCombinations * Abs (bySegment[whichBP]);
		}
		
		if (totalCombinations)
		{
			addedThisTime = 0; totalTried = 0;
			while (totalAddedCRFGroups < _randomizedCRFGroups && addedThisTime < totalCombinations && totalTried < 1000)
			{
				individual	   = {howManyABP,2};	
				for (whichBP = 0; whichBP < howManyABP; whichBP += 1)
				{
					recStructure =CRFStructure[
					(CRFGroups[(bySegment[whichBP])[Random(0,Abs(bySegment[whichBP])-0.00000001)$1]])[0]
					];
					if (Abs(recStructure) != 1)
					{
						break;
					}
					individual[whichBP][0] = (recStructure[0])[0];
					individual[whichBP][1] = (recStructure[0])[1]-1;	
				}	
				totalTried 			+= 1;
				if (whichBP  == howManyABP)
				{
					individual = encodeIndividual(individual % 1);
					//fprintf (stdout, ConvertToPartString(individual), "\n");
					
					if (Abs(initialPops[howManyABP-1])==0)
					{
						initialPops[howManyABP-1] = {};
					}
					(initialPops[howManyABP-1]) + individual;

					addedThisTime 		+= 1;
					totalAddedCRFGroups += 1;
				}
			}
		}
	}
}

//fprintf (stdout, initialPops, "\n");

//fprintf (stdout, "\nCRFGroups\n",CRFGroups,"\n\n",bySegment, "\n", CRFStructure, "\n");
//assert (0);



sequenceProbArray = {filteredDataBig.sites, 4};

for (s = 0; s<filteredDataBig.sites; s=s+1)
{
	sm = map2[s];
	GetDataInfo (charInfo, filteredDataBig, querySequenceID,sm);
	for (v = 0; v<4; v=v+1)
	{
		sequenceProbArray[s][v] = charInfo[v];
	}
}

ExecuteAFile	 ("ModelFitter.ibf");

UseModel (USE_NO_MODEL);

GetString (originalSeqName,ds, ds.species-1);
cleanedUpSeqName = "QUERY";



for (_h=1; _h<=tbc; _h=_h+1)
{
	/*refTopAVL			  = referenceTopology^0;*/
	/*if (verboseFlag)
	{
		fprintf (stdout, "[BRANCH ", (refTopAVL[_h])["Name"],"]\n");
	}*/

	treeString = Format (referenceTopology,1,0);
	Topology suppliedTree = treeString;
	suppliedTree + {"NAME": cleanedUpSeqName, 
					"PARENT": "_SPLICE_NODE_",
					"WHERE": (refTopAVL[_h])["Name"]
					};
	treeStringCache[_h]   = Format (suppliedTree,1,0);

	if (hasBannedBP[nodeIDMap[_h]] == 0)
	{
		paramSpec = {{filteredData.sites,nodeIDMap[_h]}};
		outRes    = runAModel (paramSpec,branchOptionValue);
		
		myDF  = perPartParameterCount + baseParams + 2;
		if (icOption == 0)
		{
			myAIC = -2*(outRes[0]-myDF*(baseSites/(baseSites-myDF-1)));
		}
		else
		{
			myAIC = -2*outRes[0]+myDF*Log(baseSites);	
		}
		thisSample   				 = {branchBits,1};
		decimalToBinary ("thisSample",0,branchBits, _h-1);
		ConvertToPartString			   (thisSample);
		MasterList [_ModelKeyString] = -myAIC;
	
	
		if (myAIC < bestAIC)
		{
			overallBestFound	= thisSample;
			bestAIC 			= myAIC;
			bestSBP 			= _h;
			baseLL				= outRes[0];
		}
	}
	else
	{
		/*if (verboseFlag)
		{
			fprintf (stdout, "Skipping branch ",(refTopAVL[_h])["Name"], "\n"); 
			fprintf (stdout, bannedBreakpointLocations[nodeIDMap[_h]],"\n");
		}*/
	}
	SetParameter (STATUS_BAR_STATUS_STRING, "Initial sister lineage scan ("+_h+"/"+tbc+" done)",0);

}

totalBitSize   = branchBits;

branchAttach = (refTopAVL[bestSBP])["Name"];
if (Abs(correctModel))
{
	cmp				 = Abs(correctModel);
	correctModelList = {cmp,2};
	correctModelKeys = Rows (correctModel);
	
	correctBP = {cmp,2};
	correctBP[0][1] = -1;
	for (k=0; k<cmp; k=k+1)
	{
		correctModelList[k][0] = 0+correctModelKeys[k];
		correctModelList[k][1] = nodeToID[correctModel[correctModelKeys[k]]];
		correctBP[k][0]		   = correctModelList[k][1];
		if (k < cmp-1)
		{
			correctBP[k+1][1]         = correctModelList[k][0];
		}
	}
	
	lf_MLES = runAModel (correctModelList,branchOptionValue);
	myLL  = lf_MLES[0];
	myDF  = baseParams+lf_MLES[1]+tbc-1;	

	if (adjustAICScores)
	{
		correctModelAIC = bppMap;
		bppMap = {1,baseSites}["_MATRIX_ELEMENT_COLUMN_-1"];
		myAIC = adjustAICScore (myLL, correctBP);
		bppMap = correctModelAIC;
	}
	else
	{
		if (icOption == 1)
		{
			myAIC = 2*(myLL)-myDF*Log(baseSites);	
		}
		else
		{
			myAIC = 2*(myLL-myDF*(baseSites/(baseSites-myDF-1)));
		}
	}
	correctModelAIC = -myAIC;
	if (runInMPIMode == 0)
	{
		fprintf (stdout, "Correct model IC:", correctModelAIC, "\n");
	}
	else
	{
		if (stepByStepLogging)
		{
			fprintf (MESSAGE_LOG, "Correct model IC:", correctModelAIC, "\n", correctBP, "\n");
		}	
	}
}

if (verboseFlag)
{
	fprintf 	(stdout,"Best splice at ", branchAttach, ". LogL = ", baseLL, ". IC = ", bestAIC," .\n");
}
if (stepByStepLogging)
{
	fprintf 	(MESSAGE_LOG,"\nBest splice at ", branchAttach, ". LogL = ", baseLL, ". DF = ", myDF ,". IC = ", bestAIC," .\n");
}	

currentSubtypeAssignment = _subtypeAssignmentByNode[branchAttach];

if (runInMPIMode == 0)
{
	fprintf (stdout, "Initial subtype assignment: ", currentSubtypeAssignment, "\n");
}
else
{
	returnAVL = {};
	if (stepByStepLogging)
	{
		if (Type(currentSubtypeAssignment) == "String")
		{
			bestAssignment 		= AssembleSubtypeAssignment (overallBestFound,1);
			
	
			fprintf (MESSAGE_LOG, "\nStep :",0,
							 "\nIC: ", bestAIC,
							 "\nPredicted subtype     : ", bestAssignment);
		}
	}		
}

currentPopulation  = {};
sortedScores	   = {populationSize,2};

sortedScores[0][0] 		= -bestAIC;
sortedScores[0][1] 		= 0;
crapAIC 		   		= -sortedScores[0][0];
startTimer		  		= Time (1);
MPINodeState 			= {MPI_NODE_COUNT-1,3};


currentBEST_IC 			 = crapAIC;
ibfPath 				 = "GA_CHC_2_new.ibf";

current_BPP_done_counter = 0;

lastImprovedBPC = 0;

DataSetFilter allData  = CreateFilter (ds,1,"",validTaxonNames[dataNames[speciesIndex]]);

lastIMPBPC = 0;




for (currentBPC = startBPC; currentBPC < maxBPP; currentBPC = currentBPC + 1)
{
	totalModelCounter 		 = 1;
	kf 						 = 1;
	
	for (k=1; k <= partCount; k=k+1)
	{
		totalModelCounter = totalModelCounter * (Abs(bppMap)-k+1)*tbc;
		kf 				  = kf * k;
	} 
	totalModelCounter = totalModelCounter / kf;

	current_BPP_done_counter = Abs (MasterList);
	partCount 				 = currentBPC;
	totalBitSize 			 = branchBits + partCount*bitsPerPart;
	stateVectorDimension 	 = totalBitSize;
	
	subPopulations			= {};
	subTopScores			= {};
	saveCurrentpopulation   = currentPopulation;
	saveStopCriterion		= stoppingCriterion;
	saveLocalMutInterval	= localMutationInterval;
	
	for (subpop = 0; subpop < runSubpopulations+1; subpop = subpop + 1)
	{
		if (subpop < runSubpopulations)
		{
			stoppingCriterion     = saveStopCriterion$3;
			localMutationInterval = stoppingCriterion$2;
			currentPopulation = saveCurrentpopulation;
		}
		else
		{
			stoppingCriterion     = saveStopCriterion;
			localMutationInterval = saveLocalMutInterval;
			if (verboseFlag)
			{
				fprintf (stdout, subTopScores,"\n");
			}			
		}
	
		if (subpop == runSubpopulations && runSubpopulations)
		{
			currentPopulation = {};
			populationStrings = {};
			mixingStep		  = populationSize $ runSubpopulations;
			offset			  = 0;
			for (indCounter = 0; indCounter < runSubpopulations; indCounter = indCounter + 1)
			{
				if (indCounter == runSubpopulations-1)
				{
					mixingStep = populationSize-offset;
				}
				
				for (mixingIdx = 0; mixingIdx < mixingStep; mixingIdx = mixingIdx+1)
				{
					currentPopulation[offset+mixingIdx] = (subPopulations[indCounter])[mixingIdx];
					populationStrings[currentPopulation[offset+mixingIdx]] = 1;
				}
				offset = offset + mixingStep;
			}
		}
		else
		{
			/*if (currentBPC == startBPC)
			{
				populationStrings = {};
				for (k=0; k<populationSize; k=k+1)
				{
					currentPopulation[k] = IsChildViable (SpawnRandomString(rateClassesCount));		
					populationStrings[sampleString] = 1;
				}
			}
			else*/
			{
				populationStrings = {};
				k = 0;
				/*if (subpop%2 == 0 || runSubpopulations == 0)
				{*/
				if (currentBPC > startBPC)
				{
					for (k=0; k<populationSize$3; k=k+1)
					{
						cString			 = addOne (currentPopulation[populationSize-1],Random(0,1)>0.5);
						currentPopulation[k] = IsChildViable(cString);
						populationStrings[sampleString] = 1;
						if (verboseFlag > 1)
						{
							fprintf (stdout, sampleString, "\n");
						}
					}
					
					cString			 = addOne (currentPopulation[populationSize-1],Random(0,1)>0.5);
					currentPopulation[populationSize$3] = IsChildViable(cString);
					populationStrings[sampleString] = 1;
					cString			 = addOne (currentPopulation[populationSize-1],Random(0,1)>0.5);
					currentPopulation[populationSize$3+1] = IsChildViable(cString);
					populationStrings[sampleString] = 1;
					
					k=populationSize$3+2;
				}
	
	
				availablePresets = Abs(initialPops[currentBPC]);
				h = populationSize-k;
				if (availablePresets && h )
				{
					runTo = Min(availablePresets,h);
					stencil = {1,availablePresets};
					stencil = stencil["_MATRIX_ELEMENT_COLUMN_"];
					stencil = (Random (stencil,0))[{{0,0}}][{{0,runTo-1}}];
					
					
					if (verboseFlag)
					{
						fprintf (stdout, "\nRunning CRF presets (",currentBPC," breakpoints)\n", stencil, "\n");
					}
					
					for (h2 = 0; h2 < runTo; h2=h2+1)
					{
						currentPopulation[k+h2] = (initialPops[currentBPC])[stencil[h2]];
						if (verboseFlag)
						{
							fprintf (stdout, "\n***\n", h2, "\n", ConvertToPartString (currentPopulation[k+h2]), "\n");
						}
						populationStrings[currentPopulation[k+h2]] = 1;
					}
					k = k+runTo;
				}
				/*}*/
				for (; k<populationSize; k=k+1)
				{
					currentPopulation[k] = IsChildViable (SpawnRandomString(rateClassesCount));
					populationStrings[currentPopulation[k]] = 1;
				}
			}
			
		}	
	
		children = {};
		ExecuteAFile (ibfPath);
		
		if (subpop < runSubpopulations)
		{
			subPopulations[subpop] = {};
			for (k=populationSize-1; k>=0; k=k-1)
			{
				(subPopulations[subpop])[populationSize-1-k] = currentPopulation[sortedScores[k][1]];
			}
			subTopScores  [subpop] = sortedScores[populationSize-1][0];
		}
	}
			
	kf						 = -sortedScores[populationSize-1][0];
	
	if (tryBreakpointFusion && currentBPC > 1)
	{
		saveCurrentBestAIC		 = -sortedScores[populationSize-1][0];
		if (currentBEST_IC > saveCurrentBestAIC)
		{
			if (verboseFlag)
			{
				fprintf (stdout, "\nTRYING BREAKPOINT FUSION\n");
			}		
			saveOldBestAIC			 = currentBEST_IC;
			saveCurrentpopulation    = currentPopulation;
			saveStopCriterion		 = stoppingCriterion;
			saveLocalMutInterval	 = localMutationInterval;
			saveSortedScores		 = sortedScores;
			
			
			currentBPC				 = currentBPC - 1;
			current_BPP_done_counter = Abs (MasterList);
			partCount 				 = currentBPC;
			totalBitSize 			 = branchBits + partCount*bitsPerPart;
			stateVectorDimension 	 = totalBitSize;
			children				 = {};
			fuseFromMe				 = currentPopulation[populationSize-1];
			currentBEST_IC			 = saveCurrentBestAIC;
			
			for (k=0; k<populationSize; k=k+1)
			{
				currentPopulation[k] = fuseBreakpoints(fuseFromMe);
			}
			
			stoppingCriterion     = saveStopCriterion$3;
			localMutationInterval = stoppingCriterion$2;
			ExecuteAFile (ibfPath);
			
			kf						 = -sortedScores[populationSize-1][0];
			
			stoppingCriterion	  = saveStopCriterion;
			localMutationInterval = saveLocalMutInterval;
			
			if (saveCurrentBestAIC <= kf) /* fusion failed */
			{
				currentBPC 				 = currentBPC + 1;
				partCount 				 = currentBPC-1;
				totalBitSize 			 = branchBits + partCount*bitsPerPart;
				stateVectorDimension 	 = totalBitSize;
				currentBEST_IC			 = saveOldBestAIC;
				kf						 = saveCurrentBestAIC;
				currentPopulation		 = saveCurrentpopulation;
				sortedScores			 = saveSortedScores;
			}
		}		
	}
		
	if (stepByStepLogging)
	{
		bestAssignment 		= AssembleSubtypeAssignment (currentPopulation[populationSize-1],1);
		bestPC				= GetBreakpoints (ConvertToPartString (currentPopulation[populationSize-1]));
		

		fprintf (MESSAGE_LOG, "\nStep :",currentBPC,
						 "\nIC: ", kf,
						 "\nIC-best:", kf-currentBEST_IC,  
						 "\nPredicted subtype     : ", bestAssignment, 
						 "\nPredicted breakpoints : ", bestPC);
						 
		if (Abs(correctModel))
		{
			fprintf (MESSAGE_LOG,   "Correct-this         : ", correctModelAIC-kf, "\n",
							        "Correct-best         : ", correctModelAIC-Min(kf,currentBEST_IC), "\n");
		}
	}
	

	if (currentBEST_IC > kf || currentBPC - lastIMPBPC < 2)
	{
		if (currentBEST_IC > kf)
		{
			lastImprovedBPC  = currentBPC;
			overallBestFound = currentPopulation[populationSize-1];
			lastIMPBPC	 = currentBPC;
		}
		currentBEST_IC = Min(kf,currentBEST_IC);
		stoppingCriterion 	  = stoppingCriterion + convergenceBump;
		populationSize	  	  = populationSize    + populationSizeBump;
		localMutationMultiplier = Max(localMutationMultiplier-localMutationMultiplierDecrease,localMutationMultiplierMin);
		localMutationInterval	= Min(stoppingCriterion-1,stoppingCriterion * localMutationMultiplier $ 1);
		newScores		  = {populationSize,2};
		for (k=0; k<populationSize-populationSizeBump; k=k+1)
		{
			newScores[k][0] = sortedScores[k][0];
			newScores[k][1] = sortedScores[k][1];
		}
		for (;k<populationSize; k=k+1)
		{
			currentPopulation [k] = currentPopulation[populationSize-populationSizeBump-1];
		}	
		sortedScores = newScores;
	}
	else
	{
		break;
	}

	if (Time(1) - scuealTimer > timeLimit)
	{
		prematureTermination = 1;
		break;
	}

	if (verboseFlag)
	{
		fprintf (stdout, "\nStep-up: ",
						 "\n\tPopulation size: ", populationSize,
						 "\n\tCovergence criterion: ", stoppingCriterion,
						 "\n\tRate variation option: ", rvChoice,
						 "\n\tMaximum number of breakpoints: ", noMoreBPsThan,
						 "\n\tMinimum fragment length: ", BICMinLength,
						 "\n\tLocal Mut. Interval: ", localMutationInterval,
						 "\n");
	}
}


masterKeys     		= Rows(MasterList);
totalTried	   		= Rows(masterKeys)*Columns(masterKeys);
akaikeWeights  		= {totalTried,2};

modelSupportByType	= {};
bestByModelType		= {};
bestAWModelType		= {};
totalSum	   		= 0;
recombSupport		= 0;
intraSupport		= 0;

credibleKeys		= {};

for (_mc=totalTried-1; _mc>=0; _mc=_mc-1)
{
	aKey 								= masterKeys[_mc];
	anIC 								= -MasterList[aKey];
	aw	 								= Exp((currentBEST_IC-anIC)*0.5);
	if (aw > 0.001 / totalTried)
	{	
		_cmc								= Abs (credibleKeys);
		credibleKeys[_cmc]					= aKey;
		totalSum 							= totalSum + aw;
		akaikeWeights[_cmc][0]   			= aw;
		thisModelType 						= AssembleSubtypeAssignment (aKey,0);
		modelSupportByType[thisModelType] 	= modelSupportByType[thisModelType] + aw;
		if (classType)
		{
			recombSupport = recombSupport + aw;
			if (classType == 2)
			{
				intraSupport = intraSupport + aw;
			}
		}
	
		if (aw > bestAWModelType[thisModelType])
		{
			bestAWModelType[thisModelType] = aw;
			bestByModelType[thisModelType] = _cmc;
		}		
	}	
}

if (verboseFlag)
{
	fprintf (stdout, Abs (credibleKeys), "/", totalTried, " credible models\n");
}

totalTried   = Abs (credibleKeys);
masterKeys   = credibleKeys;
credibleKeys = 0;

recombSupport = recombSupport / totalSum;
intraSupport  = intraSupport  / totalSum;

h = Abs (modelSupportByType);
v = Rows (modelSupportByType);

typesWithSupport = {};


for (k = 0; k < h; k = k + 1)
{
	if (modelSupportByType[v[k]]/totalSum >= 0.05)
	{
		typesWithSupport[v[k]] = modelSupportByType[v[k]]/totalSum;
	}
}

h 				   = Abs	(typesWithSupport);
v				   = Rows 	(typesWithSupport);

wellSupported = {h,2};
for (k = 0; k < h; k = k+1)
{
	wellSupported[k][0] = k;
	wellSupported[k][1] = typesWithSupport[v[k]];
}

wellSupported    		= wellSupported%1;
bestAssignment   		= v[wellSupported[h-1][0]];
overallBestFound 		= ConvertFromPartString(masterKeys[bestByModelType[bestAssignment]]);

if (currentBEST_IC > 1e10)
{
	fprintf (stdout,   "ERROR: \nAnalysis options and recombinants in the reference alignment created a situation where no valid models could be found\n");
	if (runInMPIMode == 1)
	{
		returnAVL = {};
	}
	return 1;
}

bestAssignment 			= AssembleSubtypeAssignment (overallBestFound,1);
bestAssignmentSimple	= AssembleSubtypeAssignmentSimple (overallBestFound,1);
runAModel			 	(PrepareSampleForARun(sortedBP, "is_banned"), branchOptionValue);
bestModelIC				= 0-MasterList[ConvertToPartString(overallBestFound)];
bestModelBL				= modelBLEstimates;
bestPC					= GetBreakpoints (ConvertToPartString (overallBestFound));
matchingSum	   			= 0;

breakPointSupport 		= {filteredData.sites,1};
otherBPSupport			= {filteredData.sites,1};
branchSupport			= {filteredData.sites,tbc};
stencil		 			= {filteredData.sites,1};

for (_mc=totalTried-1; _mc>=0; _mc=_mc-1)
{
	aKey 					= masterKeys[_mc];
	anIC 					= -MasterList[aKey];
	aw	 					= Exp((currentBEST_IC-anIC)*0.5);
	thisModelType 			= AssembleSubtypeAssignment (aKey,0);
	
	aBPList = GetBreakpoints (aKey);
	if (bestAssignment == thisModelType)
	{
		matchingSum 		= matchingSum + aw;
		for (_s2 = 0; _s2 < Rows(aBPList); _s2 = _s2+1)
		{
			breakPointSupport[aBPList[_s2]-1] = breakPointSupport[aBPList[_s2]-1] + aw; 
		}
	}
	else
	{
		for (_s2 = 0; _s2 < Rows(aBPList); _s2 = _s2+1)
		{
			otherBPSupport[aBPList[_s2]-1] = otherBPSupport[aBPList[_s2]-1] + aw; 
		}	
	}
	if (aw/totalSum > 0.1/totalTried)
	{
		splits   = splitOnRegExp (aKey[0][whereisTheSpace-1], ",");
		splits2  = splitOnRegExp (aKey[whereisTheSpace+1][Abs(aKey)-1], ",");
		bitCount = Abs(splits2);
		for (k2 = 0; k2 < bitCount; k2=k2+1)
		{
			bp    = splitOnRegExp(splits[k2],"-");
			bps   = 0+bp[0];
			bpe   = 0+bp[1];
			bn	  = splits2[k2];
						
			toAdd = stencil["(_MATRIX_ELEMENT_ROW_>=bps__)*(_MATRIX_ELEMENT_ROW_<=bpe__)"];
			bi    = nodeNameToAVL[bn]-1;
			for (idx = 0; idx < filteredData.sites; idx = idx + 1)
			{
				branchSupport[idx][bi] = branchSupport[idx][bi] + toAdd[idx] * aw/totalSum;
			}
			bnp   = bn;
		}
	}
}

otherBPSupport			 = (otherBPSupport+breakPointSupport)*(1/totalSum);
breakPointSupport   	 = breakPointSupport * (1/matchingSum);
significantBranches		 = {};

branchThresh			 = 0.05;
for (_s2 = 0; _s2 < tbc; _s2 = _s2 + 1)
{
	if (Max(branchSupport[-1][_s2],0)>=branchThresh)
	{
		significantBranches[_s2] = 1;
	}
}
stencil 			= branchSupport["significantBranches[_MATRIX_ELEMENT_COLUMN_]"];
stencil	 			= branchSupport[stencil];

_s					= Abs(significantBranches);
sigBrSupport 		= {filteredData.sites,_s+2};
labeler				= {_s+1,2};

h 		  = Rows (significantBranches);
plotRGB   = {_s,3};

//plotColors = _hyDefaultPSColors;
defColorCount = Rows(_hyDefaultPSColors);
plotColors = {_s+1,3};
branchColorMapAux = {_s,1};

for (v=0; v<_s; v=v+1)
{
	_s2 = (refTopAVL[1+h[v]])["Name"];
	labeler[v][0] = reduceSubtype(_s2);
	branchColorMapAux[v] = _s2;
	plotColors[v][0] = _hyDefaultPSColors[v%defColorCount][0];
	plotColors[v][1] = _hyDefaultPSColors[v%defColorCount][1];
	plotColors[v][2] = _hyDefaultPSColors[v%defColorCount][2];
	if (Abs(labeler[v][0])>10)
	{
		labeler[v][0] = (labeler[v][0])[0][9] + "...";
		labeler[v][1] = "";
	}
}
branchColorMap = stringMatrixToAVL("branchColorMapAux");


labeler[v][0] = "Breakpoints";
labeler[v][1] = "Impulse";
plotColors[v][0] = 0; plotColors[v][1] = 0; plotColors[v][2] = 0;

_s2 	= 0;
for (h=0; h<filteredData.sites; h=h+1)
{
	sigBrSupport[h][0] = h+1;
	for (v=1; v<=_s; v=v+1)
	{
		sigBrSupport[h][v] = 100*stencil[_s2];
		_s2 			   = _s2 + 1;
	}
	sigBrSupport[h][v] = 100*otherBPSupport[h];
}

									  
branchSupport = 	0; 	stencil = 0;


LoadFunctionLibrary ("TreeTools.ibf");

refTopAVL				=   referenceTopology^0;
mrca 		  		    = _findMRCA (refTopAVL, avlKeysToMatrix (significantBranches)%0 + 1);

totalPlotWidth			= 600;
ExecuteCommands			(ExportAMatrix(overallBestFound,bestModelBL,(refTopAVL[mrca])["Name"]));
exportedPartCount		= Rows (sortedBp);
Export					(lfExportString,lf_export);


SetDialogPrompt ("Save the results to:");

postScriptOut = ""; postScriptOut * 128;
psTranslate   = {{0,310+treeImageHeight__}};
postScriptOut * ("0 "+(treeImageHeight+10)+" translate\n");
postScriptOut * SimpleGraph("sigBrSupport", {{1,filteredData.sites},{0,100}}, "Times-Roman", {{totalPlotWidth,300,8}}, plotColors, {{"Inferred Mosaic", "Nucleotide", "Model-averaged support"}},labeler, 2);
postScriptOut * "0 300 translate\n";
								  


_mc 		  = Abs(typesWithSupport) - 1;
ibp			  = Rows(bestPC);
ibpt		  = (ibp+1)*18;

if (ibpt>20)
{
	postScriptOut * ("50 "+ibpt+" translate\n");
	psTranslate   = psTranslate + {{50,ibpt}};
}

if (runInMPIMode == 0)
{
	fprintf (stdout, "\nPredicted subtype                      : ", bestAssignment,
			 "\nPredicted simple subtype               : ",  bestAssignmentSimple,
				 "\nModel averaged support                 : ", Format(matchingSum/totalSum*100,8,4), "%",
				 "\nSupport for recombination              : ", Format(recombSupport*100,8,4), "%",
				 "\nSupport for intra-subtype recombination: ", Format(intraSupport*100,8,4), "%\n"
				 );
}

if (_mc > 0)
{
	if (runInMPIMode == 0)
	{
		fprintf (stdout, 	"\nThere are ", _mc, " other mosaic types with model-averaged support over 5%");
	}
	keys 				= Rows(typesWithSupport);
	summaryMatrix 		= {_mc+1,2};
	summaryMatrix		[0][0] = "Alternative subtype";
	summaryMatrix		[0][1] = "Model averaged support";
	
	h					= 1;
	for (_k = _mc; _k >= 0; _k=_k-1)
	{
		if (keys[wellSupported[_k][0]] != bestAssignment)
		{
			if (runInMPIMode == 0)
			{
				fprintf (stdout, "\n\tAlternative subtype        :", keys[wellSupported[_k][0]],
							 	"\n\tModel averaged support     :", Format(wellSupported[_k][1]*100,8,4), "%");
			}			 
			summaryMatrix[h][0] = keys[wellSupported[_k][0]];					 
			summaryMatrix[h][1] = Format(wellSupported[_k][1]*100,0,4) + "\\%";					 
			h = h + 1;
		}
	}
	if (runInMPIMode == 0)
	{
		fprintf (stdout, "\n");
	}
	h = 16*(_mc+1);
	if (psTranslate[0] == 0)
	{
		postScriptOut * ("50 0 translate\n");
	}
	postScriptOut * (_HYPSSetFont ("Courier",12)+"\n"+_HYPSTextTable (500,h,12,summaryMatrix,summaryMatrix["0+12*(_MATRIX_ELEMENT_ROW_==0)"]));
	postScriptOut * ("0 "+h+" translate\n");
	psTranslate = psTranslate + {{50*(psTranslate[0]==0),h}};
}
				 
summaryMatrix = {4,2};
summaryMatrix[0][0] = "Predicted subtype"; summaryMatrix[0][1] = bestAssignment;
summaryMatrix[1][0] = "Model averaged support"; summaryMatrix[1][1] = Format(matchingSum/totalSum*100,0,4) + "\\%";
summaryMatrix[2][0] = "Recombinant"; summaryMatrix[2][1] = Format(recombSupport*100,0,4) + "\\%";
summaryMatrix[3][0] = "Intra-subtype recombinant"; summaryMatrix[3][1] = Format(intraSupport*100,0,4) + "\\%";

if (psTranslate[0] == 0)
{
	postScriptOut * ("50 0 translate\n");
}

postScriptOut *( _HYPSSetFont ("Courier",12)+_HYPSTextTable (500,60,12,summaryMatrix,summaryMatrix["0"]));

if (Abs(correctModel) && runInMPIMode == 0)
{
	fprintf (stdout, "IC difference with the", 
					 "\ncorrect model         : ", correctModelAIC-currentBEST_IC, "\n"); 
}

psTranslate = psTranslate + {{50*(psTranslate[0]==0),0}};

postScriptOut * ("-50 60 translate\n");
postScriptOut * (_HYPSSetFont ("Times-Roman",18));
if (Abs(originalSeqName)>30)
{
	originalSeqName = originalSeqName[0][26] + "...";
}

postScriptOut * ("0 10 translate " + totalPlotWidth/2 + " 0 " + totalPlotWidth + "(SCUEAL subtyping report for " + (originalSeqName&&3) + ") scalecentertext\n");
psTranslate = psTranslate + {{-50,80}};

summaryMatrix = {ibp+1,1};
summaryMatrix[0][0] = "Breakpoint locations";

bpSupport = {ibp,3};

for (_mc=0; _mc<ibp; _mc=_mc+1)
{
	lb = bestPC[_mc];
	ub = bestPC[_mc];
	s  = breakPointSupport[bestPC[_mc]-1];
	while (s<0.95)
	{
		lb = lb-1;
		ub = ub+1;
		
		if (lb>=0)
		{
			s = s + breakPointSupport[lb];
		}
		if (ub<filteredData.sites)
		{
			s = s + breakPointSupport[ub];
		}
	}
	
	summaryMatrix [_mc+1] = Format(bestPC[_mc]+1,0,0) + "bp, 95\\% confidence range: " + (1+Max(lb,0)) +  "-" + (1+Min(ub,filteredData.sites-1)) + " bp.";
	
	if (runInMPIMode == 0)
	{
		fprintf (stdout, "Breakpoint ", Format(_mc+1,3,0), ": ", Format(bestPC[_mc]+1,8,0), 
						 "bp, 95% confidence range: ",	
						  1+Max(lb,0), "-", 1+Min(ub,filteredData.sites-1), " bp.\n");
	}
	else
	{
		bpSupport[_mc][0] = bestPC[_mc]+1;
		bpSupport[_mc][1] = 1+Max(lb,0);
		bpSupport[_mc][2] = 1+Min(ub,filteredData.sites-1);		
	}
}

totalPlotHeight= psTranslate[1];
if (ibpt>20)
{
	h = psTranslate[1]-treeImageHeight-325;
	postScriptOut * ("50 -" + h + " translate\n");
	psTranslate[1] = psTranslate[1] - h;
	postScriptOut * (_HYPSSetFont ("Courier",12)+_HYPSTextTable (500,ibpt,12,summaryMatrix,summaryMatrix["0+12*(_MATRIX_ELEMENT_ROW_==0)"]));
	postScriptOut * ("-50 0 translate\n");
}
h = totalPlotWidth/Abs(_psTreePlots);
postScriptOut * ("0 " + (psTranslate[0]-5) +" sub 0 " + (psTranslate[1]-20) +" sub translate\n");

for (k=0; k<Abs(_psTreePlots); k=k+1)
{
	 if (k)
	 {
		postScriptOut * (""+h +" 0 translate\n");
	 }
	 postScriptOut*_psTreePlots[k];
	 
}
postScriptOut * "\nshowpage";
postScriptOut * 0;
postScriptOut = _HYPSPageHeader (totalPlotWidth+5, totalPlotHeight+10, "SCUEAL report for " + originalSeqName&&3) + "\n" + postScriptOut;

GetDataInfo (refSequenceString, filteredData, querySequenceID);

if (runInMPIMode == 1)
{
	returnAVL["SUBTYPE"] 		= bestAssignment;
	returnAVL["SIMPLE_SUBTYPE"] = bestAssignmentSimple;
	returnAVL["SUPPORT"] 		= matchingSum/totalSum;
	returnAVL["SCORE"]			= currentBEST_IC;
	returnAVL["DELTA_IC"]		= correctModelAIC-currentBEST_IC;
	returnAVL["RECOMB"] 		= recombSupport;
	returnAVL["INTRA"] 			= intraSupport;
	returnAVL["LF"]				= lfExportString;
	returnAVL["PS"]				= postScriptOut;
	returnAVL["TERMINATOR"]		= prematureTermination;
	returnAVL["SEQUENCE"] 		= refSequenceString;
	if (Abs(_extraResult))
	{
		returnAVL["EXTRA"] = _extraResult;
	}
}
else
{
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,s);
	psOutFile = LAST_FILE_PATH;

	fprintf (psOutFile, CLEAR_FILE, postScriptOut);
	psOutFile = psOutFile + ".lf";
	fprintf (psOutFile, CLEAR_FILE, lfExportString);
}	

				 
if (runInMPIMode)
{
	returnAVL["BREAKPOINTS"] = bpSupport;
}


if (Rows(bestPC))
{
	if (runInMPIMode == 0)
	{
		fprintf (stdout, refSequenceString);
	}
	/*else
	{
		returnAVL["SEQUENCE"] = refSequenceString;
	}*/
}


/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


function SCUEALTreeAVL2String (treeAVL)
{
	rootNode = treeAVL[0];
	rootNode = rootNode["Root"];
	return SCUEALsubtreeAVLStr (rootNode,0,0);
}


/*************************************************************************************/

function SCUEALsubtreeAVLStr (nodeIndex,k,treeString)
{
	nodeInfo = treeAVL[nodeIndex];
	k = Abs(nodeInfo["Children"])-1;
	if (k>=0)
	{
		while (k>=0)
		{
			nodeInfo = treeAVL[nodeIndex];
			cNodes = nodeInfo["Children"];
			cNodes = cNodes[k];
			if (k < Abs(nodeInfo["Children"])-1)
			{
				ExecuteCommands("treeString=SCUEALsubtreeAVLStr (cNodes,k,treeString)+\",\"+treeString;");
			}
			else
			{
				ExecuteCommands("treeString=SCUEALsubtreeAVLStr (cNodes,k,treeString)+\")\";");
			}
			k=k-1;
		}
		return "("+treeString+(treeAVL[nodeIndex])["Name"];
	}
	else
	{
		callLevel = callLevel - 1;
		return nodeInfo["Name"];
	}
}

/*************************************************************************************/

function SCUEALInsertANode (theAVL&,insertAt,newNodeName,newParentName)
{
	nodeInfo 	 = theAVL[insertAt];
	newNodeNames = {};
	newNodeNames [nodeInfo["Name"]]    = 1;
	newNodeNames [newNodeName] = 1;
	newNodeNames [newParentName] = 1;
	
	needToReroot = 1;
	
	if (Abs(nodeInfo))
	{
		nparent = nodeInfo["Parent"];
		if (nparent > 0)
		{
			lastIndex = Abs(theAVL);
			myDepth = nodeInfo["Depth"];
			newParentNode = {};
			newParentNode ["Name"] = newParentName;
			newParentNode ["Parent"] = nparent;
			newParentNode ["Depth"] = myDepth;
			
			parentInfo = theAVL[nparent];
			if (parentInfo["Parent"] <= 0)
			{
				needToReroot = 0;
			}
			else
			{
				rerootNodeName = parentInfo["Name"];
			}
			
			newChildNode = {};
			newChildNode ["Name"] = newNodeName;
			newChildNode ["Parent"] = lastIndex;
			newChildNode ["Depth"] = myDepth + 1;
			
			pChildren = {};
			pChildren [0] = insertAt;
			pChildren [1] = lastIndex+1;
			newParentNode ["Children"] = pChildren;
			
			theAVL[lastIndex] = newParentNode;
			theAVL[lastIndex+1] = newChildNode;

			/* update the parent*/

			nodeInfo ["Parent"] = lastIndex;
			theAVL[insertAt] = nodeInfo;
			
			/* update the list of children at the parent node*/
			
			parentChildren = parentInfo["Children"];
			
			for (nic = Abs(parentChildren)-1; nic >= 0; nic = nic-1)
			{
				if (parentChildren[nic] == insertAt)
				{
					break;
				}
			}

			parentChildren[nic] = lastIndex;
			parentInfo["Children"] = parentChildren;
			theAVL[nparent] = parentInfo;
			
			/* now update the depths at new NodeName and all of its children */
			
			nodeCache    = {};
			nodeCache[0] = insertAt;
			cacheIndex   = 0;
			
			while (cacheIndex <= Abs(nodeCache))
			{
				nparent 			= nodeCache[cacheIndex];
				nodeInfo 			= theAVL[nparent];
				nodeInfo["Depth"] 	= nodeInfo["Depth"] + 1;
				theAVL[nparent] 	= nodeInfo;
				nodeChildren 		= nodeInfo["Children"];
				for (nic = Abs(nodeChildren)-1; nic >=0; nic = nic-1)
				{
					nodeCache [Abs(nodeCache)] = nodeChildren[nic];
				}
				cacheIndex = cacheIndex + 1;
			}
			
			nodeCache = 0;
		}
	}
	return newNodeNames;
}


/*************************************************************************************/

function	ModifyDepth    (nIndex, modAmount)
{
	nodeInfo = theAVL[nIndex];
	nodeInfo ["Depth"] = nodeInfo ["Depth"] + modAmount;
	theAVL[nIndex] = nodeInfo;

}

/*************************************************************************************/

function    fuseBreakpoints (inVector)
{
	availableParts = decodeIndividual (inVector);
	localPartCount = (Columns (inVector) * Rows (inVector) - branchBits) $ bitsPerPart + 1;
	if (localPartCount > 2) /* have at least two breakpoints */
	{
		whichPart 	   = Random (1,localPartCount)$1;
		newIndividual  = {localPartCount-1,2};
		for (_k = 0; _k < localPartCount; _k=_k+1)
		{
			if (_k!=whichPart)
			{
				newIndividual[_k-(_k>whichPart)][0] = availableParts[_k][0]-1;
				newIndividual[_k-(_k>whichPart)][1] = availableParts[_k][1];
			}
		}
		newIndividual[whichPart-1][0] = availableParts[whichPart-(Random(0,1)>0.5)][0]-1;
		return encodeIndividual (newIndividual);
	}
	return inVector;
}

/*************************************************************************************/

function addOne (inVector, moveBranch)
{
	availableParts = decodeIndividual (inVector);
	localPartCount = (Columns (inVector) * Rows (inVector) - branchBits) $ bitsPerPart + 1;
	
	outVectorDim   = Rows(inVector)+bitsPerPart;
	outVector 	   = {outVectorDim,1};
	
	/* decide which part to split */
	
	goOn = 1;
	while (goOn)
	{
		partToDuplicate = Random (0,localPartCount-0.000001)$1;
		if (partToDuplicate == 0)
		{
			startPos = 0;
			endPos	 = availableParts[1][1];
		}
		else
		{
			startPos = availableParts[partToDuplicate][1]+1;
			if (partToDuplicate < localPartCount - 1)
			{
				endPos	 =  availableParts[partToDuplicate+1][1]+1;
			}
			else
			{
				endPos	 =  bppMapSize;
			}
		}
		span = endPos-startPos-1;
		if (span >= 2)
		{
			goOn   = 0;
			newPos = Random(startPos,endPos-0.00001)$1;
			if (verboseFlag > 1)
			{
				fprintf (stdout, partToDuplicate, "(",localPartCount,"):", startPos, ":", endPos, ":", span, ":", newPos, "\n");
			}
			pc = 0;
			if (partToDuplicate)
			{
				endPos = (partToDuplicate-1)*bitsPerPart + branchBits;
				pc3    = endPos;
			}
			else
			{
				endPos = branchBits;
				pc3    = -bppSize;
			}
			
			for (pc = 0; pc < endPos; pc = pc+1)
			{
				outVector[pc] = inVector[pc];
			} 
			
			for (pc2 = bppSize-1; pc2 >= 0; pc2 = pc2 - 1)
			{
				outVector [pc+pc2] = newPos%2;
				newPos = newPos$2;
			}
			
			if (moveBranch>0)
			{
				for (pc2 = bppSize; pc2 <bitsPerPart; pc2 = pc2 + 1)
				{
					if (Random(0,1)<0.5)
					{
						outVector [pc+pc2] = 1-inVector[pc3+pc2];
					}
					else
					{
						outVector [pc+pc2] = inVector[pc3+pc2];					
					}
					newPos = newPos$2;
				}
			}
			else
			{
				for (pc2 = bppSize; pc2 <bitsPerPart; pc2 = pc2 + 1)
				{
					outVector [pc+pc2] = inVector[pc3+pc2];
					newPos = newPos$2;
				}
			}
			for (pc = endPos; pc < Rows(inVector) ; pc = pc+1)
			{
				outVector[pc+bitsPerPart] = inVector[pc];
			} 
		}
	}
	return outVector;
}
/*************************************************************************************/

function	EchoSubtypeUpdate (theDatum)
{
	newSub = AssembleSubtypeAssignment(theDatum,1);
	if (newSub != currentSubtypeAssignment)
	{
		currentSubtypeAssignment = newSub;
		fprintf (stdout, AssembleSubtypeAssignment(theDatum,1),"\n");
	}
	if (Abs(correctModel))
	{
		fprintf (stdout, "Correct model Delta IC = ", sortedScores[populationSize-1][0] + correctModelAIC, "\n");
	}
	return 0;
}

/*************************************************************************************/

function	GetBreakpoints (theDatumString)
{
	whereisTheSpace = 0;
	while (theDatumString[whereisTheSpace] != "\n")
	{
		whereisTheSpace = whereisTheSpace+1;
	}
	bpList = splitOnRegExp(theDatumString[0][whereisTheSpace-1],",");
	if (Abs(bpList) > 1)
	{
		bpOut = {Abs(bpList)-1,1};
		for (_s = 1; _s < Abs(bpList); _s=_s+1)
		{
			bpOut [_s-1] = 0+(splitOnRegExp(bpList[_s],"-"))[0];
		}
		return bpOut;
	}
	return {{}};
}

/*----------------------------------------------------------------*/

function computeSubtypeReductions ()
{
	subtypeReductions 			 = {};	
	if (Abs(_subtypeAssignmentByNode))
	{
		
		refTreeAVL					 = referenceTopology^1;
		
	
		for (node = 1; node < Abs(refTreeAVL); node = node + 1)
		{
			children = Abs((refTreeAVL[node])["Children"]);
			if (children)
			{
				parent   = (refTreeAVL[node])["Parent"];
				nodeName = (refTreeAVL[node])["Name"];
				
				nodest   = _subtypeAssignmentByNode[nodeName];
				
				if (nodest == 0) continue;
				
				if (parent > 0)
				{
					isSimple = nodest $ "/";
					if (isSimple[0] < 0) /* simple subtype */
					{
						(refTreeAVL[node])["Subtype"] = nodest;
					}		
					else
					{
						if (Abs((refTreeAVL[parent])["Subtype"]))
						{
							(refTreeAVL[node])["Subtype"] = (refTreeAVL[parent])["Subtype"];
							subtypeReductions[nodeName]   = (refTreeAVL[parent])["Subtype"];
						}
					}
				}
			}
		}
	}
	
	return subtypeReductions;
}


/*----------------------------------------------------------------*/

function splitOnRegExp (string, splitter)
{
	matched = string || splitter;
	splitBits = {};
	if (matched [0] < 0)
	{
		splitBits[0] = string;
	}
	else
	{
		mc = 0;
		if (matched[0] == 0)
		{
			fromPos = matched[1]+1;
			mc = 2;
		}
		else
		{
			fromPos = 0;
			toPos	= 0;
		}
		for (; mc < Rows (matched); mc = mc+2)
		{
			toPos = matched[mc]-1;
			splitBits [Abs(splitBits)] = string[fromPos][toPos];
			fromPos    = matched[mc+1]+1;
		}
		splitBits [Abs(splitBits)] = string[fromPos][Abs(string)-1];
	}
	return splitBits;
}

