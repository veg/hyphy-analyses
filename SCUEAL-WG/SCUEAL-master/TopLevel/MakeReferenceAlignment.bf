ExecuteAFile ("../HBF/utils.ibf");
ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
								   + "Utility" + DIRECTORY_SEPARATOR + 
								   "GrabBag.bf");
								   
SetDialogPrompt ("Reference sequences?");

ChoiceList (alignmentType, "Codons or Nucleotides", 1, SKIP_NONE, "Codon", 		"In-frame (universal code) codon alignment",
																  "Nucleotide", "Nucleotide alignment");

if (alignmentType < 0)
{
	return 0;
}

onlyFilterSequenceNames = alignmentType;

DataSet ds 			   = ReadDataFile (PROMPT_FOR_FILE);
inFileToAlign 		   = LAST_FILE_PATH;
inFileInfo			   = splitFilePath (inFileToAlign);

alignmentOptions 	   = {};
if (alignmentType == 0)
{
	alignmentOptions ["0"] = "HIV 25%";
	alignmentOptions ["1"] = "No penalty";
	alignmentOptions ["2"] = "First in file";
	alignmentOptions ["3"] = "No";
	alignmentOptions ["4"] = inFileToAlign;
	alignmentOptions ["41"] = "No";
	alignmentOptions ["5"] = "No";
	alignmentOptions ["6"] = "Universal";
	alignmentOptions ["7"] = inFileInfo["DIRECTORY"] + inFileInfo["FILENAME"] + 	
							".aln";
}
else
{
	alignmentOptions ["0"]  = "No penalty";
	alignmentOptions ["2"]  = "First in file";
	alignmentOptions ["3"]  = "No";
	alignmentOptions ["4"]  = inFileToAlign;
	alignmentOptions ["41"] = "No";
	alignmentOptions ["7"]  = inFileInfo["DIRECTORY"] + inFileInfo["FILENAME"] + 	
							".aln";

}
						
ChoiceList (doAlign, "Align the sequences?", 1, SKIP_NONE, "Yes", "Align reference sequences",
													 	   "No", "Assume that sequences are already aligned");
				
cleanNames = {};
cleanNames ["2"] = "Yes/No";


hadTree = 0;
if (doAlign == 0)
{
	fprintf (stdout, "\n[PHASE 1]. Sequence alignment (wrt the first input sequence)\n");			
					
	if (alignmentType==0)
	{
		cleanNames ["0"] = "Universal";
		ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
									   + "SeqAlignmentCodon.bf", alignmentOptions);
	}
	else
	{
		ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
									   + "SeqAlignmentNuc.bf", alignmentOptions);
	}
	
	fprintf (stdout, "\n[PHASE 2]. Cleaning up sequence names and removing duplicates (if any)\n");			
	cleanNames ["1"] = alignmentOptions ["7"] + ".nuc";
}
else
{
	if (alignmentType==0)
	{
		cleanNames ["0"] = "Universal";
	}
	cleanNames ["1"] = inFileToAlign;
}
cleanNames ["3"] = inFileToAlign;

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
								   	   + "CleanStopCodons.bf", cleanNames);
								   	   
DataSet reloadMe = ReadDataFile (inFileToAlign);

useFileTree = Abs(DATAFILE_TREE);
if (useFileTree)
{
	ChoiceList (useFileTree, "Use file tree?", 1, SKIP_NONE, "Yes", "Use the tree already present in the file",
				  	      "No", "Use a TN93 NJ tree");
	if (useFileTree < 0)
	{
		return -1;
	}
	useFileTree = 1- useFileTree;
}

if (useFileTree == 0)
{
	fprintf (stdout, "\n[PHASE 3]. Making a TN93 NJ tree");			
}

njOptions = {};
njOptions ["0"] = "Distance formulae";
njOptions ["1"] = "Nucleotide/Protein";
njOptions ["2"] = cleanNames ["3"];
njOptions ["3"] = "Force Zero";
njOptions ["4"] = "TN93";
njOptions ["5"] = "y";
njOptions ["6"] = inFileInfo["DIRECTORY"] + inFileInfo["FILENAME"] + "_nj.tree";

if (useFileTree == 0)
{
	ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
									   + "NeighborJoining.bf", njOptions);
}
else
{
	outFile = njOptions ["6"];
	fprintf (outFile, CLEAR_FILE, DATAFILE_TREE);
}
skipSwap = 1;

ChoiceList (skipSwap, "Branch Swapping", 1, SKIP_NONE, "Yes", "Do branch swapping",
												     "No", "Skip branch swapping");
												     
if (skipSwap == 0)
{						
	fprintf (stdout, "\n[PHASE 4]. Doing Branch Swapping");		
	
	swappingOptions = {};
	swappingOptions ["0"] = "Greedy";
	swappingOptions ["1"] = "Nucleotide/Protein";
	swappingOptions ["2"] = njOptions ["2"];
	swappingOptions	["3"] = "GRM";
	swappingOptions	["4"] = "Global w/variation";
	swappingOptions	["5"] = "General Discrete";
	swappingOptions	["6"] = "3";
	swappingOptions	["7"] = njOptions ["6"];
	swappingOptions	["8"] = "Quick and Dirty";
	swappingOptions	["9"] = njOptions ["6"];
	
	VERBOSITY_LEVEL = 1;

	ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
								   + "BranchSwap.bf", swappingOptions);

	VERBOSITY_LEVEL = 0;
}

fprintf (stdout, "\n[PHASE 5]. Fitting a model to estimate the ancestor");	

VERBOSITY_LEVEL = 1;
OPTIMIZATION_PROGRESS_QUANTUM = 0.5;
OPTIMIZATION_PROGRESS_STATUS  = "OPTIMIZING THE LIKELIHOOD FUNCTION";
OPTIMIZATION_PROGRESS_TEMPLATE = "$1 Log(L) = $2 ($3% done) Time elapsed: $4 LF evals/second: $5 CPU load: $6";

if (alignmentType == 0)
{
	codonOptions = {};
	codonOptions ["0"] = "Universal";
	codonOptions ["1"] = cleanNames ["3"];
	codonOptions ["2"] = "MG94CUSTOM";
	codonOptions ["3"] = "Global";
	codonOptions ["4"] = "012345";
	if (IS_TREE_PRESENT_IN_DATA)
	{
		codonOptions ["5"] = "n";
	}
	codonOptions ["6"] = njOptions ["6"];
	codonOptions ["7"] = "Estimate";
	
	ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
									   + "AnalyzeCodonData.bf", codonOptions);
}

if (alignmentType == 1)
{
	codonOptions = {};
	codonOptions ["0"] = cleanNames ["3"];
	codonOptions ["1"] = "GRM";
	codonOptions ["2"] = "Global w/variation";
	codonOptions ["3"] = "General Discrete";
	codonOptions ["4"] = "3";
	if (IS_TREE_PRESENT_IN_DATA)
	{
		codonOptions ["5"] = "n";
	}
	codonOptions ["6"] = njOptions ["6"];
	codonOptions ["7"] = "Don't Display";
	
	ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR 
									   + "AnalyzeNucProtData.bf", codonOptions);
}

VERBOSITY_LEVEL = 0;

fprintf (stdout, "\n[PHASE 6]. Reconstructing ancestors.");		

DataSet			ancestralSequences = ReconstructAncestors (lf);
DataSet			jointDS			   = Combine (ds,ancestralSequences);
DataSetFilter	referenceFilter	   = CreateFilter (jointDS,1,"",speciesIndex <= filteredData.species);

IS_TREE_PRESENT_IN_DATA 		   = 1;
DATAFILE_TREE					   = Format(givenTree,1,1);
DATA_FILE_PRINT_FORMAT			   = 6;

SetDialogPrompt ("Save reference alignment and tree to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,referenceFilter);
outFilter = LAST_FILE_PATH + ".labels";

fprintf (stdout, "\n[PHASE 7]. Label sequences.\n");
sequenceLabels = {};

tc = TipCount (givenTree);

/*
for (k=0; k < tc; k = k+1)
{
	nodeName 	  = TipName (givenTree,k);
	fprintf 		(stdout, "Label for tip ", nodeName, " : ");
	fscanf			(stdin, "String", nodeLabel);
	sequenceLabels [nodeName&&1] = nodeLabel;
}*/

LoadFunctionLibrary ("ReadDelimitedFiles");

for (k=0; k < tc; k = k+1)
{
	nodeName 	  = TipName (givenTree,k);
	subexp = extractAllExpressions (nodeName, "[^_]+", "");
	sequenceLabels [nodeName&&1] = subexp[0];
	
	fprintf (stdout, nodeName, "->", sequenceLabels [nodeName&&1], "\n");
}


fprintf (stdout, "Auto-generating internal node labels"); 
tc = BranchCount (givenTree);
for (k=0; k < tc; k = k+1)
{
	nodeName 	  = BranchName (givenTree,k);
	subtreeAVL	  = givenTree[nodeName];
	nodeNames	  = Rows(subtreeAVL);
	crfType			= {};
	pureType		= {};
	
	for (k2 = 0; k2 < Abs (subtreeAVL); k2 = k2+1)
	{
		if (((nodeNames[k2]&&1)$"^NODE")[0] < 0)
		{
			thisLabel = sequenceLabels[nodeNames[k2]&&1];
			if (thisLabel/"CRF*")
			{
				crfType [thisLabel[3][Abs(thisLabel)-1]] = 1;
			}
			else
			{
				pureType [thisLabel] = 1;
			}
		}
	}
	
	_bufferString = ""; _bufferString * 128;
	/* reduce CRF labels here */
	_bufferCount = 0;
	if (Abs(pureType) > 1)
	{
		pureType["buildINodeLabel"][""];
		if (Abs (crfType))
		{
			_bufferCount = 0;
			_bufferString * "/CRF";
			crfType["buildINodeLabel"][""];
		}
	}
	else
	{
		if (Abs(pureType) == 1)
		{
			nodeLabel = (Rows(pureType))[0];
		}
		else
		{
			_bufferString * "CRF";
			crfType["buildINodeLabel"][""];
		}
	}
	_bufferString * 0;
	if (Abs(_bufferString))
	{
		sequenceLabels [nodeName&&1] = _bufferString;
	}
	else
	{
		sequenceLabels [nodeName&&1] = nodeLabel;
	}
}

fprintf (outFilter, CLEAR_FILE, "_subtypeAssignmentByNode = ", sequenceLabels, ";");

