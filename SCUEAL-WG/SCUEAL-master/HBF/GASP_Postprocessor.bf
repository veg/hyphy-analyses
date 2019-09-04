ExecuteAFile (HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "ReadDelimitedFiles.bf");
RequireVersion ("0.9920060821");

treeString = Format(referenceTopology,1,1);

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


		
SetDialogPrompt ("Save SCUEAL detailed resutls for later use here:");		
		
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, _subtypeAssignmentByNode, 
						"screen_results={};\n", 
					    "_hyphyAssociativeArray[\"Branch Lengths\"] = \"",
					    bestModelBL,
					    "\";",
					    "_hyphyAssociativeArray[\"IC value\"] = \"",
					    bestModelIC,
					    "\";",
					    "_hyphyAssociativeArray[\"Tree\"] = \"",
					    Format(baselineTree,1,1),
					    "\";",
					    "_hyphyAssociativeArray[\"Models\"] = ",
					    models,
					    ";",
					   "_hyphyAssociativeArray[\"Scores\"] = ",
					   scores,
					   ";",
					   "\nscreen_results[\"",
					   sName,
					   "\"] = _hyphyAssociativeArray;");

lastFilePath		 = LAST_FILE_PATH + ".fit";



modelAIC = {};
modelBP	 = {};
minAIC	 = 1e100;

Topology	refTop  = treeString;

bNames				= BranchName (refTop,-1);
bCount  			= Columns (bNames)-1;

siteCount			= ds.sites;
modelCount			= Abs(MasterList);
bpSupport			= {siteCount,2};

supportByBranch 	= {};
modelSupport		= {modelCount,1};
branchSupport		= {siteCount, bCount};
nameToIndex			= {};
for (k=0; k< bCount; k=k+1)
{
	nameToIndex[bNames[k]] = k;
}

resMatrix0 = Rows(MasterList);

bestAIC = 1e100;
for (k=0; k<modelCount; k=k+1)
{
	modelSupport[k] = 0 - (MasterList[resMatrix0[k]]);
	bestAIC				= Min (bestAIC, modelSupport[k]);
}

modelSupport = modelSupport["Exp(bestAIC-_MATRIX_ELEMENT_VALUE_)/2"];
bestAIC		 = (Transpose(modelSupport["1"])*modelSupport)[0];
modelSupport = modelSupport * (1/bestAIC);

stencil		 = {siteCount,1};

for (k=0; k<modelCount; k=k+1)
{
	if (modelSupport[k] > 0.00001)
	{
		currentKey = resMatrix0[k];
		whereisTheSpace = 0;
		while (currentKey[whereisTheSpace] != "\n")
		{
			whereisTheSpace = whereisTheSpace +1;
		}
		splits  = splitOnRegExp (currentKey[0][whereisTheSpace-1],   ",");
		splits2 = splitOnRegExp (currentKey[whereisTheSpace+1][Abs(currentKey)-1], ",");
		bitCount = Abs(splits);
		for (k2 = 0; k2 < bitCount; k2=k2+1)
		{
			bp    = splitOnRegExp(splits[k2],"-");
			bps   = 0+bp[0];
			bpe   = 0+bp[1];
			bn	  = splits2[k2];
			
			if (k2)
			{
				if (bn!=bnp)
				{
					bpSupport [bps-1][0] = bpSupport [bps-1][0] + modelSupport[k];
				}
				else
				{
					bpSupport [bps-1][1] = bpSupport [bps-1][1] + modelSupport[k];
				}
			}
			
			toAdd = stencil["(_MATRIX_ELEMENT_ROW_>=bps__)*(_MATRIX_ELEMENT_ROW_<=bpe__)"];
			bi    = nameToIndex[bn];
			for (idx = 0; idx < siteCount; idx = idx + 1)
			{
				branchSupport[idx][bi] = branchSupport[idx][bi] + toAdd[idx] * modelSupport[k];
			}
			bnp   = bn;
		}
	}
}

thresh   = 0.01;
displayB = {};

for	(k=0; k<bCount; k=k+1)
{
	if(((branchSupport[-1][k])%0)[siteCount-1] >= thresh)
	{
		displayB[Abs(displayB)] = k;
	}
}

selectedBranches = Abs(displayB);
labels = {1,selectedBranches};

labelsToPlot = "";

for	(k=0; k<selectedBranches; k=k+1)
{
	labels [k] = _subtypeAssignmentByNode[bNames[displayB[k]]];
	if (k)
	{
		labelsToPlot = labelsToPlot + ";";
	}
	labelsToPlot = labelsToPlot + labels [k];
}

dataMatrix = {siteCount, selectedBranches};

for (k2 = 0; k2 < selectedBranches; k2=k2+1)
{
	dbi = displayB[k2];
	for (k=0; k<siteCount; k=k+1)
	{
		dataMatrix[k][k2] = branchSupport[k][dbi];
	}
}


OpenWindow (CHARTWINDOW,{{"GASP Clustering Results"}
		{"labels"}
		{"dataMatrix"}
		{"Step Plot"}
		{"Index"}
		{labelsToPlot}
		{"Nucleotide"}
		{""}
		{"Support"}
		{"3"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Hoefler Text:14:0;Hoefler Text:14:0;Hoefler Text:14:1"}
		{"0;0;16777215;5000268;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"500;500;70;70");

columnHeaders = {{"Hard","Soft"}};
bpCounter     = {1,Rows (bpSupport)} ["1"] * bpSupport;

hbpc = (bpCounter[0]+0.5)$1;
sbpc = (bpCounter[1]+0.5)$1;

fprintf (stdout, "\n", hbpc, " hard breakpoint(s).\n", sbpc, " soft breakpoint(s).\n");

OpenWindow (CHARTWINDOW,{{"Breakpoint Placement"}
		{"columnHeaders"}
		{"bpSupport"}
		{"Bar Chart"}
		{"Index"}
		{"Hard;Soft"}
		{""}
		{""}
		{""}
		{"3"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Hoefler Text:14:0;Hoefler Text:14:0;Hoefler Text:14:1"}
		{"0;0;16777215;5000268;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"500;500;300;300");

