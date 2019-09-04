fscanf (stdin, "String", sName);
fscanf (stdin, "String", sData);
fscanf (stdin, "Number", alignmentType);

referenceFile = "data/reference.nex";
rfp		      = "../Configs/settings.ibf";
if (!rfp)
{
	ExecuteAFile ("../Configs/settings.ibf");
	referenceFile = "data/"+referenceAlignmentFileName;
	HAVE_PRESET_VALUES = 1;
}

inOptions2 = {};
inOptions2["0"] = "Universal";
inOptions2["1"] = "../"+referenceFile;
inOptions2["2"] = sName;
inOptions2["3"] = sData;

ExecuteAFile ("../" + referenceFile + ".labels");
ExecuteAFile ("I_am_the_aligner.bf", inOptions2);

if (frameSelectByAlignment == 1)
{
	bestFrameReturned = bestFrame;
}
else
{
	bestFrameReturned = 0;
}

if (stepByStepLogging)
{
	fprintf (MESSAGE_LOG,"\nProcessing sequence ", sName, "\n");
}

if (Abs(outputAlignment))
{
	if (Abs(saveAlignmentToFile))
	{
		alignmentLogFile = saveAlignmentToFile + sName; 
		fprintf (alignmentLogFile, CLEAR_FILE, outputAlignment);
	}

	runInMPIMode = 1;
	ExecuteAFile ("SingleSequenceScan2.bf");
	if (stepByStepLogging)
	{
		fprintf (MESSAGE_LOG,"\nSCUEAL results: ", returnAVL, "\n");
	}
	if (Type (returnAVL) == "AssociativeList")
	{
		returnAVL ["FRAME"] = bestFrameReturned;
	}
	else
	{
		return {};
	}
	return returnAVL;
}
else
{
	return {};
}

