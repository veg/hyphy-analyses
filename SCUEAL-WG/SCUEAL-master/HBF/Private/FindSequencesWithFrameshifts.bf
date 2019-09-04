alignmentType = 2;

SetDialogPrompt ("POL sequence file:");
DataSet ds_in = ReadDataFile (PROMPT_FOR_FILE);


referenceFile = "../data/reference.nex";
rfp		      = "../../Configs/settings.ibf";
if (!rfp)
{
	ExecuteAFile ("../../Configs/settings.ibf");
	referenceFile = "../data/"+referenceAlignmentFileName;
}

for (whichSeq = 0; whichSeq < ds_in.species; whichSeq += 1)
{
	inOptions2 = {};
	inOptions2["1"] = referenceFile;
	inOptions2["0"] = "Universal";
	GetString (sName, ds_in, whichSeq);
	inOptions2["2"] = sName;
	DataSetFilter ds_fil = CreateFilter (ds_in,1);
	GetDataInfo (sData, ds_fil, whichSeq);
	inOptions2["3"] = sData;
	
	fprintf (stdout, "Screening ", sName, "\n");
	
	ExecuteAFile ("../" + referenceFile + ".labels");
	ExecuteAFile ("../I_am_the_aligner.bf", inOptions2);

}
