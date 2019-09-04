ChoiceList (alignmentType, "Codons or Nucleotides", 1, SKIP_NONE, "Codon", 		"In-frame (universal code) codon alignment",
																  "Nucleotide", "Nucleotide alignment",
																  "Codon-direct", "Directly align codons (good for frameshift detection)");

if (alignmentType < 0)
{
	return 0;
}

SetDialogPrompt ("POL sequence file:");
DataSet ds_in = ReadDataFile (PROMPT_FOR_FILE);

whichSeq = -1;
if (ds_in.species > 1)
{
	ChoiceList (whichSeq, "Sequence to screen", 1, SKIP_NONE, ds_in);
}
else
{
	if (ds_in.species == 1)
	{
		whichSeq = 0;
	}
}

referenceFile = "data/reference.nex";
rfp		      = "../Configs/settings.ibf";
if (!rfp)
{
	ExecuteAFile ("../Configs/settings.ibf");
	referenceFile = "data/"+referenceAlignmentFileName;
}

if (whichSeq >= 0)
{
	inOptions2 = {};
	inOptions2["1"] = PATH_TO_CURRENT_BF + "../" + referenceFile;
	inOptions2["0"] = "Universal";
	GetString (sName, ds_in, whichSeq);
	inOptions2["2"] = sName;
	DataSetFilter ds_fil = CreateFilter (ds_in,1);
	GetDataInfo (sData, ds_fil, whichSeq);
	inOptions2["3"] = sData;
	
	fprintf (stdout, "Screening ", sName, "\n");
	
	ExecuteAFile ("../" + referenceFile + ".labels");
	ExecuteAFile ("../HBF/I_am_the_aligner.bf", inOptions2);

	if (Abs(outputAlignment))
	{
		ExecuteAFile ("../HBF/SingleSequenceScan2.bf");
		/*ExecuteAFile ("../HBF/GASP_Postprocessor.bf");*/
	}
}
