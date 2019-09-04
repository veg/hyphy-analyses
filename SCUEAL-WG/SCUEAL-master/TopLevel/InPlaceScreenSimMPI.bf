SetDialogPrompt ("A sequence file:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

whichSeq = -1;
if (ds.species > 3)
{
	ChoiceList (whichSeq, "Sequence to screen", 1, SKIP_NONE, ds);
}
else
{
	fprintf (stdout, "[ERROR:] At least 3 sequences are required\n");
	return 0;
}

fprintf (stdout, "Reference tree:");
fscanf 	(stdin,  "String", DATAFILE_TREE);
fprintf (stdout, "Correct model file?:");
fscanf 	(stdin,  "String", correctModelFile);

if (whichSeq >= 0)
{
	EXISTING_ALIGNMENT_RUN = 1;
	runInMPIMode 		   = 1;
	ExecuteAFile ("../HBF/SingleSequenceScan2.bf");
	ExecuteAFile ("../HBF/GASP_Postprocessor.bf");
}

return returnAVL;