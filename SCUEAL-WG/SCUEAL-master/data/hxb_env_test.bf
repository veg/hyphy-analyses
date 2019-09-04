ExecuteAFile				 (HYPHY_BASE_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"HXB2Mapper.bf");
ExecuteAFile				 (HYPHY_BASE_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"ReadDelimitedFiles.bf");

DataSet 				ds  	= ReadDataFile ("/Users/sergei/Desktop/SCUEAL/data/HIV1_ALL_2007_ENV_DNA.fasta.txt");
DataSetFilter			dsf		= CreateFilter (ds,1);
c2gp41							= {{"c2","v3","c3","v4","c4","v5","post","gp41ecto","gp41endo"}};


seqIDs							= {};
fastaOut						= ""; fastaOut * 256;
for (s = 0; s < dsf.species; s=s+1)
{
	GetDataInfo			  			  (seqInfo, dsf, s);
	seqInfo							= seqInfo ^ {{"[^A-Z]",""}};
	GetString 						  (seqID, dsf, s);
	fprintf (stdout, "Working on ", seqID, "\n");
	fastaOut						* (">" + normalizeSequenceID (seqID, "seqIDs") + "\n" + selectHXB2subsequence 			(seqInfo, c2gp41)+"\n");
}

fastaOut * 0;
SetDialogPrompt ("Save FASTA to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, fastaOut);
