SetDialogPrompt (".labels file");
ExecuteAFile (PROMPT_FOR_FILE);

haveSubtype    = {};
alreadyDone = {};

_subtypeAssignmentByNode["grabCRFs"][""];

fprintf (stdout, haveSubtype);

function grabCRFs (key,value)
{
	if ((key $ "^NODE")[0] < 0)
	{
		lKey = key ^ {{"\\_CRF\\_[0-9]+",""}};
		if (alreadyDone[lKey] == 0)
		{
			alreadyDone[lKey] = 1;
			haveSubtype[value] = haveSubtype[value]+1;
			
		}
	}
	return 0;
}
