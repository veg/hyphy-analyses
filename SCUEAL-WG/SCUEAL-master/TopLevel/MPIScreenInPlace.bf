if (MPI_NODE_COUNT <= 1)
{
	fprintf (stdout, "[ERROR] This script requires MPI \n");
	return 0;
}

SetDialogPrompt ("Which configuration file:?");
ExecuteAFile 	(PROMPT_FOR_FILE);

MPI_NODE_STATUS = {MPI_NODE_COUNT-1,1}; /* sequence indices being processed */




fprintf 						(stdout, "\nProcessing ", 
											inputFileCount, 
											" files from prefix",
											fileNamePrefix, "\n");
											
SetDialogPrompt					("Write a tab-separated file result file to:");
fprintf 						(PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"Index\tFilepath\tType\tSupport\tIC Difference\tBreakpoints\tSequence");
resultsFile				= 		LAST_FILE_PATH;

fprintf (stdout, "\n");
jobsFinished = 0;

for (seqID = startScreeningFrom; seqID < inputFileCount; seqID = seqID + 1)
{
	SendAJob (seqID);
}

/* clean up MPI jobs */

howManyPending = 0;
for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
{
	if (MPI_NODE_STATUS[mpiNode])
	{
		howManyPending = howManyPending + 1;
	}
}

for (; howManyPending; howManyPending = howManyPending-1)
{
	ReceiveAJob (0);
}

fprintf 						(resultsFile, CLOSE_FILE);



/*------------------------------------------------------------------------*/

function SendAJob (sequenceID)
{
	inOptions 	   = {};
	inOptions["0"] = fileNamePrefix + sequenceID + fileNameSuffix;
	inOptions["1"] = recombinantSequenceID;
	inOptions["2"] = referenceTreeStringIn;
	inOptions["3"] = correctModelPath;
	
	for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
	{
		if (MPI_NODE_STATUS[mpiNode] == 0) /* free node */
		{
			break;
		}
	}
	if (mpiNode == MPI_NODE_COUNT-1) /* all busy */
	{
		mpiNode = ReceiveAJob (0);
	}
	fprintf (stdout, "[SEND] File ", sequenceID+1, " to MPI node ", mpiNode + 1, "\n");
	MPI_NODE_STATUS [mpiNode] = sequenceID+1;
	MPISend (mpiNode+1, "InPlaceScreenSimMPI.bf", inOptions);
	return 0;
}

/*------------------------------------------------------------------------*/

function ReceiveAJob (dummy)
{
	MPIReceive 		(-1, whichNode, returnValue);
	whichNode	  = whichNode-1;
	processedID   =	MPI_NODE_STATUS [whichNode]-1; 
	processedName = fileNamePrefix + processedID + fileNameSuffix;
	
	MPI_NODE_STATUS [whichNode] = 0;
	ExecuteCommands	("returnAVL="+returnValue);
	jobsFinished    = jobsFinished + 1;
	fprintf (stdout, "[RECEIVE] File ", processedID+1, " from node ", whichNode + 1, " (", (inputFileCount-startScreeningFrom-jobsFinished), " alignments remaining)");
	
	subtypeFound = returnAVL["SUBTYPE"];
	if (Abs(subtypeFound) == 0) /* error */
	{
		fprintf (stdout, ": Error/ alignment failed\n");
		fprintf (resultsFile, "\n", processedID+1, "\t", processedName, "\tError: alignment failed");
	}
	else
	{
		fprintf (stdout, ": ", subtypeFound, "\n");
		supp 	= returnAVL["SUPPORT"];
		bps	 	= returnAVL["BREAKPOINTS"];
		icfound = returnAVL["DELTA_IC"];
		seq		= returnAVL["SEQUENCE"];
		
		fprintf (resultsFile, "\n", processedID+1, "\t", processedName, "\t", subtypeFound, "\t", supp, "\t", icfound);
		if (Abs(seq))
		{
			fprintf (resultsFile, "\t");
			for (bpID = 0; bpID < Rows (bps); bpID = bpID + 1)
			{
				if (bpID > 0)
				{
					fprintf (resultsFile, ";");
				}
				fprintf (resultsFile, bps[bpID][0], "(", bps[bpID][1], "-", bps[bpID][2], ")");
			}
			fprintf (resultsFile, "\t", seq);
		}
	}
	return 			whichNode;
}
