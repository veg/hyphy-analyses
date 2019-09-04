RequireVersion ("2.5.0");

if (MPI_NODE_COUNT <= 1) {
	fprintf (stdout, "[ERROR] This script requires MPI \n");
	return 0;
}


ChoiceList (alignmentType, "Codons or Nucleotides", 1, SKIP_NONE, "Codon", 		"In-frame (universal code) codon alignment",
																  "Nucleotide", "Nucleotide alignment",
																  "Codon-direct", "Directly align codons (good for frameshift detection)");

if (alignmentType < 0)
{
	return 0;
}



MPI_NODE_STATUS = {MPI_NODE_COUNT-1,1}; /* sequence indices being processed */

ExecuteAFile 			("../Configs/settings.ibf");
	
SetDialogPrompt 		("A sequence file to screen:");
DataSet ds_in 			= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter ds_fil 	= CreateFilter (ds_in,1);
GetString				(sequenceNames, ds_fil, -1);


fprintf 						(stdout, "\nRead ", ds_in.species, " sequences\n");
SetDialogPrompt					("Write a tab-separated file to:");

/* check if the output file already exists */
fprintf 						(PROMPT_FOR_FILE,"");
resultsFile				= 		LAST_FILE_PATH;

headerString			= 		"Index\tName\tSubtype\tSimplified Subtype\tSupport\tRecombination Support\tIntra-subtype Support\tReading Frame\tBreakpoints\tSequence";

for (k = 0; k < Abs(_extraOutputColumns); k+=1)
{	
	headerString += "\t"+_extraOutputColumns[k];
}


fscanf							(resultsFile, "Lines", testMe);

alreadyDoneIDs = {};
fprintf 						(resultsFile,CLEAR_FILE,headerString);

if (Rows(testMe))
{
	ExecuteAFile			(HYPHY_BASE_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"ReadDelimitedFiles.bf");
	if (testMe[0] == headerString)
	{
		previouslyProcessedIDs = 0;
		for (lineID = 1; lineID < Columns (testMe); lineID = lineID + 1)
		{
			thisLine = splitStringByTab(testMe[lineID]);
			if (Abs (thisLine) == 9)
			{
				previouslyProcessedIDs = previouslyProcessedIDs + 1;
				fprintf (resultsFile, "\n", testMe[lineID]);
				alreadyDoneIDs [thisLine[1] && 1] = 1;
			}
		}
		
		fprintf (stdout, "Read ", previouslyProcessedIDs, " previously processed seqeuences\n");
	}
}

fprintf (stdout, "\n");

ChoiceList (resultType, "Result output option", 1, SKIP_NONE,  "Summary", 		"Only write out a summary for the run",
															   "Summary+Detail", "Generate a PostScript report and a likelihood function with the best fitting model for each sequence (2 files per query sequence).");

if (resultType < 0)
{
	return 0;
}

detailedResults = 0;

if (resultType)
{
	fprintf 					(stdout, "Write detailed results to this path:");
	fscanf 						(stdin,"String", detailedResults);	
}


jobsFinished = 0;

for (seqID = 0; seqID < ds_in.species; seqID = seqID + 1)
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



/*------------------------------------------------------------------------*/

function SendAJob (sequenceID)
{
	if (alreadyDoneIDs[sequenceNames[sequenceID]&&1])
	{
		fprintf (stdout, "[SKIP] Sequence ", sequenceNames[sequenceID], " has been subtyped previously\n");
		jobsFinished = jobsFinished + 1;
		return 0;
	}

	inOptions 	   = {};
	inOptions["0"] = sequenceNames[sequenceID];
	GetDataInfo (theSeq, ds_fil, sequenceID);
	inOptions["1"] = theSeq;
	inOptions["2"] = ""+alignmentType;
	
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
	fprintf (stdout, "[SEND] Sequence ", inOptions["0"], " to MPI node ", mpiNode + 1, "\n");
	MPI_NODE_STATUS [mpiNode] = sequenceID+1;
	MPISend (mpiNode+1, "../HBF/MPI_Wrapper.bf", inOptions);
	return 0;
}

/*------------------------------------------------------------------------*/

function ReceiveAJob (dummy)
{
	MPIReceive 		(-1, whichNode, returnValue);
	whichNode	  = whichNode-1;
	processedID   =	MPI_NODE_STATUS [whichNode]-1; 
	processedName = sequenceNames[processedID];
	
	MPI_NODE_STATUS [whichNode] = 0;
	ExecuteCommands	("returnAVL = " + returnValue);
	jobsFinished    = jobsFinished + 1;
	fprintf (stdout, "[RECEIVE] Sequence ", processedName, " from node ", whichNode + 1, " (", (ds_in.species-jobsFinished), " alignments remaining)");
	subtypeFound = returnAVL["SUBTYPE"];
	simpleSubtype = returnAVL["SIMPLE_SUBTYPE"];
	if (Abs(subtypeFound) == 0) /* error */
	{
		fprintf (stdout, ": Error/ alignment failed\n");
		fprintf (resultsFile, "\n", processedID+1, "\t", processedName, "\tError: alignment failed");
	}
	else
	{
		fprintf (stdout, ": ", subtypeFound, " (simple subtype ",simpleSubtype,")\n");
		supp 			= returnAVL["SUPPORT"];
		seq 		    = returnAVL["SEQUENCE"];
		bps	 			= returnAVL["BREAKPOINTS"];
		recs 			= returnAVL["RECOMB"];
		recsi 			= returnAVL["INTRA"];
		readFrame       = returnAVL["FRAME"];
		
		fprintf (resultsFile, "\n", processedID+1, "\t", processedName, "\t", subtypeFound, "\t", simpleSubtype, "\t", supp, "\t", recs, "\t", recsi,"\t", readFrame);
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
		else
		{
			fprintf (resultsFile, "\t");		
		}
		
		extra			= returnAVL["EXTRA"];
		
		if (Abs(extra))
		{
			for (k = 0; k < Abs(_extraOutputColumns); k+=1)
			{	
				fprintf (resultsFile, "\t", extra[_extraOutputColumns[k]]);
			}
		}
		
		if (resultType)
		{
			outPath = detailedResults + (processedID+1) + ".ps";
			ps 		= returnAVL["PS"];
			fprintf (outPath,CLEAR_FILE,ps);
			outPath = detailedResults + (processedID+1) + ".lf";
			ps 		= returnAVL["LF"];
			fprintf (outPath,CLEAR_FILE,ps);
		}
	}
	return 			whichNode;
}
