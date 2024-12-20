RequireVersion ("2.5.60");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Read an alignment and a tree and rename all the sequences to conform with HyPhy naming requirements.
                            The result will be written as a combined alignment, with the format specified by the DATA_FILE_PRINT_FORMAT HyPhy environment variable 
                            ",
                            terms.io.version :          "0.0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA and a tree."
                          };


io.DisplayAnalysisBanner (filter.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

KeywordArgument                     ("msa", "The MSA to process");
SetDialogPrompt                     ("The MSA to process");

info = alignments.ReadNucleotideDataSet ("msa", TRUE);

DataSetFilter msa.filter = CreateFilter (msa, 1);


if (None != info[terms.data.name_mapping]) {
    name_mapping = {};
    
    for (k, v; in; info[terms.data.name_mapping]) {
        name_mapping[v] = k;
    }    
} else {
    name_mapping = {};
    
    for (v; in; alignments.GetSequenceNames ("msa")) {
        name_mapping[v] = v;
    }
}



KeywordArgument                     ("tree", "The tree");
SetDialogPrompt                     ("The tree");

fscanf (PROMPT_FOR_FILE, "Raw", match.tree_string);
match.tree  = trees.LoadAnnotatedTopology (match.tree_string);
Topology T = match.tree[terms.trees.newick_with_lengths];


//ExecuteCommands ("Topology T = " + T);

io.CheckAssertion("info[terms.data.sequences]==TipCount(T)", "The number of tips in the tree does not match the number of sequences in the alignment");
match.existing = match.tree [terms.trees.model_map];

KeywordArgument  ("regexp", "Use the following pattern to match sequences and tree labels","^(.+)$");
match.regexp    = io.PromptUserForString ("Use the following regular expression to select a subset of leaves");
match.sequences = {};

for (original,renamed; in;name_mapping) {

    match.me = regexp.FindSubexpressions (original,match.regexp);
    io.CheckAssertion("None!=match.me", "Could not match pattern for sequence name `original`");
    io.CheckAssertion("Abs (match.me) == 2", "Could not extract subexpression match pattern for sequence name `original`");
    io.CheckAssertion("match.sequences/match.me[1] == 0", "Duplicate match pattern for sequence `original`");
    match.sequences[match.me[1]] = renamed;
}

theAVL = T^0;
_ost = "";
_ost * 256;


lastLevel = 0;
treeSize  = Abs(theAVL);
treeInfo  = theAVL[0];
rootIndex = treeInfo["Root"];
lastDepth = 0;


for (nodeIndex = 1; nodeIndex < treeSize; nodeIndex += 1) {
    nodeInfo = theAVL[nodeIndex];
    myDepth = nodeInfo["Depth"];
    if (lastDepth < myDepth) {
        if (lastDepth) {
            _ost * ",";
        }
        for (pidx = lastDepth; pidx < myDepth; pidx += 1) {
            _ost * "(";
        }
    } else {
        if (lastDepth > myDepth) {
            for (pidx = myDepth; pidx < lastDepth; pidx += 1) {
                _ost * ")";
            }
        } else {
            _ost * ",";
        }
    }
    
    
    if (Abs(nodeInfo["Children"]) == 0) {
        match.me = regexp.FindSubexpressions (nodeInfo["Name"],match.regexp);
        io.CheckAssertion("None!=match.me", "Could not match pattern for tree tip `nodeInfo['Name']`");
        io.CheckAssertion("Abs (match.me) == 2", "Could not extract subexpression match pattern for tree tip `nodeInfo['Name']`");
        match.me = match.me[1];
        io.CheckAssertion('match.sequences / match.me', "Tip `nodeInfo["Name"]` has no match in the alignment");
        _ost * match.sequences[match.me];
        match.sequences - match.me;
    }
    
    if (nodeIndex < treeSize - 1) {
        if (Abs (match.existing[nodeInfo["Name"]]) > 0) {
            _ost * '{';
            _ost * match.existing[nodeInfo["Name"]];
            _ost * '}';
            
        }
        _ost * ":";
        _ost * (""+nodeInfo ["Length"]);
        
    }
    lastDepth = myDepth;
}

_ost * 0;

Tree T = _ost;

utility.SetEnvVariable ("DATAFILE_TREE", _ost);
utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA", TRUE);

KeywordArgument ("output", "Write cleaned MSA+Tree to");
path = io.PromptUserForFilePath ("Write cleaned MSA+Tree to");
fprintf (path, CLEAR_FILE, msa.filter);
