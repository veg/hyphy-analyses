RequireVersion ("2.5.0");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Read an alignment (and, optionally, a tree) remove duplicate sequences, and prune the tree accordingly.
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA and, optionally, a tree"
                          };


io.DisplayAnalysisBanner (filter.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);

KeywordArgument                     ("msa", "The MSA to remove duplicate sequences from");
SetDialogPrompt                     ("The MSA to remove duplicate sequences from");

DataSet filter.dataset              = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filter.input          = CreateFilter (filter.dataset,1);

console.log ("> Loaded an alignment with `filter.input.species` sequences and `filter.input.sites` sites from `LAST_FILE_PATH`");
filter.unique_count = alignments.CompressDuplicateSequences ("filter.input", "filter.datafilter.unique", FALSE);
console.log ("\nThere are **`filter.unique_count`** unique sequences in alignment ");

if (filter.unique_count == filter.input.species) {
    console.log ("### No duplicate sequences found");
}

KeywordArgument     ("tree", "An optional tree file to trim", "None");
filter.tree = io.PromptUserForString ("An optional tree file to trim");

if (filter.tree != "None") {
    fscanf (filter.tree, "Raw", filter.tree_string);
    filter.tree  = trees.LoadAnnotatedTopology (filter.tree_string);

    Topology T = filter.tree[terms.trees.newick_with_lengths];

    filter.valid_names = {};
    for (n; in; alignments.GetSequenceNames ("filter.datafilter.unique")) {
        filter.valid_names [n] = TRUE;
    }


    filter.delete_leaves = {};
    for (k, s; in; filter.tree[terms.trees.partitioned]) {
        if (s == terms.tree_attributes.leaf) {
            if ( filter.valid_names / k == FALSE) {
                filter.delete_leaves [k] = TRUE;
            }
        }
    }

    T - utility.Keys(filter.delete_leaves);
    utility.SetEnvVariable ("DATAFILE_TREE", Format (T,1,1));
    utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA", TRUE);
}


KeywordArgument ("output", "Write de-duplicated MSA to");
filter.path = io.PromptUserForFilePath ("Write de-duplicated MSA to");
fprintf (filter.path, CLEAR_FILE, filter.datafilter.unique);
