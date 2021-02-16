
LoadFunctionLibrary ("libv3/all-terms.bf");

LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");


function IgSCUEAL.set_up (settings) {
    IgSCUEAL.settings = settings;
    io.CheckAssertion ("Type (IgSCUEAL.settings) == 'AssociativeList'",
                        "**Missing {Dict} IgSCUEAL.settings**");

    IgSCUEAL.universal_code = alignments.LoadGeneticCode ("Universal");

    IgSCUEAL.references     = utility.Map (IgSCUEAL.settings["segments"], "_ref_",
                                       "IgSCUEAL.load_reference_alignment (_ref_[1], _ref_[2], _ref_[0])");

    IgSCUEAL.alignment_options = IgSCUEAL.define_alignment_settings (IgSCUEAL.universal_code);
    IgSCUEAL.alignment_options["code"] = IgSCUEAL.universal_code;

    if (IgSCUEAL.settings / "d") {
        IgSCUEAL.settings ["d"] = IgSCUEAL.load_d_region (IgSCUEAL.settings["d"]);
    } else {
        IgSCUEAL.settings ["d"] = None;
    }

    IgSCUEAL.precomputed_likelihoods = {};

}

function IgSCUEAL.screen_a_read (read) {
    sanitized_read = IgSCUEAL.sanitize_sequence(read);
    mapped_read     = IgSCUEAL.handle_sequence_alignment (sanitized_read, IgSCUEAL.references, IgSCUEAL.alignment_options);
    annotated_read = None;
    if (None != mapped_read) {
        annotated_read  = IgSCUEAL.annotate (mapped_read, IgSCUEAL.references, IgSCUEAL.universal_code, IgSCUEAL.settings ["d"]);
        extracted_data = IgSCUEAL.extracted_for_phylomap (mapped_read, IgSCUEAL.references);
        utility.ForEach (mapped_read["SEGMENTS"], "_segment_", "IgSCUEAL.phylo.prepare_alignment (_segment_['index'], ((IgSCUEAL.settings ['segments'])[_segment_['index']])[1], IgSCUEAL.precomputed_likelihoods)");
        utility.Extend (annotated_read, IgSCUEAL.phylo.place_sequence (extracted_data, IgSCUEAL.precomputed_likelihoods));
        annotated_read = utility.Extend(annotated_read, {"input-sequence" : sanitized_read});
    }

    return annotated_read;
}

// ---------------------------------------------------------------------------------------------------------
// Phylogenetic placement functions
// ---------------------------------------------------------------------------------------------------------

namespace IgSCUEAL.phylo {

    lfunction place_sequence (phylomap_data, reference_data) {

        placement_result = {};

        numerical_filter = {"FILTER_NAMES" : {{"Clade1","Clade2","Query"}},
                            "FILTER_ARRAYS" : {}};

        mapped_segments     = utility.Keys (phylomap_data);
        segment_count       = utility.Array1D (mapped_segments);

        for (segment_id = 0; segment_id < segment_count; segment_id += 1) {

            segment_name = mapped_segments[segment_id];

            branch_names = utility.Keys ((reference_data[segment_name])["conditionals"]);
            branch_count = utility.Array1D (branch_names) - 1; // subtract "Node0"

            ExecuteCommands ((reference_data[segment_name])["model"]);
            Tree screening_tree = (Clade1, Clade2, Query);

            assert (Abs ((phylomap_data[segment_name])["mapped_qry"]) > 0 && Type ((phylomap_data[segment_name])["mapped_qry"]) == "String" ,
                    "Internal error in place_sequence ");

            DataSet qry           = ReadFromString (">qry\n" + (phylomap_data[segment_name])["mapped_qry"]);


            DataSetFilter   qry_f = CreateFilter (qry, 1);


            unique_sites = utility.Array1D (qry_f.site_freqs);
            resolutions  = {};
            GetDataInfo (map_to_unique, qry_f);
            for (p = 0; p < unique_sites; p+=1) {
                GetDataInfo (character_resolution, qry_f, 0, p);
                resolutions[p] = character_resolution;
            }

            //console.log ((phylomap_data[segment_name])["mapped_qry"]);
            numeric_query = {qry.sites, 4};
            for (p = 0; p < qry.sites; p+=1) {
                pattern = resolutions[map_to_unique[p]];
                //console.log ("Site " + (p+1) + " => " + pattern + "\n");
                for (s = 0; s < 4; s+=1) {
                    numeric_query[p][s] = pattern[s];
                }
            }

            this_branch = branch_names [0];
            reference_sites = Rows ((((reference_data[segment_name])["conditionals"])[this_branch])["up"]);

            (numerical_filter["FILTER_ARRAYS"])[2] = numeric_query;

            filter   = (phylomap_data[segment_name])["filter_def"];
            selector = Transpose(utility.Map (utility.Keys (filter), "_site_", "0 + _site_"));
            numerical_filter["FILTER_FREQS"] = {1,qry.sites}["1"];


            logL_by_branch = {};
            utility.SetEnvVariable ("OPTIMIZATION_METHOD", 4);
            utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
            utility.SetEnvVariable ("USE_LAST_VALUES", 0);
            utility.SetEnvVariable ("OPTIMIZATION_PRECISION", 0.001);

            best_AIC = 1e100;

            for (branch = 0; branch < branch_count; branch += 1) {
                this_branch = branch_names [branch];
                (numerical_filter["FILTER_ARRAYS"])[0] = ((((reference_data[segment_name])["conditionals"])[this_branch])["up"])[selector];
                (numerical_filter["FILTER_ARRAYS"])[1] = ((((reference_data[segment_name])["conditionals"])[this_branch])["down"])[selector];

                //console.log (numerical_filter);

    			DataSetFilter numeric_filter_object = CreateFilter (numerical_filter);
                LikelihoodFunction placement_lf    = (numeric_filter_object, screening_tree);
                Optimize					  (res_lf,placement_lf);

                //console.log (res_lf);

                //assert (0);

                logL_by_branch [this_branch] = -2*res_lf[1][0];
                if (logL_by_branch [this_branch] < best_AIC) {
                    best_AIC = logL_by_branch [this_branch];
                    divergence = BranchLength (screening_tree, "Query") + ((reference_data[segment_name])["lengths"]) [this_branch];

                }
            }

            labels = ((reference_data[segment_name])["labels"]);
            logL_by_branch = utility.Map (logL_by_branch, "_logL_", "Exp(0.5*(^'`&best_AIC`'-_logL_))");
            sum = +logL_by_branch;
            raw_names_supported = utility.Map (logL_by_branch, "_logL_", "_logL_/^'`&sum`'");

            supported = {};
            utility.ForEachPair (raw_names_supported, "_branch_", "_support_",
                                                      "(^'`&supported`')[(^'`&labels`')[_branch_]] += _support_");



            best_type = Max (supported);

            supported = utility.Filter (supported, "_support_", "_support_ >= 0.01");
            if (utility.Array1D (supported) == 0) {
                // ensure that there's at least one supported assignment
                supported [best_type["key"]] = best_type["value"];
            }

            placement_result [segment_name] = {"best" : {"type": best_type["key"], "support" : best_type ["value"]},
                                               "credible" : supported,
                                               "divergence" : divergence };

        }

        rearrangements = {};

        options = utility.CatersianProduct (utility.Map (placement_result, "_segment_", "_segment_['credible']"));

        for (k = 0; k < Abs (options); k+=1) {
            weight = 1;
            for (w = 0; w < segment_count; w += 1) {
                weight = weight * ((placement_result[mapped_segments[w]])['credible'])[(options[k])[w]];
            }
            rearrangements [Join (",", options[k])] = weight;
        }

        best_rearrangement = Max (rearrangements);
        rearrangements = utility.Filter (rearrangements, "_support_", "_support_ >= 0.01");



       return {"PHYLO-PLACEMENT" : placement_result, "REARRANGEMENTS" :
                                               {"best" : {"type": best_rearrangement["key"], "support" : best_rearrangement ["value"]},
                                               "credible" : rearrangements }};
    }

    lfunction standardize_label (label) {
        parts = regexp.Split (label, "_");

        assert (Abs (parts) >= 2, "Expected at least two parts for each reference sequence label (_ delimited): `label` is non-conforming");

        new_label = regexp.Replace (parts[0], "^IGH", "");
        for (k = 1; k < Abs (parts) - 1; k += 1) {
            letter_prefix = regexp.Find (parts[k], "^OR+");
            if (None != letter_prefix) {
                new_label += "/" + parts[k];
            } else {
                new_label += "-" + parts[k];
            }
        }
        new_label += "*" + parts[Abs (parts)-1];
        return new_label;
    }

    lfunction merge_labels (current_parent_label, my_label) {
        shared = 0;


        current_parent_label_parts = regexp.Split (current_parent_label, "[\*\-]");
        my_label_parts = regexp.Split (my_label, "[\*\-]");

        upto = Min (Abs(current_parent_label_parts), Abs (my_label_parts));

        while (current_parent_label_parts[shared] == my_label_parts[shared]) {
            shared += 1;
            if (shared >= upto) {
                break;
            }
        }

        shared_label = "";

        if (shared == 0) {
            upto = Min (Abs(current_parent_label_parts[0]), Abs (my_label_parts[0]));
            shared = 0;
            //console.log ("`current_parent_label_parts[0]` `my_label_parts[0]`");
            while ((current_parent_label_parts[0])[shared] == (my_label_parts[0])[shared]) {
                shared += 1;
                if (shared >= upto) {
                    break;
                }
            }
            shared_label = (my_label_parts[0])[0][shared-1];

        } else {
            shared_label = "";
            current_index = 0;
            for (k = 0; k < shared; k+=1) {
                shared_label += current_parent_label_parts[k];
                current_index += Abs(current_parent_label_parts[k]) + (k>0);
                if (k < shared - 1) {
                    shared_label  += current_parent_label[current_index];
                }
            }

        }

        //console.log ("`current_parent_label` x `my_label` => `shared_label`");
        return shared_label;
    }

    lfunction generate_labels (tree_name) {
        tip_avl = (^tree_name) ^ 0;
        labels  = {};
        zero_branch_lock = {};
        lengths = {};

        for (node = 1; node < Abs (tip_avl) - 1; node += 1) {
            child_count = Abs ((tip_avl[node])["Children"]);
            node_name = (tip_avl[node])["Name"];
            if (child_count == 0) {
                my_label = standardize_label (node_name);
                labels [(tip_avl[node])["Name"]] = my_label;
                (tip_avl[node])["Total tips below"] = 1;
                (tip_avl[node])["Total length below"] =  (tip_avl[node])["Length"];
            } else {
                my_label = labels [node_name];
                // compute average branch lengths
                for (child = 0; child < child_count; child += 1) {
                    (tip_avl[node])["Total tips below"] += (tip_avl[((tip_avl[node])["Children"])[child]])["Total tips below"];
                    (tip_avl[node])["Total length below"] += (tip_avl[((tip_avl[node])["Children"])[child]])["Total length below"];
                }
            }

            lengths [node_name] = (tip_avl[node])["Total length below"] / (tip_avl[node])["Total tips below"];

            parent_name = (tip_avl[(tip_avl[node])["Parent"]])["Name"];

            if ((tip_avl[node])["Length"] <= 1e-8) {

                if (zero_branch_lock / parent_name == TRUE && labels [parent_name] != my_label) {
                    labels[parent_name] = merge_labels (labels[parent_name], my_label);
                } else {
                    labels[parent_name] = my_label;
                }



                zero_branch_lock [parent_name] = (tip_avl[node])["Name"];
            }

            if (labels / parent_name) { // need to merge
                if ( zero_branch_lock / parent_name  == FALSE) {

                    current_parent_label = labels[parent_name];
                    shared_label = merge_labels (labels[parent_name], my_label);
                    // split the label

                    assert (Abs(shared_label) > 0, "Reference error: incompatible labels in reference tree `tree_name`; they all MUST share a prefix of at least one letter: `current_parent_label` and `my_label`");

                    labels [parent_name] = shared_label;

                    //console.log ("`my_label` and `current_parent_label` => `shared_label`");
                }
            } else {
                //console.log ("Locking `parent_name` to `my_label`");
                labels [parent_name] = my_label;
            }
        }

        return {"labels" : labels,
                "lengths" : lengths};
    }

    lfunction prepare_alignment (index, path, cache) {
        if (cache / index) {
            return TRUE;
        } else {
            data_path = path + ".hyphy-data";
            if (io.FileExists (data_path) == FALSE) {
                fit_path = path + ".hyphy";
                //GetString (existing_likelihood_functions, LikelihoodFunction,-1);
                lf_count = Rows(LikelihoodFunction);
                if (io.FileExists (fit_path) == FALSE) {
                    namespaced = &alignment;
                    info = alignments.ReadNucleotideDataSet (namespaced, path);
                    tree = utility.GetEnvVariable ("DATAFILE_TREE");
                    assert (Type (tree) == "String" && Abs (tree) > 0, 'Missing a required tree string in `path`');
                    names = alignments.GetSequenceNames (namespaced);
                    DataSetFilter not_mrca = CreateFilter (alignment, 1, "", names[speciesIndex] != 'MRCA');
                    utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);
                    result = estimators.FitGTR_Ext (&not_mrca, trees.LoadAnnotatedTopology(tree), None, {"retain-lf-object": TRUE});

                    utility.ForEachPair (result["global"], "_key_", "_varname_", "parameters.SetConstraint (_varname_['ID'], _varname_['MLE'], '');");
                    Export (result, ^(result["LF"]));
                    fprintf (Min(fit_path,1), CLEAR_FILE, result);
                } else {
                    ExecuteAFile (fit_path);
                }

                GetString   (fitted, LikelihoodFunction, lf_count);
                GetString   (lf_info, ^fitted, -1);
                ConstructCategoryMatrix (data_matrix, ^((lf_info["Trees"])[0]));


                node_count = utility.Array1D (data_matrix["Nodes"]);
                GetDataInfo (duplicate_map, ^((lf_info["Datafilters"])[0]));


                Export (model_string, ^((lf_info["Models"])[0]));


                branch_level_conditionals = {"conditionals" : {},
                                             "model" : model_string};

                utility.Extend (branch_level_conditionals, generate_labels ((lf_info["Trees"])[0]));

                site_count                = utility.Array1D (duplicate_map);

                global_array_index_down  := node_count * pattern_index * 8 +  4 * n;
                global_array_index_up    := global_array_index_down + 4 * node_count;

                for (n = 0; n < node_count; n += 1) {
                    branch_name = (data_matrix["Nodes"])[n];


                    down = {site_count,4};
                    up   = {site_count,4};

                    for (s = 0; s < site_count; s+=1) {
                        pattern_index = duplicate_map[s];
                        for (i = 0; i < 4; i+=1) {
                            down [s][i] = (data_matrix["Values"])[global_array_index_down + i];
                            up   [s][i] = (data_matrix["Values"])[global_array_index_up + i];
                        }
                   }

                   (branch_level_conditionals["conditionals"]) [branch_name] = {"up" : up, "down" : down};
                }


                fprintf (Min(data_path,1), CLEAR_FILE, branch_level_conditionals);
                DeleteObject (^fitted);
            } else {
                fscanf (data_path, "Raw", branch_level_conditionals);
                branch_level_conditionals = Eval (branch_level_conditionals);
            }
        }
        cache [index] = branch_level_conditionals;
        return FALSE;
    }
}

// ---------------------------------------------------------------------------------------------------------
// Reference alignment functions
// ---------------------------------------------------------------------------------------------------------

namespace IgSCUEAL {

     lfunction check_alignment (path, filter, qry) {

        namespace = (&alignment);
        d = alignments.ReadNucleotideDataSet (namespace, path);
        DataSetFilter df = CreateFilter (alignment, 1, filter[siteIndex] > 0);
        fprintf (stdout, df, "\n", qry, "\n");

     }

     lfunction extracted_for_phylomap (mapped_data, ref) {
        extracted = {};
        current_offset = 0;
        up_to = 0;
        qry_string = mapped_data["QRY"];

         for (k = 0; k < Abs (ref); k+=1) {
           if (mapped_data["SEGMENTS"] / k) {
                up_to += (ref[k])["SITES"];
                filter_def = {};
                reduced_qry = ""; reduced_qry * Abs (qry_string);
                last_index = 0;
                for (s = current_offset; s < up_to; s += 1) {
                    if ((mapped_data["reduced"])[s]) {
                        filter_def [s-current_offset] = TRUE;
                        reduced_qry * (qry_string[(mapped_data["mapping"])[s]]);
                    }
                }


                reduced_qry * 0;

                if (Abs (reduced_qry)) {
                    extracted[k] = {"filter_def" : Eval (filter_def) , "mapped_qry" : reduced_qry};
                }

                current_offset = up_to;
           }
        }
        return extracted;
    }

     lfunction annotate (mapped_data, ref, code, d_region) {

        if (None != mapped_data) {
            current_offset = 0;
            annotation = {};
            locations  = {};

            for (k = 0; k < Abs (ref); k+=1) {
               if (mapped_data["SEGMENTS"] / k) {
                   fragments = ((ref[k])["FRAGMENTS"]);
                   if (Type (fragments) == "AssociativeList"){
                        fragment_names = utility.Keys (fragments);
                        fragment_count = Abs (fragments);

                        for (f = 0; f < fragment_count; f += 1) {
                            this_fragment = fragment_names[f];
                            this_fragment_span = fragments [this_fragment];

                            from = (mapped_data["mapping"])[current_offset + this_fragment_span[0]];
                            to = (mapped_data["mapping"])[current_offset + this_fragment_span[1]];

                            locations[this_fragment] = {{from__,to__}};

                            if (from < to) {
                                annotation[this_fragment] = strip_in_frame_indels ((mapped_data["QRY"])[from][to]);
                            } else {
                                annotation[this_fragment] = "";
                            }

                            annotation[this_fragment + "_AA"] = alignments.TranslateCodonsToAminoAcids (annotation[this_fragment], 0, code);
                        }
                   }
                   current_offset += (ref[k])["SITES"];
                }
            }

            if (locations / "JUNCTION_CDR3" && locations / "JUNCTION_J") {
                junction = (mapped_data["QRY"])[(locations["JUNCTION_CDR3"])[0]][(locations["JUNCTION_J"])[0]+2];
                annotation ["CDR3"] = strip_in_frame_indels(junction [3][Abs (junction)-4]);
                annotation ["CDR3_AA"] = alignments.TranslateCodonsToAminoAcids (annotation ["CDR3"], 0, code);
                annotation ["JUNCTION"] = strip_in_frame_indels(junction);
                annotation ["JUNCTION_AA"] = alignments.TranslateCodonsToAminoAcids (annotation["JUNCTION"], 0, code);

                if (None != d_region) {
                    _dAlignOptions = {};
                    _dAlignOptions ["SEQ_ALIGN_CHARACTER_MAP"]="ACGT";


                    _dAlignOptions ["SEQ_ALIGN_SCORE_MATRIX"]   = 	{
                        {5,-4,-4,-4}
                        {-4,5,-4,-4}
                        {-4,-4,5,-4}
                        {-4,-4,-4,5}
                    };

                    _dAlignOptions ["SEQ_ALIGN_GAP_OPEN"]		= 	10;
                    _dAlignOptions ["SEQ_ALIGN_GAP_OPEN2"]	    = 	10;
                    _dAlignOptions ["SEQ_ALIGN_GAP_EXTEND"]	    = 	1;
                    _dAlignOptions ["SEQ_ALIGN_GAP_EXTEND2"]	= 	1;
                    _dAlignOptions ["SEQ_ALIGN_AFFINE"]		    =   1;
                    _dAlignOptions ["SEQ_ALIGN_NO_TP"]		    =   1;
                    _dAlignOptions ["code"] = code;
                    _dAlignOptions ["E"]    = 0.1;
                    _dAlignOptions ["REPORT_VARIANTS" ] = TRUE;

                    d = align_sequence_to_reference_set (junction, d_region, _dAlignOptions);
                     if (None != d) {
                        annotation ["D_ALLELE"] = d ["REFNAME"];
                    }
                }
            }

            utility.Extend (mapped_data, {"ANNOTATION" : annotation});
        }

        return mapped_data;
     }

     lfunction handle_sequence_alignment (seq, ref, options) {

        direct = handle_sequence_alignment_1 (seq, ref, options);
        reverse = handle_sequence_alignment_1 (compute_rc (seq), ref, options);

        if (None != direct) {
            if (None != reverse) {
                if (direct["SCORE"] < reverse["SCORE"]) {
                    return reverse;
                }
            }
            return direct
        } else {
            if (None != reverse) {
                return reverse;
            }
        }

        return None;
    }

   lfunction handle_sequence_alignment_1 (seq, ref, options) {
        seq_processed = strip_gaps (seq);
        reference_count = Abs (ref);
        segments = {};
        current_sequence = seq_processed['stripped'];
        mapped_segments = {};
        for (i = 0; i < reference_count; i+=1) {
            segments [i] = align_sequence_to_reference_set (current_sequence, ref[i], options);
            if (None != segments[i]) {

                mapped_reference_sequence = (ref[i])[(segments[i]) ["REFNAME"]];


                mapped_segments [i] =  {"index" : i,
                                   "reference" : (segments[i]) ["REFNAME"],
                                   "span" : (segments[i]) ["SPAN"]};

                if ((segments[i])["SUFFIX"]) {
                    current_sequence = (segments[i])["SUFFIX"];
                    continue;
                }
                break;
            }
        }
        if (Abs (mapped_segments)) {
            // paste all the mapped segments together
            ref_aligned  = Join ("",utility.Map (mapped_segments, "_value_", "(((`&ref`)[_value_['index']])[_value_['reference']])['aligned']"));
            ref_segments = Join ("",utility.Map (mapped_segments, "_value_", "(((`&ref`)[_value_['index']])[_value_['reference']])['stripped']"));
            joint_references = align_sequence_to_reference_set (seq_processed['stripped'], {"ref" : {"stripped" : ref_segments}}, options);

            if (None != joint_references) {
                utility.Extend (joint_references, {"SEGMENTS": mapped_segments});
                utility.Extend (joint_references, alignments.MapAlignmentToReferenceCoordinates (ref_aligned, joint_references["REF"], joint_references["QRY"], joint_references["OFFSET"]));
                return joint_references;
            }

        }

        return None;
    }

    // -------------------------------------------------------------------------- //

    lfunction  compute_rc (sequence) {

        nucleotide_rc = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "M": "K",
            "R": "Y",
            "W": "W",
            "S": "S",
            "Y": "R",
            "K": "M",
            "B": "V",  /* not A */
            "D": "H",  /* not C */
            "H": "D",  /* not G */
            "V": "B",  /* not T */
            "N": "N"
        };

        _seqOut = ""; _seqOut*128;
        _seqL   = Abs(sequence);
        for (_r = _seqL-1; _r >=0 ; _r = _r-1) {
            _seqOut * (nucleotide_rc[sequence[_r]]);
        }
        _seqOut*0;
        return _seqOut;
    }

    // -------------------------------------------------------------------------- //

    lfunction sanitize_sequence (sequence) {
        return strip_non_letters (sequence && 1)['stripped'];
    }

    lfunction strip_in_frame_indels (sequence) {
        return sequence ^ {{"-{3}",""}};
    }

    lfunction  strip_gaps (sequence) {
        return {'aligned': sequence, 'stripped' : sequence ^ {{"\-",""}}};
    }

    lfunction  strip_non_letters (sequence) {
        return {'aligned': sequence, 'stripped' : sequence ^ {{"[^A-Z,a-z]",""}}};
    }

    lfunction ensure_required_partitions (list, partitions) {
         assert (Type (partitions) == "Matrix", "Missing required NEXUS partitions");
         loaded_partitions = {};
         for (k = 0; k < Rows (partitions); k+=1) {
            loaded_partitions [partitions[k][0]] = regexp.Split (partitions[k][1], "-");
         }
         req_count = utility.Array1D (list);

         for (k = 0; k < req_count; k+=1) {
            assert (loaded_partitions / list[k], "Missing required partition '`list[k]`'");
            assert (Abs (loaded_partitions [list[k]]) == 2, "Invalid required partition specification '`loaded_partitions [list[k]]`'");
         }

         return {"FRAGMENTS" : Eval (loaded_partitions)};
    }

    // -------------------------------------------------------------------------- //

    lfunction  load_reference_alignment (file_path, filter, label) {
       namespace = (&alignment);

       data  = alignments.ReadCodonDataSetFromPathGivenCode (namespace, file_path,
                                                           (^"IgSCUEAL.universal_code")["code"], (^"IgSCUEAL.universal_code")["stops"]);




       align_with_these = utility.Filter (alignments.GetSequenceNames (namespace), "_name_", "regexp.Find (_name_, `&filter`)");
       alignments.GetSequenceByName (namespace, None);

       result = {'SITES' : data["sites"]};

       utility.ForEach (align_with_these, "_value_",
        "`&result`[_value_] = IgSCUEAL.strip_gaps (alignments.GetSequenceByName ('`namespace`', _value_));"
        );



        if (label == "V") {
            utility.Extend (result, ensure_required_partitions ({{"FW1","CDR1","FW2","CDR2","FW3","CDR3","JUNCTION_CDR3"}}, data["partitions"]));
        }
        if (label == "J") {
            utility.Extend (result, ensure_required_partitions ({{"J","JUNCTION_J"}}, data["partitions"]));
        }
        if (label == "C") {
            utility.Extend (result, ensure_required_partitions ({{"CH"}}, data["partitions"]));
        }
       return  result;
    }

     // -------------------------------------------------------------------------- //

    lfunction  load_d_region (file_path) {
       namespace = (&alignment);
       d = alignments.ReadNucleotideDataSet (namespace, file_path);
       d = utility.Map (alignments.GetSequenceNames (namespace), "_name_", "_name_");
       result = {};

       alignments.GetSequenceByName (namespace, None);

       utility.ForEach (d, "_value_",
        "`&result`[_value_] = IgSCUEAL.strip_non_letters (alignments.GetSequenceByName ('`namespace`', _value_));"
        );

       return  result;
    }

    // ---------------------------------------------------------------------------------------------------------

    lfunction align_sequence_to_reference_set (seq, references, alignment_settings, seq_name) {

        ref_names = utility.Keys (references);
        ref_count = utility.Array1D (references);

        input     = {1,2};
        input [1] = seq;

        overall = {'SCORE' : -1e100, 'REFNAME' : ''};

        for (k = 0; k < ref_count; k+=1) {
            if (Type (references[ref_names[k]]) == "AssociativeList") {
                if (references[ref_names[k]] / "stripped") {
                    input [0] = (references[ref_names[k]])["stripped"];

                    AlignSequences (result, input, alignment_settings);

                    result = result[0];

                    if (result [0] >= overall['SCORE']) {

                        if (alignment_settings ["REPORT_VARIANTS"] && result [0] == overall['SCORE'] && Abs (overall['REFNAME'])) {
                            overall['REFNAME'] += "|" + ref_names[k];
                        } else {
                            overall['REFNAME'] = ref_names[k];
                        }
                        overall['SCORE'] = result[0];
                        overall['RAW-REF'] = result[1];
                        overall['RAW-QRY'] = result[2];
                    }
                }
            }
        }


        if (Type (overall['RAW-REF']) == "String") {

            computed_score = (overall["SCORE"] - 30 * alignment_settings["MATCH"] * Exp (-Abs(seq)/3) ) / Abs (seq) * 3 ;

            if (alignment_settings["E"] <= computed_score) {
                utility.Extend (overall, correctReadUsingCodonAlignedData (overall['RAW-REF'], overall['RAW-QRY'], alignment_settings["code"]));


                if (alignment_settings["SEQ_ALIGN_CODON_ALIGN"] == TRUE && overall["SPAN"] <= 3) {
                    //assert (0, "Internal error in align_sequence_to_reference_set" + overall + "\nComputed score `computed_score`; expected score " + alignment_settings["E"] + "; match score " +  alignment_settings["MATCH"] + "\nInput sequence: `seq`");
                    return None;
                }
                return overall;
            }
        }


        return None;
     }

    // -------------------------------------------------------------------------- //

lfunction define_alignment_settings (code) {

        igh_matrix = {
            {                 6,                -3,                -4,                -4,                -2,                -2,                -2,                -1,                -3,                -3,                -3,                -2,                -2,                -4,                -2,                 0,                -1,                -5,                -3,                -1,                -4,                -2,                -2,                -7}
            {                -3,                 8,                -2,                -4,                -6,                 0,                -2,                -5,                -2,                -6,                -4,                 1,                -3,                -5,                -4,                -2,                -3,                -5,                -3,                -5,                -3,                -1,                -2,                -7}
            {                -4,                -2,                 8,                 0,                -5,                -1,                -2,                -2,                 0,                -6,                -6,                -1,                -4,                -5,                -4,                 0,                -1,                -7,                -4,                -5,                 6,                -2,                -2,                -7}
            {                -4,                -4,                 0,                 8,                -6,                -2,                 0,                -3,                -3,                -5,                -6,                -2,                -6,                -6,                -3,                -1,                -3,                -7,                -6,                -6,                 6,                 0,                -3,                -7}
            {                -2,                -6,                -5,                -6,                10,                -5,                -7,                -5,                -5,                -3,                -3,                -6,                -3,                -5,                -5,                -2,                -2,                -4,                -4,                -2,                -5,                -6,                -4,                -7}
            {                -2,                 0,                -1,                -2,                -5,                 8,                 1,                -4,                 0,                -6,                -4,                 0,                -1,                -6,                -3,                -1,                -2,                -3,                -3,                -4,                -1,                 6,                -2,                -7}
            {                -2,                -2,                -2,                 0,                -7,                 1,                 7,                -4,                -1,                -6,                -5,                 0,                -4,                -6,                -3,                -1,                -2,                -5,                -4,                -4,                 0,                 6,                -2,                -7}
            {                -1,                -5,                -2,                -3,                -5,                -4,                -4,                 7,                -4,                -7,                -6,                -3,                -5,                -5,                -4,                -2,                -4,                -4,                -5,                -6,                -2,                -4,                -4,                -7}
            {                -3,                -2,                 0,                -3,                -5,                 0,                -1,                -4,                10,                -6,                -5,                -2,                -3,                -3,                -4,                -2,                -4,                -5,                 0,                -6,                -1,                -1,                -3,                -7}
            {                -3,                -6,                -6,                -5,                -3,                -6,                -6,                -7,                -6,                 6,                 0,                -5,                 0,                -1,                -5,                -5,                -2,                -5,                -3,                 2,                -5,                -6,                -2,                -7}
            {                -3,                -4,                -6,                -6,                -3,                -4,                -5,                -6,                -5,                 0,                 6,                -5,                 1,                -1,                -5,                -5,                -3,                -3,                -3,                 0,                -6,                -5,                -2,                -7}
            {                -2,                 1,                -1,                -2,                -6,                 0,                 0,                -3,                -2,                -5,                -5,                 7,                -3,                -6,                -2,                -1,                -2,                -5,                -3,                -4,                -2,                 0,                -2,                -7}
            {                -2,                -3,                -4,                -6,                -3,                -1,                -4,                -5,                -3,                 0,                 1,                -3,                 9,                -1,                -5,                -3,                -2,                -3,                -3,                 0,                -5,                -2,                -1,                -7}
            {                -4,                -5,                -5,                -6,                -5,                -6,                -6,                -5,                -3,                -1,                -1,                -6,                -1,                 8,                -6,                -4,                -4,                 0,                 1,                -3,                -6,                -6,                -3,                -7}
            {                -2,                -4,                -4,                -3,                -5,                -3,                -3,                -4,                -4,                -5,                -5,                -2,                -5,                -6,                 9,                -2,                -3,                -6,                -5,                -4,                -4,                -3,                -4,                -7}
            {                 0,                -2,                 0,                -1,                -2,                -1,                -1,                -2,                -2,                -5,                -5,                -1,                -3,                -4,                -2,                 7,                 0,                -5,                -3,                -4,                -1,                -1,                -2,                -7}
            {                -1,                -3,                -1,                -3,                -2,                -2,                -2,                -4,                -4,                -2,                -3,                -2,                -2,                -4,                -3,                 0,                 7,                -4,                -3,                -1,                -2,                -2,                -2,                -7}
            {                -5,                -5,                -7,                -7,                -4,                -3,                -5,                -4,                -5,                -5,                -3,                -5,                -3,                 0,                -6,                -5,                -4,                12,                 0,                -6,                -7,                -4,                -4,                -7}
            {                -3,                -3,                -4,                -6,                -4,                -3,                -4,                -5,                 0,                -3,                -3,                -3,                -3,                 1,                -5,                -3,                -3,                 0,                 9,                -3,                -5,                -3,                -2,                -7}
            {                -1,                -5,                -5,                -6,                -2,                -4,                -4,                -6,                -6,                 2,                 0,                -4,                 0,                -3,                -4,                -4,                -1,                -6,                -3,                 6,                -6,                -4,                -2,                -7}
            {                -4,                -3,                 6,                 6,                -5,                -1,                 0,                -2,                -1,                -5,                -6,                -2,                -5,                -6,                -4,                -1,                -2,                -7,                -5,                -6,                 7,                -1,                -3,                -7}
            {                -2,                -1,                -2,                 0,                -6,                 6,                 6,                -4,                -1,                -6,                -5,                 0,                -2,                -6,                -3,                -1,                -2,                -4,                -3,                -4,                -1,                 7,                -2,                -7}
            {                -2,                -2,                -2,                -3,                -4,                -2,                -2,                -4,                -3,                -2,                -2,                -2,                -1,                -3,                -4,                -2,                -2,                -4,                -2,                -2,                -3,                -2,                -2,                -7}
            {                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                 1}
        };

        max_score = Max (igh_matrix,0);
        penalty = Max (max_score, -Min (igh_matrix,0));


        options = {
            "SEQ_ALIGN_CHARACTER_MAP" : "ACGTN",
            "SEQ_ALIGN_GAP_OPEN"	 : 	1.5*penalty,
            "SEQ_ALIGN_GAP_OPEN2"	 :	1.5*penalty,
            "SEQ_ALIGN_GAP_EXTEND"   : 	0.1*penalty,
            "SEQ_ALIGN_GAP_EXTEND2"  : 	0.1*penalty,
            "SEQ_ALIGN_FRAMESHIFT"   :  2*penalty,
            "SEQ_ALIGN_NO_TP"        :  1,
            "SEQ_ALIGN_AFFINE"       :  1,
            "SEQ_ALIGN_CODON_ALIGN"  :  1
        };

       base_frequencies = {
                                {   0.069352301}
                                {   0.021000333}
                                {   0.049862283}
                                {   0.026563029}
                                {   0.033026253}
                                {     0.1054545}
                                {   0.008970554}
                                {   0.036648952}
                                {   0.036329907}
                                {   0.065164842}
                                {   0.021805663}
                                {   0.032992805}
                                {     0.0348093}
                                {   0.036818766}
                                {   0.054168098}
                                {    0.12149678}
                                {   0.082464001}
                                {   0.053564744}
                                {   0.038383113}
                                {    0.07112377}
                            };



        options["SEQ_ALIGN_SCORE_MATRIX"] = pSM2cSM(igh_matrix, "ACDEFGHIKLMNPQRSTVWY", code["code"], code["ordering"]);
        shift_penalty = computeExpectedPerBaseScore (.5,igh_matrix,base_frequencies);

        _cdnaln_partialScoreMatrices = cSM2partialSMs(options["SEQ_ALIGN_SCORE_MATRIX"],
                {{shift_penalty__*1.5,shift_penalty__,shift_penalty__,shift_penalty*1.5}});


        options ["SEQ_ALIGN_PARTIAL_3x1_SCORES"] = _cdnaln_partialScoreMatrices["3x1"];
        options ["SEQ_ALIGN_PARTIAL_3x2_SCORES"] = _cdnaln_partialScoreMatrices["3x2"];
        options ["SEQ_ALIGN_PARTIAL_3x4_SCORES"] = _cdnaln_partialScoreMatrices["3x4"];
        options ["SEQ_ALIGN_PARTIAL_3x5_SCORES"] = _cdnaln_partialScoreMatrices["3x5"];

        options ["MATCH"] = computeExpectedPerBaseScore (1,igh_matrix,base_frequencies);
        options ["E"] = Max (0.1, computeExpectedPerBaseScore (.4,igh_matrix,base_frequencies));

        return options;
}

lfunction cSM2partialSMs(_cdnScoreMatrix, penalties) {

        m3x5  =  { 126, 1250 };
        m3x4  =  { 126, 500 };
        m3x2  =  { 126,  75 };
        m3x1  =  { 126,  15 };


        if (utility.Array1D (penalties) == 4) {
            p3x1 = penalties [0];
            p3x2 = penalties [1];
            p3x4 = penalties [2];
            p3x5 = penalties [3];
        } else {
            p3x5 = 0;
            p3x4 = 0;
            p3x2 = 0;
            p3x1 = 0;
        }

        for ( thisCodon = 0; thisCodon < 64; thisCodon += 1 ) {
            for ( d1 = 0; d1 < 5; d1 += 1 ) {
                max100 = -1e100;
                max010 = -1e100;
                max001 = -1e100;

                for ( d2 = 0; d2 < 5; d2 += 1 ) {
                    partialCodon = 5 * d1 + d2;
                    max110 = -1e100;
                    max101 = -1e100;
                    max011 = -1e100;

                    for ( d3 = 0; d3 < 5; d3 += 1 ) {
                        thisCodon2 = 5 * partialCodon + d3;
                        thisScore = _cdnScoreMatrix[ thisCodon ][ thisCodon2 ];

                        // this is the trivial and stupid way of doing it, but it should work
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 0 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 1 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 2 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 3 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 4 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 5 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 6 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 7 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 8 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 9 ] = thisScore - p3x5;

                        m3x4[ thisCodon ][ 4 * thisCodon2 + 0 ] = thisScore - p3x4;
                        m3x4[ thisCodon ][ 4 * thisCodon2 + 1 ] = thisScore - p3x4;
                        m3x4[ thisCodon ][ 4 * thisCodon2 + 2 ] = thisScore - p3x4;
                        m3x4[ thisCodon ][ 4 * thisCodon2 + 3 ] = thisScore - p3x4;

                        // d1 is 1
                        max100 = Max( max100, _cdnScoreMatrix[ thisCodon ][ 25 * d1 + 5 * d2 + d3 ] );
                        max010 = Max( max010, _cdnScoreMatrix[ thisCodon ][ 25 * d2 + 5 * d1 + d3 ] );
                        max001 = Max( max001, _cdnScoreMatrix[ thisCodon ][ 25 * d2 + 5 * d3 + d1 ] );

                        // d1 and d2 are 1
                        max110 = Max( max110, _cdnScoreMatrix[ thisCodon ][ 25 * d1 + 5 * d2 + d3 ] );
                        max101 = Max( max101, _cdnScoreMatrix[ thisCodon ][ 25 * d1 + 5 * d3 + d2 ] );
                        max011 = Max( max011, _cdnScoreMatrix[ thisCodon ][ 25 * d3 + 5 * d1 + d2 ] );
                    }

                    m3x2[ thisCodon ][ 3 * partialCodon + 0 ] = max110 - p3x2;
                    m3x2[ thisCodon ][ 3 * partialCodon + 1 ] = max101 - p3x2;
                    m3x2[ thisCodon ][ 3 * partialCodon + 2 ] = max011 - p3x2;
                }

                m3x1[ thisCodon ][ 3 * d1 + 0 ] = max100 - p3x1;
                m3x1[ thisCodon ][ 3 * d1 + 1 ] = max010 - p3x1;
                m3x1[ thisCodon ][ 3 * d1 + 2 ] = max001 - p3x1;
            }
        }


        return { "3x1": m3x1, "3x2": m3x2, "3x4": m3x4, "3x5": m3x5 };
    }



    lfunction _digits (index) {
        return {{index__$25, index__ % 25 $ 5, index__ % 5}};
    };

    lfunction _hazN (digits) {
        return digits [0] == 4 || digits [1] == 4 || digits [2] == 4;
    };

    lfunction _map_to_nuc (digits) {
        return digits [0] * 16 + digits [1] * 4 + digits[2];
    };

    lfunction _generate_resolutions (digits) {
        resolutions = {};
        for (k = 0; k < 4; k += 1) {
            try = digits;
            if (digits[0] == 4) {
                try [0] = k;
            }
            for (k2 = 0; k2 < 4; k2 += 1) {
                if (digits[1] == 4) {
                    try [1] = k;
                }
                for (k3 = 0; k3 < 4; k3 += 1) {
                    if (digits[2] == 4) {
                        try [2] = k;
                    }
                    resolutions[_map_to_nuc (try)] = 1;
                }
            }
        }
        return Rows(resolutions);
    };

    // -------------------------------------------------------------------------- //

    lfunction pSM2cSM (_scorematrix, _letters, code, ordering) {

        _cdnScoreMatrix  = { 126 , 126 };
        _mapping      = utility.MapStrings ( ordering, _letters );


        for ( _k = 0; _k < 125; _k += 1 ) {

            letters1 = _digits (_k);

            if (_hazN (letters1)) {
                 letters1 = _generate_resolutions (letters1);
                  for ( _k2 = _k; _k2 < 125 ; _k2 += 1 ) {
                    letters2 = _digits (_k2);
                    if (_hazN (letters2) == 0) {
                        codon2 = _map_to_nuc(letters2);
                        _mappedK2 = _mapping[ code[ codon2 ] ];
                        _aScore = -1e4;
                        if (_mappedK2 >= 0) {

                            res_count = utility.Array1D (letters1);
                            for (r = 0; r < res_count; r += 1) {
                                resolution_codon = 0 + letters1[r];
                                resolution_aa = _mapping[ code[ resolution_codon  ] ];
                                if (resolution_aa >= 0) {
                                    try_score = _scorematrix[ resolution_aa ][ _mappedK2 ] - 1;
                                    if (resolution_aa == _mappedK2 && codon2 != resolution_codon) {
                                        try_score = try_score - 1;
                                    }
                                    _aScore = Max (_aScore, try_score);
                                }
                            }
                        }


                        _cdnScoreMatrix[ _k ][ _k2 ] = _aScore;
                        _cdnScoreMatrix[ _k2 ][ _k ] = _aScore;
                    }
                 }
            } else {
                codon1 = _map_to_nuc(letters1);

                _mappedK = _mapping[ code[ codon1 ] ];

                if ( _mappedK >= 0) {
                    for ( _k2 = _k; _k2 < 125 ; _k2 += 1 ) {
                        letters2 = _digits (_k2);

                        if (_hazN (letters2)) {
                           _aScore = -1e4;
                           letters2 = _generate_resolutions (letters2);
                            res_count = utility.Array1D (letters2);
                            for (r = 0; r < res_count; r += 1) {
                                resolution_codon = 0 + letters2[r];
                                resolution_aa = _mapping[ code[ resolution_codon  ] ];
                                if (resolution_aa >= 0) {
                                    try_score = _scorematrix[ _mappedK ][ resolution_aa ] - 1;
                                    if (resolution_aa == _mappedK && codon1 != resolution_codon) {
                                        try_score = try_score - 1;
                                    }
                                    _aScore = Max (_aScore, try_score);
                                }
                            }
                            _cdnScoreMatrix[ _k ][ _k2 ] = _aScore;
                            _cdnScoreMatrix[ _k2 ][ _k ] = _aScore;
                            continue;
                        }

                        codon2 = _map_to_nuc(letters2);

                        _mappedK2 = _mapping[ code[ codon2 ] ];
                        if ( _mappedK2 >= 0 ) {
                            _aScore = _scorematrix[ _mappedK ][ _mappedK2 ];
                            if ( _mappedK == _mappedK2 && _k2 > _k ) {
                                _aScore = _aScore - 1; // synonymous match
                            }
                        } else {
                            // stop codons don't match anything
                            _aScore = -1e4;
                        }
                        _cdnScoreMatrix[ _k ][ _k2 ] = _aScore;
                        _cdnScoreMatrix[ _k2 ][ _k ] = _aScore;
                    }
                } else { // stop codons here
                    for ( _k2 = _k; _k2 < 125; _k2 += 1 ) {

                        letters2 = _digits (_k2);

                        if (_hazN (letters2)) {
                            continue;
                        }

                        codon2 = _map_to_nuc(letters2);

                        _mappedK2 = _mapping[ code[ codon2 ] ];

                        if ( _mappedK2 < 0 ) {
                            // don't penalize stop codons matching themselves
                            _cdnScoreMatrix[ _k ][ _k2 ] = 0;
                            _cdnScoreMatrix[ _k2 ][ _k ] = 0;
                        } else {
                            _cdnScoreMatrix[ _k ][ _k2 ] = -1e4;
                            _cdnScoreMatrix[ _k2 ][ _k ] = -1e4;
                        }
                    }
                }
            }
        }

        return _cdnScoreMatrix;
    }

    // -------------------------------------------------------------------------- //

    lfunction computeExpectedPerBaseScore( _expectedIdentity, _cdnaln_scorematrix, _cdnaln_base_freqs ) {
        meanScore = 0;

        for (_aa1 = 0; _aa1 < 20; _aa1 += 1) {
            for (_aa2 = 0; _aa2 < 20; _aa2 += 1) {
                if ( _aa1 != _aa2 ) {
                    meanScore += ( 1 - _expectedIdentity ) * _cdnaln_scorematrix[_aa1][_aa2] * _cdnaln_base_freqs[_aa1] * _cdnaln_base_freqs[_aa2];
                } else {
                    meanScore += _expectedIdentity * _cdnaln_scorematrix[_aa1][_aa1] * _cdnaln_base_freqs[_aa1] * _cdnaln_base_freqs[_aa1];
                }
            }
        }

        return meanScore;
    }


    // -------------------------------------------------------------------------- //

    lfunction	computeCorrection (str) {
        result = {1,2};

        result[0]	 = (str$"^\\-+")[1]+1;
        result[1]	 = (str$"\\-+$")[0];

        if (result[1] >= 0) {
            result[1] = Abs(str)-result[1];
        }
        else {
            result[1] = 0;
        }
        return result;
    }

    // -------------------------------------------------------------------------- //

    lfunction correctReadUsingCodonAlignedData (ref, qry, code) {
        reference_shifts = computeCorrection(ref);

        /*reference_shifts is the starting,ending nucleotide on the reference relative to the read. if reference is longer than the read, then both are 0*/

        offsetFrom = (qry$"^\\-+")[1]+1;
        offsetTo   = (qry$"\\-+$")[0]-1;

        /* the $ looks for the regular expression in bestAl[2] and returns a 2x1 array with the starting and ending 0-based positions of the regular expression. in this case multiple indels, -. returns -1 for both if the regular expression is not found.
            i.e. 0-based index leading indels start at (bestAl[2]$"^\\-+")[0] and end at (bestAl[2]$"^\\-+")[1]; trailing indels start at (bestAl[2]$"\\-+$")[0] and end at (bestAl[2]$"\\-+$")[0];

            so offSetFrom to offSetTo will return the reference sequence co-ordinates overlapping with the read.
        */


        if (offsetTo < 0) {
            offsetTo = Abs(qry)-1; /*if no trailing indels then to end of read*/
        }

        // check to see if the prefix in REF has some out-of-frame indels so that we can start offsetFrom in the correct frame */

        if (offsetFrom > 0) {
            frame_skips = 0;
            for (i = 0; i < offsetFrom; i += 1) {
               if ((ref[i] && 1) != ref[i]) {
                    frame_skips += 1;
               }
            }
            offsetFrom += -frame_skips;
        }

        seqOffset  = offsetFrom;          /*set the offset of the read relative to the reference. ie the number of indels needed on the read to align to the reference */
        offsetFrom +=  reference_shifts[0];           /*if the read starts before the reference then shift to start of reference ie. by reference_shifts[0] */
        offsetTo    =  offsetTo	- reference_shifts[1];           /*if the read extends beyond the reference then shift to end of reference ie. by reference_shifts[1] */

        theSeq     = qry;
        theSeq	   = qry[reference_shifts[0]][Abs(qry)-reference_shifts[1]-1]; /*the nucleotide sequence of the read that overlaps with the reference sequence */

        nucSeq	   = qry[offsetFrom][offsetTo]; /*read sequence pruned to exactly overlapping region*/
        nucSeqRef  = ref[offsetFrom][offsetTo]; /*reference sequence pruned to exactly overlapping region*/

        extra_keys = {};

        if (reference_shifts[0] > 0) { // 'qry' has a prefix that's not aligned to the reference
            extra_keys['PREFIX'] = qry[0][reference_shifts[0]-1];
        }
        if (reference_shifts[1] > 0) {
            l = Abs (qry);
            extra_keys['SUFFIX'] = qry[l-reference_shifts[1]][l-1];
        }


        return utility.Extend (igg_alignment_cleanup (nucSeqRef, nucSeq,seqOffset, code), extra_keys);
    }

    // -------------------------------------------------------------------------- //

    lfunction igg_alignment_cleanup (reference, query, offset_nuc, code) {

        too_short = 0;
        too_long  = 0;
        span      = 0; // how many nucleotides in the reference were covered by non-gaps
        _seqL     = Abs (reference);

        //console.log("\nSEQL");
        //console.log(_seqL);

        ref_cleaned = ""; ref_cleaned * 128;
        qry_cleaned = ""; qry_cleaned * 128;

        _codon_in_reference = 0;

        for ( _rcidx = 0; _rcidx < _seqL; _rcidx += 1 ) {
            _del1 = reference [_rcidx] != (reference [_rcidx]&&1);
            if (_del1) {
                too_short += 1;
                _codon_in_reference += 1;
                ref_cleaned * (reference [_rcidx]&&1);
                qry_cleaned * (query [_rcidx]&&1);
            } else {
                _del1 = query [_rcidx] != (query [_rcidx]&&1);
                if (_del1) {
                    if (_seqL-_rcidx < 3 && _codon_in_reference % 3 == 0) {
                        break;
                    }
                    too_long += 1;
                } else {
                    ref_cleaned * (reference [_rcidx]&&1);
                    qry_cleaned * (query [_rcidx]&&1);
                    span += 1;
                    _codon_in_reference +=1;
                }
            }
        }
        ref_cleaned * 0; qry_cleaned * 0;

        return {"REF": ref_cleaned, "QRY": qry_cleaned, "TOO_SHORT" : too_short, "TOO_LONG": too_long, "SPAN": span, "OFFSET_AA" :  offset_nuc$3 + (offset_nuc % 3 > 0),"OFFSET" :  offset_nuc, "AA" : alignments.TranslateCodonsToAminoAcids (qry_cleaned, (3-offset_nuc%3)%3, code), "AA_REF" : alignments.TranslateCodonsToAminoAcids (ref_cleaned, (3-offset_nuc%3)%3, code)};
    }
}
