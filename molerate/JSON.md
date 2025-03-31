# Molerate JSON Schema Description


The top-level JSON object contains the following keys:

* **`analysis`** (Object): Contains metadata about the analysis performed.
    * `authors` (String): Name of the author(s).
    * `contact` (String): Contact email address.
    * `info` (String): Brief description of the analysis.
    * `labeling strategy` (String): The strategy used for labeling branches (e.g., "all-descendants").
    * `model` (String): The substitution model used (e.g., "WAG").
    * `rate variation` (String): The model used for rate variation across sites (e.g., "GDD").
    * `requirements` (String): Description of the input data required.
    * `version` (String): Version of the analysis tool.

* **`branch attributes`** (Object): Contains attributes for phylogenetic branches, grouped by partition (e.g., "0") and model.
    * `<partition_id>` (Object, e.g., "0"): Represents a data partition.
        * `<branch_name>` (Object, e.g., "HLacaPus1", "Node1"): Represents a specific branch in the tree.
            * `Proportional` (Number): Estimated branch length for the Proportional model.
            * `Proportional Partitioned` (Number): Estimated branch length for the Proportional Partitioned model.
            * `Reference` (Number): Branch length from the reference tree.
            * `Unconstrained Test` (Number): Estimated branch length for the Unconstrained Test model.
            * *( Potentially other model names as keys with numeric values )*
    * `attributes` (Object): Defines metadata for the models listed under branch names.
        * `<model_name>` (Object, e.g., "Proportional"):
            * `attribute type` (String): Type of attribute (e.g., "branch length").
            * `display order` (Number): Suggested order for displaying this model's results.

* **`branch level analysis`** (Object): Contains detailed likelihood ratio test results for individual branches designated as 'test'.
    * `<branch_name>` (Object, e.g., "HLamaAes1"): Represents a specific 'test' branch.
        * `LogL` (Number): Log-likelihood for the `Proportional+1` model at this branch.
        * `alternative` (Object): `Proportional+1` vs `Unconstrained Test`
            * `Corrected P-value` (Number): Corrected (Holm-Bonferroni) p-value.
            * `LRT` (Number): Likelihood Ratio Test statistic.
            * `p-value` (Number): Uncorrected p-value.
        * `null` (Object):  `Proportional` vs `Proportional+1`
            * `Corrected P-value` (Number): orrected (Holm-Bonferroni) p-value.
            * `LRT` (Number): Likelihood Ratio Test statistic.
            * `p-value` (Number): Uncorrected p-value.

* **`fits`** (Object): Contains details about the statistical fit of each evolutionary model tested.
    * `<model_name>` (Object, e.g., "Proportional", "Unconstrained Test"):
        * `AIC-c` (Number): Corrected Akaike Information Criterion value.
        * `Log Likelihood` (Number): Log-likelihood value for the model fit.
        * `Rate Distributions` (Object): Parameters related to site rate variation.
            * *( Keys vary based on the rate variation model, e.g., GDD category rates, mixture weights, branch scalers )* (Number)
        * `display order` (Number): Suggested order for displaying this model's fit results.
        * `estimated parameters` (Number): Number of parameters estimated for this model.

* **`input`** (Object): Describes the input data used for the analysis.
    * `file name` (String): Path to the input alignment file.
    * `number of sequences` (Number): Count of sequences in the alignment.
    * `number of sites` (Number): Count of sites (columns) in the alignment.
    * `partition count` (Number): Number of data partitions (usually 0 or 1 if not partitioned).
    * `trees` (Object): Contains the input tree structure(s).
        * `<partition_id>` (String, e.g., "0"): Tree structure in Newick format string.

* **`runtime`** (String):  Version of `HyPhy` used to execute the analysis.

* **`test results`** (Object): Contains results of Likelihood Ratio Tests comparing pairs of nested models.
    * `<model_comparison>` (Object, e.g., "Proportional Partitioned:Unconstrained Test"): Represents a comparison between two models (Null:Alternative).
        * `Corrected P-value` (Number): Corrected (Holm-Bonferroni) p-value for the comparison.
        * `LRT` (Number): Likelihood Ratio Test statistic.
        * `Uncorrected P-value` (Number): Uncorrected p-value.

* **`tested`** (Object): Maps each branch name from the input tree to its classification in the analysis.
    * `<branch_name>` (String, e.g., "HLacaPus1"): Value is either "test" or "background".

* **`timers`** (Object): Records the time duration for specific parts of the analysis.
    * `<component_name>` (Object, e.g., "Overall", "Proportional"):
        * `order` (Number): Order index for the component.
        * `timer` (Number): Time duration in seconds.