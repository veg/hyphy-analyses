RequireVersion ("2.5.25");


LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

namespace mss {
    neutral_reference = "NEUTRAL";
};

terms.model.mss.codon_classes = "classes for synonymous codons";


lfunction mss.LoadClasses (file) {

    SetDialogPrompt ("A TSV file with three columns (AA, Codon, Class) which is used to partition synonymous substitutions into groups");
    classes = io.ReadDelimitedFile (file, "\t", TRUE);
    headers = utility.Array1D(classes[^'terms.io.header']);
    io.CheckAssertion("`&headers`==3", "Expected a TSV file with 3 columns");
    codons_by_class = {};
    for (_record_; in; classes [^"terms.io.rows"]) {
        codons_by_class[_record_[1]] = _record_[2];
    }
    
    classes = utility.UniqueValues(codons_by_class);
    class_count = utility.Array1D(classes);
    io.CheckAssertion("`&class_count`>=2", "Expected at least 2 codon classes");

    choices = {class_count,2};
    for (i = 0; i < class_count; i += 1) {
        choices[i][0] = classes[i];
        choices[i][1] = "Codon class " + classes[i];
    }

    ^"mss.neutral_reference" = io.SelectAnOption  (choices, "Select the codon class which will serve as the neutral rate reference (relative rate = 1)");
    
    return codons_by_class;
}

// requires the definition of a global mss.codons_by_class variable for now


lfunction models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff(fromChar, toChar);
    
    if (None != _GenerateRate.diff) {
        _GenerateRate.p[model_type] = {};
        _GenerateRate.p[utility.getGlobalValue("terms.global")] = {};

        if (_GenerateRate.diff[utility.getGlobalValue("terms.diff.from")] > _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")]) {
            nuc_rate = "theta_" + _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")] + _GenerateRate.diff[utility.getGlobalValue("terms.diff.from")];
        } else {
            nuc_rate = "theta_" + _GenerateRate.diff[utility.getGlobalValue("terms.diff.from")] + _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")];
        }

        nuc_rate = parameters.ApplyNameSpace(nuc_rate, namespace);
        (_GenerateRate.p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible(_GenerateRate.diff[utility.getGlobalValue("terms.diff.from")], _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")])] = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (_GenerateRate.p[model_type])[omega_term] = aa_rate;
            } else {
                aa_rate = beta;
                (_GenerateRate.p[model_type])[beta_term] = aa_rate;
            }
            _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + aa_rate;
        } else {
            class_from = (^"mss.codons_by_class")[fromChar];
            class_to   = (^"mss.codons_by_class")[toChar];
            if (class_from == class_to) {
                if (class_from == ^"mss.neutral_reference") {
                   if (model_type == utility.getGlobalValue("terms.local")) {
                        (_GenerateRate.p[model_type])[alpha_term] = alpha;
                        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + alpha;
                    } else {
                        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
                    }
                } else {
                    codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from, namespace);
                    (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_within" + class_from] = codon_rate;
                    _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + codon_rate;
                }
            } else {
                if (class_from > class_to) {
                    codon_rate = class_to;
                    class_to = class_from;
                    class_from = codon_rate;
                }
                if (model_type == utility.getGlobalValue("terms.local")) {
                    codon_rate = alpha + "_" + class_from + "_" + class_to;
                } else {
                    codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from + "_" + class_to, namespace);
                }
                (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_between" + class_from + " and "  + class_to] = codon_rate;
                _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + codon_rate;
            }        
        }
   }
    

    return _GenerateRate.p;
}