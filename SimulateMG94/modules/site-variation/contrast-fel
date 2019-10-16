LoadFunctionLibrary("libv3/all-terms.bf"); 
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/models/parameters.bf");
LoadFunctionLibrary("libv3/convenience/random.bf");

global simulator.site.alpha = 1;
simulator.site_betas = {};

lfunction simulator.prepare_site_distribution (model, site_count, tree, tree_info) {
 
    rates = (model[^"terms.parameters"])[^"terms.local"];
        
    classes = utility.UniqueValues (tree_info[^"terms.trees.model_map"]);
    class_count = utility.Array1D (classes);
    
    assert (class_count > 1, "Contrast-FEL requires at least two partitions in the tree");
    
    class_count = 0;
    
    beta_vars = {};
    
    utility.ForEach (classes, "_class_name_", '
        simulator.site_betas  [_class_name_] = "simulator.omega.class" + ^"`&class_count`";
        ^"`&class_count`" += 1;
        parameters.DeclareGlobalWithRanges (simulator.site_betas  [_class_name_], 1, 0, 10000);
        `&beta_vars` + simulator.site_betas  [_class_name_];
         
    '); 
    
    
    utility.ForEach (^tree, "_node_", '
         simulator.branch_class = ((^"`&tree_info`")[^"terms.trees.model_map"])[_node_];
         parameters.SetConstraint ("`tree`."+ _node_ + ".`rates[^'terms.parameters.nonsynonymous_rate']`",
                                    simulator.site_betas [simulator.branch_class] + "*`tree`."+ _node_ + ".`rates[^'terms.parameters.synonymous_rate']`", "");
         parameters.SetConstraint ("`tree`."+ _node_ + ".`rates[^'terms.parameters.synonymous_rate']`",
                                    "simulator.site.alpha*`tree`."+ _node_ + ".`rates[^'terms.parameters.synonymous_rate']`__", "");
    ');
    
    
    KeywordArgument   ("gamma-shape-srv", "The shape of the gamma parameter for synonymous rates (unit mean)", "1.");            
    alpha_srv = io.PromptUser ("The shape of the gamma parameter for SRV (unit mean)", 1.0, 1e-6, 20, FALSE);
    KeywordArgument   ("gamma-shape", "The shape of the gamma parameter for omegas (rate = beta = 1)", "0.5");            
    alpha = io.PromptUser ("The shape of the gamma parameter for omegas (rate = beta = 1)", 1.0, 1e-6, 20, FALSE);
    KeywordArgument   ("fraction-same", "The fraction of sites where omegas are the same", "0.8");            
    fraction_same = io.PromptUser ("The fraction of sites where omegas are the same ", 0.8, 0, 1, FALSE);

    site_profile = {site_count, 1 + class_count};
    
    for (i = 0; i < site_count; i+=1) {
          site_profile[i][0] = (random.gamma (alpha_srv,alpha_srv) * 100 $ 1)/100;  
          if (Random (0,1) < fraction_same) {
              shared_omega = (random.gamma_fast (alpha) * 100 $ 1)/100 ;
              for (j = 1; j <= class_count; j+=1) {
                site_profile[i][j] =  shared_omega;
              }
         
          } else {
              for (j = 1; j <= class_count; j+=1) {
                site_profile[i][j] = (random.gamma_fast (alpha) * 100 $ 1)/100;  
              }
          }
    }
    
    return {
        "variables" : beta_vars,
        "rates" : site_profile
    } ;
    
} 

lfunction simulator.apply_site_distribution (model, site, tree) {
    site = Eval (site);
    ^"simulator.site.alpha" = site[0];
    for (k = 0; k < Abs (^"simulator.site_betas"); k+=1) {
        parameters.SetValue (((^"simulator.site_profile")["variables"])[k], site[k+1]);
    }
}

lfunction simulator.set_site_omega (model, site_id, branch_info) {
/**
    no site-to-site rate variation
*/

    return ((^"simulator.site_profile")["rates"])[site_id][-1];

}