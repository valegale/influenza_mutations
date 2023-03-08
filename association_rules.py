# Function that generate association rules.
# Modified version of the association_rules function 
# from the package mlxtend by Sebastian Raschka

from itertools import combinations
import numpy as np
import pandas as pd
from antigenic_sites import contain_antigenic_site, contain_transition_site

def generate_association_rules(df, metric="confidence", min_threshold=0.8, support_only=False, only_antigenic_sites = False, only_transition_sites = False):
    """ Original function at: http://rasbt.github.io/mlxtend/user_guide/frequent_patterns/association_rules/
	
	Modified here, as there are two additional parameters: 
	 - only_antigenic_sites (default = False)
	 - only_transition_sites default = False)
	 
	When only_antigenic_sites is True, only rules with at least one position in the 
	antigenic site are generated.
	
	Alternatively, with only_transition_sites = True, only rules are generated that include 
	the seven positions that have been shown to cause major antigenic changes 
	References: https://pubmed.ncbi.nlm.nih.gov/24264991/
	 
	Check script antigenic_sites.py for the specific positions.
    """
    if not df.shape[0]:
        raise ValueError(
            "The input DataFrame `df` containing " "the frequent itemsets is empty."
        )

    # check for mandatory columns
    if not all(col in df.columns for col in ["support", "itemsets"]):
        raise ValueError(
            "Dataframe needs to contain the\
                         columns 'support' and 'itemsets'"
        )

    def conviction_helper(sAC, sA, sC):
        confidence = sAC / sA
        conviction = np.empty(confidence.shape, dtype=float)
        if not len(conviction.shape):
            conviction = conviction[np.newaxis]
            confidence = confidence[np.newaxis]
            sAC = sAC[np.newaxis]
            sA = sA[np.newaxis]
            sC = sC[np.newaxis]
        conviction[:] = np.inf
        conviction[confidence < 1.0] = (1.0 - sC[confidence < 1.0]) / (
            1.0 - confidence[confidence < 1.0]
        )

        return conviction

    def zhangs_metric_helper(sAC, sA, sC):
        denominator = np.maximum(sAC * (1 - sA), sA * (sC - sAC))
        numerator = metric_dict["leverage"](sAC, sA, sC)

        with np.errstate(divide="ignore", invalid="ignore"):
            # ignoring the divide by 0 warning since it is addressed in the below np.where
            zhangs_metric = np.where(denominator == 0, 0, numerator / denominator)

        return zhangs_metric

    # metrics for association rules
    metric_dict = {
        "antecedent support": lambda _, sA, __: sA,
        "consequent support": lambda _, __, sC: sC,
        "support": lambda sAC, _, __: sAC,
        "confidence": lambda sAC, sA, _: sAC / sA,
        "lift": lambda sAC, sA, sC: metric_dict["confidence"](sAC, sA, sC) / sC,
        "leverage": lambda sAC, sA, sC: metric_dict["support"](sAC, sA, sC) - sA * sC,
        "conviction": lambda sAC, sA, sC: conviction_helper(sAC, sA, sC),
        "zhangs_metric": lambda sAC, sA, sC: zhangs_metric_helper(sAC, sA, sC),
    }

    columns_ordered = [
        "antecedent support",
        "consequent support",
        "support",
        "confidence",
        "lift",
        "leverage",
        "conviction",
        "zhangs_metric",
    ]

    # check for metric compliance
    if support_only:
        metric = "support"
    else:
        if metric not in metric_dict.keys():
            raise ValueError(
                "Metric must be 'confidence' or 'lift', got '{}'".format(metric)
            )

    # get dict of {frequent itemset} -> support
    keys = df["itemsets"].values
    values = df["support"].values
    frozenset_vect = np.vectorize(lambda x: frozenset(x))
    frequent_items_dict = dict(zip(frozenset_vect(keys), values))

    # prepare buckets to collect frequent rules
    rule_antecedents = []
    rule_consequents = []
    rule_supports = []
    i = 0
    # iterate over all frequent itemsets
    for k in frequent_items_dict.keys():
		
        if only_antigenic_sites and not contain_antigenic_site(k): #modified here
            i += 1
            pass
		
        elif only_transition_sites and not contain_transition_site(k):
            i += 1
            pass

        else:
			
	        sAC = frequent_items_dict[k]

	        # to find all possible combinations
	        for idx in range(len(k) - 1, 0, -1):
	            # of antecedent and consequent
	            for c in combinations(k, r=idx):
	                antecedent = frozenset(c)
	                consequent = k.difference(antecedent)
	
	                if support_only:
	                    # support doesn't need these,
	                    # hence, placeholders should suffice
	                    sA = None
	                    sC = None
	
	                else:
	                    try:
	                        sA = frequent_items_dict[antecedent]
	                        sC = frequent_items_dict[consequent]
	                    except KeyError as e:
	                        s = (
	                            str(e) + "You are likely getting this error"
	                            " because the DataFrame is missing "
	                            " antecedent and/or consequent "
	                            " information."
	                            " You can try using the "
	                            " `support_only=True` option"
	                        )
	                        raise KeyError(s)
	                    # check for the threshold
	
	                score = metric_dict[metric](sAC, sA, sC)
	                if score >= min_threshold:
	                    rule_antecedents.append(antecedent)
	                    rule_consequents.append(consequent)
	                    rule_supports.append([sAC, sA, sC])
        
            
			
    #print ("removed sequences: " + str(i))

    # check if frequent rule was generated
    if not rule_supports:
        return pd.DataFrame(columns=["antecedents", "consequents"] + columns_ordered)

    else:
        # generate metrics
        rule_supports = np.array(rule_supports).T.astype(float)
        df_res = pd.DataFrame(
            data=list(zip(rule_antecedents, rule_consequents)),
            columns=["antecedents", "consequents"],
        )

        if support_only:
            sAC = rule_supports[0]
            for m in columns_ordered:
                df_res[m] = np.nan
            df_res["support"] = sAC

        else:
            sAC = rule_supports[0]
            sA = rule_supports[1]
            sC = rule_supports[2]
            for m in columns_ordered:
                df_res[m] = metric_dict[m](sAC, sA, sC)

        return df_res