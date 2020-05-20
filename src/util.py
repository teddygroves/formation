import numpy as np
import pandas as pd
from math import isnan

def get_stoichiometry(tecrdb):
    out = pd.DataFrame()
    formula_col = "Reaction formula in python dictionary"
    rxn_stoichiometries = tecrdb.groupby('reaction_id')[formula_col].first()
    for rid, stoichiometry in rxn_stoichiometries.iteritems():
        if type(stoichiometry) is not str:
            continue
        for cid, coef in eval(stoichiometry).items():
            out = out.append({
                "compound_id": cid, "reaction_id": rid, "coefficient": coef
            }, ignore_index=True)
    return out
