"""Make Stan input file(s) based on the Du et al data."""

import numpy as np
import os
import pandas as pd
from typing import Any, List, Dict, Union, Iterable
from cmdstanpy.utils import jsondump, rdump
from util import get_stoichiometry

INPUT_PATH_MICROSPECIES = "../data/processed/microspecies.csv"
INPUT_PATH_TECRDB = "../data/processed/tecrdb.csv"
INPUT_PATH_COMPOUNDS = "../data/processed/compounds.csv"
OUTPUT_PATH_JSON = "../data/model_input/stan_model_input.json"
OUTPUT_PATH_R = "../data/model_input/stan_model_input.r"
OUTPUT_PATH_STAN_CODES = "../data/model_input/stan_codes.json"
OUTPUT_PATH_TECRDB = "../data/model_input/tecrdb.csv"
OUTPUT_PATH_COMPOUNDS = "../data/model_input/compounds.csv"
OUTPUT_PATH_PKAS = "../data/model_input/pkas.csv"


def filter_tecrdb(t: pd.DataFrame) -> pd.Series:
    bad_evals = ["D", "E"]
    return (
        (t["exclude"] == "no")
        & (t["Ionic strength"].notnull())
        & (t["pH"].notnull())
        & (t["pMg"].notnull())
        & (t["T(K)"].notnull())
        & (~t["Evaluation"].isin(bad_evals))
    )
    

def stan_codify(i: Iterable[str]) -> Dict[str, int]:
    uniques = sorted(list(set(i)))
    return dict(zip(uniques, range(1, len(uniques) + 1)))


def stan_encode(i: Iterable[str], codes: Dict[str, int]) -> Iterable[int]:
    return [codes[s] for s in i]


def get_group_first_ix(df: pd.DataFrame, groupcol: str) -> List[int]:
    return df.groupby(groupcol).apply(lambda d: d.index[0] + 1).tolist()


def get_stan_codes(stoichiometry: pd.DataFrame) -> Dict[str, Dict[str, int]]:
    return {
        "reaction_id": stan_codify(stoichiometry["reaction_id"]),
        "compound_id": stan_codify(stoichiometry["compound_id"])
    }


def get_stan_input(
    tecrdb: pd.DataFrame,
    stoichiometry: pd.DataFrame,
    microspecies: pd.DataFrame,
    compounds: pd.DataFrame,
    pkas: pd.DataFrame,
    pkmgs: pd.DataFrame,
    stan_codes: Dict[str, Dict[str, int]], 

) -> Dict[str, Any]:
    cpd_codes, rxn_codes = stan_codes["compound_id"], stan_codes["reaction_id"]
    stoichiometry["reaction_id_stan"] = stan_encode(
        stoichiometry["reaction_id"], rxn_codes
    )
    tecrdb["reaction_id_stan"] = stan_encode(tecrdb["reaction_id"], rxn_codes)
    stoichiometry["compound_id_stan"] = stan_encode(
        stoichiometry["compound_id"], cpd_codes
    )
    microspecies["compound_id_stan"] = stan_encode(
        microspecies["compound_id"], cpd_codes
    )
    compounds["compound_id_stan"] = stan_encode(compounds["compound_id"], cpd_codes)
    pkas["compound_id_stan"] = stan_encode(pkas["compound_id"], cpd_codes)
    pkmgs["compound_id_stan"] = stan_encode(pkmgs["compound_id"], cpd_codes)
    n_pka, n_pkmg = (
        microspecies.groupby("compound_id_stan")[n]
        .transform(lambda s: s.rank(method='dense') - 1)
        .astype(int)
        for n in ["nh", "nmg"]
    )
    n_pka_cpd, n_pkmg_cpd = (
        pks.groupby("compound_id_stan").size().sort_index()
        for pks in [pkas, pkmgs] 
    )
    compound_to_first_pka, compound_to_first_pkmg = (
        dict(zip(
            p["compound_id_stan"].unique(),
            get_group_first_ix(p, 'compound_id_stan')
        ))
        for p in [pkas, pkmgs]
    )
    first_pka, first_pkmg = (
        microspecies["compound_id_stan"]
        .map(d)
        .fillna(0)
        .astype(int)
        .tolist()
        for d in [compound_to_first_pka, compound_to_first_pkmg]
    )
    compounds = compounds.sort_values("compound_id_stan")
    compounds["prior_loc_dgf"] = np.where(
        compounds["dgf_obs"].notnull(),
        compounds["dgf_obs"],
        0
    )
    compounds["prior_regime_dgf"] = np.where(compounds["dgf_obs"].notnull(), 1, 2)

    stan_input = {
        "N_measurement_kpr": len(tecrdb),
        "N_measurement_pka": len(pkas),
        "N_measurement_pkmg": len(pkmgs),
        "N_microspecies": len(microspecies),
        "N_pka": len(pkas),
        "N_pkmg": len(pkmgs),
        "N_compound": stoichiometry["compound_id"].nunique(),
        "N_stoic": len(stoichiometry),
        "N_reaction":  stoichiometry["reaction_id"].nunique(),
        "N_Hable_compound": pkas["compound_id"].nunique(),
        "N_Mgable_compound": pkmgs["compound_id"].nunique(),
        "n_cpd": stoichiometry.groupby("reaction_id_stan").size().tolist(),
        "charge": microspecies["charge"].astype(float).tolist(),
        "nH": microspecies["nh"].astype(int).tolist(),
        "nMg": microspecies["nmg"].astype(int).tolist(),
        "stoic_cpd": stoichiometry["compound_id_stan"].tolist(),
        "ms_cpd": microspecies["compound_id_stan"].tolist(),
        "stoic_coef": stoichiometry["coefficient"].tolist(),
        "fst_stoic": get_group_first_ix(stoichiometry, "reaction_id_stan"),
        "n_pka": n_pka.tolist(),
        "n_pka_cpd": n_pka_cpd.tolist(),
        "fst_pka": first_pka,
        "n_pkmg": n_pkmg.tolist(),
        "n_pkmg_cpd": n_pkmg_cpd.tolist(),
        "fst_pkmg": first_pkmg,
        "n_ms": microspecies.groupby("compound_id_stan").size().tolist(),
        "fst_ms": get_group_first_ix(microspecies.reset_index(), "compound_id_stan"),
        "pka_obs_ix": list(pkas.reset_index().index + 1),
        "pka_obs": pkas["pka"].tolist(),
        "pkmg_obs_ix": list(pkmgs.reset_index().index + 1),
        "pkmg_obs": pkmgs["pkmg"].tolist(),
        "kpr_obs": tecrdb["K'"].tolist(),
        "rxn": tecrdb["reaction_id_stan"].tolist(),
        "temperature": tecrdb["T(K)"].tolist(),
        "I": tecrdb["Ionic strength"].tolist(),
        "pH": tecrdb["pH"].tolist(),
        "pMg": tecrdb["pMg"].tolist(),
        "prior_loc_dgf": compounds["prior_loc_dgf"].tolist(),
        "prior_regime_dgf": compounds["prior_regime_dgf"].tolist(),
    }
    return stan_input


def main():
    here = os.path.dirname(os.path.realpath(__file__))
    input_path_tecrdb = os.path.join(here, INPUT_PATH_TECRDB)
    input_path_compounds = os.path.join(here, INPUT_PATH_COMPOUNDS)
    input_path_microspecies = os.path.join(here, INPUT_PATH_MICROSPECIES)
    tecrdb = pd.read_csv(input_path_tecrdb, index_col=0).loc[filter_tecrdb]
    stoichiometry = get_stoichiometry(tecrdb)
    cids = stoichiometry['compound_id'].unique()
    rids = stoichiometry['reaction_id'].unique()
    compounds = (
        pd.read_csv(input_path_compounds, index_col=0)
        .loc[lambda df: df['compound_id'].isin(cids)]
    )
    microspecies = (
        pd.read_csv(input_path_microspecies, index_col=0)
        .loc[lambda df: df['compound_id'].isin(cids)]
    )
    pkas, pkmgs = (
        microspecies
        .groupby(['compound_id', n])[[pk, "name"]].first()
        .dropna()
        .reset_index()
        for n, pk in [('nh', 'pka'), ('nmg', 'pkmg')]
    )
    stan_codes = get_stan_codes(stoichiometry)
    stan_input = get_stan_input(
        tecrdb, stoichiometry, microspecies, compounds, pkas, pkmgs, stan_codes
    )
    pkas.to_csv(os.path.join(here, OUTPUT_PATH_PKAS))
    compounds.to_csv(os.path.join(here, OUTPUT_PATH_COMPOUNDS))
    tecrdb.to_csv(os.path.join(here, OUTPUT_PATH_TECRDB))
    jsondump(os.path.join(here, OUTPUT_PATH_STAN_CODES), stan_codes)
    jsondump(os.path.join(here, OUTPUT_PATH_JSON), stan_input)
    rdump(os.path.join(here, OUTPUT_PATH_R), stan_input)


if __name__ == "__main__":
    main()
