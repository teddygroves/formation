"""Write some tables based on the spreadsheet, namely

   - microspecies
   - TECRDB keq measurements
   - compounds

"""

import math
import numpy as np
import os
import pandas as pd
from typing import Dict
from util import get_stoichiometry

SHEET_PATH = "../data/raw/mmc2.xlsx"

OUTPUT_PATH_MICROSPECIES = "../data/processed/microspecies.csv"
OUTPUT_PATH_TECRDB = "../data/processed/tecrdb.csv"
OUTPUT_PATH_COMPOUNDS = "../data/processed/compounds.csv"


def process_dgf(dgf_in: pd.DataFrame) -> pd.DataFrame:
    out = dgf_in.copy()
    out["compound_id"], _ = zip(*out["Aqueous species id"].str.rsplit("_", 1))
    out["dgf_obs"] = out["ΔfG° (kJ/mol)"].copy()
    return out.dropna(subset=["dgf_obs"])


def process_tecrdb(tecrdb_in: pd.DataFrame) -> pd.DataFrame:
    out = tecrdb_in.copy()
    formula_col = "Reaction formula in python dictionary"
    formulas = out[formula_col].unique()
    rxn_ids = ["rxn_" + str(i) for i, _ in enumerate(formulas)]
    formula_to_rxn_id = dict(zip(formulas, rxn_ids))
    out['reaction_id'] = out[formula_col].map(formula_to_rxn_id.get)
    for col in ["pH", "pMg", "Ionic strength", "K'"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def get_protonation_states(pka_in: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame()
    charge_col = "Charges for different protonation states"
    pka_col = "pKas"
    nH_col = "Number of hydrogen atoms for different protonated states"
    for _, row in pka_in.iterrows():
        for i, (charge, nh) in enumerate(zip(
            eval(row[charge_col]), eval(row[nH_col])
        )):
            pka = eval(row[pka_col])[i-1] if i > 0 else np.nan
            out = out.append({
                "compound_id": row["Compound id"],
                "charge": int(charge),
                "nh": nh,
                "pka": pka
            }, ignore_index="True")
    out["compound_id"] = out["compound_id"].astype("string")
    return out

            
def get_magnesium_binding_states(pkmg_in: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame()
    id_col = "Aqueous species id"
    pkmg1_col = "pKMg"
    pkmg2_col = "pKMg2 (binding the second Mg)"
    data_col = "Data or predicted"
    for i, row in pkmg_in.iterrows():
        cid, charge = row[id_col].rsplit("_", 1)
        # if there are two or more recorded values, keep only the first
        if not charge.lstrip('-').isnumeric():
            if charge[-1] == 'a':
                charge = charge[:-1]
            else:
                continue
        pkmgs = [np.nan, row[pkmg1_col], row[pkmg2_col]]
        no_second_pkmg = math.isnan(row[pkmg2_col])
        nMgs = [0, 1] if no_second_pkmg else [0, 1, 2]
        for nMg, pkmg in zip(nMgs, pkmgs):
            predicted = row[data_col] == "predicted" if nMg > 0 else np.nan
            out = out.append({
                "nmg": nMg,
                "compound_id": cid,
                "charge": int(charge),
                "pkmg_is_predicted": predicted,
                "pkmg": pkmg
            }, ignore_index=True)
    out["compound_id"] = out["compound_id"].astype("string")
    return out


def get_microspecies(pka: pd.DataFrame, pkmg: pd.DataFrame) -> pd.DataFrame:
    H = get_protonation_states(pka)
    Mg = get_magnesium_binding_states(pkmg)
    compound_id_to_name = (
        pka.groupby("Compound id")["Compound name in TECRDB"].first().to_dict()
    )
    out = H.merge(Mg, on=["compound_id", "charge"], how='left')
    out["nmg"] = out["nmg"].fillna(0)
    out["microspecies_id"] = out.index.map(lambda i: "ms_" + str(i))
    out["name"] = out["compound_id"].map(compound_id_to_name.get)
    return out


def get_compounds(dgf_in, tecrdb, microspecies, pka):
    dgf = dgf_in.copy().rename(columns={"Charge": "charge"})
    stoichiometry = get_stoichiometry(tecrdb)
    cids = stoichiometry["compound_id"].unique()
    least_protonated_charge = (
        microspecies
        .loc[lambda df:
             (df["pka"].isnull()) & (df["compound_id"].isin(cids))
        ]
        .groupby("compound_id")["charge"].first()
    )
    compound_id_to_name = (
        pka.groupby("Compound id")["Compound name in TECRDB"].first().to_dict()
    )
    out = least_protonated_charge.loc[cids].reset_index()
    dgf["compound_id"], _ = zip(*dgf_in["Aqueous species id"].str.rsplit("_", 1))
    dgf["dgf_obs"] = dgf["ΔfG° (kJ/mol)"]
    out["name"] = out["compound_id"].map(compound_id_to_name.get)
    return out.merge(
        dgf[["compound_id", "charge", "dgf_obs"]],
        on=["compound_id", "charge"],
        how="left"
    )


def main():
    here = os.path.dirname(os.path.realpath(__file__))
    sheet_path = os.path.join(here, SHEET_PATH)
    dgf_in = pd.read_excel(sheet_path, sheet_name="Table S3. Compound thermo data")
    pka_in = pd.read_excel(sheet_path, sheet_name="Table S4. pKa data")
    tecrdb_in = pd.read_excel(sheet_path, sheet_name="Table S1. TECRDB Keqs")
    pkmg_in = pd.read_excel(sheet_path, sheet_name="Table S5. Mg binding data")

    microspecies = get_microspecies(pka_in, pkmg_in) 
    tecrdb = process_tecrdb(tecrdb_in)
    compounds = get_compounds(dgf_in, tecrdb, microspecies, pka_in)

    microspecies.to_csv(os.path.join(here, OUTPUT_PATH_MICROSPECIES))
    tecrdb.to_csv(os.path.join(here, OUTPUT_PATH_TECRDB))
    compounds.to_csv(os.path.join(here, OUTPUT_PATH_COMPOUNDS))


if __name__ == "__main__":
    main()
