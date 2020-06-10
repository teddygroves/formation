"""Functions for adjusting thermodynamic measurements given conditions."""

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from typing import Dict, List


def dh(temperature: float, ionic_strength: float) -> float:
    """Get the Debeye-Hukel number for a temperature and ionic strength.

    Source for the hardcoded numbers: 
    https://gitlab.com/equilibrator/equilibrator-cache/-/blob/develop/src/equilibrator_cache/thermodynamic_constants.py

    """
    a1 = 9.20483e-3
    a2 = 1.284668e-5
    a3 = 4.95199e-8
    alpha = a1 * temperature + a2 * temperature ** 2 + a3 * temperature ** 3
    sqrtI = np.sqrt(ionic_strength)
    return alpha * sqrtI / (1 + 1.6 * sqrtI)

def lt(
    dgf_std,
    pkas,
    pkmgs,
    nH,
    pH,
    pMg,
    I,
    charge,
    t,
    R
) -> float:
    """Get the condition-specific formation energy of a microspecies using a
    Legendre transform.

    source for the hardcoded numbers: 
    https://gitlab.com/equilibrator/equilibrator-cache/-/blob/develop/src/equilibrator_cache/thermodynamic_constants.py
    """
    log10 = np.log(10)
    RT = R * t
    dgfmg = -455.3
    dhfmg = -467.0
    dgfmg_prime = (t / 298.15) * dgfmg + (1.0 - t / 298.15) * dhfmg
    nMg = len(pkmgs) + 1;
    ddg_over_rt = (
        nH * RT * log10 * pH
        + nMg * (RT * log10 * pMg - dgfmg_prime)
        + dh(t, I) * (nH + 4 * nMg - charge ** 2)
    )
    if len(pkas) > 0:
      ddg_over_rt -= sum(log10 * np.array(pkas))
    if len(pkmgs) > 0:
      ddg_over_rt += len(pkmgs) * dgfmg/RT - sum(log10 * np.array(pkmgs))
    return dgf_std + ddg_over_rt * RT


def get_ms_pks(full_pk_list, pk_order):
    """Get a list of pks describing the binding reactions separating a microspecies
    from its least-bound state.

    :param pk_dict: dictionary mapping microspecies ids to their pk

    :param ms_id: id of the microspecies

    :param microspecies: dataframe of microspecies. Columns for "lowest_nh",
    "lowest_nmg", "nh" and "nmg". Indexed by microspecies ids.

    :param ncol: string indicating the number of bound ions. "nh" for pkas and
    "nmg" for pkmgs

    """
    n = microspecies.loc[ms_id, ncol]
    relevant_ms_ids = (
        microspecies
        .sort_values(ncol)
        .loc[lambda df: (df[ncol] <= n) & ~df["lowest_" + ncol]]
        .index
    )
    return [pk_dict[ms_id] for ms_id in relevant_ms_ids]


def get_dgr_prime(
    I: float,
    t: float,
    pH: float,
    pMg: float,
    R: float,
    microspecies: pd.DataFrame,
    pka: Dict[str, List[float]],
    pkmg: Dict[str, List[float]],
    reaction_stoichiometry: Dict[str, float],
    dgf: Dict[str, float]
) -> float:
    """Get the Gibbs energy change of a reaction in given conditions."""
    out = 0;
    RT = R * t;
    for compound_id, coef in reaction_stoichiometry.items():
        cpd_microspecies = (
            microspecies.loc[lambda df: df["compound_id"] == compound_id].copy()
        )
        ms_dgf_primes = []
        full_pka_list = pka[compound_id] if compound_id in pka.keys() else []
        full_pkmg_list = pkmg[compound_id] if compound_id in pkmg.keys() else []
        cpd_microspecies["pka_order"], cpd_microspecies["pkmg_order"] = (
            cpd_microspecies[n].rank(method="dense").subtract(1)
            for n in ["nh", "nmg"]
        )
        for ms_id, ms in cpd_microspecies.iterrows():
            ms_pkas = full_pka_list[:int(ms["pka_order"])]
            ms_pkmgs = full_pkmg_list[:int(ms["pkmg_order"])]
            ms_dgf_prime = lt(
                dgf[compound_id], ms_pkas, ms_pkmgs, ms["nh"],
                pH, pMg, I, ms["charge"], t, R
            )
            ms_dgf_primes.append(ms_dgf_prime)
        dgf_prime = -RT * logsumexp(-1 * np.array(ms_dgf_primes) / RT);
        out += coef * dgf_prime
    return out
