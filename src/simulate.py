"""Generate some simulated data about the adk reaction"""

from legendre import get_dgr_prime
import numpy as np
import pandas as pd
from typing import Dict, List

STOICHIOMETRY = pd.DataFrame({
    "reaction_id": ["rxn_61", "rxn_61", "rxn_61"],
    "compound_id": ["CHB_16761", "CHB_15422", "CHB_16027"],
    "coef": [-2, 1, 1],
})
MICROSPECIES = pd.DataFrame({
    "microspecies_id": [
        'ms_1845', 'ms_1846', 'ms_1847', 'ms_1848', 'ms_1849', 'ms_1850',
        'ms_1851', 'ms_1852', 'ms_1853', 'ms_1854', 'ms_1948', 'ms_1949',
        'ms_1950', 'ms_1951', 'ms_1952', 'ms_1953', 'ms_1954', 'ms_1955',
        'ms_1956', 'ms_1998', 'ms_1999', 'ms_2000', 'ms_2001', 'ms_2002',
        'ms_2003', 'ms_2004', 'ms_2005', 'ms_2006', 'ms_2007', 'ms_2008',
        'ms_2009', 'ms_2010']
    ,
    "compound_id": [
        'CHB_16761', 'CHB_16761', 'CHB_16761', 'CHB_16761', 'CHB_16761',
        'CHB_16761', 'CHB_16761', 'CHB_16761', 'CHB_16761', 'CHB_16761',
        'CHB_16027', 'CHB_16027', 'CHB_16027', 'CHB_16027', 'CHB_16027',
        'CHB_16027', 'CHB_16027', 'CHB_16027', 'CHB_16027', 'CHB_15422',
        'CHB_15422', 'CHB_15422', 'CHB_15422', 'CHB_15422', 'CHB_15422',
        'CHB_15422', 'CHB_15422', 'CHB_15422', 'CHB_15422', 'CHB_15422',
        'CHB_15422', 'CHB_15422'
    ],
    "charge": [
        -4.0, -3.0, -3.0, -3.0, -2.0, -2.0, -1.0, -1.0, 0.0, 1.0, -3.0, -3.0,
        -2.0, -2.0, -1.0, -1.0, 0.0, 0.0, 1.0, -5.0, -4.0, -4.0, -4.0, -3.0,
        -3.0, -2.0, -2.0, -1.0, -1.0, 0.0, 0.0, 1.0
    ],
    "nh": [
        11.0, 12.0, 12.0, 12.0, 13.0, 13.0, 14.0, 14.0, 15.0, 16.0, 11.0,
        11.0, 12.0, 12.0, 13.0, 13.0, 14.0, 14.0, 15.0, 11.0, 12.0, 12.0,
        12.0, 13.0, 13.0, 14.0, 14.0, 15.0, 15.0, 16.0, 16.0, 17.0
    ],
    "nmg": [
        0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0,
        1.0, 0.0, 1.0, 0.0, 1.0, 0.0
    ],
})
TRUE_PARAMS = {
    "dgf": {
        "CHB_16761": -1893.88,
        "CHB_15422": -2749.05,
        "CHB_16027": -1009.73,
    },
    "pka": {
        "CHB_16761": [12.46, 7.05, 4.37, 2.05, 1.3],
        "CHB_15422": [12.6, 7.47, 4.67, 2.56, 1.54, 0.79],
        "CHB_16027": [12.46, 6.64, 3.37, 1.73],

    },
    "pkmg": {
        "CHB_16027": [3.49021041183173],
    }
}
DG_MEASUREMENT_CONDITIONS = pd.DataFrame({
    "reaction_id": [
        "rxn_61", "rxn_61", "rxn_61", "rxn_61", "rxn_61", "rxn_61", "rxn_61",
        "rxn_61", "rxn_61", "rxn_61", "rxn_61", "rxn_61", "rxn_61", "rxn_61",
        "rxn_61"
    ],
    "ph": [
        6.04, 6.67, 6.77, 6.83, 6.88, 7.04, 7.41, 7.45, 7.79, 7.93, 7.98, 8.13,
        8.16, 8.28, 8.30
    ],
    "pmg": [
        6.34, 6.84, 6.97, 3.51, 7.05, 3.54, 7.22, 7.24, 7.13, 7.15, 6.84, 6.84,
        6.78, 4.7, 5.49
    ],
    "I": [
        0.026, 0.037, 0.034, 0.249, 0.038, 0.245, 0.044, 0.048, 0.052, 0.061,
        0.071, 0.078, 0.072, 0.061, 0.069
    ],
    "t": [
        298.15, 298.15, 298.15, 311.15, 298.15, 311.15, 298.15, 298.15, 298.15,
        298.15, 310.25, 310.25, 304.15, 298.15, 298.15
    ],
})
PK_MEASUREMENT_POSITIONS = {
    "pka": {
        "CHB_16761": [0, 1, 2, 3, 4],
        "CHB_15422": [0, 1, 2, 3, 4, 5],
        "CHB_16027": [0, 1, 2, 3],
    },
    "pkmg": {
        "CHB_16027": [0]
    }
}
OUTDIR = "../data/simulated"
R = 8.314e-3
ERROR_SDS = {
    'pka': 0.2,
    'pkmg': 0.2,
    'dg': 0.2,
}


def simulate_pk_measurements(
    measurement_positions: Dict[str, List[int]],
    true_pks: Dict[str, List[float]],
    pk_type: str,
    error_sd: float
):
    out = pd.DataFrame()
    for cpd_id, cpd_true_pks in true_pks.items():
        cpd_measurement_positions = measurement_positions[cpd_id]
        pk_hat = [cpd_true_pks[pos] for pos in cpd_measurement_positions]
        pk_meas = np.random.normal(pk_hat, error_sd)
        cpd_measurements = pd.DataFrame({
            "compound_id": cpd_id,
            pk_type + "_order": cpd_measurement_positions,
            "true_" + pk_type: pk_hat,
            "measured_" + pk_type: pk_meas
        })
        out = pd.concat([out, cpd_measurements])
    return out

def simulate_dg_measurements(
    dg_measurement_conditions: pd.DataFrame,
    true_dgfs: Dict[str, float],
    true_pkas: Dict[str, float],
    true_pkmgs: Dict[str, float],
    microspecies: pd.DataFrame,
    stoichiometry: pd.DataFrame,
    error_sd: float
):
    dg_measurements = dg_measurement_conditions.copy()
    for i, row in dg_measurements.iterrows():
        rxn_stoichiometry = (
            stoichiometry
            .loc[lambda df: df["reaction_id"] == row["reaction_id"]]
            .set_index("compound_id")
            ["coef"]
            .to_dict()
        )
        dg_hat = get_dgr_prime(
            row["I"],
            row["t"],
            row["ph"],
            row["pmg"],
            R,
            microspecies,
            true_pkas,
            true_pkmgs,
            rxn_stoichiometry,
            true_dgfs
        )
        dg_meas = np.random.normal(dg_hat, error_sd)
        dg_measurements.loc[i, "true_dgr_prime"] = dg_hat
        dg_measurements.loc[i, "measured_dgr_prime"] = dg_meas
    return dg_measurements


def main():
    pka_measurements, pkmg_measurements = (
        simulate_pk_measurements(
            PK_MEASUREMENT_POSITIONS[p], TRUE_PARAMS[p], p, ERROR_SDS[p]
        ) for p in ["pka", "pkmg"]
    )
    dg_measurements = simulate_dg_measurements(
        DG_MEASUREMENT_CONDITIONS,
        TRUE_PARAMS["dgf"],
        TRUE_PARAMS["pka"],
        TRUE_PARAMS["pkmg"],
        MICROSPECIES,
        STOICHIOMETRY,
        ERROR_SDS["dg"]
    )
    print(dg_measurements)

if __name__ == "__main__":
    main()
