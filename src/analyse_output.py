import arviz as az
import json
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
plt.style.use('sparse.mplstyle')
import numpy as np
import os
import pandas as pd

SAMPLE_FOLDER = "samples"
R = 8.314e-3


def plot_predictions(d):
    y = d["dgr_obs"].rank()
    obs = d["dgr_obs"]
    hi, low = d[0.1], d[0.9]
    f, ax = plt.subplots(figsize=[12.5, 7])
    ax.scatter(obs, y, c=d["pH"], cmap="viridis")
    ax.hlines(y, hi, low, color='tab:orange', zorder=0)
    for i, row in d.iterrows():
        ax.text(obs.loc[i], y.loc[i], row["Reaction"], fontsize=8)
    plt.savefig("img/pred.png")


def get_infd(csvs, stan_codes, tecrdb):
    return az.from_cmdstan(
        csvs,
        coords={
            "compound": list(stan_codes["compound_id"].keys()),
            "measurement": list(tecrdb.index)
        },
        dims={"dgf": ["compound"], "dgr_hat": ["measurement"]}
    )

def main():
    stan_codes = json.load(open("artefacts/stan_codes.json", "r"))
    tecrdb = pd.read_csv("artefacts/tecrdb.csv", index_col=0)
    pkas = pd.read_csv("artefacts/pkas.csv", index_col=0)
    compounds = pd.read_csv("artefacts/compounds.csv", index_col=0)
    csvs = [
        os.path.join(SAMPLE_FOLDER, f)
        for f in os.listdir(SAMPLE_FOLDER)
        if f.endswith("csv")
    ]
    infd = get_infd(csvs, stan_codes, tecrdb)
    samples = {
        p: infd.posterior[p].to_series().unstack() 
        for p in ["dgf", "dgr_hat", "pka", "pkmg"]
    }
    tecrdb["dgr_obs"] = -R * tecrdb["T(K)"] * np.log(tecrdb["K'"])
    tecrdb = tecrdb.join(samples["dgr_hat"].quantile([0.1, 0.5, 0.9]).T)
    compounds = compounds.join(samples["dgf"].quantile([0.1, 0.5, 0.9]).T, on="compound_id")
    pkas = pkas.join(samples["pka"].quantile([0.1, 0.5, 0.9]).T)
    plot_predictions(tecrdb)
    plt.close("all")



if __name__ == "__main__":
    main()
