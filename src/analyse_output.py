import arviz as az
import json
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import numpy as np
import os
import pandas as pd

RELPATHS = {
    "samples_folder": "../data/model_output/samples",
    "stan_codes": "../data/model_input/stan_codes.json",
    "tecrdb": "../data/model_input/tecrdb.csv",
    "pkas": "../data/model_input/pkas.csv",
    "compounds": "../data/model_input/compounds.csv",
    "img_folder": "../data/analysis/img",
}
R = 8.314e-3
HERE = os.path.dirname(os.path.realpath(__file__))



def plot_predictions(d):
    y = d["K'"].rank(method="first")
    obs = d["K'"]
    lo, hi = d[0.1], d[0.9]
    f, ax = plt.subplots(figsize=[12.5, 10])
    f.suptitle("K' Observations vs simulated measurements")
    ax.get_yaxis().set_visible(False)
    ax.set(xlabel="K'")
    ax.scatter(obs, y, label="Observed value")
    ax.hlines(y, hi, lo, color='tab:orange', zorder=0, label="10%-90% credible region")
    for i, row in d.iterrows():
        ax.text(hi.loc[i], y.loc[i]-0.2, row["Reaction"], fontsize=8)
    img = os.path.join(HERE, RELPATHS["img_folder"])
    plt.semilogx()
    ax.legend(frameon=False, loc="upper left")
    plt.savefig(os.path.join(img, "pred.png"))


def plot_formation_energy(d):
    med = d[0.5]
    obs = d["dgf_obs"].dropna()
    y = med.rank()
    lo, hi = d[0.1], d[0.9]
    f, ax = plt.subplots(figsize=[12.5, 7])
    f.suptitle("compound formation energy estimates")
    ax.get_yaxis().set_visible(False)
    ax.set(
        xlabel="formation energy (kj/mol)",
        ylabel="compound"
    )
    ax.scatter(obs, y.loc[obs.index], label="observed value")
    ax.hlines(y, hi, lo, color='tab:orange', zorder=0, label="10%-90% credible region")
    for i, row in d.iterrows():
        ax.text(hi.loc[i], y.loc[i]-0.2, row["name"], fontsize=8)
    img = os.path.join(HERE, RELPATHS["img_folder"])
    ax.legend(frameon=False, loc="upper left")
    plt.savefig(os.path.join(img, "form.png"))


def plot_pka(d):
    med = d[0.5]
    obs = d["pka"].dropna()
    y = med.rank()
    lo, hi = d[0.1], d[0.9]
    labs = d["name"].str.cat(d["nh"].astype(int).astype(str), sep=", nH ")
    f, ax = plt.subplots(figsize=[12.5, 15])
    f.suptitle("pKa estimates")
    ax.get_yaxis().set_visible(False)
    ax.set(xlabel="pKa")
    ax.scatter(obs, y.loc[obs.index], label="observed value")
    ax.hlines(y, hi, lo, color='tab:orange', zorder=0, label="10%-90% credible region")
    for i in d.index:
        ax.text(hi.loc[i], y.loc[i]-0.2, labs.loc[i], fontsize=8)
    img = os.path.join(HERE, RELPATHS["img_folder"])
    ax.legend(frameon=False, loc="upper left")
    plt.savefig(os.path.join(img, "pka.png"))


def get_infd(csvs, stan_codes, tecrdb):
    return az.from_cmdstan(
        csvs,
        coords={
            "compound": list(stan_codes["compound_id"].keys()),
            "measurement": list(tecrdb.index)
        },
        dims={
            "dgf": ["compound"],
            "dgr_hat": ["measurement"],
            "kpr_rep": ["measurement"]
        }
    )

def main():
    plt.style.use(os.path.join(HERE, 'sparse.mplstyle'))
    sample_folder_path = os.path.join(HERE, RELPATHS["samples_folder"])
    stan_codes = json.load(open(os.path.join(HERE, RELPATHS["stan_codes"]), "r"))
    tecrdb = pd.read_csv(os.path.join(HERE, RELPATHS["tecrdb"]), index_col=0)
    pkas = pd.read_csv(os.path.join(HERE, RELPATHS["pkas"]), index_col=0)
    compounds = pd.read_csv(os.path.join(HERE, RELPATHS["compounds"]), index_col=0)
    csvs = [
        os.path.join(sample_folder_path, f)
        for f in os.listdir(sample_folder_path)
        if f.endswith("csv")
    ]
    infd = get_infd(csvs, stan_codes, tecrdb)
    samples = {
        p: infd.posterior[p].to_series().unstack() 
        for p in ["dgf", "dgr_hat", "pka", "pkmg", "kpr_rep"]
    }
    tecrdb = tecrdb.join(samples["kpr_rep"].quantile([0.1, 0.5, 0.9]).T)
    compounds = compounds.join(samples["dgf"].quantile([0.1, 0.5, 0.9]).T, on="compound_id")
    pkas = pkas.join(samples["pka"].quantile([0.1, 0.5, 0.9]).T)
    plot_predictions(tecrdb)
    plot_formation_energy(compounds)
    plot_pka(pkas)
    plt.close("all")



if __name__ == "__main__":
    main()
