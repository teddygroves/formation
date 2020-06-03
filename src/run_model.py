from cmdstanpy import CmdStanModel
import os

RELATIVE_PATHS = {
    "model": "stan/model.stan",
    "model_input": "../data/model_input/stan_model_input.json",
    "model_output_dir": "../data/model_output/samples"
}

SAMPLE_ARGS = {
    "chains": 2,
    "cores": 2,
    "iter_warmup": 400,
    "save_warmup": True,
    "iter_sampling": 400,
    "max_treedepth": 10
}

def standardise_csv_filename(path_in: str, output_dir: str) -> str:
    directory, filename = os.path.split(path_in)
    _, chain_num, _ = filename.rsplit('-', 2)
    new_filename = '-'.join(["samples", str(chain_num)]) + ".csv"
    return os.path.join(output_dir, new_filename)


def main():
    here = os.path.dirname(os.path.realpath(__file__))
    model_path = os.path.join(here, RELATIVE_PATHS["model"])
    input_path = os.path.join(here, RELATIVE_PATHS["model_input"])
    output_dir = os.path.join(here, RELATIVE_PATHS["model_output_dir"])
    model = CmdStanModel(stan_file=model_path)
    mcmc = model.sample(data=input_path, **SAMPLE_ARGS)
    for f in mcmc.runset.csv_files:
        os.replace(f, standardise_csv_filename(f, output_dir))



if __name__ == "__main__":
    main()
