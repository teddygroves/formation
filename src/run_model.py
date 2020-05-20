from cmdstanpy import CmdStanModel
import os

HERE = os.path.dirname(os.path.realpath(__file__))
MODEL_PATH = os.path.join(HERE, "stan/model.stan")
INPUT_PATH = os.path.join(HERE, "../data/model_input/stan_model_input.json")
OUTPUT_DIR = os.path.join(HERE, "../data/model_output/samples")


SAMPLE_ARGS = {
    "chains": 2,
    "cores": 2,
    "iter_warmup": 500,
    "save_warmup": True,
    "iter_sampling": 500,
    "output_dir": OUTPUT_DIR,
    "max_treedepth": 10
}

def main():
    model = CmdStanModel(stan_file=MODEL_PATH)
    mcmc = model.sample(data=INPUT_PATH, **SAMPLE_ARGS)
    mcmc.diagnose()


if __name__ == "__main__":
    main()
