.PHONY: clean clean-stan clean-output output

RAW_DATA = data/raw/mmc2.xlsx
PROCESSED_DATA = data/processed/microspecies.csv \
	data/processed/tecrdb.csv \
	data/processed/compounds.csv
MODEL_INPUT = data/model_input/compounds.csv \
	data/model_input/pkas.csv \
	data/model_input/stan_codes.json \
	data/model_input/tecrdb.csv \
	data/model_input/stan_model_input.r \
	data/model_input/stan_model_input.json
SAMPLE_DIR = data/model_output/samples/
ANALYSIS = data/analysis/img/pred.png
STAN_FILES = src/stan/model.stan \
	src/stan/legendre.stan \
	src/stan/ordered_ragged_array.stan


all: output	

$(RAW_DATA): src/fetch_data.py
	python src/fetch_data.py

$(PROCESSED_DATA): src/process_data.py $(RAW_DATA)
	python $<

$(MODEL_INPUT): src/prepare_stan_input.py $(PROCESSED_DATA)
	python $<

output: src/run_model.py $(MODEL_INPUT) $(STAN_FILES)
	python $<

clean-stan: AUXS=$(shell find src/stan -type f -not -name "*.stan")
clean-stan:
	$(RM) $(AUXS)

clean-output: CSVS=$(shell find $(SAMPLE_DIR) -type f -name '*.csv')
clean-output: TXTS=$(shell find $(SAMPLE_DIR) -type f -name '*.txt')
clean-output:
	$(RM) $(CSVS) $(TXTS)

clean: clean-stan clean-output
	$(RM) $(PROCESSED_DATA) $(MODEL_INPUT)
