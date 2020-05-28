.PHONY: clean clean-stan clean-output clean-analysis output analysis

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
ANALYSIS = data/analysis/img/pred.png \
	data/analysis/img/form.png \
	data/analysis/img/pka.png
STAN_FILES = src/stan/model.stan \
	src/stan/legendre.stan \
	src/stan/ordered_ragged_array.stan
SAMPLE_FILES = $(wildcard data/model_output/samples/*.csv)


all: analysis

$(RAW_DATA): src/fetch_data.py
	python src/fetch_data.py

$(PROCESSED_DATA): src/process_data.py $(RAW_DATA)
	python $<

$(MODEL_INPUT): src/prepare_stan_input.py $(PROCESSED_DATA)
	python $<

output: src/run_model.py $(MODEL_INPUT) $(STAN_FILES)
	python $<

analysis: src/analyse_output.py $(SAMPLE_FILES)
	python $<

clean-stan: AUXS=$(shell find src/stan -type f -not -name "*.stan")
clean-stan:
	$(RM) $(AUXS)

clean-output: CSVS=$(shell find $(SAMPLE_DIR) -type f -name '*.csv')
clean-output: TXTS=$(shell find $(SAMPLE_DIR) -type f -name '*.txt')
clean-output:
	$(RM) $(CSVS) $(TXTS)

clean-analysis:
	$(RM) $(ANALYSIS)

clean: clean-stan clean-output clean-analysis
	$(RM) $(PROCESSED_DATA) $(MODEL_INPUT)
