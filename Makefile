.PHONY: clean clean-stan clean-samples clean-analysis samples analysis paper

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
ANALYSIS = analysis/img/pred.png \
	analysis/img/form.png \
	analysis/img/pka.png
STAN_FILES = src/stan/model.stan \
	src/stan/legendre.stan \
	src/stan/ordered_ragged_array.stan
SAMPLE_FILES = data/model_output/samples/samples-1.csv \
	data/model_output/samples/samples-2.csv
PAPER = paper.pdf
BIBLIOGRAPHY = bibliography.bib
PANDOCFLAGS =                          \
  --from=org                           \
  --highlight-style=pygments           \
  --pdf-engine=xelatex                 \
  --bibliography=$(BIBLIOGRAPHY)

all: paper

samples: $(SAMPLE_FILES)

analysis: $(ANALYSIS)

paper: $(PAPER)

$(RAW_DATA): src/fetch_data.py
	python src/fetch_data.py

$(PROCESSED_DATA): src/process_data.py $(RAW_DATA)
	python $<

$(MODEL_INPUT): src/prepare_stan_input.py $(PROCESSED_DATA)
	python $<

$(SAMPLE_FILES): src/run_model.py $(MODEL_INPUT) $(STAN_FILES)
	python $<

$(ANALYSIS): src/analyse_output.py $(SAMPLE_FILES)
	python $<

$(PAPER): paper.org $(ANALYSIS) $(BIBLIOGRAPHY)
	pandoc $< -o $@ $(PANDOCFLAGS)

clean-stan: AUXS=$(shell find src/stan -type f -not -name "*.stan")
clean-stan:
	$(RM) $(AUXS)

clean-samples:
	$(RM) $(SAMPLE_FILES)

clean-paper:
	$(RM) $(PAPER)

clean-analysis:
	$(RM) $(ANALYSIS)

clean: clean-stan clean-samples clean-analysis
	$(RM) $(PROCESSED_DATA) $(MODEL_INPUT)
