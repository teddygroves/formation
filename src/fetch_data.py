"""Fetch Du et al's spreadsheet from the ncbi website."""


import os
from urllib.request import urlretrieve

URL = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129446/bin/mmc2.xlsx"
OUTPATH = "../data/raw/mmc2.xlsx"

def main(url, outpath):
    urlretrieve(url, outpath)

if __name__ == "__main__":
    here = os.path.dirname(os.path.realpath(__file__))
    main(URL, os.path.join(here, OUTPATH))
