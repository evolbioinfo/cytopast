import logging
import os

import pandas as pd

from cytopast import read_tree
from pypastml.acr import acr
from pypastml.ml import MPPA, F81

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'examples', 'Ebola', 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Makona_1610_cds_ig.MCC.tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'Makona_1610_metadata_2016-06-23.csv')
HTML = os.path.join(DATA_DIR, 'tree_{}_{}_{}.html')
HTML_MAP = os.path.join(DATA_DIR, 'map_{}_{}_{}.html')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%H:%M:%S", filename=None)

    column = 'location'
    df = pd.read_csv(STATES_INPUT, index_col=1, header=0)[[column]]
    df[df == '?'] = None

    method = MPPA
    model = F81
    acr(read_tree(TREE_NWK), df, prediction_method=method, model=model)

