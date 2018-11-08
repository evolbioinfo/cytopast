import logging
import os

import pandas as pd

from cytopast import read_tree
from pypastml import acr, JOINT, MPPA, MAP

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'examples', 'Albania', 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')
HTML = os.path.join(DATA_DIR, 'tree_{}.html')
HTML_MAP = os.path.join(DATA_DIR, 'map_{}.html')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%H:%M:%S", filename=None)

    df = pd.read_csv(STATES_INPUT, index_col=0, header=0)
    column = 'Country'

    for method in (MPPA, MAP, JOINT):
        acr(read_tree(TREE_NWK), df[[column]], html=HTML.format(method), html_compressed=HTML_MAP.format(method),
            prediction_method=method)
