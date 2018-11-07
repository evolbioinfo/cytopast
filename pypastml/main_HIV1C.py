import logging
import os

import pandas as pd

from cytopast import read_tree
from pypastml import acr

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'examples', 'HIV1C', 'data')
TREE_NWK = os.path.join(DATA_DIR, 'best', 'pastml_fast_tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'metadata.tab')
HTML = os.path.join(DATA_DIR, 'tree.html')
MAP = os.path.join(DATA_DIR, 'map.html')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%H:%M:%S", filename=None)

    tree = read_tree(TREE_NWK)
    df = pd.read_table(STATES_INPUT, index_col=0, header=0)
    columns = ['RT:K103N', 'RT:M184V']

    acr(tree, df[columns], html=HTML, html_compressed=MAP)
