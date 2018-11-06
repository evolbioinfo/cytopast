import logging
import os

import pandas as pd

from cytopast import read_tree
from pypastml import reconstruct_ancestral_states


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'examples', 'Albania', 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')
HTML = os.path.join(DATA_DIR, 'tree.html')
MAP = os.path.join(DATA_DIR, 'map.html')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%H:%M:%S", filename=None)

    tree = read_tree(TREE_NWK)
    df = pd.read_csv(STATES_INPUT, index_col=0, header=0)
    column = 'Country'

    for _ in tree:
        if _.name in df.index:
            value = df.loc[_.name, column]
            if not pd.isna(value):
                _.add_feature(column, value)
    reconstruct_ancestral_states(tree, column, [_ for _ in df[column].unique() if not pd.isnull(_)],
                                 html=HTML, html_compressed=MAP)
