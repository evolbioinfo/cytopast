import logging
import os

import pandas as pd

from cytopast import read_tree
from pypastml.acr import acr, COPY
from pypastml.ml import MPPA, MAP, JOINT, JC, F81, EFT
from pypastml.parsimony import DOWNPASS, ACCTRAN, DELTRAN

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'examples', 'Albania', 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')
HTML = os.path.join(DATA_DIR, 'tree_{}_{}.html')
HTML_MAP = os.path.join(DATA_DIR, 'map_{}_{}.html')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%H:%M:%S", filename=None)

    df = pd.read_csv(STATES_INPUT, index_col=0, header=0)
    column = 'Country'

    # method = MPPA
    # model = F81
    # acr(read_tree(TREE_NWK), df[[column]], html=HTML.format(method, model),
    #     html_compressed=HTML_MAP.format(method, model),
    #     prediction_method=method, model=model)

    acr(read_tree(TREE_NWK), df[[column]], prediction_method=COPY)

    for method in (MPPA, MAP, JOINT):
        for model in (JC, F81, EFT):
            acr(read_tree(TREE_NWK), df[[column]], prediction_method=method, model=model)

    for method in (DOWNPASS, ACCTRAN, DELTRAN):
            acr(read_tree(TREE_NWK), df[[column]], prediction_method=method, model=None)
