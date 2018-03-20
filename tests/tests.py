import logging
import os
from pastml import *

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'tree.nwk')
DATA_TAB = os.path.join(DATA_DIR, 'data.tab')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    for model in (F81, JC):
        for pred_method in (MARGINAL, MARGINAL_APPROXIMATION, MAX_POSTERIORI, JOINT):
            pastml_pipeline(data=DATA_TAB, tree=TREE_NWK, columns=['Location'],
                            html_compressed=os.path.join(DATA_DIR, 'map_{}_{}.html'.
                                                         format(model, pred_method)),
                            html=os.path.join(DATA_DIR, 'tree_{}_{}.html'.
                                              format(model, pred_method)),
                            model=model, prediction_method=pred_method, work_dir=os.path.join(DATA_DIR, 'pastml'),
                            verbose=True)
