import os

from cytopast.pastml_analyser import pastml_pipeline
from pastml import JOINT, MARGINAL_APPROXIMATION, MAX_POSTERIORI, MARGINAL, JC, F81

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.txt')

if '__main__' == __name__:
    for model in (F81, JC):
        for method in (JOINT, MAX_POSTERIORI, MARGINAL_APPROXIMATION, MARGINAL):
            pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                            html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}_{}.html'.format(model, method)),
                            html=os.path.join(DATA_DIR, 'trees', 'tree_{}_{}.html'.format(model, method)),
                            data_sep=',', model=model, verbose=True, prediction_method=method,
                            work_dir=os.path.join(DATA_DIR, 'pastml'))
