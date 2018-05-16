import os

from cytopast.pastml_analyser import pastml_pipeline
from pastml import JOINT, MARGINAL_APPROXIMATION, MAX_POSTERIORI, MARGINAL, JC, F81, DOWNPASS, ACCTRAN, DELTRAN

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')

if '__main__' == __name__:
    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                    html=os.path.join(DATA_DIR, 'trees', 'Albanian_tree_initial.html'),
                    html_compressed=os.path.join(DATA_DIR, 'maps', 'Albanian_map_initial.html'),
                    data_sep=',', verbose=True, copy_columns=['Country'], work_dir=os.path.join(DATA_DIR, 'pastml'))
    # for model in (F81, JC):
    for model in (F81,):
        for method in (JOINT, MAX_POSTERIORI, MARGINAL_APPROXIMATION, MARGINAL, DOWNPASS, ACCTRAN, DELTRAN):
            pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                            html_compressed=os.path.join(DATA_DIR, 'maps', 'Albanian_map_{}_{}.html'.format(model, method)),
                            html=os.path.join(DATA_DIR, 'trees', 'Albanian_tree_{}_{}.html'.format(model, method)),
                            data_sep=',', model=model, verbose=True, prediction_method=method,
                            work_dir=os.path.join(DATA_DIR, 'pastml'))
