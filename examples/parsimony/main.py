import os

from cytopast.pastml_analyser import pastml_pipeline
from pastml import DOWNPASS, ACCTRAN, DELTRAN

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data.tab')

if '__main__' == __name__:
    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                    html=os.path.join(DATA_DIR, 'trees', 'tree_initial.html'),
                    verbose=True, copy_columns=['Symbol'], work_dir=os.path.join(DATA_DIR, 'pastml'))
    for method in (DOWNPASS, ACCTRAN, DELTRAN):
            pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                            html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}.html'.format(method)),
                            html=os.path.join(DATA_DIR, 'trees', 'tree_{}.html'.format(method)),
                            verbose=True, prediction_method=method,
                            work_dir=os.path.join(DATA_DIR, 'pastml'))
