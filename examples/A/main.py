import os

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data2.txt')

if '__main__' == __name__:
    for model in ('JC',):
        pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                        html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}.html'.format(model)),
                        html=os.path.join(DATA_DIR, 'trees', 'tree_{}.html'.format(model)),
                        model=model, verbose=True, columns=['Problematic Sequence', 'Infection Country'])
