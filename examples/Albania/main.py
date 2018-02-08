import os

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.txt')

if '__main__' == __name__:
    for model in ('F81', 'JC'):
        pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                        html_compressed=os.path.join(DATA_DIR, 'map_{}.html'.format(model)),
                        html=os.path.join(DATA_DIR, 'tree_{}.html'.format(model)),
                        data_sep=',', model=model)
