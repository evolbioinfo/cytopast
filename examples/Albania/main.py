import logging
import os

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.txt')

if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    for model in ('F81', 'JC'):
        pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                        html_compressed=os.path.join(DATA_DIR, 'albania_compressed_{}.html'.format(model)),
                        data_sep=',', work_dir=DATA_DIR, model=model)
