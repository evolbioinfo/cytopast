import csv
import logging
import os

import cytopast as cy
from cytopast.cytoscape_manager import save_as_cytoscape_html

import pandas as pd

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'sequence-homo-align.tree')
STATES_IN = os.path.join(DATA_DIR, 'sequence-homo1.txt')
STATES_OUTPUT = os.path.join(DATA_DIR, 'sequence-homo_states.txt')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    df = pd.read_csv(STATES_IN, header=0, index_col=0)
    df['country'] = df['country'].apply(lambda s: s.replace(' ', '_'))
    df.to_csv(STATES_OUTPUT, columns=['country'], index=True, header=False, quoting=csv.QUOTE_ALL, quotechar='"')
    res_annotations = cy.apply_pastml(annotation_file=STATES_OUTPUT, tree_file=TREE_NWK)

    res_html_comp = os.path.join(DATA_DIR, 'tree_compressed.html')
    res_html = os.path.join(DATA_DIR, 'tree.html')
    tree, categories = cy.annotate_tree_with_metadata(TREE_NWK, res_annotations)
    save_as_cytoscape_html(tree, res_html, categories, graph_name='Zika')
    tree = cy.compress_tree(tree, categories, bin=False)
    save_as_cytoscape_html(tree, res_html_comp, categories, graph_name='Zika (compressed)', layout='cose')
