import csv
import logging
import os

import cytopast
from cytopast.cytoscape_manager import save_as_cytoscape_html

import pandas as pd
import re

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
# TREE_NWK = os.path.join(DATA_DIR, 'sequence-homo-align.tree')
TREE_NWK = os.path.join(DATA_DIR, 'booster.tree')
# STATES_IN = os.path.join(DATA_DIR, 'sequence-homo1.txt')
STATES_IN = os.path.join(DATA_DIR, 'Zika_annotations_Sep2017-homo-asian-outgroup3.csv')
# STATES_OUTPUT = os.path.join(DATA_DIR, 'sequence-homo_states.txt')
STATES_OUTPUT = os.path.join(DATA_DIR, 'zika_states.txt')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    df = pd.read_csv(STATES_IN, header=0, index_col=0)
    df.index = df.index.map(lambda s: re.sub(r'\.\d*', '', s))
    # df['Location'] = df['Location'].apply(lambda s: s.replace(' ', '_').replace('/', '_'))
    df.to_csv(STATES_OUTPUT, columns=['Location'], index=True, header=False, quoting=csv.QUOTE_ALL, quotechar='"')
    res_annotations = cytopast.apply_pastml(annotation_file=STATES_OUTPUT, tree_file=TREE_NWK)

    res_html_comp = os.path.join(DATA_DIR, 'tree_compressed.html')
    res_html = os.path.join(DATA_DIR, 'tree.html')
    tree, categories = cytopast.annotate_tree_with_metadata(TREE_NWK, res_annotations)
    save_as_cytoscape_html(tree, res_html, categories, graph_name='Zika')
    tree = cytopast.compress_tree(tree, categories, bin=False)
    save_as_cytoscape_html(tree, res_html_comp, categories, graph_name='Zika (compressed)', add_fake_nodes=False)
