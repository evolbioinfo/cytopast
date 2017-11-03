import logging
import os

import cytopast as cy
from cytopast.cytoscape_manager import save_as_cytoscape_html

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Result_treeIDs.3036.taxa.14.states.tre')
STATES_OUTPUT = os.path.join(DATA_DIR, 'Result_states_probs.FULL.3036.taxa.14.states.txt')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)
    res_tree, res_annotations = TREE_NWK, STATES_OUTPUT
    res_html_comp = os.path.join(DATA_DIR, 'C_tree_compressed.html')
    res_html = os.path.join(DATA_DIR, 'C_tree.html')
    tree, categories = cy.annotate_tree_with_metadata(res_tree, res_annotations, sep=', ')
    save_as_cytoscape_html(tree, res_html, categories, graph_name='C')
    tree = cy.compress_tree(tree, categories)
    save_as_cytoscape_html(tree, res_html_comp, categories, graph_name='C (compressed)', add_fake_nodes=False)
