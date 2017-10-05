import logging
import os

import cytopast as cy
from cytopast.cytoscape_manager import save_as_cytoscape_html

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.txt')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)
    res_tree = TREE_NWK
    res_annotations = cy.apply_pastml(annotation_file=STATES_INPUT, tree_file=res_tree)
    res_html = os.path.join(DATA_DIR, 'albania.html')
    res_html_comp = os.path.join(DATA_DIR, 'albania_compressed.html')

    tree, categories = cy.annotate_tree_with_metadata(res_tree, res_annotations)
    save_as_cytoscape_html(tree, res_html, categories, graph_name='Albania')
    tree = cy.compress_tree(tree, categories)
    save_as_cytoscape_html(tree, res_html_comp, categories, graph_name='Albania (compressed)',
                           layout='cose')
