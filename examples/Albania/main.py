import logging
import os

import cytopast as cy
from cytopast.cytoscape_manager import save_as_cytoscape_html
import pandas as pd
import numpy as np

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.txt')
STATES_INPUT_PERCENT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.randomized_{}_percent.txt')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    res_tree = TREE_NWK
    res_annotations = cy.apply_pastml(annotation_file=STATES_INPUT, tree_file=res_tree)
    res_html = os.path.join(DATA_DIR, 'albania.html')
    res_html_comp = os.path.join(DATA_DIR, 'albania_compressed.html')

    tree, categories = cy.annotate_tree_with_metadata(res_tree, res_annotations)
    save_as_cytoscape_html(tree, res_html, categories, graph_name='Albania')
    tree = cy.compress_tree(tree)
    save_as_cytoscape_html(tree, res_html_comp, categories, graph_name='Albania (compressed)', add_fake_nodes=False)

    for percent in (5, 10, 25):
        df = pd.read_csv(STATES_INPUT, index_col=0, header=None)
        random_index = df.sample(frac=percent / 100).index
        new_values = np.random.choice(df[1].unique(), size=len(random_index))
        df.loc[random_index, 1] = new_values
        randomized_state_file = STATES_INPUT_PERCENT.format(percent)
        df.to_csv(randomized_state_file, header=False, index=True)

        res_annotations = cy.apply_pastml(annotation_file=randomized_state_file, tree_file=res_tree)
        res_html = os.path.join(DATA_DIR, 'albania_randomized_{}_percent.html'.format(percent))
        res_html_comp = os.path.join(DATA_DIR, 'albania_compressed_randomized_{}_percent.html'.format(percent))

        tree, categories = cy.annotate_tree_with_metadata(res_tree, res_annotations)
        save_as_cytoscape_html(tree, res_html, categories, graph_name='Albania {} precent'.format(percent))
        tree = cy.compress_tree(tree)
        save_as_cytoscape_html(tree, res_html_comp, categories,
                               graph_name='Albania (compressed) {} precent'.format(percent), add_fake_nodes=False)

