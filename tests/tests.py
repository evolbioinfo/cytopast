import glob
import logging
import os

import pandas as pd
from ete3 import Tree

from cytopast import pasml_annotations2cytoscape_annotation, annotate_tree_with_cyto_metadata, compress_tree
from cytopast.cytoscape_manager import save_as_cytoscape_html
from cytopast.pastml_analyser import get_enough_colours

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'JTTvsJC', 'rate0.1')
TREE_NWK = os.path.join(DATA_DIR, '1000.taxa.b.0.1.1.ID.tre')
STATES_TRUE = os.path.join(DATA_DIR, 'TrueScenario.txt')
STATES_OPT = os.path.join(DATA_DIR, 'OptScenario.txt')

HTML_TRUE = os.path.join(DATA_DIR, 'TrueScenario.html')
HTML_TRUE_FULL = os.path.join(DATA_DIR, 'TrueScenario_full.html')
HTML_OPT = os.path.join(DATA_DIR, 'OptScenario.html')
HTML_OPT_FULL = os.path.join(DATA_DIR, 'OptScenario_full.html')


def visualise(tree, res_data, html_compressed, html=None):
    new_res_data = res_data + '.pastml.csv'
    pd.read_table(res_data, sep=', ', header=0, index_col=0).astype(bool).to_csv(new_res_data, sep=',',
                                                                                 index=True)
    res_data = new_res_data

    columns = ['State']
    col2annotation_files = {columns[0]: res_data}
    res_annotations = os.path.join(DATA_DIR, 'combined_annotations_{}.tab'.format('_'.join(columns)))

    logging.info('Combining the data from different columns...')
    pasml_annotations2cytoscape_annotation(col2annotation_files, res_annotations, sep='\t')

    if len(col2annotation_files) == 1:
        df = pd.read_table(list(col2annotation_files.values())[0], sep=',', index_col=0, header=0)
        comb_df = pd.read_table(res_annotations, sep='\t', index_col=0, header=0)
        if set(df.columns) & set(comb_df.columns):
            df = df.join(comb_df, lsuffix='category', rsuffix='')
        else:
            df = df.join(comb_df)
        df.to_csv(res_annotations, sep='\t')

    tree, categories = annotate_tree_with_cyto_metadata(tree, res_annotations, sep='\t',
                                                        one_state=len(columns) == 1)

    name_column = columns[0]
    categories.remove(name_column)

    name2colour = {}
    colours = get_enough_colours(len(categories))
    for cat, col in zip(categories, colours):
        name2colour['{}_{}'.format(cat, True)] = col
    n2tooltip = lambda n, categories: ', '.join('{}'.format(cat) for cat in categories
                                                if hasattr(n, cat) and bool(getattr(n, cat, False)))

    if html:
        save_as_cytoscape_html(tree, html, categories=categories, graph_name='Tree', name2colour=name2colour,
                               name_feature=name_column, n2tooltip=n2tooltip)
    tree = compress_tree(tree, categories=(categories + [name_column]) if name_column else categories,
                         name_feature=name_column)
    save_as_cytoscape_html(tree, html_compressed, categories, graph_name='Summary map',
                           name2colour=name2colour, add_fake_nodes=False, n2tooltip=n2tooltip)


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    for work_dir in glob.glob(os.path.join(os.path.dirname(__file__), '*', '*')):
        tree_file = glob.glob(os.path.join(work_dir, '*.ID.tre'))[0]
        states_true = glob.glob(os.path.join(work_dir, 'TrueScenario.txt'))[0]
        states_opt = glob.glob(os.path.join(work_dir, 'OptScenario.txt'))[0]
        
        html_true = os.path.join(work_dir, 'TrueScenario_map.html')
        html_opt = os.path.join(work_dir, 'OptScenario_map.html')
        html_true_tree = os.path.join(work_dir, 'TrueScenario_tree.html')
        html_opt_tree = os.path.join(work_dir, 'OptScenario_tree.html')

        logging.info('Processing {}'.format(tree_file))

        new_tree = tree_file + '.pastml.tre'

        with open(tree_file, 'r') as tf:
            nwk = tf.readline()
            nwk = nwk.replace(';', ':0;')

        tree = Tree(nwk, 3)
        tree.name = 'ROOT'
        tree.write(outfile=new_tree, format=3, format_root_node=True)

        visualise(new_tree, states_true, html_compressed=html_true, html=html_true_tree)
        visualise(new_tree, states_opt, html_compressed=html_opt, html=html_opt_tree)
