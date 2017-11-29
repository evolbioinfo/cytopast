import logging
import os
import random
import shutil
import tempfile
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

from cytopast import apply_pastml, compress_tree, STATES_TAB_PASTML_OUTPUT, read_tree, \
    pasml_annotations2cytoscape_annotation, annotate_tree_with_cyto_metadata
from cytopast.cytoscape_manager import save_as_cytoscape_html

COLOURS = ['#a6dba0', '#a50026', '#fdae61', '#313695', '#d73027',
           '#fee090', '#4575b4', '#f46d43', '#abd9e9', '#ffffbf',
           '#74add1', '#e0f3f8']
WHITE = '#ffffff'


def random_hex_color():
    r = lambda: random.randint(0, 255)
    return '#%02X%02X%02X' % (r(), r(), r())


def work(args):
    tree, df, work_dir, column = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)
    rep_dir = os.path.join(work_dir, category)
    unique_states = df.unique()
    n_tips = len(read_tree(tree).get_leaves())
    res_file = os.path.join(rep_dir, STATES_TAB_PASTML_OUTPUT).format(tips=n_tips,
                                                                      states=len(unique_states))
    if os.path.exists(res_file):
        return category, res_file
    os.makedirs(rep_dir, exist_ok=True)
    state_file = os.path.join(rep_dir, 'state_{}.csv'.format(category))

    # For binary states make sure that 0 or false is not mistaken by missing data
    if len(unique_states) == 2 and np.any(pd.isnull(unique_states)):
        state = unique_states[~pd.isnull(unique_states)][0]
        other_state = not state if isinstance(state, bool) \
            else (0 if isinstance(state, (int, float, complex)) and state != 0
                  else 'other' if state != 'other' else 'unknown')
        df.replace(np.nan, other_state, inplace=True)

    df.to_csv(state_file, index=True, header=False)
    return category, apply_pastml(annotation_file=state_file, tree_file=tree)


def col_name2cat(column):
    column_string = ''.join(s for s in column if s.isalnum())
    return column_string


def infer_ancestral_states(tree, data, work_dir, res_annotations):
    with ThreadPool() as pool:
        col2annotation_files = \
            pool.map(func=work, iterable=((tree, data[column], work_dir, column) for column in data.columns))
    logging.info('Combining the data from different columns...')
    pasml_annotations2cytoscape_annotation(dict(col2annotation_files), res_annotations)


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    import argparse

    parser = argparse.ArgumentParser(description="Processes data files.")

    parser.add_argument('--tree', help="the input tree in newick format.", type=str, required=True)
    parser.add_argument('--data', required=True, type=str,
                        help="the annotation file in tab/csv format with the first row containing the column names.")  
    parser.add_argument('--data_sep', required=False, type=str, default='\t',
                        help="the column separator for the data table. By default is set to tab, i.e. for tab file. " \
                             "Set it to ',' if your file is csv.")
    parser.add_argument('--id_index', required=False, type=int, default=0,
                        help="the index of the column in the data table that contains the tree tip names, "
                             "indices start from zero (by default is set to 0).")
    parser.add_argument('--columns', nargs='*',
                        help="names of the data table columns that contain states to be analysed with PASTML,"
                             "if not specified all columns will be considered.",
                        type=str)
    parser.add_argument('--name_column', type=str, default=None,
                        help="name of the data table column that should be used for node names in the visualisation"
                             "(should be one of those specified in columns, if columns are specified)."
                             "If the data table contains only one column it will be used by default.")
    parser.add_argument('--for_names_only', action='store_true', 
                        help="If specified, and if we are to analyse multiple state (specified in columns),"
                             "and the name_column is specified,"
                             "then the name_column won't be assigned a coloured section on the nodes, "
                             "but will only be shown as node names.")
    parser.add_argument('--work_dir', required=False, default=None, type=str,
                        help="the working dir for PASTML (if not specified a temporary dir will be created).")
    parser.add_argument('--html', required=False, default=None, type=str,
                        help="the output tree visualisation file (html).")
    parser.add_argument('--html_compressed', required=True, default=None, type=str,
                        help="the output summary map visualisation file (html).")
    params = parser.parse_args()

    using_temp_dir = False

    if not params.work_dir:
        using_temp_dir = True
        params.work_dir = tempfile.mkdtemp()

    df = pd.read_table(params.data, sep=params.data_sep, index_col=params.id_index, header=0)

    if not params.columns:
        params.columns = df.columns

    res_annotations = \
        os.path.join(params.work_dir, 'combined_annotations_{}.tab'.format('_'.join(params.columns)))

    if not os.path.isfile(res_annotations):
        infer_ancestral_states(tree=params.tree, work_dir=params.work_dir,
                               res_annotations=res_annotations, data=df[params.columns])

    tree, categories = annotate_tree_with_cyto_metadata(params.tree, res_annotations)

    if not params.name_column and len(params.columns) == 1:
        params.name_column = params.columns[0]

    if params.name_column:
        params.name_column = col_name2cat(params.name_column)
        if params.for_names_only and len(params.columns) > 1:
            categories.remove(params.name_column)

    df = pd.read_csv(res_annotations, index_col=0, header=0)
    name2colour = {}
    for cat in categories:
        unique_values = df[cat].unique()
        unique_values = sorted(unique_values[~pd.isnull(unique_values)].astype(str))
        num_unique_values = len(unique_values)
        colours = list(COLOURS)
        if len(colours) < num_unique_values:
            for _ in range(len(colours), num_unique_values):
                colours.append(random_hex_color())
        colours = colours[::int(len(COLOURS) / num_unique_values)]
        for value, col in zip(sorted(unique_values), colours):
            name2colour['{}_{}'.format(cat, value)] = col
        # let ambiguous values be white
        name2colour['{}_'.format(cat)] = WHITE

    if params.html:
        save_as_cytoscape_html(tree, params.html, categories=categories, graph_name='Tree', name2colour=name2colour)
    tree = compress_tree(tree, categories=categories, name_feature=params.name_column)
    save_as_cytoscape_html(tree, params.html_compressed, categories, graph_name='Summary map',
                           name2colour=name2colour, add_fake_nodes=False)

    if using_temp_dir:
        shutil.rmtree(params.work_dir)
