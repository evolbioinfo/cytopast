import logging
import os
import random
import shutil
import tempfile
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd
from ete3 import Tree

from cytopast import apply_pastml, compress_tree, STATES_TAB_PASTML_OUTPUT, read_tree, \
    pasml_annotations2cytoscape_annotation, annotate_tree_with_cyto_metadata, name_tree
from cytopast.cytoscape_manager import save_as_cytoscape_html

NUM2COLOURS = {
    1: ['#fdc086'],
    2: ['#b2df8a', '#1f78b4'],
    3: ['#8dd3c7', '#ffffb3', '#bebada'],
    4: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3'],
    5: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00'],
    6: ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f'],
    7: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f'],
    8: ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5'],
    9: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6'],
    10: ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd'],
    11: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99'],
    12: ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f']
}
WHITE = '#ffffff'


def random_hex_color():
    r = lambda: random.randint(100, 255)
    return '#%02X%02X%02X' % (r(), r(), r())


def work(args):
    tree, df, work_dir, column = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)
    rep_dir = os.path.join(work_dir, category)
    unique_states = df.unique()
    n_tips = len(Tree(tree, 3).get_leaves())
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


def infer_ancestral_states(tree, data, work_dir, res_annotations, sep='\t'):
    with ThreadPool() as pool:
        col2annotation_files = \
            pool.map(func=work, iterable=((tree, data[column], work_dir, column) for column in data.columns))
    logging.info('Combining the data from different columns...')
    col2annotation_files = dict(col2annotation_files)
    pasml_annotations2cytoscape_annotation(col2annotation_files, res_annotations, sep=sep)
    if len(col2annotation_files) == 1:
        df = pd.read_table(list(col2annotation_files.values())[0], sep=',', index_col=0, header=0)
        comb_df = pd.read_table(res_annotations, sep=sep, index_col=0, header=0)
        if set(df.columns) & set(comb_df.columns):
            df = df.join(comb_df, lsuffix='category', rsuffix='')
        else:
            df = df.join(comb_df)
        df.to_csv(res_annotations, sep=sep)


def pastml(tree, data, html_compressed, html=None, data_sep='\t', id_index=0, columns=None, name_column=None,
           for_names_only=False, work_dir=None, all=False):
    using_temp_dir = False

    if not work_dir:
        using_temp_dir = True
        work_dir = tempfile.mkdtemp()
    else:
        os.makedirs(work_dir, exist_ok=True)

    df = pd.read_table(data, sep=data_sep, index_col=id_index, header=0)

    if not columns:
        columns = df.columns

    res_annotations = \
        os.path.join(work_dir, 'combined_annotations_{}.tab'.format('_'.join(columns)))

    new_tree = os.path.join(work_dir, os.path.basename(tree) + '.pastml.nwk')
    if not os.path.exists(new_tree):
        root = read_tree(tree)
        names = df.index.astype(np.str)
        nodes = [n for n in root.iter_leaves() if n.name in names]
        n_tips = len(nodes)
        need_to_prune = len(root.get_leaves()) > n_tips
        if need_to_prune:
            logging.info('Pruning...')
            root.prune(nodes, preserve_branch_length=True)
        name_tree(root)
        root.write(outfile=new_tree, format=3, format_root_node=True)
    tree = new_tree

    # if not os.path.isfile(res_annotations):
    infer_ancestral_states(tree=tree, work_dir=work_dir,
                           res_annotations=res_annotations, data=df[columns], sep=data_sep)

    tree, categories = annotate_tree_with_cyto_metadata(tree, res_annotations, sep=data_sep,
                                                        one_state=len(columns) == 1)

    if not name_column and len(columns) == 1:
        name_column = columns[0]

    if name_column:
        name_column = col_name2cat(name_column)
        if for_names_only or len(columns) == 1:
            categories.remove(name_column)

    df = pd.read_table(res_annotations, index_col=0, header=0, sep=data_sep)
    name2colour = {}
    if len(columns) > 1:
        for cat in categories:
            unique_values = df[cat].unique()
            unique_values = sorted(unique_values[~pd.isnull(unique_values)].astype(str))
            num_unique_values = len(unique_values)
            colours = get_enough_colours(num_unique_values)
            for value, col in zip(sorted(unique_values), colours):
                name2colour['{}_{}'.format(cat, value)] = col
            # let ambiguous values be white
            name2colour['{}_'.format(cat)] = WHITE
    else:
        colours = get_enough_colours(len(categories))
        for cat, col in zip(categories, colours):
            name2colour['{}_{}'.format(cat, True)] = col

    if len(columns) == 1:
        n2tooltip = lambda n, categories: ', '.join('{}'.format(cat) for cat in categories
                                                    if hasattr(n, cat) and bool(getattr(n, cat, False)))
    else:
        n2tooltip = lambda n, categories: \
            ', '.join('{}:{}'.format(cat, getattr(n, cat)) for cat in categories if hasattr(n, cat))

    if html:
        save_as_cytoscape_html(tree, html, categories=categories, graph_name='Tree', name2colour=name2colour,
                               name_feature=name_column, n2tooltip=n2tooltip)
    tree = compress_tree(tree, categories=(categories + [name_column]) if name_column else categories,
                         name_feature=name_column, cut=not all)
    save_as_cytoscape_html(tree, html_compressed, categories, graph_name='Summary map',
                           name2colour=name2colour, add_fake_nodes=False, n2tooltip=n2tooltip)

    if using_temp_dir:
        shutil.rmtree(work_dir)


def get_enough_colours(num_unique_values):
    if num_unique_values <= 12:
        colours = NUM2COLOURS[num_unique_values]
    else:
        colours = NUM2COLOURS[12]
        for _ in range(len(colours), num_unique_values):
            colours.append(random_hex_color())
    return colours


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    import argparse

    parser = argparse.ArgumentParser(description="Processes data files.")

    parser.add_argument('--tree', help="the input tree in newick format.", type=str, required=True)
    parser.add_argument('--data', required=True, type=str,
                        help="the annotation file in tab/csv format with the first row containing the column names.")
    parser.add_argument('--html_compressed', required=True, default=None, type=str,
                        help="the output summary map visualisation file (html).")
    parser.add_argument('--html', required=False, default=None, type=str,
                        help="the output tree visualisation file (html).")
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
    parser.add_argument('--all', action='store_true', help="Keep all the nodes in the map, even the minor ones.")
    params = parser.parse_args()

    pastml(**vars(params))
