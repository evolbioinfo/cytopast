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
    11: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
         '#ffff99'],
    12: ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd',
         '#ccebc5', '#ffed6f']
}
WHITE = '#ffffff'


def random_hex_color():
    """
    Generates and returns a random HEX colour.
    :return: str, HEX colour
    """
    r = lambda: random.randint(100, 255)
    return '#%02X%02X%02X' % (r(), r(), r())


def _work(args):
    tree, df, work_dir, column, model = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)
    rep_dir = os.path.join(work_dir, category, model)
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
    return category, apply_pastml(annotation_file=state_file, tree_file=tree, model=model)


def _do_nothing(args):
    tree, df, work_dir, column, model = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)

    states = [s for s in df.unique() if not pd.isnull(s)]
    n_states = len(states)
    logging.info('States are {}'.format(states))

    tree = Tree(tree, 3)
    n_tips = len(tree.get_leaves())

    rep_dir = os.path.join(work_dir, category, model)
    res_file = os.path.join(rep_dir, STATES_TAB_PASTML_OUTPUT).format(tips=n_tips, states=n_states)
    if os.path.exists(res_file):
        return category, res_file
    os.makedirs(rep_dir, exist_ok=True)

    res_df = pd.DataFrame(index=[str(n.name) for n in tree.traverse()], columns=states)
    res_df.fillna(value=True, inplace=True)
    df.index = df.index.map(str)
    df = df.filter(res_df.index, axis=0)

    logging.info(df.head(3))

    for state in states:
        res_df.loc[df[df == state].index.tolist(), res_df.columns[res_df.columns != state]] = False

    res_df.astype(bool).to_csv(res_file, sep=',', index=True)

    return category, res_file


def col_name2cat(column):
    """
    Reformats the column string to make sure it contains only numerical or letter characters.
    :param column: str, column name to be reformatted
    :return: str, the column name with illegal characters removed
    """
    column_string = ''.join(s for s in column if s.isalnum())
    return column_string


def apply_pastml_to_all_columns(tree, data, work_dir, res_annotations, sep='\t', model='JC', copy_data=None):
    """
    Applies PASTML as many times as there are categories (columns) in the data,
    infers ancestor states and reformats them into an output tab/csv file.
    :param model: str (optional, default is 'JC'), model to be used by PASTML.
    :param tree: str, path to the tree in newick format.
    :param data: pandas.DataFrame containing tree tip names as indices and categories as columns.
    :param work_dir: str, path to the working dir where PASTML can place its temporary files.
    :param res_annotations: str, path to the file where the output table will be created.
    :param sep: char (optional, by default '\t'), the column separator for the output table.
    By default is set to tab, i.e. for tab file. Set it to ',' for csv.
    :return: void
    """
    with ThreadPool() as pool:
        col2annotation_files = \
            pool.map(func=_work, iterable=((tree, data[column], work_dir, column, model)
                                           for column in data.columns))
    if copy_data is not None:
        with ThreadPool() as pool:
            col2annotation_files.extend(
                pool.map(func=_do_nothing, iterable=((tree, copy_data[column], work_dir, column, model)
                                                     for column in copy_data.columns)))

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


def pastml_pipeline(tree, data, html_compressed, html=None, data_sep='\t', id_index=0, columns=None, name_column=None,
                    for_names_only=False, work_dir=None, all=False, model='JC', copy_columns=None):
    """
    Applies PASTML to the given tree with the specified states and visualizes the result (as html maps).
    :param tree: str, path to the input tree in newick format.
    :param data: str, path to the annotation file in tab/csv format with the first row containing the column names.
    :param html_compressed: str, path where the output summary map visualisation file (html) will be created.
    :param html: str (optional), path where the output tree visualisation file (html) will be created.
    :param data_sep: char (optional, by default '\t'), the column separator for the data table.
    By default is set to tab, i.e. for tab file. Set it to ',' if your file is csv.
    :param id_index: int (optional, by default is 0) the index of the column in the data table that contains the tree tip names,
    indices start from zero.
    :param columns: array of str (optional), names of the data table columns that contain states to be analysed with PASTML,
    if not specified all columns will be considered.
    :param name_column: str (optional), name of the data table column to be used for node names in the visualisation
    (must be one of those specified in columns, if columns are specified). If the data table contains only one column,
    it will be used by default.
    :param for_names_only: bool (optional, by default is False) If is set to True, and we are to analyse multiple states
    (specified in columns),and the name_column is specified,
    then the name_column won't be assigned a coloured section on the nodes, but will only be shown as node names.
    :param work_dir: str (optional), path to the working dir for PASTML (if not specified a temporary dir will be created).
    :param all: bool (optional, by default is False), if to keep all the nodes in the map, even the minor ones.
    :param model: str (optional, default is 'JC'), model to be used by PASTML.
    :return: void
    """
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
        os.path.join(work_dir, 'combined_annotations_{}_{}.tab'.format('_'.join(columns), model))

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
    else:
        root = read_tree(new_tree)

    apply_pastml_to_all_columns(tree=new_tree, data=df[columns], work_dir=work_dir, res_annotations=res_annotations,
                                sep=data_sep, model=model, copy_data=df[copy_columns] if copy_columns else None)

    past_vis(root, res_annotations, html_compressed, html, data_sep=data_sep,
             columns=(columns + copy_columns) if copy_columns else columns, name_column=name_column,
             for_names_only=for_names_only, all=all)

    if using_temp_dir:
        shutil.rmtree(work_dir)


def past_vis(tree, res_annotations, html_compressed, html=None, data_sep='\t', columns=None, name_column=None,
             for_names_only=False, all=False):
    """
    Applies PASTML to the given tree with the specified states and visualizes the result (as html maps).
    :param tree: ete3.Tree, the tree of interest.
    :param html_compressed: str, path where the output summary map visualisation file (html) will be created.
    :param html: str (optional), path where the output tree visualisation file (html) will be created.
    :param data_sep: char (optional, by default '\t'), the column separator for the annotation table.
    By default is set to tab, i.e. for tab file. Set it to ',' if your file is csv.
    :param columns: array of str, names of the data table columns that contain states to be analysed with PASTML.
    :param name_column: str (optional), name of the data table column to be used for node names in the visualisation
    (must be one of those specified in columns, if columns are specified). If columns contain only one column,
    it will be used by default.
    :param for_names_only: bool (optional, by default is False) If is set to True, and we are to analyse multiple states
    (specified in columns),and the name_column is specified,
    then the name_column won't be assigned a coloured section on the nodes, but will only be shown as node names.
    :param all: bool (optional, by default is False), if to keep all the nodes in the map, even the minor ones.
    :return: void
    """
    one_column = len(columns) == 1
    tree, categories = annotate_tree_with_cyto_metadata(tree, res_annotations, sep=data_sep, one_state=one_column)

    if not name_column and one_column:
        name_column = columns[0]

    if name_column:
        name_column = col_name2cat(name_column)
        if for_names_only or one_column:
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

    if one_column:
        n2tooltip = lambda n, categories: ', '.join('{}'.format(cat) for cat in categories
                                                    if hasattr(n, cat) and bool(getattr(n, cat, False)))
    else:
        n2tooltip = lambda n, categories: \
            ', '.join('{}:{}'.format(cat, getattr(n, cat)) for cat in categories
                      if hasattr(n, cat) and bool(getattr(n, cat, False)))

    if html:
        save_as_cytoscape_html(tree, html, categories=categories, name2colour=name2colour,
                               name_feature=name_column, n2tooltip=n2tooltip)
    tree = compress_tree(tree,
                         categories=([name_column] + categories) if name_column and not one_column else categories,
                         name_feature=name_column, cut=not all)
    save_as_cytoscape_html(tree, html_compressed, categories,
                           name2colour=name2colour, add_fake_nodes=False, n2tooltip=n2tooltip)


def get_enough_colours(num_unique_values):
    """
    Generates and returns an array of `num_unique_values` HEX colours.
    :param num_unique_values: int, number of colours to be generated.
    :return: array of str, containing colours in HEX format.
    """
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
                        help="names of the data table columns that contain states to be analysed with PASTML, "
                             "if not specified all columns will be considered.",
                        type=str)
    parser.add_argument('--copy_columns', nargs='*',
                        help="names of the data table columns that contain tip states to be used for visualisation "
                             "without applying PASTML (ancestral states would stay unresolved).",
                        type=str)
    parser.add_argument('--name_column', type=str, default=None,
                        help="name of the data table column to be used for node names in the visualisation"
                             "(must be one of those specified in columns, if columns are specified)."
                             "If the data table contains only one column it will be used by default.")
    parser.add_argument('--for_names_only', action='store_true',
                        help="If specified, and if we are to analyse multiple states (specified in columns),"
                             "and the name_column is specified,"
                             "then the name_column won't be assigned a coloured section on the nodes, "
                             "but will only be shown as node names.")
    parser.add_argument('--model', required=False, default='JC', type=str,
                        help="the model to be used by PASTML (can be JC or F81).")
    parser.add_argument('--work_dir', required=False, default=None, type=str,
                        help="the working dir for PASTML (if not specified a temporary dir will be created).")
    parser.add_argument('--all', action='store_true', help="Keep all the nodes in the map, even the minor ones.")
    params = parser.parse_args()

    pastml_pipeline(**vars(params))
