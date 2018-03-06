import logging
import os
import pastml
import shutil
import tempfile
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

from cytopast import compress_tree, read_tree, \
    pasml_annotations2cytoscape_annotation, annotate_tree_with_cyto_metadata, name_tree
from cytopast.colour_generator import get_enough_colours, WHITE
from cytopast.cytoscape_manager import save_as_cytoscape_html

STATES_TAB_PASTML_OUTPUT = 'Result.tree_{tree}.category_{category}.txt'


def _work(args):
    tree, df, work_dir, column, model, cache = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)
    res_dir = os.path.abspath(os.path.join(work_dir, category, model))

    # For binary states make sure that 0 or false is not mistaken by missing data
    if len(df.unique()) == 2 and np.any(pd.isnull(df.unique())):
        state = df.unique()[~pd.isnull(df.unique())][0]
        other_state = not state if isinstance(state, bool) \
            else (0 if isinstance(state, (int, float, complex)) and state != 0
                  else 'other' if state != 'other' else 'unknown')
        df.replace(np.nan, other_state, inplace=True)

    unique_states = [s for s in df.unique() if not pd.isnull(s)]
    logging.info('States are {}'.format(unique_states))
    tree_name = os.path.splitext(os.path.basename(tree))[0]
    res_file = os.path.join(res_dir, STATES_TAB_PASTML_OUTPUT).format(tree=tree_name, category=category)
    if cache and os.path.exists(res_file):
        return column, res_file
    os.makedirs(res_dir, exist_ok=True)
    state_file = os.path.join(res_dir, 'state_{}.csv'.format(category))

    # Keep only records corresponding to the tips in the state file
    names = df.index.astype(np.str)
    df = df[np.in1d(names, [n.name for n in read_tree(tree).iter_leaves()])]
    # Prepare the state file for PASTML
    df.to_csv(state_file, index=True, header=False)

    res_tree = os.path.join(res_dir, 'temp.tree.pastml.{tree}.{category}.nwk').format(tree=tree_name, category=category)

    hide_warnings = logging.getLogger().getEffectiveLevel() >= logging.ERROR
    try:
        pastml.infer_ancestral_states(state_file, os.path.abspath(tree), res_file, res_tree, model,
                                      1 if hide_warnings else 0)
    except:
        logging.error("PASTML could not infer states for {}, so we'll keep the tip states only.".format(category))
        return column, None

    # Try to remove a useless tree produced by PASTML
    try:
        os.remove(res_tree)
    except:
        pass
    return column, res_file


def _do_nothing(args):
    tree, df, work_dir, column, model, cache = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)

    states = [s for s in df.unique() if not pd.isnull(s)]
    logging.info('States are {}'.format(states))

    tree_name = os.path.splitext(os.path.basename(tree))[0]
    res_dir = os.path.abspath(os.path.join(work_dir, category, model))
    res_file = os.path.join(res_dir, STATES_TAB_PASTML_OUTPUT).format(tree=tree_name, category=category)

    if cache and os.path.exists(res_file):
        return column, res_file
    os.makedirs(res_dir, exist_ok=True)

    res_df = pd.DataFrame(index=[str(n.name) for n in read_tree(tree).traverse()], columns=states)
    res_df.fillna(value=1, inplace=True)
    df.index = df.index.map(str)
    df = df.filter(res_df.index, axis=0)

    for state in states:
        res_df.loc[df[df == state].index.tolist(), res_df.columns[res_df.columns != state]] = 0

    res_df.to_csv(res_file, sep=',', index=True)
    return column, res_file


def col_name2cat(column):
    """
    Reformats the column string to make sure it contains only numerical or letter characters.
    :param column: str, column name to be reformatted
    :return: str, the column name with illegal characters removed
    """
    column_string = ''.join(s for s in column if s.isalnum())
    return column_string


def get_ancestral_states_for_all_columns(tree, df, columns, work_dir, res_annotations, sep='\t', model='JC',
                                         copy_columns=None, cache=False):
    """
    Creates a table with node states by applying PASTML to infer ancestral states for categories specified in columns,
    and copying the states from copy_columns.
    :param cache: bool, if True the results of previous PASTML runs on this data will be reused when possible.
    :param columns: list of columns to be processed with PASTML.
    :param copy_columns: list of columns to be copied as-is.
    :param model: str (optional, default is 'JC'), model to be used by PASTML.
    :param tree: str, path to the tree in newick format.
    :param df: pandas.DataFrame containing tree tip names as indices and categories as columns.
    :param work_dir: str, path to the working dir where PASTML can place its temporary files.
    :param res_annotations: str, path to the file where the output table will be created.
    :param sep: char (optional, by default '\t'), the column separator for the output table.
    By default is set to tab, i.e. for tab file. Set it to ',' for csv.
    :return: void
    """
    col2annotation_files = []
    if columns is not None and len(columns):
        with ThreadPool() as pool:
            col2annotation_files = \
                dict(pool.map(func=_work, iterable=((tree, df[column], work_dir, column, model, cache)
                                                    for column in columns)))
    if copy_columns is None:
        copy_columns = []
    for col, af in col2annotation_files.items():
        if af is None:
            copy_columns.append(col)

    if copy_columns:
        with ThreadPool() as pool:
            col2annotation_files.update(
                dict(pool.map(func=_do_nothing, iterable=((tree, df[column], work_dir, column, model, cache)
                                                          for column in copy_columns))))

    col2annotation_files = {col_name2cat(col): af for (col, af) in col2annotation_files.items()}

    logging.info('Combining the data from different columns...')
    pasml_annotations2cytoscape_annotation(col2annotation_files, res_annotations, sep=sep)
    if len(col2annotation_files) == 1:
        df = pd.read_table(list(col2annotation_files.values())[0], sep=',', index_col=0, header=0).astype(bool)
        comb_df = pd.read_table(res_annotations, sep=sep, index_col=0, header=0)
        if set(df.columns) & set(comb_df.columns):
            df = df.join(comb_df, lsuffix='category', rsuffix='')
        else:
            df = df.join(comb_df)
        df.to_csv(res_annotations, sep=sep)


def quote(str_list):
    return ', '.join('"{}"'.format(_) for _ in str_list) if str_list is not None else ''


def pastml_pipeline(tree, data, out_data=None, html_compressed=None, html=None, data_sep='\t', id_index=0, columns=None,
                    name_column=None, work_dir=None, all=False, model='JC', copy_columns=None,
                    verbose=False, cache=False):
    """
    Applies PASTML to the given tree with the specified states and visualizes the result (as html maps).
    :param copy_columns: list of str (optional), names of the data table columns that contain states to be copied as-is,
    without applying PASTML (the missing states will stay unresolved).
    :param cache: bool, if True the results of previous PASTML runs on this data will be reused when possible.
    :param verbose: bool, print information on the progress of the analysis.
    :param out_data: str, path to the output annotation file with the states inferred by PASTML.
    :param tree: str, path to the input tree in newick format.
    :param data: str, path to the annotation file in tab/csv format with the first row containing the column names.
    :param html_compressed: str, path where the output summary map visualisation file (html) will be created.
    :param html: str (optional), path where the output tree visualisation file (html) will be created.
    :param data_sep: char (optional, by default '\t'), the column separator for the data table.
    By default is set to tab, i.e. for tab file. Set it to ',' if your file is csv.
    :param id_index: int (optional, by default is 0) the index of the column in the data table
    that contains the tree tip names, indices start from zero.
    :param columns: list of str (optional), names of the data table columns that contain states
    to be analysed with PASTML, if not specified all columns will be considered.
    :param name_column: str (optional), name of the data table column to be used for node names in the visualisation
    (must be one of those specified in columns, if columns are specified). If the data table contains only one column,
    it will be used by default.
    :param work_dir: str (optional), path to the working dir for PASTML
    (if not specified a temporary dir will be created).
    :param all: bool (optional, by default is False), if to keep all the nodes in the map, even the minor ones.
    :param model: str (optional, default is 'JC'), model to be used by PASTML.
    :return: void
    """
    logging.basicConfig(level=logging.INFO if verbose else logging.ERROR,
                        format='%(asctime)s: %(message)s', datefmt="%H:%M:%S", filename=None)

    using_temp_dir = False

    if not work_dir:
        using_temp_dir = True
        work_dir = tempfile.mkdtemp()
    else:
        os.makedirs(work_dir, exist_ok=True)

    df = pd.read_table(data, sep=data_sep, index_col=id_index, header=0)

    unknown_columns = ((set(columns) if columns else set())
                       | (set(copy_columns) if copy_columns else set())) - set(df.columns)
    if unknown_columns:
        raise ValueError('{} of the specified columns ({}) {} not found among the annotation columns: {}.'
                         .format('One' if len(unknown_columns) == 1 else 'Some',
                                 quote(unknown_columns),
                                 'is' if len(unknown_columns) == 1 else 'are',
                                 quote(df.columns)))

    if not columns and not copy_columns:
        columns = df.columns

    if name_column and (not columns or name_column not in columns) \
            and (not copy_columns or name_column not in copy_columns):
        raise ValueError('The name column ({}) should be one of those specified as columns ({}) or copy_columns ({}).'
                         .format(quote([name_column]), quote(columns), quote(copy_columns)))

    tree_name = os.path.basename(tree)

    res_annotations = out_data if out_data else \
        os.path.join(work_dir, 'combined_annotations_{}_{}_{}.tab'.format('_'.join(columns), model, tree_name))

    new_tree = os.path.join(work_dir, tree_name + '.pastml.nwk')
    if not cache or not os.path.exists(new_tree):
        root = read_tree(tree)
        names = df.index.astype(np.str)
        node_names = [n.name for n in root.iter_leaves() if n.name not in names]
        if node_names:
            missing_value_df = pd.DataFrame(index=node_names, columns=df.columns)
            df = df.append(missing_value_df)
        name_tree(root)
        root.write(outfile=new_tree, format=3, format_root_node=True)
    else:
        root = read_tree(new_tree)

    get_ancestral_states_for_all_columns(tree=new_tree, df=df, columns=columns, work_dir=work_dir,
                                         res_annotations=res_annotations,
                                         sep=data_sep, model=model, copy_columns=copy_columns, cache=cache)

    _past_vis(root, res_annotations, html_compressed, html, data_sep=data_sep,
              columns=(columns + copy_columns) if copy_columns else columns, name_column=name_column, all=all)

    if using_temp_dir:
        shutil.rmtree(work_dir)


def _past_vis(tree, res_annotations, html_compressed=None, html=None, data_sep='\t', columns=None, name_column=None,
              all=False):
    one_column = len(columns) == 1
    tree, categories = annotate_tree_with_cyto_metadata(tree, res_annotations, sep=data_sep, one_state=one_column)

    if not name_column and one_column:
        name_column = columns[0]

    if name_column:
        name_column = col_name2cat(name_column)
        if one_column:
            categories.remove(name_column)

    df = pd.read_table(res_annotations, index_col=0, header=0, sep=data_sep)
    name2colour = {}
    if len(columns) > 1:
        for cat in categories:
            unique_values = df[cat].unique()
            unique_values = sorted(unique_values[~pd.isnull(unique_values)].astype(str))
            num_unique_values = len(unique_values)
            colours = get_enough_colours(num_unique_values)
            for value, col in zip(unique_values, colours):
                name2colour['{}_{}'.format(cat, value)] = col
            # let ambiguous values be white
            name2colour['{}_'.format(cat)] = WHITE
    else:
        colours = get_enough_colours(len(categories))
        for cat, col in zip(categories, colours):
            name2colour['{}_{}'.format(cat, True)] = col

    if one_column:
        n2tooltip = lambda n, cats: ', '.join('{}'.format(_) for _ in cats
                                              if hasattr(n, _) and getattr(n, _, '') != '')
    else:
        n2tooltip = lambda n, cats: \
            ', '.join('{}:{}'.format(_, getattr(n, _)) for _ in cats if hasattr(n, _) and getattr(n, _, '') != '')

    if html:
        save_as_cytoscape_html(tree, html, categories=categories, name2colour=name2colour, n2tooltip=n2tooltip,
                               name_feature='name')

    if html_compressed:
        tree = compress_tree(tree,
                             categories=([name_column] + categories) if name_column and not one_column else categories,
                             name_feature=name_column, cut=not all)
        save_as_cytoscape_html(tree, html_compressed, categories,
                               name2colour=name2colour, add_fake_nodes=False, n2tooltip=n2tooltip)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Visualisation of annotated phylogenetic trees (as html maps).")

    annotation_group = parser.add_argument_group('annotation-related arguments')
    annotation_group.add_argument('-d', '--data', required=True, type=str,
                                  help="the annotation file in tab/csv format with the first row "
                                       "containing the column names.")
    annotation_group.add_argument('-s', '--data_sep', required=False, type=str, default='\t',
                                  help="the column separator for the data table. "
                                       "By default is set to tab, i.e. for tab file. " 
                                       "Set it to ',' if your file is csv.")
    annotation_group.add_argument('-i', '--id_index', required=False, type=int, default=0,
                                  help="the index of the column in the data table that contains the tree tip names, "
                                       "indices start from zero (by default is set to 0).")
    annotation_group.add_argument('-c', '--columns', nargs='*',
                                  help="names of the data table columns that contain states "
                                       "to be analysed with PASTML. "
                                       "If neither columns nor copy_columns are specified, "
                                       "then all columns will be considered for PASTMl analysis.",
                                  type=str)
    annotation_group.add_argument('--copy_columns', nargs='*',
                                  help="names of the data table columns that contain states to be copied as-is, "
                                       "without applying PASTML (the missing states will stay unresolved).",
                                  type=str)

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('-t', '--tree', help="the input tree in newick format.", type=str, required=True)

    pastml_group = parser.add_argument_group('ancestral-state inference-related arguments')
    pastml_group.add_argument('-m', '--model', required=False, default='JC', type=str,
                              help="the evolutionary model to be used by PASTML (can be JC or F81).")
    pastml_group.add_argument('--work_dir', required=False, default=None, type=str,
                              help="the working dir for PASTML to put intermediate files into "
                                   "(if not specified a temporary dir will be created).")
    pastml_group.add_argument('--cache', action='store_true',
                              help="if set, the results of previous PASTML runs on this data will be reused "
                                   "when possible")

    vis_group = parser.add_argument_group('visualisation-related arguments')
    vis_group.add_argument('-n', '--name_column', type=str, default=None,
                           help="name of the data table column to be used for node names "
                                "in the compressed map visualisation"
                                "(must be one of those specified in columns or copy_columns if they are specified)."
                                "If the data table contains only one column it will be used by default.")
    vis_group.add_argument('-a', '--all', action='store_true',
                           help="Keep all the nodes in the compressed map visualisation, "
                                "even the minor ones.")

    out_group = parser.add_argument_group('output-related arguments')
    out_group.add_argument('-o', '--out_data', required=False, type=str,
                           help="the output annotation file with the states inferred by PASTML.")
    out_group.add_argument('-p', '--html_compressed', required=False, default=None, type=str,
                           help="the output summary map visualisation file (html).")
    out_group.add_argument('-l', '--html', required=False, default=None, type=str,
                           help="the output tree visualisation file (html).")

    parser.add_argument('-v', '--verbose', action='store_true',
                        help="print information on the progress of the analysis")
    params = parser.parse_args()

    pastml_pipeline(**vars(params))


if '__main__' == __name__:
    main()
