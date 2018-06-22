import logging
import os
import shutil
import tempfile
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd
import pastml
from pastml import JOINT, MARGINAL, MARGINAL_APPROXIMATION, MAX_POSTERIORI, JC, F81, DOWNPASS, ACCTRAN, DELTRAN

from cytopast import compress_tree, read_tree, \
    pasml_annotations2cytoscape_annotation, annotate_tree_with_cyto_metadata, name_tree, collapse_zero_branches, \
    col_name2cat, REASONABLE_NUMBER_OF_TIPS, date_tips, DATE
from cytopast.colour_generator import get_enough_colours, WHITE
from cytopast.cytoscape_manager import save_as_cytoscape_html

STATES_TAB_PASTML_OUTPUT = 'Result.tree_{tree}.category_{category}.csv'
PARAM_TAB_PASTML_OUTPUT = 'Result.tree_{tree}.category_{category}.params.csv'
PARAM_TAB_PASTML_INPUT = 'Result.tree_{tree}.category_{category}.input.params.csv'


def _work(args):
    tree, df, work_dir, column, model, prob_method, params = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)
    res_dir = os.path.abspath(os.path.join(work_dir, category, model, prob_method))

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
    os.makedirs(res_dir, exist_ok=True)

    res_file = os.path.join(res_dir, STATES_TAB_PASTML_OUTPUT).format(tree=tree_name, category=category)
    in_param_file = ''
    if params:
        in_param_file = parse_parameters(category, params, res_dir, tree_name, unique_states)
    out_param_file = os.path.join(res_dir, PARAM_TAB_PASTML_OUTPUT).format(tree=tree_name, category=category)
    state_file = os.path.join(res_dir, 'state_{}.csv'.format(category))

    percentage_unknown = df.isnull().sum() / len(df)
    if percentage_unknown > .9:
        raise ValueError('%.1f%% of tip annotations for %s are unknown, not enough data to infer ancestral states. '
                         'Check your annotation file and if its id column corresponds to the tree tip names.'
                         % (percentage_unknown * 100, column))

    percentage_unique = df.nunique() / df.count()
    if df.count() > 100 and percentage_unique > .5:
        raise ValueError('The column {} seem to contain non-categorical data: {} of values are unique.'
                         'PASTML cannot infer ancestral states for a tree with too many tip states.'
                         .format(column, int(100 * percentage_unique)))
    # Prepare the state file for PASTML
    df.to_csv(state_file, index=True, header=False, encoding='ascii')

    res_tree = os.path.join(res_dir, 'temp.tree.pastml.{tree}.{category}.nwk').format(tree=tree_name, category=category)

    hide_warnings = logging.getLogger().getEffectiveLevel() >= logging.ERROR
    try:
        pastml.infer_ancestral_states(state_file, os.path.abspath(tree), in_param_file, res_file, res_tree,
                                      out_param_file, model, prob_method, 1 if hide_warnings else 0)
    except:
        logging.error("PASTML could not infer states for {}, so we'll keep the tip states only.".format(category))
        return column, None

    # Try to remove a useless tree produced by PASTML
    try:
        os.remove(res_tree)
    except:
        pass
    return column, res_file


def parse_parameters(category, params, res_dir, tree_name, unique_states):
    if isinstance(params, dict):
        param_df = pd.Series(data=list(params.values()), index=params.keys())
        param_df = param_df[np.in1d(param_df.index, unique_states + ['epsilon', 'scaling factor'])]
        freq_df = param_df[np.in1d(param_df.index, unique_states)]
        if len(freq_df):
            try:
                freq_df = freq_df.astype(np.float)
            except:
                logging.error('Specified frequencies ({}) are not float,'
                              'ignoring them.'.format(freq_df))
                param_df = param_df[np.in1d(param_df.index, ['epsilon', 'scaling factor'])]
            if len(freq_df) != len(unique_states):
                logging.error('Frequency parameters are specified ({}), but not for all of the states ({}), '
                              'ignoring them.'.format(freq_df.columns, unique_states))
                param_df = param_df[np.in1d(param_df.index, ['epsilon', 'scaling factor'])]
            elif sum(freq_df) != 1:
                logging.error('Specified frequencies ({}) do not sum up to one,'
                              'ignoring them.'.format(freq_df))
                param_df = param_df[np.in1d(param_df.index, ['epsilon', 'scaling factor'])]
            elif any(freq_df < 0):
                logging.error('Specified frequencies ({}) must not be negative,'
                              'ignoring them.'.format(freq_df))
                param_df = param_df[np.in1d(param_df.index, ['epsilon', 'scaling factor'])]
        if 'epsilon' in param_df.index:
            try:
                epsilon = float(param_df['epsilon'])
                if epsilon < 0:
                    logging.error('Epsilon ({}) cannot be negative, ignoring it.'.format(param_df['epsilon']))
                    param_df = param_df[np.in1d(param_df.index, ['epsilon'])]
            except:
                logging.error('Epsilon ({}) is not float, ignoring it.'.format(param_df['epsilon']))
                param_df = param_df[np.in1d(param_df.index, ['epsilon'])]
        if 'scaling factor' in param_df.index:
            try:
                epsilon = float(param_df['scaling factor'])
                if epsilon < 0:
                    logging.error(
                        'Scaling factor ({}) cannot be negative, ignoring it.'.format(param_df['scaling factor']))
                    param_df = param_df[np.in1d(param_df.index, ['scaling factor'])]
            except:
                logging.error('Scaling factor ({}) is not float, ignoring it.'.format(param_df['scaling factor']))
                param_df = param_df[np.in1d(param_df.index, ['scaling factor'])]
        if len(param_df):
            in_param_file = os.path.join(res_dir, PARAM_TAB_PASTML_INPUT).format(tree=tree_name, category=category)
            param_df.to_csv(in_param_file)
    elif isinstance(params, str):
        if not os.path.exists(params):
            raise ValueError('You have specified some parameters ({}) but such a file does not exist!'
                             .format(params))
        in_param_file = params
    else:
        raise ValueError('Parameters must be specified either as a dict or as a path to a csv file, not as {}!'
                         .format(type(params)))
    return in_param_file


def _do_nothing(args):
    tree, df, work_dir, column = args
    logging.info('Processing {}'.format(column))
    category = col_name2cat(column)

    # For binary states make sure that 0 or false is not mistaken by missing data
    if len(df.unique()) == 2 and np.any(pd.isnull(df.unique())):
        state = df.unique()[~pd.isnull(df.unique())][0]
        other_state = not state if isinstance(state, bool) \
            else (0 if isinstance(state, (int, float, complex)) and state != 0
                  else 'other' if state != 'other' else 'unknown')
        df.replace(np.nan, other_state, inplace=True)

    states = [s for s in df.unique() if not pd.isnull(s)]
    logging.info('States are {}'.format(states))

    tree_name = os.path.splitext(os.path.basename(tree))[0]
    res_dir = os.path.abspath(os.path.join(work_dir, category, 'copy'))
    res_file = os.path.join(res_dir, STATES_TAB_PASTML_OUTPUT).format(tree=tree_name, category=category)
    os.makedirs(res_dir, exist_ok=True)

    res_df = pd.DataFrame(index=[str(n.name) for n in read_tree(tree).traverse()], columns=states)
    res_df.fillna(value=0, inplace=True)
    df.index = df.index.map(str)
    df = df.filter(res_df.index, axis=0)

    for state in states:
        res_df.loc[df[df == state].index.tolist(), state] = 1

    res_df.to_csv(res_file, index=True)
    return column, res_file


def get_ancestral_states_for_all_columns(tree, df, tip_df, columns, work_dir, res_annotations, sep='\t', model=JC,
                                         copy_columns=None, prediction_method=MARGINAL_APPROXIMATION,
                                         column2parameters=None):
    """
    Creates a table with node states by applying PASTML to infer ancestral states for categories specified in columns,
    and copying the states from copy_columns.
    :param prediction_method: str (optional, default is pastml.MARGINAL_APPROXIMATION),
    ancestral state prediction method to be used by PASTML.
    :param columns: list of columns to be processed with PASTML.
    :param copy_columns: list of columns to be copied as-is.
    :param model: str (optional, default is pastml.JC), model to be used by PASTML.
    :param tree: str, path to the tree in newick format.
    :param df: pandas.DataFrame containing tree node names as indices and categories as columns.
    :param tip_df: pandas.DataFrame containing only tree tip names as indices and categories as columns.
    :param work_dir: str, path to the working dir where PASTML can place its temporary files.
    :param res_annotations: str, path to the file where the output table will be created.
    :param sep: char (optional, by default '\t'), the column separator for the output table.
    By default is set to tab, i.e. for tab file. Set it to ',' for csv.
    :return: void
    """
    col2annotation_files = {}
    if columns is not None and len(columns):
        with ThreadPool() as pool:
            col2annotation_files = \
                dict(pool.map(func=_work, iterable=((tree, tip_df[column], work_dir, column, model, prediction_method,
                                                     column2parameters[column] if column2parameters
                                                                                  and column in column2parameters
                                                     else None) for column in columns)))
    if copy_columns is None:
        copy_columns = []
    for col, af in col2annotation_files.items():
        if af is None:
            copy_columns.append(col)

    if copy_columns:
        with ThreadPool() as pool:
            col2annotation_files.update(
                dict(pool.map(func=_do_nothing, iterable=((tree, df[column], work_dir, column)
                                                          for column in copy_columns))))

    col2annotation_files = {col_name2cat(col): af for (col, af) in col2annotation_files.items()}

    pasml_annotations2cytoscape_annotation(col2annotation_files, res_annotations, sep=sep)


def quote(str_list):
    return ', '.join('"{}"'.format(_) for _ in str_list) if str_list is not None else ''


def pastml_pipeline(tree, data, out_data=None, html_compressed=None, html=None, data_sep='\t', id_index=0, columns=None,
                    name_column=None, work_dir=None, tip_size_threshold=REASONABLE_NUMBER_OF_TIPS,
                    model=F81, prediction_method=MARGINAL_APPROXIMATION,
                    copy_columns=None, verbose=False, date_column=None, column2parameters=None):
    """
    Applies PASTML to the given tree with the specified states and visualizes the result (as html maps).

    :param date_column: str (optional), name of the data table column that contains tip dates.
    (following pandas specification, by default is inferred).
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
    :param copy_columns: list of str (optional), names of the data table columns that contain states to be copied as-is,
    without applying PASTML (the missing states will stay unresolved).
    :param name_column: str (optional), name of the data table column to be used for node names in the visualisation
    (must be one of those specified in columns, if columns are specified). If the data table contains only one column,
    it will be used by default.
    :param work_dir: str (optional), path to the working dir for PASTML
    (if not specified a temporary dir will be created).
    :param tip_size_threshold: int (optional, by default is 25), remove the tips of size less than threshold-th
    from the compressed map (set to inf to keep all).
    :param model: str (optional, default is pastml.F81), model to be used by PASTML.
    :param prediction_method: str (optional, default is pastml.MARGINAL_APPROXIMATION),
    ancestral state prediction method to be used by PASTML.
    :param verbose: bool, print information on the progress of the analysis.
    :param column2parameters: dict, an optional way to fix some parameters, must be in a form {column: {param: value}},
    where param can be a state (then the value should specify its frequency between 0 and 1),
    "scaling factor" (then the value should be the scaling factor for three branches,
    e.g. set to 1 to keep the original branches),
    or "epsilon" (the values specifies a min tip branch length to use for smoothing)
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

    if not copy_columns:
        copy_columns = []
    if not columns:
        columns = []

    df = pd.read_table(data, sep=data_sep, index_col=id_index, header=0)
    df.index = df.index.map(str)

    if date_column and date_column not in df.columns:
        raise ValueError('The date column {} not found among the annotation columns: {}.'
                         .format(date_column, quote(df.columns)))
    if date_column:
        if df[date_column].dtype == float:
            df[date_column] = pd.to_datetime(df[date_column], format='%Y.0')
        else:
            df[date_column] = pd.to_datetime(df[date_column], infer_datetime_format=True)

    unknown_columns = (set(columns) | set(copy_columns)) - set(df.columns)
    if unknown_columns:
        raise ValueError('{} of the specified columns ({}) {} not found among the annotation columns: {}.'
                         .format('One' if len(unknown_columns) == 1 else 'Some',
                                 quote(unknown_columns),
                                 'is' if len(unknown_columns) == 1 else 'are',
                                 quote(df.columns)))

    if not columns and not copy_columns:
        columns = list(df.columns)

    if not columns and not copy_columns:
        raise ValueError('Could not find any states in the annotation file {}. '
                         'Make sure that the file separator (--data_sep {}) is set correctly.'
                         .format(data, data_sep))

    if name_column and (name_column not in columns) and (name_column not in copy_columns):
        raise ValueError('The name column ({}) should be one of those specified as columns ({}) or copy_columns ({}).'
                         .format(quote([name_column]), quote(columns), quote(copy_columns)))

    tree_name = os.path.basename(tree)

    res_annotations = out_data if out_data else \
        os.path.join(work_dir, 'combined_annotations_{}_{}_{}_{}.tab'
                     .format('_'.join(columns + copy_columns), model, prediction_method, tree_name))
    new_tree = os.path.join(work_dir, tree_name + '.pastml.nwk')
    root = read_tree(tree)
    name_tree(root)
    collapse_zero_branches(root)
    # this is expected by PASTML
    root.name = 'ROOT'
    root.dist = 0
    root.write(outfile=new_tree, format=3, format_root_node=True)

    # append null-valued rows for missing tips (needed for PASTML)
    tip_names = np.array([n.name for n in root])
    missing_tip_names = np.setdiff1d(tip_names, df.index.astype(np.str), assume_unique=True)
    if len(missing_tip_names):
        missing_value_df = pd.DataFrame(index=missing_tip_names, columns=df.columns)
        df = df.append(missing_value_df)

    get_ancestral_states_for_all_columns(tree=new_tree, df=df, tip_df=df[np.in1d(df.index.astype(np.str), tip_names)],
                                         columns=columns, work_dir=work_dir,
                                         res_annotations=res_annotations,
                                         sep=data_sep, model=model, copy_columns=copy_columns,
                                         prediction_method=prediction_method, column2parameters=column2parameters)

    if html or html_compressed:
        if date_column:
            min_date, max_date = date_tips(root, df[date_column])
        else:
            min_date, max_date = 0, 0
        logging.info("Dates vary between {} and {}.".format(min_date, max_date))
        root = _past_vis(root, res_annotations, html_compressed, html, data_sep=data_sep,
                         columns=(columns + copy_columns) if copy_columns else columns, name_column=name_column,
                         tip_size_threshold=tip_size_threshold, min_date=min_date, max_date=max_date)

    if using_temp_dir:
        shutil.rmtree(work_dir)
    return root


def _past_vis(tree, res_annotations, html_compressed=None, html=None, data_sep='\t', columns=None, name_column=None,
              tip_size_threshold=REASONABLE_NUMBER_OF_TIPS, min_date=0, max_date=0):
    one_column = len(columns) == 1
    tree, categories = annotate_tree_with_cyto_metadata(tree, res_annotations, columns=columns, sep=data_sep)

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
            logging.info('Mapped the values to colours as following: {}, {}'.format(unique_values, colours))
            # let ambiguous values be white
            name2colour['{}_'.format(cat)] = WHITE
    else:
        colours = get_enough_colours(len(categories))
        for cat, col in zip(categories, colours):
            name2colour['{}_{}'.format(cat, True)] = col
        logging.info('Mapped the values to colours as following: {}, {}'.format(categories, colours))

    # set internal node dates to min of its tips' dates
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            if not hasattr(n, DATE):
                n.add_feature(DATE, 0)
        else:
            n.add_feature(DATE, min(getattr(_, DATE) for _ in n))

    def get_category_str(n):
        if one_column:
            return '{}: {}'.format(columns[0], ' or '.join('{}'.format(_) for _ in categories if hasattr(n, _)
                                                           and getattr(n, _, '') != ''))
        return '<br>'.join('{}: {}'.format(_, getattr(n, _)) for _ in categories if hasattr(n, _)
                           and getattr(n, _, '') != '')

    if html:
        save_as_cytoscape_html(tree, html, categories=categories, name2colour=name2colour,
                               n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               name_feature='name', min_date=min_date, max_date=max_date)

    if html_compressed:
        tree = compress_tree(tree, categories=categories, tip_size_threshold=tip_size_threshold)
        save_as_cytoscape_html(tree, html_compressed, categories=categories,
                               name2colour=name2colour, add_fake_nodes=False,
                               n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               min_date=min_date, max_date=max_date, name_feature=name_column)
    return tree


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
    annotation_group.add_argument('--date_column', required=False, default=None,
                                  help="name of the data table column that contains tip dates.",
                                  type=str)

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('-t', '--tree', help="the input tree in newick format.", type=str, required=True)

    pastml_group = parser.add_argument_group('ancestral-state inference-related arguments')
    pastml_group.add_argument('-m', '--model', required=False, default=F81, choices=[JC, F81], type=str,
                              help='the evolutionary model to be used by PASTML, by default {}.'.format(F81))
    pastml_group.add_argument('--prediction_method', required=False, default=MARGINAL_APPROXIMATION,
                              choices=[MARGINAL_APPROXIMATION, MARGINAL, MAX_POSTERIORI, JOINT,
                                       DOWNPASS, ACCTRAN, DELTRAN], type=str,
                              help='the ancestral state prediction method to be used by PASTML, '
                                   'by default {}.'.format(MARGINAL_APPROXIMATION))
    pastml_group.add_argument('--column2parameters', required=False, default=None, type=dict,
                              help='optional way to fix some parameters, must be in a form {column: {param: value}}, '
                                   'where param can be a state (then the value should specify its frequency between 0 and 1),'
                                   '"scaling factor" (then the value should be the scaling factor for three branches, '
                                   'e.g. set to 1 to keep the original branches), '
                                   'or "epsilon" (the values specifies a min tip branch length to use for smoothing)')
    pastml_group.add_argument('--work_dir', required=False, default=None, type=str,
                              help="the working dir for PASTML to put intermediate files into "
                                   "(if not specified a temporary dir will be created).")

    vis_group = parser.add_argument_group('visualisation-related arguments')
    vis_group.add_argument('-n', '--name_column', type=str, default=None,
                           help="name of the data table column to be used for node names "
                                "in the compressed map visualisation"
                                "(must be one of those specified in columns or copy_columns if they are specified)."
                                "If the data table contains only one column it will be used by default.")
    vis_group.add_argument('--tip_size_threshold', type=int, default=REASONABLE_NUMBER_OF_TIPS,
                           help="Remove the tips of size less than the threshold-th from the compressed map "
                                "(set to inf to keep all tips).")

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
