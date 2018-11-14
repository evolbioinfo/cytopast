import logging
import os

import numpy as np
import pandas as pd

from cytopast import col_name2cat
from pypastml.ml import is_ml, is_marginal

CYTOPAST_ANCESTRAL_STATES_CSV_ML = 'combined_ancestral_states.states_{states}.method_{method}.model_{model}.csv'
CYTOPAST_ANCESTRAL_STATES_CSV_MP = 'combined_ancestral_states.states_{states}.method_{method}.csv'
CYTOPAST_NAMED_TREE_NWK = 'named.tree_{tree}'

PASTML_ANCESTRAL_STATES_CSV_ML = 'ancestral_states.state_{state}.method_{method}.model_{model}.csv'
PASTML_ANCESTRAL_STATES_CSV_MP = 'ancestral_states.state_{state}.method_{method}.csv'
PASTML_TIP_STATES_CSV = 'input_tip_states.state_{state}.csv'
PASTML_ML_PARAMS_CSV = 'params.state_{state}.method_{method}.model_{model}.csv'
PASTML_MP_PARAMS_CSV = 'params.state_{state}.method_{method}.csv'
PASTML_MARGINAL_PROBS_CSV = 'marginal_probabilities.state_{state}.model_{model}.csv'
PASTML_PARAMS_CSV_INPUT = 'input_params.state_{state}.method_{method}.model_{model}.csv'


def get_pastml_parameter_file(method, model, column, is_input=False):
    """
    Get the filename where the PastML parameters are saved
    (for non-ML methods and input parameters will be None, as they have no parameters).
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param method: str, the ancestral state prediction method used by PASTML.
    :param model: str, the state evolution model used by PASTML.
    :param column: str, the column for which ancestral states are reconstructed with PASTML.
    :param is_input: bool, whether this is an input or an output file for PASTML.
    :return: str, filename or None for non-ML methods
    """
    ml = is_ml(method)
    if not ml and is_input:
        return None
    template = PASTML_PARAMS_CSV_INPUT if is_input else (PASTML_ML_PARAMS_CSV if ml else PASTML_MP_PARAMS_CSV)
    return template.format(state=col_name2cat(column), method=method, model=model)


def get_pastml_ancestral_state_file(method, model, column):
    """
    Get the filename where the PastML ancestral states are saved.
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param method: str, the ancestral state prediction method used by PASTML.
    :param model: str, the state evolution model used by PASTML.
    :param column: str, the column for which ancestral states are reconstructed with PASTML.
    :return: str, filename
    """
    template = PASTML_ANCESTRAL_STATES_CSV_ML if is_ml(method) else PASTML_ANCESTRAL_STATES_CSV_MP
    return template.format(state=col_name2cat(column), method=method, model=model)


def get_combined_ancestral_state_file(method, model, columns):
    """
    Get the filename where the combined ancestral states are saved (for one or several columns).
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param method: str, the ancestral state prediction method used by PASTML.
    :param model: str, the state evolution model used by PASTML.
    :param columns: list of str, the column(s) for which ancestral states are reconstructed with PASTML.
    :return: str, filename
    """
    template = CYTOPAST_ANCESTRAL_STATES_CSV_ML if is_ml(method) else CYTOPAST_ANCESTRAL_STATES_CSV_MP
    return template.format(states='_'.join(col_name2cat(column) for column in columns),
                           method=method, model=model)


def get_named_tree_file(tree):
    """
    Get the filename where the PastML tree (input tree but named and with collapsed zero branches) is saved.
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param tree: str, the input tree in newick format.
    :return: str, filename
    """
    return CYTOPAST_NAMED_TREE_NWK.format(tree=os.path.basename(tree))


def get_pastml_tip_state_file(column):
    """
    Get the filename where the PastML input tip states are saved.
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param column: str, the column for which ancestral states are reconstructed with PASTML.
    :return: str, filename
    """
    return PASTML_TIP_STATES_CSV.format(state=col_name2cat(column))


def get_pastml_marginal_prob_file(method, model, column):
    """
    Get the filename where the PastML marginal probabilities of node states are saved (will be None for non-marginal methods).
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param method: str, the ancestral state prediction method used by PASTML.
    :param model: str, the state evolution model used by PASTML.
    :param column: str, the column for which ancestral states are reconstructed with PASTML.
    :return: str, filename or None if the method is not marginal.
    """
    if not is_marginal(method):
        return None
    return PASTML_MARGINAL_PROBS_CSV.format(state=col_name2cat(column), model=model)


def parse_parameters(params, unique_states, param_file):
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
            param_df.to_csv(param_file)
    elif isinstance(params, str):
        if not os.path.exists(params):
            raise ValueError('You have specified some parameters ({}) but such a file does not exist!'
                             .format(params))
        param_file = params
    else:
        raise ValueError('Parameters must be specified either as a dict or as a path to a csv file, not as {}!'
                         .format(type(params)))
    return param_file
