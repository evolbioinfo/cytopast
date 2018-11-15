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


def get_pastml_parameter_file(method, model, column):
    """
    Get the filename where the PastML parameters are saved
    (for non-ML methods and input parameters will be None, as they have no parameters).
    This file is inside the work_dir that can be specified for the pastml_pipeline method.
    :param method: str, the ancestral state prediction method used by PASTML.
    :param model: str, the state evolution model used by PASTML.
    :param column: str, the column for which ancestral states are reconstructed with PASTML.
    :return: str, filename or None for non-ML methods
    """
    ml = is_ml(method)
    template = PASTML_ML_PARAMS_CSV if ml else PASTML_MP_PARAMS_CSV
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
