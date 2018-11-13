import logging
from collections import namedtuple
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

from cytopast import collapse_zero_branches, name_tree
from pypastml import rescale_tree, value2list
from pypastml.annotation import preannotate_tree, get_tree_stats
from pypastml.ml import MPPA, ml_acr, F81, is_ml, MAP, JOINT
from pypastml.parsimony import parsimonious_acr, is_parsimonious, ACCTRAN, DELTRAN, DOWNPASS
from pypastml.visualization import visualize

COPY = 'copy'

ACRCopyResult = namedtuple('ACRCopyResult', field_names=['character', 'states', 'method'])


def reconstruct_ancestral_states(tree, feature, states, avg_br_len, prediction_method=MPPA, model=None):
    logging.info('ACR settings for {}:\n\tMethod:\t{}{}.\n'.format(feature, prediction_method,
                                                                   '\n\tModel:\t{}'.format(model)
                                                                   if model and is_ml(prediction_method) else ''))
    if is_ml(prediction_method):
        return ml_acr(tree, feature, prediction_method, model, states, avg_br_len)

    if is_parsimonious(prediction_method):
        return parsimonious_acr(tree, feature, prediction_method, states)

    if COPY == prediction_method:
        return ACRCopyResult(character=feature, states=states, method=prediction_method)

    raise ValueError('Method {} is unknown, should be one of ML ({}, {}, {}), one of MP ({}, {}, {}) or {}'
                     .format(prediction_method, MPPA, MAP, JOINT, ACCTRAN, DELTRAN, DOWNPASS, COPY))


def acr(tree, df, prediction_method=MPPA, model=F81, html=None, html_compressed=None):
    collapse_zero_branches(tree)
    name_tree(tree)
    columns = preannotate_tree(df, tree)

    avg_br_len = get_tree_stats(tree)

    logging.info('\n=============RECONSTRUCTING ANCESTRAL STATES=============\n')
    rescale_tree(tree, 1. / avg_br_len)

    def _work(args):
        return reconstruct_ancestral_states(*args)

    prediction_methods = value2list(len(columns), prediction_method, MPPA)
    models = value2list(len(columns), model, F81)

    with ThreadPool() as pool:
        acr_results = \
            pool.map(func=_work, iterable=((tree, column, np.sort([_ for _ in df[column].unique()
                                                                   if pd.notnull(_) and _ != '']),
                                            avg_br_len, method, model)
                                           for (column, method, model) in zip(columns, prediction_methods, models)))

    if html or html_compressed:
        logging.info('\n=============VISUALIZING=============\n')
        visualize(tree, columns, html=html, html_compressed=html_compressed)

    return acr_results
