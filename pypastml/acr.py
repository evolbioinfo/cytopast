import logging
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

from cytopast import collapse_zero_branches
from pypastml import rescale_tree, value2list
from pypastml.annotation import preannotate_tree, get_tree_stats
from pypastml.ml import MPPA, ml_acr, F81, is_ml
from pypastml.parsimony import parsimonious_acr
from pypastml.visualization import visualize


def reconstruct_ancestral_states(tree, feature, states, avg_br_len, prediction_method=MPPA, model=None):
    logging.info('ACR settings for {}:\n\tMethod:\t{}{}.\n'.format(feature, prediction_method,
                                                                   '\n\tModel:\t{}'.format(model) if model else ''))
    if is_ml(prediction_method):
        return ml_acr(tree, feature, prediction_method, model, states, avg_br_len)

    return parsimonious_acr(tree, feature, prediction_method, states)


def acr(tree, df, prediction_method=MPPA, model=F81, html=None, html_compressed=None):
    collapse_zero_branches(tree)
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
