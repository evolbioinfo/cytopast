import logging
import os

import pandas as pd
from ete3 import Tree
from pastml import F81, MARGINAL_APPROXIMATION
from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, '4', 'phyml_tree.rooted.collapsed_support_0.5.collapsed_dist_0.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data.tab')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    tree = Tree(TREE_NWK, format=0)

    df = pd.read_table(STATES_INPUT, header=0, index_col=0)
    df.index = df.index.map(str)
    df = df[df.index.isin({_.name for _ in tree})]
    mutations = sorted([_ for _ in df.columns if _.startswith('RT:') or _.startswith('PR:')],
                       key=lambda drm: len(df[df[drm].isnull()]))[:10]
    logging.info('10 top mutations are {}'.format(mutations))

    model = F81
    prediction_method = MARGINAL_APPROXIMATION
    mutation = 'RT:D67N'

    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                    html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation], prediction_method=prediction_method,
                    work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')

    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                    html_compressed=os.path.join(DATA_DIR, 'maps',
                                                 'map_Location_{}.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation, 'Location'],
                    prediction_method=prediction_method,
                    name_column='Location',
                    work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')

    # pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
    #                 html_compressed=os.path.join(DATA_DIR, 'maps',
    #                                              'map_Location.html'),
    #                 model=model, verbose=False, columns=['Location'],
    #                 prediction_method=prediction_method,
    #                 work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')
    #
    # mutations = ['RT:K103N', 'RT:D67N']
    # pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
    #                 html_compressed=os.path.join(DATA_DIR, 'maps',
    #                                              'map_{}.html'.format("-".join(mutations))),
    #                 model=model, verbose=False, columns=mutations,
    #                 prediction_method=prediction_method,
    #                 work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')
