import logging
import os
from collections import defaultdict, Counter
from pastml import F81, JC, ACCTRAN, JOINT, MARGINAL_APPROXIMATION
from random import sample

import pandas as pd
import numpy as np
from ete3 import Tree

from cytopast import TIPS_INSIDE, name_tree
from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, '4', 'phyml_tree.rooted.collapsed_support_0.5.collapsed_dist_0.nwk')
TREE_NWK_NAMED = os.path.join(DATA_DIR, '4', 'phyml_tree.rooted.collapsed_support_0.5.collapsed_dist_0.named.nwk')
TREE_NWK_DATE = os.path.join(DATA_DIR, '4', 'phyml_tree.rooted.collapsed_support_0.5.collapsed_dist_0_%d.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data_loc.tab')
DATES_TAB = os.path.join(DATA_DIR, 'dates.tab')


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def remove_certain_leaves(tr, to_remove=lambda node: False):
    """
    Removes all the branches leading to naive leaves from the given tree.
    :param tr: the tree of interest (ete3 Tree)
    [(state_1, 0), (state_2, time_of_transition_from_state_1_to_2), ...]. Branch removals will be added as '!'.
    :param to_remove: a method to check is a leaf should be removed.
    :return: the tree with naive branches removed (ete3 Tree) or None is all the leaves were naive in the initial tree.
    """

    def _merge_node_with_its_child(nd, child=None):
        if not child:
            child = nd.children[0]
        child.dist += nd.dist
        if nd.is_root():
            child.up = None
        else:
            parent = nd.up
            parent.remove_child(nd)
            parent.add_child(child)
        return child

    for node in tr.traverse("postorder"):
        # If this node has only one child branch
        # it means that the other child branch used to lead to a naive leaf and was removed.
        # We can merge this node with its child
        # (the child was already processed and either is a leaf or has 2 children).
        if len(node.children) == 1:
            merged_node = _merge_node_with_its_child(node)
            if merged_node.is_root():
                tr = merged_node
                break
        elif node.is_leaf() and to_remove(node):
            if node.is_root():
                return None
            node.up.remove_child(node)
    return tr


def parameter_bootstrap(year2nwk, mutations):
    for mutation in mutations:
        year2params = defaultdict(list)
        for year, nwk in year2nwk.items():
            print("Mutation {}, year {}".format(mutation, year))

            for rep in range(10):
                tree = Tree(nwk, format=0)
                random_leaves = set(sample(tree.get_leaves(), int(len(tree) * .8)))
                tree = remove_certain_leaves(tree, to_remove=lambda _: _ not in random_leaves)
                tree_nwk = '{}_{}.nwk'.format(nwk, rep)
                tree.write(outfile=tree_nwk, format=0)

                wd = os.path.join(DATA_DIR, 'pastml_{}'.format(rep))
                pastml_pipeline(data=STATES_INPUT, tree=tree_nwk,
                                # html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}_{}.html'.format(mutation, year)),
                                model=model, verbose=False, columns=[mutation], prediction_method=prediction_method,
                                work_dir=wd)
                mutation_cat = mutation.replace(':', '')
                year2params[year].append(os.path.join(wd, mutation_cat, model, prediction_method,
                                                      'Result.tree_{}.pastml.category_{}.params.csv'
                                                      .format(os.path.split(tree_nwk)[1], mutation_cat)))

        for year, params in year2params.items():
            print('\nYEAR {}======================='.format(year))
            freq_true = []
            for par in params:
                df = pd.read_csv(par, index_col=0, header=0)
                print(df)
                freq_true.append(df.loc['True', 'value'])
            print('\n{}: {}\n'.format(freq_true, np.mean(freq_true)))


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    tree = Tree(TREE_NWK, format=0)

    dates_df = pd.read_table(DATES_TAB, skiprows=[0], header=None, index_col=0)
    dates_df.columns = ['Date']
    dates_df.index = dates_df.index.map(str)
    dates_df = dates_df[dates_df.Date.apply(isfloat)]
    dates_df['Date'] = dates_df['Date'].astype(np.float)

    for n in tree:
        if n.name in dates_df.index:
            n.add_feature('date', dates_df.loc[n.name, 'Date'])

    name_tree(tree)
    tree.write(outfile=TREE_NWK_NAMED, format=3)
    min_date, max_date = min(_.date for _ in tree), max(_.date for _ in tree)
    logging.info("The tips are sampled between {} and {}".format(min_date, max_date))

    df = pd.read_table(STATES_INPUT, header=0, index_col=0)
    df.index = df.index.map(str)
    df = df[df.index.isin({_.name for _ in tree})]
    mutations = sorted([_ for _ in df.columns if _.startswith('RT:') or _.startswith('PR:')],
                       key=lambda drm: len(df[df[drm].isnull()]))[:5]
    logging.info('5 top mutations are {}'.format(mutations))
    # year2nwk = {}
    # for year in np.arange(np.ceil(max_date), np.ceil(min_date), step=-1):
    #     year = int(year)
    #     tree = remove_certain_leaves(tree, to_remove=lambda _: not hasattr(_, 'date') or getattr(_, 'date') > year)
    #     tree_year = TREE_NWK_DATE % year
    #     tree.write(outfile=tree_year, format=0)
    #     year2nwk[year] = tree_year
    #     logging.info('Found {} tips sampled before year {}'.format(len(tree), year))

    # parameter_bootstrap(year2nwk, mutations)

    model = F81
    prediction_method = MARGINAL_APPROXIMATION

    # for year, nwk in year2nwk.items():
    #     pastml_pipeline(data=STATES_INPUT, tree=nwk,
    #                     html_compressed=os.path.join(DATA_DIR, 'maps',
    #                                                  'map_Location_{}.html'.format(year)),
    #                     model=model, verbose=False, columns=['Location'], prediction_method=prediction_method,
    #                     work_dir=os.path.join(DATA_DIR, 'pastml'))

    mutation = 'RT:D67N'
    # exit()

    # tree = pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
    #                        html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}.html'.format(mutation)),
    #                        model=model, verbose=False, columns=[mutation], prediction_method=prediction_method,
    #                        work_dir=os.path.join(DATA_DIR, 'pastml'))
    # for _ in tree.traverse():
    #     tips = getattr(_, TIPS_INSIDE, [])
    #     print(_.state, tips)
    #     if tips and isinstance(tips[0], str):
    #         print(Counter(df.loc[tips, 'Location']))
    #     else:
    #         print([Counter(df.loc[_, 'Location']) for _ in tips])
    #
    #     print()

    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK_NAMED,
                    html_compressed=os.path.join(DATA_DIR, 'maps',
                                                 'map_Loc_{}.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation, 'Loc'],
                    prediction_method=prediction_method,
                    name_column='Loc',
                    work_dir=os.path.join(DATA_DIR, 'pastml'))
    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK_NAMED,
                    html_compressed=os.path.join(DATA_DIR, 'maps',
                                                 'map_Loc.html'),
                    model=model, verbose=False, columns=['Loc'],
                    prediction_method=prediction_method,
                    work_dir=os.path.join(DATA_DIR, 'pastml'))
    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK_NAMED,
                    html_compressed=os.path.join(DATA_DIR, 'maps',
                                                 'map_{}.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation],
                    prediction_method=prediction_method,
                    work_dir=os.path.join(DATA_DIR, 'pastml'))

    exit()

    # for mutation in mutations:
    #     for year, nwk in year2nwk.items():
    #         print("Mutation {}, year {}".format(mutation, year))
    #         pastml_pipeline(data=STATES_INPUT, tree=nwk,
    #                         html_compressed=os.path.join(DATA_DIR, 'maps',
    #                                                      'map_Location_{}_{}.html'.format(mutation, year)),
    #                         model=model, verbose=False, columns=[mutation, 'Location'],
    #                         prediction_method=prediction_method,
    #                         name_column='Location',
    #                         work_dir=os.path.join(DATA_DIR, 'pastml'))
