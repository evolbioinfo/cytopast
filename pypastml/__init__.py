import logging
from collections import defaultdict
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from cytopast import col_name2cat, REASONABLE_NUMBER_OF_TIPS, compress_tree, DATE, collapse_zero_branches
from cytopast.colour_generator import get_enough_colours, WHITE
from cytopast.cytoscape_manager import save_as_cytoscape_html

BU_LH = 'BOTTOM_UP_LIKELIHOOD'
TD_LH = 'TOP_DOWN_LIKELIHOOD'
LH = 'LIKELIHOOD'
LH_SF = 'LIKELIHOOD_SF'
BU_LH_SF = 'BOTTOM_UP_LIKELIHOOD_SF'
BU_LH_JOINT_STATES = 'BOTTOM_UP_LIKELIHOOD_JOINT_STATES'
TD_LH_SF = 'TOP_DOWM_LIKELIHOOD_SF'


JOINT = 'joint'
MPPA = 'MPPA'
MAP = 'MAP'

JC = 'JC'
F81 = 'F81'
EFT = 'EFT'

DOWNPASS = 'downpass'
ACCTRAN = 'AccTran'
DELTRAN = 'DelTran'


def is_marginal(method):
    """
    Checks if the method is marginal, i.e. is either marginal itself, or MAP, or MPPA.
    :param method: str, the ancestral state prediction method used by PASTML.
    :return: bool
    """
    return method in {MPPA, MAP}


def is_ml(method):
    """
    Checks if the method is max likelihood, i.e. is either joint or one of the marginal ones
    (marginal itself, or MAP, or MPPA).
    :param method: str, the ancestral state prediction method used by PASTML.
    :return: bool
    """
    return method == JOINT or is_marginal(method)


MIN_VALUE = np.log10(np.finfo(np.float64).eps)
MAX_VALUE = np.log10(np.finfo(np.float64).max)


def get_mu(frequencies):
    """
    Calculates the mutation rate for F81 (and JC that is a simplification of it),
    as \mu = 1 / (1 - sum_i \pi_i^2). This way the overall rate of mutation -\mu trace(\Pi Q) is 1.
    See [Gascuel "Mathematics of Evolution and Phylogeny" 2005] for further details.

    :param frequencies: numpy array of state frequencies \pi_i
    :return: mutation rate \mu = 1 / (1 - sum_i \pi_i^2)
    """
    return 1. / (1. - frequencies.dot(frequencies))


def get_pij(frequencies, mu, t, sf):
    """
    Calculate the probability of substitution i->j over time t, given the mutation rate mu:
    For K81 (and JC which is a simpler version of it)
    Pij(t) = \pi_j (1 - exp(-mu t)) + exp(-mu t), if i == j, \pi_j (1 - exp(-mu t)), otherwise
    [Gascuel "Mathematics of Evolution and Phylogeny" 2005].

    :param frequencies: numpy array of state frequencies \pi_i
    :param mu: float, mutation rate: \mu = 1 / (1 - sum_i \pi_i^2)
    :param t: float, time t
    :param sf: float, scaling factor by which t should be multiplied.
    :return: numpy matrix Pij(t) = \pi_j (1 - exp(-mu t)) + exp(-mu t), if i == j, \pi_j (1 - exp(-mu t)), otherwise
    """
    # if mu == inf (e.g. just one state) and t == 0, we should prioritise mu
    exp_mu_t = 0. if (mu == np.inf) else np.exp(-mu * t * sf)
    return (1 - exp_mu_t) * frequencies + np.eye(len(frequencies)) * exp_mu_t


def get_bottom_up_likelihood(tree, feature, frequencies, sf, is_marginal=True):
    """
    Calculates the bottom-up likelihood for the given tree.
    The likelihood for each node is stored in the corresponding feature,
    given by get_personalised_feature_name(feature, BU_LH).

    :param is_marginal: bool, whether the likelihood reconstruction is marginal (true) or joint (false)
    :param tree: ete3.Tree tree
    :param feature: str, character for which the likelihood is calculated
    :param frequencies: numpy array of state frequencies \pi_i
    :param sf: float, scaling factor
    :return: float, the log likelihood
    """
    lh_sf_feature = get_personalized_feature_name(feature, BU_LH_SF)
    lh_feature = get_personalized_feature_name(feature, BU_LH)
    lh_joint_state_feature = get_personalized_feature_name(feature, BU_LH_JOINT_STATES)

    mu = get_mu(frequencies)
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            node.add_feature(lh_sf_feature, 0)
            continue

        likelihood_array = np.ones(len(frequencies), dtype=np.float64)

        for child in node.children:
            child_likelihoods = get_pij(frequencies, mu, child.dist, sf) * getattr(child, lh_feature)

            if is_marginal:
                likelihood_array *= child_likelihoods.sum(axis=1)
            else:
                likelihood_array *= child_likelihoods.max(axis=1)
                child_states = child_likelihoods.argmax(axis=1)
                child.add_feature(lh_joint_state_feature, child_states)

        if np.all(likelihood_array == 0):
            return -np.inf

        factors = rescale(likelihood_array, node, going_up=True)
        node.add_feature(lh_feature, likelihood_array)
        node.add_feature(lh_sf_feature, factors + sum(getattr(_, lh_sf_feature) for _ in node.children))
    return np.log(getattr(tree, lh_feature).dot(frequencies)) - getattr(tree, lh_sf_feature) * np.log(10)


def rescale(likelihood_array, node, going_up=True):
    """
    Rescales the likelihood array if it gets too small/large, by multiplying it by a factor of 10.
    :param likelihood_array: numpy array containing the likelihood to be rescaled
    :param node: ete3.TreeNode who's likelihood is about to be rescaled.
    :param going_up: bool, whether we are going bottom-up (true) or top-down (false) in the tree with this likelihood calculation.
    :return: float, factor of 10 by which the likelihood array has been multiplies.
    """

    min_lh_value = np.log10(np.min(likelihood_array[np.nonzero(likelihood_array)]))
    max_lh_value = np.log10(np.max(likelihood_array[np.nonzero(likelihood_array)]))

    num_siblings = 2
    # if we go up, this likelihood will be multiplied by those of our siblings at the next step
    if going_up and not node.is_root():
        num_siblings = len(node.up.children)
    # if we go down, this likelihood will be multiplied by those of our children at the next step
    elif not going_up and not node.is_leaf():
        num_siblings = len(node.children)

    factors = 0

    if max_lh_value > MAX_VALUE / num_siblings:
        factors = MAX_VALUE / num_siblings - max_lh_value
        likelihood_array *= np.power(10, factors)
    elif min_lh_value < MIN_VALUE / num_siblings:
        factors = min(-min_lh_value, MAX_VALUE / num_siblings - max_lh_value)
        likelihood_array *= np.power(10, factors)
    return factors


def initialize_tip_bottom_up_likelihoods(tree, feature, state2index):
    """
    Initializes the bottom-up likelihood arrays for tips based on their states given by the feature.
    :param tree: ete3.Tree, tree for which the tip likelihoods are to be initialized
    :param feature: str, feature in which the tip states are stored (the value could be None for a missing state)
    :param state2index: dict, mapping state to its index in the likelihood array.
    :return: void, adds the get_personalised_feature_name(feature, BU_LH) feature to tree tips.
    """
    lh_feature = get_personalized_feature_name(feature, BU_LH)

    # tips likelihood arrays won't be modified so we might as well just share them
    all_ones = np.ones(len(state2index), np.float64)
    state2array = {}
    for state, index in state2index.items():
        likelihood_array = np.zeros(len(state2index), np.float64)
        likelihood_array[index] = 1.
        state2array[state] = likelihood_array
    state2array[None] = all_ones
    state2array[''] = all_ones

    for tip in tree:
        tip.add_feature(lh_feature, state2array[getattr(tip, feature, '')])


def optimize_likelihood_params(tree, feature, frequencies, sf, optimise_sf=True, optimise_frequencies=True):
    """
    Optimizes the likelihood parameters (state frequencies and scaling factor) for the given tree.
    :param tree: ete3.Tree, tree of interest
    :param feature: str, character for which the likelihood is optimised
    :param frequencies: numpy array of initial state frequencies
    :param sf: float, initial scaling factor
    :param optimise_sf: bool, whether the scaling factor needs to be optimised
    :param optimise_frequencies: bool, whether the state frequencies need to be optimised
    :return: tuple (frequencies, scaling_factor) with the optimized parameters
    """
    bounds = []
    if optimise_frequencies:
        bounds += [np.array([1e-6, 10e6], np.float64)] * (len(frequencies) - 1)
    if optimise_sf:
        bounds += [np.array([0.001, 10.])]
    bounds = np.array(bounds, np.float64)

    def get_freq_sf_from_params(ps):
        freqs = frequencies
        if optimise_frequencies:
            freqs = np.hstack((ps[: (len(frequencies) - 1)], [1.]))
            freqs /= freqs.sum()
        sf_val = ps[(len(frequencies) - 1) if optimise_frequencies else 0] if optimise_sf else sf
        return freqs, sf_val

    def get_v(ps):
        if np.any(pd.isnull(ps)):
            return np.nan
        freqs, sf_val = get_freq_sf_from_params(ps)
        res = get_bottom_up_likelihood(tree, feature, freqs, sf_val, True)
        # logging.info('{}\t{}\t->\t{}'.format(frequencies if optimise_frequencies else '',
        #                                      sf_val if optimise_sf else '', res))
        return np.inf if pd.isnull(res) else -res

    params = None
    optimum = None
    for i in range(5):
        if i == 0:
            vs = np.hstack((frequencies[:-1] / frequencies[-1] if optimise_frequencies else [],
                            [sf] if optimise_sf else []))
        elif i == 1 and optimise_frequencies:
            vs = np.hstack((np.ones(len(frequencies) - 1, np.float64),
                            [sf] if optimise_sf else []))
        else:
            vs = np.random.uniform(bounds[:, 0], bounds[:, 1])
        fres = minimize(get_v, x0=vs, method='L-BFGS-B', bounds=bounds)
        if fres.success and not np.any(np.isnan(fres.x)):
            if optimum is None or fres.fun < optimum:
                params = fres.x
                optimum = fres.fun
                # logging.info('Found a new optimum: {}'.format(-fres.fun))
    if params is None:
        return None, None
    return get_freq_sf_from_params(params)


def rescale_tree(tree, sf):
    """
    Multiplies all the tree branches by the given scaling factor.

    :param tree: ete3.Tree, the tree to be rescaled
    :param sf: float, scaling factor
    :return: void, the initial tree is modified
    """
    for n in tree.traverse():
        n.dist *= sf


def calculate_top_down_likelihood(tree, feature, frequencies, sf):
    """
    Calculates the top-down likelihood for the given tree.
    The likelihood for each node is stored in the corresponding feature,
    given by get_personalised_feature_name(feature, TD_LH).

    To calculate the top-down likelihood of a node, we assume that the tree is rooted in this node
    and combine the likelihoods of the “up-subtrees”,
    e.g. to calculate the top-down likelihood of a node N1 being in a state i,
    given that its parent node is P and its brother node is N2, we imagine that the tree is re-rooted in N1,
    therefore P becoming the child of N1, and N2 its grandchild.
    We then calculate the bottom-up likelihood from the P subtree:
    L_top_down(N1, i) = \sum_j P(i -> j, dist(N1, P)) * L_top_down(P) * \sum_k P(j -> k, dist(N2, P)) * L_bottom_up (N2).

    For the root node we assume its top-down likelihood to be 1 for all the states.

    :param tree: ete3.Tree, the tree of interest (with bottom-up likelihood precalculated)
    :param frequencies: numpy array of state frequencies
    :return: void, stores the node top-down likelihoods in the get_personalised_feature_name(feature, TD_LH) feature.
    """

    lh_feature = get_personalized_feature_name(feature, TD_LH)
    lh_sf_feature = get_personalized_feature_name(feature, TD_LH_SF)
    bu_lh_feature = get_personalized_feature_name(feature, BU_LH)
    bu_lh_sf_feature = get_personalized_feature_name(feature, BU_LH_SF)

    mu = get_mu(frequencies)

    for node in tree.traverse('preorder'):
        if node.is_root():
            node.add_feature(lh_feature, np.ones(len(frequencies), np.float64))
            node.add_feature(lh_sf_feature, 0)
            continue

        parent = node.up
        parent_bu_likelihood = getattr(parent, bu_lh_feature)

        node_pjis = np.transpose(get_pij(frequencies, mu, node.dist, sf))
        node_contribution = getattr(node, bu_lh_feature).dot(node_pjis)

        parent_likelihood = getattr(parent, lh_feature) * parent_bu_likelihood
        parent_likelihood[np.nonzero(parent_likelihood)] /= node_contribution[np.nonzero(parent_likelihood)]
        factors = getattr(parent, lh_sf_feature) + getattr(parent, bu_lh_sf_feature) - getattr(node, bu_lh_sf_feature)

        td_likelihood = parent_likelihood.dot(node_pjis)
        factors += rescale(td_likelihood, node, going_up=False)
        
        node.add_feature(lh_feature, td_likelihood)
        node.add_feature(lh_sf_feature, factors)


def alter_zero_tip_likelihoods(tree, feature):
    """
    Alters the bottom-up likelihood arrays for zero-distance tips
    to make sure they do not contradict with other zero-distance tip siblings.

    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the likelihood is altered
    :return: void, modifies the get_personalised_feature_name(feature, BU_LH) feature to zero-distance tips.
    """
    zero_parent2tips = defaultdict(list)

    lh_feature = get_personalized_feature_name(feature, BU_LH)

    for tip in tree:
        if tip.dist == 0:
            state = getattr(tip, feature, None)
            if state is not None and state != '':
                zero_parent2tips[tip.up].append(tip)

    # adjust zero tips to contain all the zero tip options as states
    for parent, zero_tips in zero_parent2tips.items():
        likelihood_array = None
        for tip in zero_tips:
            if likelihood_array is None:
                likelihood_array = getattr(tip, lh_feature).copy()
            else:
                tip_likelihood_array = getattr(tip, lh_feature)
                likelihood_array[np.nonzero(tip_likelihood_array)] = 1.
            tip.add_feature(lh_feature, likelihood_array)


def unalter_zero_tip_likelihoods(tree, feature, state2index):
    """
    Unalters the bottom-up likelihood arrays for zero-distance tips
    to contain ones only in their states.

    :param state2index: dict, mapping between states and their indices in the likelihood array
    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the likelihood was altered
    :return: void, modifies the get_personalised_feature_name(feature, BU_LH) feature to zero-distance tips.
    """
    lh_feature = get_personalized_feature_name(feature, BU_LH)
    for tip in tree:
        if tip.dist > 0:
            continue
        state = getattr(tip, feature, None)
        if state is not None and state != '':
            likelihood_array = np.zeros(len(state2index), np.float64)
            likelihood_array[state2index[state]] = 1.
            tip.add_feature(lh_feature, likelihood_array)


def calculate_marginal_likelihoods(tree, feature, frequencies):
    """
    Calculates marginal likelihoods for each tree node
    by multiplying state frequencies with their bottom-up and top-down likelihoods.

    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the likelihood is calculated
    :param frequencies: numpy array of state frequencies
    :return: void, stores the node marginal likelihoods in the get_personalised_feature_name(feature, LH) feature.
    """
    bu_lh_feature = get_personalized_feature_name(feature, BU_LH)
    bu_lh_sf_feature = get_personalized_feature_name(feature, BU_LH_SF)
    td_lh_feature = get_personalized_feature_name(feature, TD_LH)
    td_lh_sf_feature = get_personalized_feature_name(feature, TD_LH_SF)
    lh_feature = get_personalized_feature_name(feature, LH)
    lh_sf_feature = get_personalized_feature_name(feature, LH_SF)

    for node in tree.traverse('preorder'):
        likelihood = getattr(node, bu_lh_feature) * getattr(node, td_lh_feature) * frequencies
        node.add_feature(lh_feature, likelihood)
        node.add_feature(lh_sf_feature, getattr(node, td_lh_sf_feature) + getattr(node, bu_lh_sf_feature))

        node.del_feature(bu_lh_feature)
        node.del_feature(bu_lh_sf_feature)
        node.del_feature(td_lh_feature)
        node.del_feature(td_lh_sf_feature)


def check_marginal_likelihoods(tree, feature):
    """
    Tests if marginal likelihoods were correctly calculated
    by comparing the likelihoods of all the nodes (should be all the same).

    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the likelihood is to be checked
    :return: void, assertion fails if there is a problem with the likelihood.
    """
    lh_feature = get_personalized_feature_name(feature, LH)
    lh_sf_feature = get_personalized_feature_name(feature, LH_SF)

    for node in tree.traverse('preorder'):
        if not node.is_root() and not (node.is_leaf() and node.dist == 0):
            node_loglh = np.log10(getattr(node, lh_feature).sum()) - getattr(node, lh_sf_feature)
            parent_loglh = np.log10(getattr(node.up, lh_feature).sum()) - getattr(node.up, lh_sf_feature)
            assert np.round(node_loglh, 2) == np.round(parent_loglh, 2)


def convert_likelihoods_to_probabilities(tree, feature):
    """
    Normalizes each node maginal likelihoods to convert them to marginal probabilities.

    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the probabilities are calculated
    :return: void, modifies the node get_personalised_feature_name(feature, LH) feature to store the probabilities.
    """
    lh_feature = get_personalized_feature_name(feature, LH)
    for node in tree.traverse():
        lh = getattr(node, lh_feature)
        lh /= lh.sum()


def get_personalized_feature_name(character, feature):
    """
    Precedes the feature name by the character name
    (useful when likelihoods for different characters are calculated in parallel).

    :param character: str, character name
    :param feature: str, feature to be personalized
    :return: str, the personalized feature
    """
    return '{}_{}'.format(character, feature)


def choose_ancestral_states_mppa(tree, feature, states):
    """
    Chooses node ancestral states based on their marginal probabilities using MPPA method.

    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the ancestral states are to be chosen
    :param states: numpy.array of possible character states in order corresponding to the probabilities array
    :return: void, add ancestral states as the `feature` feature to each node
    (as a list if multiple states are possible or as a string if only one state is chosen)
    """
    lh_feature = get_personalized_feature_name(feature, LH)
    n = len(states)
    for node in tree.traverse():
        likelihood = getattr(node, lh_feature)
        sorted_likelihood = np.sort(likelihood)
        best_k = n
        best_correstion = np.inf
        for k in range(1, n):
            correction = np.hstack((np.zeros(n - k), np.ones(k) / k)) - sorted_likelihood
            correction = correction.dot(correction)
            if correction < best_correstion:
                best_correstion = correction
                best_k = k

        possible_states = states[sorted(range(n), key=lambda _: -likelihood[_])[:best_k]].tolist()
        node.add_feature(feature, possible_states[0] if best_k == 1 else possible_states)


def choose_ancestral_states_map(tree, feature, states):
    """
    Chooses node ancestral states based on their marginal probabilities using MAP method.

    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the ancestral states are to be chosen
    :param states: numpy.array of possible character states in order corresponding to the probabilities array
    :return: void, add ancestral states as the `feature` feature to each node
    """
    lh_feature = get_personalized_feature_name(feature, LH)
    for node in tree.traverse():
        likelihood = getattr(node, lh_feature)
        best_state = states[likelihood.argmax()]
        node.add_feature(feature, best_state)


def choose_ancestral_states_joint(tree, feature, states, frequencies):
    """
    Chooses node ancestral states based on their marginal probabilities using joint method.

    :param frequencies: numpy array of state frequencies
    :param tree: ete3.Tree, the tree of interest
    :param feature: str, character for which the ancestral states are to be chosen
    :param states: numpy.array of possible character states in order corresponding to the probabilities array
    :return: void, add ancestral states as the `feature` feature to each node
    """
    lh_feature = get_personalized_feature_name(feature, BU_LH)
    lh_state_feature = get_personalized_feature_name(feature, BU_LH_JOINT_STATES)

    def chose_consistent_state(node, state_index):
        best_state = states[state_index]
        node.add_feature(feature, best_state)
        for child in node.children:
            chose_consistent_state(child, getattr(child, lh_state_feature)[state_index])

    chose_consistent_state(tree, (getattr(tree, lh_feature) * frequencies).argmax())


def visualize(tree, columns, name_column=None, html=None, html_compressed=None,
              tip_size_threshold=REASONABLE_NUMBER_OF_TIPS, min_date=0, max_date=0):
    one_column = len(columns) == 1

    column2values = {}
    for feature in columns:
        column2values[feature] = annotate(tree, feature, unique=one_column)

    if not name_column and one_column:
        name_column = columns[0]

    name2colour = {}
    for cat in columns:
        unique_values = column2values[cat]
        num_unique_values = len(unique_values)
        colours = get_enough_colours(num_unique_values)
        for value, col in zip(unique_values, colours):
            name2colour['{}_{}'.format(value, True) if one_column else '{}_{}'.format(cat, value)] = col
        logging.info('Mapped states to colours for {} as following: {} -> {}'.format(cat, unique_values, colours))
        # let ambiguous values be white
        if not one_column:
            name2colour['{}_'.format(cat)] = WHITE

    # set internal node dates to min of its tips' dates
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            if not hasattr(n, DATE):
                n.add_feature(DATE, 0)
        else:
            n.add_feature(DATE, min(getattr(_, DATE) for _ in n))

    def get_category_str(n):
        if one_column:
            return '{}: {}'.format(columns[0], ' or '.join('{}'.format(_)
                                                           for _ in column2values[columns[0]]
                                                           if hasattr(n, _) and getattr(n, _, '') != ''))
        return '<br>'.join('{}: {}'.format(_, getattr(n, _))
                           for _ in columns if hasattr(n, _) and getattr(n, _, '') != '')

    if html:
        save_as_cytoscape_html(tree, html, categories=column2values[columns[0]] if one_column else columns,
                               name2colour=name2colour,
                               n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               name_feature='name', min_date=min_date, max_date=max_date, is_compressed=False)

    if html_compressed:
        tree = compress_tree(tree, categories=column2values[columns[0]] if one_column else columns,
                             tip_size_threshold=tip_size_threshold)
        save_as_cytoscape_html(tree, html_compressed, categories=column2values[columns[0]] if one_column else columns,
                               name2colour=name2colour, n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               min_date=min_date, max_date=max_date, name_feature=name_column, is_compressed=True)
    return tree


def annotate(tree, feature, unique=True):
    all_states = set()
    for node in tree.traverse():
        possible_states = getattr(node, feature)
        if isinstance(possible_states, list):
            node.add_feature(feature, '')
        else:
            all_states.add(possible_states)
        if unique:
            for state in (possible_states if isinstance(possible_states, list) else [possible_states]):
                node.add_feature(state, True)
                all_states.add(state)
    return sorted(all_states)


def reconstruct_ancestral_states(tree, feature, states, avg_br_len, prediction_method=MPPA):
    state2index = dict(zip(states, range(len(states))))

    frequencies = np.zeros(len(state2index), np.float64)
    missing_frequency = 0.
    for _ in tree:
        state = getattr(_, feature, None)
        if state is not None and state != '':
            frequencies[state2index[state]] += 1
        else:
            missing_frequency += 1
    total_count = frequencies.sum() + missing_frequency
    frequencies /= total_count

    initialize_tip_bottom_up_likelihoods(tree, feature, state2index)
    alter_zero_tip_likelihoods(tree, feature)
    sf = 1.
    likelihood = get_bottom_up_likelihood(tree, feature, frequencies, sf, True)
    logging.info('Initial {} values:{}{}{}.\n'
                 .format(feature,
                         ''.join('\n\tfrequency of {}:\t{:.3f}'.format(state, frequencies[state2index[state]])
                                 for state in states),
                         '\n\tfrequency of missing data:\t{:.3f}'.format(missing_frequency / total_count)
                                          if missing_frequency else '',
                         '\n\tlog likelihood:\t{:.3f}'.format(likelihood)))

    frequencies, sf = optimize_likelihood_params(tree, feature, frequencies, sf, optimise_frequencies=False, optimise_sf=True)
    likelihood = get_bottom_up_likelihood(tree, feature, frequencies, sf, True)
    logging.info('Optimised SF for {}:\n'
                 '\tSF:\t{:.3f}, i.e. {:.3f} changes per avg branch\n'
                 '\tlog likelihood:\t{:.3f}.\n'
                 .format(feature, sf / avg_br_len, sf, likelihood))

    frequencies, sf = optimize_likelihood_params(tree, feature, frequencies, sf, optimise_frequencies=True, optimise_sf=False)
    likelihood = get_bottom_up_likelihood(tree, feature, frequencies, sf, True)
    logging.info('Optimised frequencies for {}:{}\n'
                 '\tlog likelihood:\t{:.3f}.\n'
                 .format(feature,
                         ''.join('\n\t{}:\t{:.3f}'.format(state, frequencies[state2index[state]]) for state in states),
                         likelihood))

    if is_marginal(prediction_method):
        calculate_top_down_likelihood(tree, feature, frequencies, sf)
        unalter_zero_tip_likelihoods(tree, feature, state2index)
        calculate_marginal_likelihoods(tree, feature, frequencies)
        check_marginal_likelihoods(tree, feature)
        convert_likelihoods_to_probabilities(tree, feature)
        if MPPA == prediction_method:
            logging.info('Choosing MPPA ancestral states for {}.\n'.format(feature))
            choose_ancestral_states_mppa(tree, feature, states)
        elif MAP == prediction_method:
            logging.info('Choosing MAP ancestral states for {}.\n'.format(feature))
            choose_ancestral_states_map(tree, feature, states)
    # joint
    else:
        logging.info('Choosing joint ancestral states for {}.\n'.format(feature))
        get_bottom_up_likelihood(tree, feature, frequencies, sf, False)
        choose_ancestral_states_joint(tree, feature, states, frequencies)

    return likelihood


def acr(tree, df, prediction_method=MPPA, html=None, html_compressed=None):
    collapse_zero_branches(tree)
    columns = preannotate_tree(df, tree)

    avg_br_len = get_tree_stats(tree)

    logging.info('\n=============RECONSTRUCTING ANCESTRAL STATES=============\n')
    rescale_tree(tree, 1. / avg_br_len)

    def _work(args):
        reconstruct_ancestral_states(*args)

    # create a prediction_method variable for all the columns by either
    # specifying the default value if nothing was chosen
    if prediction_method is None:
        prediction_method = MPPA
    # and propagating the chosen value to all columns
    if not isinstance(prediction_method, list):
        prediction_method = [prediction_method] * len(columns)
    # or making sure that the default value if chosen for the columns for which the method was not specified
    else:
        prediction_method += [MPPA] * (len(columns) - len(prediction_method))

    with ThreadPool() as pool:
        pool.map(func=_work,
                 iterable=((tree, column, np.sort([_ for _ in df[column].unique() if pd.notnull(_)]), avg_br_len, method)
                           for (column, method) in zip(columns, prediction_method)))

    if html or html_compressed:
        logging.info('\n=============VISUALIZING=============\n')
        visualize(tree, columns, html=html, html_compressed=html_compressed)


def get_tree_stats(tree):
    len_sum = 0
    num_zero_nodes = 0
    max_polynomy = 0
    max_len = 0
    num_tips = 0
    num_nodes = 0
    num_zero_tips = 0
    tip_len_sum = 0

    for node in tree.traverse():
        num_nodes += 1
        len_sum += node.dist
        max_polynomy = max(len(node.children), max_polynomy)
        max_len = max(max_len, node.dist)
        if not node.dist:
            num_zero_nodes += 1

        if node.is_leaf():
            num_tips += 1
            tip_len_sum += node.dist
            if not node.dist:
                num_zero_tips += 1

    avg_br_len = len_sum / (num_nodes - num_zero_nodes)
    logging.info('\n=============CHECKING THE TREE=============\n'
                 '\tnumber of tips:\t{}\n'
                 '\tnumber of zero-branch tips:\t{}\n'
                 '\tnumber of internal nodes:\t{}\n'
                 '\tmax number of children per node:\t{}\n'
                 '\tmax branch length:\t{:.5f}\n'
                 '\tavg non-zero branch length:\t{:.5f}\n'
                 .format(num_tips,
                         num_zero_tips,
                         num_nodes - num_tips,
                         max_polynomy,
                         max_len,
                         avg_br_len))
    return avg_br_len


def preannotate_tree(df, tree):
    df.columns = [col_name2cat(col) for col in df.columns]
    df.index = df.index.map(str)
    df.fillna('', inplace=True)
    for _ in tree:
        if _.name in df.index:
            _.add_features(**df.loc[_.name, :].to_dict())
    return df.columns


