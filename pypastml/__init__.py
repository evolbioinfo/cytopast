import logging

import numpy as np
import pandas as pd
from cytopast.cytoscape_manager import save_as_cytoscape_html

from cytopast.colour_generator import get_enough_colours
from scipy.optimize import minimize

from cytopast import read_tree, col_name2cat, REASONABLE_NUMBER_OF_TIPS, compress_tree

BU_LH = 'BOTTOM_UP_LIKELIHOOD'
TD_LH = 'TOP_DOWN_LIKELIHOOD'
LH = 'LIKELIHOOD'
BU_LH_SF = 'BOTTOM_UP_LIKELIHOOD_SF'
TD_LH_SF = 'TOP_DOWM_LIKELIHOOD_SF'

MIN_VALUE = np.log10(np.finfo(np.float64).eps)
MAX_VALUE = np.log10(np.finfo(np.float64).max)


def get_mu(frequencies):
    """
    Calculates the mutation rate for F81 (and JC that is a simplification of it),
    as \mu = 1 / (1 - sum_i \pi_i^2). This way the overall rate of mutation -\mu trace(\Pi Q) is 1.
    See [Gascuel "Mathematics of Evolution and Phylogeny" 2005] for further details.
    """
    return 1. / (1. - frequencies.dot(frequencies))


def get_pi(frequencies, mu, t, i, sf):
    """
    Calculate the probability of substitution i->j over time t, given the mutation rate mu:
    For K81 (and JC which is a simpler version of it)
    Pxy(t) = \pi_y (1 - exp(-mu t)) + exp(-mu t), if x == y, \pi_y (1 - exp(-mu t)), otherwise
    [Gascuel "Mathematics of Evolution and Phylogeny" 2005].
    """
    # if mu == inf (e.g. just one state) and t == 0, we should prioritise mu
    exp_mu_t = 0. if (mu == np.inf) else np.exp(-mu * t * sf)

    res = frequencies * (1. - exp_mu_t)
    res[i] += exp_mu_t
    return res


def get_bottom_up_likelihood(tree, frequencies, sf):
    mu = get_mu(frequencies)
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            node.add_feature(BU_LH_SF, 0)
            continue
        likelihood_array = getattr(node, BU_LH, np.zeros(len(frequencies), np.float64))
        for i in range(likelihood_array.shape[0]):
            likelihood_array[i] = np.prod([get_pi(frequencies, mu, child.dist, i, sf).dot(getattr(child, BU_LH))
                                           for child in node.children])

        if np.all(likelihood_array == 0):
            return -np.inf

        factors = rescale(likelihood_array, node)
        node.add_feature(BU_LH, likelihood_array)
        node.add_feature(BU_LH_SF, factors + sum(getattr(_, BU_LH_SF) for _ in node.children))
    return np.log(getattr(tree, BU_LH).dot(frequencies)) - getattr(tree, BU_LH_SF) * np.log(10)


def rescale(likelihood_array, node):
    factors = 0
    if not node.is_root():
        min_lh_value = np.log10(np.min(likelihood_array[np.nonzero(likelihood_array)]))
        max_lh_value = np.log10(np.max(likelihood_array[np.nonzero(likelihood_array)]))
        num_siblings = len(node.up.children)
        if min_lh_value < MIN_VALUE / num_siblings and max_lh_value < MAX_VALUE / num_siblings:
            factors = min(-min_lh_value, MAX_VALUE / num_siblings - max_lh_value)
            likelihood_array *= np.power(10, factors)
    return factors


def initialize_tip_states(tree, feature, state2index):
    zero_parents = set()

    for tip in tree:
        state = getattr(tip, feature, None)
        if state is not None and state != '':
            likelihood_array = np.zeros(len(state2index), np.float64)
            likelihood_array[state2index[state]] = 1.
        else:
            likelihood_array = np.ones(len(state2index), np.float64)
        tip.add_feature(BU_LH, likelihood_array)

        if tip.dist == 0 and state is not None and state != '':
            zero_parents.add(tip.up)

    # adjust zero tips
    for parent in zero_parents:
        states = set()
        zero_tips = []
        for tip in parent.children:
            if not tip.is_leaf() or tip.dist > 0.:
                continue
            state = getattr(tip, feature, None)
            if state is not None and state != '':
                zero_tips.append(tip)
                states.add(state)
        if len(states) > 1.:
            likelihood_array = np.zeros(len(state2index), np.float64)
            for state in states:
                likelihood_array[state2index[state]] = 1.
            for tip in zero_tips:
                tip.add_feature(BU_LH, likelihood_array)


def minimize_params(tree, frequencies, sf, optimise_sf=True, optimise_frequencies=True):
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
        res = get_bottom_up_likelihood(tree, freqs, sf_val)
        return np.inf if pd.isnull(res) else -res

    params = None
    optimum = None
    for i in range(5):
        res = np.inf
        while res == np.inf:
            if i == 0:
                vs = np.hstack((frequencies[:-1] / frequencies[-1] if optimise_frequencies else [],
                                [sf] if optimise_sf else []))
            elif i == 1 and optimise_frequencies:
                vs = np.hstack((np.ones(len(frequencies) - 1, np.float64),
                                [sf] if optimise_sf else []))
            else:
                vs = np.random.uniform(bounds[:, 0], bounds[:, 1])
            # print(vs)
            res = get_v(vs)
        fres = minimize(get_v, x0=vs, method='L-BFGS-B', bounds=bounds)
        if fres.success and not np.any(np.isnan(fres.x)):
            if optimum is None or fres.fun < optimum:
                params = fres.x
                optimum = fres.fun
                logging.info('Optimisation took {} iterations and finished due to {}.'
                             .format(fres.nit, str(fres.message)))
                if optimise_sf and not optimise_frequencies:
                    break
                else:
                    logging.info('Found a new optimum: {}'.format(-fres.fun))
    if params is None:
        return None, None
    return get_freq_sf_from_params(params)


def rescale_tree(tree, sf):
    for n in tree.traverse():
        n.dist *= sf


def calculate_top_down_likelihood(tree, frequencies, sf):
    """
    To calculate it we assume that the tree is rooted in this node and combine the likelihoods of the “up-subtrees”,
    e.g. to calculate the top-down likelihood of a node N1 being in a state i,
    given that its parent node is P and its brother node is N2, we imagine that the tree is re-rooted in N1,
    therefore P becoming the child of N1, and N2 its grandchild,
    and then calculate the bottom-up likelihood from the P subtree:
    L_top_down(N1, i) = \sum_j P(i -> j, dist(N1, P)) * L_top_down(P) * \sum_k P(j -> k, dist(N2, P)) * L_bottom_up (N2).

    For the root node we assume its top-down likelihood to be 1 for all the states.
    :param tree:
    :param frequencies:
    :return:
    """

    mu = get_mu(frequencies)

    for node in tree.traverse('preorder'):
        if node.is_root():
            node.add_feature(TD_LH, np.ones(len(frequencies), np.float64))
            node.add_feature(TD_LH_SF, 0)
            continue

        parent = node.up
        parent_bu_likelihood = getattr(parent, BU_LH)

        node_pjis = np.transpose(np.array([get_pi(frequencies, mu, node.dist, i, sf) for i in range(len(frequencies))]))
        node_contribution = np.transpose(getattr(node, BU_LH).dot(node_pjis))

        parent_likelihood = getattr(parent, TD_LH) * parent_bu_likelihood
        parent_likelihood[np.nonzero(parent_likelihood)] /= node_contribution[np.nonzero(parent_likelihood)]
        factors = getattr(parent, TD_LH_SF) + getattr(parent, BU_LH_SF) - getattr(node, BU_LH_SF)

        td_likelihood = np.transpose(parent_likelihood.dot(node_pjis))
        factors += rescale(td_likelihood, node)
        
        node.add_feature(TD_LH, td_likelihood)
        node.add_feature(TD_LH_SF, factors)


def unalter_problematic_tip_states(tree, feature, state2index):
    for tip in tree:
        if tip.dist > 0:
            continue
        state = getattr(tip, feature, None)
        if state is not None and state != '':
            likelihood_array = np.zeros(len(state2index), np.float64)
            likelihood_array[state2index[state]] = 1.
            tip.add_feature(BU_LH, likelihood_array)


def calculate_marginal_probabilities(tree, frequencies):
    for node in tree.traverse('preorder'):
        likelihood = getattr(node, BU_LH) * getattr(node, TD_LH) * frequencies
        node.add_feature(LH, likelihood)


def check_marginal_probabilities(tree):
    for node in tree.traverse('preorder'):
        if not node.is_root() and not (node.is_leaf() and node.dist == 0):
            node_loglh = np.log10(getattr(node, LH).sum()) - (getattr(node, TD_LH_SF) + getattr(node, BU_LH_SF))
            parent_loglh = np.log10(getattr(node.up, LH).sum()) \
                           - (getattr(node.up, TD_LH_SF) + getattr(node.up, BU_LH_SF))
            assert np.round(node_loglh, 2) == np.round(parent_loglh, 2)


def normalize_result_probabilities(tree):
    for node in tree.traverse():
        lh = getattr(node, LH)
        lh /= lh.sum()


def choose_ancestral_states(tree, n):
    for node in tree.traverse():
        likelihood = getattr(node, LH)
        sorted_likelihood = np.sort(likelihood)
        best_k = n
        best_correstion = np.inf
        for k in range(1, n):
            correction = np.hstack((np.zeros(n - k), np.ones(k))) - sorted_likelihood
            correction = correction.dot(correction)
            if correction < best_correstion:
                best_correstion = correction
                best_k = k
        zero_indices = sorted(range(n), key=lambda _: likelihood[_])[:-best_k]
        likelihood[zero_indices] = 0.
        likelihood[np.nonzero(likelihood)] = 1 / best_k


def visualize(tree, columns, categories, name_column=None, html=None, html_compressed=None,
              tip_size_threshold=REASONABLE_NUMBER_OF_TIPS, min_date=0, max_date=0):
    one_column = len(columns) == 1

    if not name_column and one_column:
        name_column = columns[0]

    if name_column:
        name_column = col_name2cat(name_column)

    name2colour = {}
    # if len(columns) > 1:
    #     for cat in columns:
    #         unique_values = df[cat].unique()
    #         unique_values = sorted(unique_values[~pd.isnull(unique_values)].astype(str))
    #         num_unique_values = len(unique_values)
    #         colours = get_enough_colours(num_unique_values)
    #         for value, col in zip(unique_values, colours):
    #             name2colour['{}_{}'.format(cat, value)] = col
    #         logging.info('Mapped the values to colours as following: {}, {}'.format(unique_values, colours))
    #         # let ambiguous values be white
    #         name2colour['{}_'.format(cat)] = WHITE
    # else:
    colours = get_enough_colours(len(categories))
    for cat, col in zip(categories, colours):
        name2colour['{}_{}'.format(cat, True)] = col
    logging.info('Mapped the values to colours as following: {}, {}'.format(categories, colours))

    # # set internal node dates to min of its tips' dates
    # for n in tree.traverse('postorder'):
    #     if n.is_leaf():
    #         if not hasattr(n, DATE):
    #             n.add_feature(DATE, 0)
    #     else:
    #         n.add_feature(DATE, min(getattr(_, DATE) for _ in n))

    def get_category_str(n):
        if one_column:
            return '{}: {}'.format([0], ' or '.join('{}'.format(_) for _ in categories if hasattr(n, _)
                                                           and getattr(n, _, '') != ''))
        return '<br>'.join('{}: {}'.format(_, getattr(n, _)) for _ in categories if hasattr(n, _)
                           and getattr(n, _, '') != '')

    if html:
        save_as_cytoscape_html(tree, html, categories=categories, name2colour=name2colour,
                               n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               name_feature='name', min_date=min_date, max_date=max_date, is_compressed=False)

    if html_compressed:
        tree = compress_tree(tree, categories=categories, tip_size_threshold=tip_size_threshold)
        save_as_cytoscape_html(tree, html_compressed, categories=categories,
                               name2colour=name2colour, n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               min_date=min_date, max_date=max_date, name_feature=name_column, is_compressed=True)
    return tree


def annotate(tree, state2index, feature, unique=True):
    category = col_name2cat(feature)
    states = np.array(sorted(state2index.keys(), key=lambda _: state2index[_]))
    for node in tree.traverse():
        likelihood = getattr(node, LH)
        possible_states = states[np.nonzero(likelihood)]
        for state in possible_states:
            node.add_feature(state, True)
        node.add_feature(category, possible_states[0] if len(possible_states) == 1 else '')


def reconstruct_ancestral_states(tree, feature, states, html=None, html_compressed=None):
    states = sorted(states)
    state2index = dict(zip(states, range(len(states))))
    initialize_tip_states(tree, feature, state2index)
    avg_br_len = np.mean([n.dist for n in tree.traverse() if n.dist])
    max_br_len = np.max([n.dist for n in tree.traverse()])
    logging.info('Max branch is {}'.format(max_br_len))
    logging.info('Avg branch is {}'.format(avg_br_len))

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
    rescale_tree(tree, 1. / avg_br_len)
    sf = 1.

    logging.info('Initial frequencies are:\n{}'
                 .format(''.join('\n\t{}:\t\t{}'.format(state, frequencies[state2index[state]]) for state in states)))
    logging.info('\tMissing data:\t\t{}'.format(missing_frequency / total_count))

    logging.info('Initial log likelihood is {}'.format(get_bottom_up_likelihood(tree, frequencies, sf)))

    frequencies, sf = minimize_params(tree, frequencies, sf, optimise_frequencies=False, optimise_sf=True)

    logging.info('Optimised SF is {}, i.e. {} changes per avg branch'.format(sf / avg_br_len, sf))
    logging.info('Optimised log likelihood so far is {}'.format(get_bottom_up_likelihood(tree, frequencies, sf)))

    frequencies, sf = minimize_params(tree, frequencies, sf, optimise_frequencies=True, optimise_sf=False)

    logging.info('Optimised frequencies are:{}'
                 .format(''.join('\n\t{}:\t\t{}'.format(state, frequencies[state2index[state]]) for state in states)))
    logging.info('Optimised log likelihood is {}'.format(get_bottom_up_likelihood(tree, frequencies, sf)))

    calculate_top_down_likelihood(tree, frequencies, sf)
    
    unalter_problematic_tip_states(tree, feature, state2index)
    
    calculate_marginal_probabilities(tree, frequencies)
    check_marginal_probabilities(tree)
    normalize_result_probabilities(tree)

    choose_ancestral_states(tree, len(frequencies))

    annotate(tree, state2index, feature, unique=True)

    visualize(tree, [feature], categories=states, html=html, html_compressed=html_compressed)


