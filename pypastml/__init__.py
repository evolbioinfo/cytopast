import numpy as np

ALLOWED_STATES = 'ALLOWED_STATES'


def rescale_tree(tree, sf):
    """
    Multiplies all the tree branches by the given scaling factor.

    :param tree: ete3.Tree, the tree to be rescaled
    :param sf: float, scaling factor
    :return: void, the initial tree is modified
    """
    for n in tree.traverse():
        n.dist *= sf


def get_personalized_feature_name(character, feature):
    """
    Precedes the feature name by the character name
    (useful when likelihoods for different characters are calculated in parallel).

    :param character: str, character name
    :param feature: str, feature to be personalized
    :return: str, the personalized feature
    """
    return '{}_{}'.format(character, feature)


def value2list(n, value, default_value):
    # create a variable for n columns
    # specifying the default value if nothing was chosen
    if value is None:
        value = default_value
    # and propagating the chosen value to all columns
    if not isinstance(value, list):
        value = [value] * n
    # or making sure that the default value if chosen for the columns for which the value was not specified
    else:
        value += [default_value] * (n - len(value))
    return value


def initialize_allowed_states(tree, feature, states):
    """
    Initializes the allowed state arrays for tips based on their states given by the feature.
    :param tree: ete3.Tree, tree for which the tip likelihoods are to be initialized
    :param feature: str, feature in which the tip states are stored (the value could be None for a missing state)
    :param states: numpy array of ordered states.
    :return: void, adds the get_personalised_feature_name(feature, ALLOWED_STATES) feature to tree tips.
    """
    allowed_states_feature = get_personalized_feature_name(feature, ALLOWED_STATES)

    all_ones, state2array = get_state2allowed_states(states)

    for node in tree.traverse():
        if node.is_leaf():
            node.add_feature(allowed_states_feature, state2array[getattr(node, feature, '')])
        else:
            node.add_feature(allowed_states_feature, all_ones)


def get_state2allowed_states(states, by_name=True):
    # tips allowed state arrays won't be modified so we might as well just share them
    n = len(states)
    all_ones = np.ones(n, np.int)
    state2array = {}
    for index, state in enumerate(states):
        allowed_state_array = np.zeros(n, np.int)
        allowed_state_array[index] = 1
        state2array[state if by_name else index] = allowed_state_array
    if by_name:
        state2array[None] = all_ones
        state2array[''] = all_ones
    return all_ones, state2array

