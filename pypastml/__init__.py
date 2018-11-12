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

