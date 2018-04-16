import logging
from collections import defaultdict, Counter
from functools import reduce
from queue import Queue

import numpy as np
import pandas as pd
from ete3 import Tree, TreeNode

REASONABLE_NUMBER_OF_TIPS = 25

CATEGORIES = 'categories'

SIZE = 'size'
MIN_NUM_TIPS_INSIDE = 'min_size'
MAX_NUM_TIPS_INSIDE = 'max_size'

TIPS_INSIDE = 'tips'

MIN_NUM_TIPS_BELOW = 'min_num_tips'
MAX_NUM_TIPS_BELOW = 'max_num_tips'

EDGE_SIZE = 'edge_size'
METACHILD = 'metachild'
FONT_SIZE = 'fontsize'


def name_tree(tree):
    """
    Names all the tree nodes that are not named, with unique names.
    :param tree: ete3.Tree, the tree to be named
    :return: void, modifies the original tree
    """
    existing_names = Counter((_.name for _ in tree.traverse() if _.name))
    i = 0
    for node in tree.traverse('postorder'):
        if node.is_root():
            node.name = 'ROOT'
            existing_names['ROOT'] += 1
        if not node.name or existing_names[node.name] > 1:
            name = None
            while not name or name in existing_names:
                name = '{}_{}'.format('tip' if node.is_leaf() else 'node', i)
                i += 1
            node.name = name


def collapse_zero_branches(tree):
    num_collapsed = 0
    for n in list(tree.traverse('postorder')):
        for child in list(n.children):
            if not child.is_leaf() and child.dist == 0:
                n.remove_child(child)
                for grandchild in child.children:
                    n.add_child(grandchild)
                num_collapsed += 1
    logging.info('Collapsed {} zero branches.'.format(num_collapsed))


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


def pasml_annotations2cytoscape_annotation(cat2file, output, sep='\t'):
    logging.info('Combining the data from different columns into {}'.format(output))

    def get_state(name, df):
        row = df.loc[name, :]
        states = df.columns[row]
        return None if len(states) != 1 else states[0]

    cat2df = {cat: pd.read_table(data_path, sep=',', index_col=0, header=0).astype(bool) for (cat, data_path) in
              cat2file.items()}
    df = pd.DataFrame(index=next(iter(cat2df.values())).index, columns=cat2df.keys())
    for cat, cat_df in cat2df.items():
        cat_df.index = cat_df.index.map(str)
        df[cat] = df.index.map(lambda name: get_state(name, cat_df))

    if len(cat2df) == 1:
        cat_df = next(iter(cat2df.values()))
        rsuffix = 'category' if set(df.columns) & set(cat_df.columns) else ''
        df = df.join(cat_df, rsuffix=rsuffix, lsuffix='')

    df.to_csv(output, sep=sep, index_label='Node')


def annotate_tree_with_cyto_metadata(tree, data_path, columns, sep='\t'):
    df = pd.read_table(data_path, sep=sep, index_col=0, header=0)
    df.index = df.index.map(str)
    df.fillna('', inplace=True)
    tree = read_tree(tree) if not isinstance(tree, TreeNode) else tree

    for n in tree.traverse():
        if len(columns) != 1:
            n.add_features(**df.loc[n.name, :].to_dict())
        else:
            category = col_name2cat(columns[0])
            data = df.loc[n.name, :]
            n.add_features(**data[data != False].to_dict())
            # In case category's value was also False (and therefore was not added), let's add it again
            n.add_feature(category, data[category])
    return tree, sorted(df.columns)


def col_name2cat(column):
    """
    Reformats the column string to make sure it contains only numerical or letter characters.
    :param column: str, column name to be reformatted
    :return: str, the column name with illegal characters removed
    """
    column_string = ''.join(s for s in column if s.isalnum())
    return column_string


def read_tree(tree_path):
    for format in (3, 2, 5, 1, 0, 3, 4, 6, 7, 8, 9):
        try:
            return Tree(tree_path, format=format)
        except:
            continue
    raise ValueError('Could not read the tree {}. Is it a valid newick?'.format(tree_path))


def get_states(n, categories):
    return {cat: getattr(n, cat) for cat in categories if hasattr(n, cat)}


def compress_tree(tree, categories, can_merge_diff_sizes=True, tip_size_threshold=REASONABLE_NUMBER_OF_TIPS,
                  name_feature=None):
    for n in tree.traverse():
        if n.is_leaf():
            n.add_feature(MIN_NUM_TIPS_INSIDE, 1)
            n.add_feature(MAX_NUM_TIPS_INSIDE, 1)
            n.add_feature(TIPS_INSIDE, [n.name])
        else:
            num_tips = len(n.get_leaves())
            n.add_feature(MIN_NUM_TIPS_BELOW, num_tips)
            n.add_feature(MAX_NUM_TIPS_BELOW, num_tips)

    collapse_vertically(tree, lambda _: get_states(_, categories))
    remove_mediators(tree, lambda _: get_states(_, categories))

    tip_sizes = set(getattr(l, MAX_NUM_TIPS_INSIDE, 0) for l in tree.iter_leaves())
    merge_different_sizes = len(tip_sizes) > 10 and can_merge_diff_sizes
    if merge_different_sizes:
        tips2bin = lambda _: int(np.log10(max(1, _)))
    else:
        tips2bin = lambda _: _

    logging.info('Gonna collapse horizontally, {}merging nodes of different sizes'
                 .format('' if merge_different_sizes else 'not '))
    collapse_horizontally(tips2bin, tree, lambda _: get_states(_, categories))

    if not merge_different_sizes and can_merge_diff_sizes and len(tree.get_leaves()) > REASONABLE_NUMBER_OF_TIPS:
        merge_different_sizes = True
        tips2bin = lambda _: int(np.log10(max(1, _)))

        logging.info('Gonna re-collapse horizontally, merging nodes of different sizes')
        collapse_horizontally(tips2bin, tree, lambda _: get_states(_, categories))

    tip_sizes = [getattr(_, MAX_NUM_TIPS_INSIDE, 0) * getattr(_, EDGE_SIZE, 1) for _ in tree.get_leaves()]
    if len(tip_sizes) > tip_size_threshold:
        threshold = sorted(tip_sizes)[-tip_size_threshold]
        logging.info('Removing tips of size less than {}'.format(threshold))
        remove_small_tips(tree,
                          to_be_removed=lambda _: getattr(_, MAX_NUM_TIPS_INSIDE, 0)
                                                  * getattr(_, EDGE_SIZE, 1) <= threshold)
        remove_mediators(tree, lambda _: get_states(_, categories))

    logging.info('Gonna collapse horizontally, {}merging nodes of different sizes'
                 .format('' if merge_different_sizes else 'not '))
    collapse_horizontally(tips2bin, tree, lambda _: get_states(_, categories))

    tip_sizes = [getattr(n, MAX_NUM_TIPS_INSIDE, 0) for n in tree.traverse() if getattr(n, MAX_NUM_TIPS_INSIDE, 0)]
    max_size = max(tip_sizes)
    min_size = min(tip_sizes)
    need_log = max_size / min_size > 100
    logging.info('Max vertical cluster size is {}, min is {}: {}need log'.format(max_size, min_size,
                                                                                 '' if need_log else 'do not '))

    e_szs = [getattr(n, EDGE_SIZE, 1) for n in tree.traverse()]
    max_e_size = max(e_szs)
    min_e_size = min(e_szs)
    logging.info('Max horizontal cluster size is {}, min is {}'.format(max_e_size, min_e_size))
    need_e_log = max_e_size / min_e_size > 100

    transform_size = lambda _: np.power(np.log10(_ + 9) if need_log else _, 1 / 2)
    transform_e_size = lambda _: np.log10(_) if need_e_log else _
    size_scaling = get_scaling_function(y_m=30, y_M=30 * min(8, int(max_size / min_size)),
                                        x_m=transform_size(min_size), x_M=transform_size(max_size))
    font_scaling = get_scaling_function(y_m=10, y_M=10 * min(3, int(max_size / min_size)),
                                        x_m=transform_size(min_size), x_M=transform_size(max_size))
    e_size_scaling = get_scaling_function(y_m=5, y_M=5 * min(5, int(max_e_size / min_e_size)),
                                          x_m=transform_e_size(min_e_size), x_M=transform_e_size(max_e_size))
    for n in tree.traverse():
        state = getattr(n, name_feature, '') if name_feature is not None else ''

        min_n_tips = getattr(n, MIN_NUM_TIPS_INSIDE, 0)
        max_n_tips = getattr(n, MAX_NUM_TIPS_INSIDE, 0)

        min_n_tips_below = getattr(n, MIN_NUM_TIPS_BELOW, 0)
        max_n_tips_below = getattr(n, MAX_NUM_TIPS_BELOW, 0)

        edge_size = getattr(n, EDGE_SIZE, 1)

        n.state = '{} {} {}'.format(state,
                                    ('{}-{}'.format(min_n_tips, max_n_tips) if min_n_tips != max_n_tips else min_n_tips)
                                    if max_n_tips > 0 else '',
                                    '({})'.format(
                                        '{}-{}'.format(min_n_tips_below,
                                                       max_n_tips_below) if min_n_tips_below != max_n_tips_below
                                        else min_n_tips_below)
                                    if max_n_tips_below > 0 else '')

        n.add_feature(SIZE, 20 if max_n_tips == 0 else size_scaling(transform_size(max_n_tips)))
        n.add_feature(FONT_SIZE, 10 if max_n_tips == 0 else font_scaling(transform_size(max_n_tips)))

        n.add_feature('edge_name', str(edge_size) if edge_size > 1 else '')
        n.add_feature(EDGE_SIZE, e_size_scaling(transform_e_size(edge_size)))

    return tree


def collapse_horizontally(tips2bin, tree, get_states):
    parents = [tree]
    while parents:
        _collapse_horizontally(get_states, parents=parents, tips2bin=tips2bin)
        parents = reduce(lambda l1, l2: l1 + l2, (p.children for p in parents))


def get_scaling_function(y_m, y_M, x_m, x_M):
    # calculate a linear function y = k x + b, where y \in [m, M]
    if x_M <= x_m:
        return lambda _: y_m
    k = (y_M - y_m) / (x_M - x_m)
    b = y_m - k * x_m
    return lambda _: int(k * _ + b)


def _collapse_horizontally(get_states, parents, tips2bin=lambda _: _):
    def get_sorted_states(n, add_edge_size=True):
        return tuple(sorted('{}:{}'.format(k, v) for (k, v) in get_states(n).items())), \
               tips2bin(getattr(n, MAX_NUM_TIPS_INSIDE, 0)), \
               (getattr(n, EDGE_SIZE, 1) if add_edge_size else -1)

    def get_configuration(n):
        queue = Queue()
        queue.put((0, n), block=False)
        config = [(0, get_sorted_states(n, False))]
        while not queue.empty():
            level, n = queue.get(block=False)
            for (child, states) in sorted(((_, get_sorted_states(_)) for _ in n.children), key=lambda _: _[1]):
                config.append((level + 1, states))
                queue.put((level + 1, child))
        return tuple(config)

    for p in parents:
        state2children = defaultdict(list)
        for c in p.children:
            state2children[get_configuration(c)].append(c)
        for children in (_ for _ in state2children.values() if len(_) > 1):
            edge_size = sum(getattr(c, EDGE_SIZE, 1) for c in children)
            child = children[0]
            for c in children[1:]:
                p.remove_child(c)
            child.add_feature(EDGE_SIZE, edge_size)
            child.add_feature(METACHILD, True)
            merge_features(child, children, (MIN_NUM_TIPS_INSIDE,), min, default_value=0)
            merge_features(child, children, (MAX_NUM_TIPS_INSIDE,), max, default_value=0)

            merge_features(child, children, (MIN_NUM_TIPS_BELOW,), min, default_value=0)
            merge_features(child, children, (MAX_NUM_TIPS_BELOW,), max, default_value=0)

            child.add_feature(TIPS_INSIDE, [getattr(_, TIPS_INSIDE, []) for _ in children])


def remove_small_tips(tree, to_be_removed):
    changed = True
    while changed:
        changed = False
        for l in tree.get_leaves():
            if l.up and to_be_removed(l):
                l.up.remove_child(l)
                changed = True


def merge_features(main_node, nodes, features, op, default_value=0):
    for feature in features:
        main_node.add_feature(feature, op([getattr(node, feature, default_value) for node in nodes]))


def collapse_vertically(tree, get_states):
    """
    Collapses a child node into its parent if they are in the same state.
    :param get_states: a function that returns a set of node states
    :param tree: ete3.Tree
    :return: void, modifies the input tree
    """
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            continue

        states = set(get_states(n).items())
        children = list(n.children)
        for child in children:
            # merge the child into this node if their states are the same
            if set(get_states(child).items()) == states:
                old_max_tips = getattr(n, MAX_NUM_TIPS_INSIDE, 0)
                old_min_tips = getattr(n, MIN_NUM_TIPS_INSIDE, 0)

                merge_features(n, (n, child), (MAX_NUM_TIPS_INSIDE, MIN_NUM_TIPS_INSIDE), sum, default_value=0)
                n.add_feature(TIPS_INSIDE, getattr(n, TIPS_INSIDE, []) + getattr(child, TIPS_INSIDE, []))

                max_tips = getattr(n, MAX_NUM_TIPS_INSIDE, 0)
                min_tips = getattr(n, MIN_NUM_TIPS_INSIDE, 0)

                n.add_feature(MAX_NUM_TIPS_BELOW, getattr(n, MAX_NUM_TIPS_BELOW) - (max_tips - old_max_tips))
                n.add_feature(MIN_NUM_TIPS_BELOW, getattr(n, MIN_NUM_TIPS_BELOW) - (min_tips - old_min_tips))

                n.remove_child(child)
                for grand_child in child.children:
                    n.add_child(grand_child)


def remove_mediators(tree, get_states):
    """
    Removes intermediate nodes that are just mediators between their parent and child states.
    :param get_states: a function that returns a set of node states
    :param tree: ete3.Tree
    :return: void, modifies the input tree
    """
    for n in tree.traverse('postorder'):
        if getattr(n, METACHILD, False):
            continue
        states = get_states(n)
        if n.is_leaf() or len(n.children) > 1:
            continue
        parent = n.up
        if not parent:
            continue
        parent_states = get_states(parent)
        child = n.children[0]
        child_states = get_states(child)
        if set(states.keys()) == set(parent_states.keys()) == set(child_states.keys()):
            compatible = next((False for (key, state) in states.items()
                               if state and (state != parent_states[key] or state != child_states[key])), True)
        else:
            compatible = set(states.items()) == set(parent_states.items()) | set(child_states.items())
        if compatible:
            old_max_tips = getattr(parent, MAX_NUM_TIPS_INSIDE, 0)
            old_min_tips = getattr(parent, MIN_NUM_TIPS_INSIDE, 0)

            merge_features(parent, (parent, n), (MAX_NUM_TIPS_INSIDE, MIN_NUM_TIPS_INSIDE), sum, default_value=0)
            parent.add_feature(TIPS_INSIDE, getattr(parent, TIPS_INSIDE, []) + getattr(n, TIPS_INSIDE, []))

            max_tips = getattr(parent, MAX_NUM_TIPS_INSIDE, 0)
            min_tips = getattr(parent, MIN_NUM_TIPS_INSIDE, 0)

            parent.add_feature(MAX_NUM_TIPS_BELOW, getattr(parent, MAX_NUM_TIPS_BELOW) - (max_tips - old_max_tips))
            parent.add_feature(MIN_NUM_TIPS_BELOW, getattr(parent, MIN_NUM_TIPS_BELOW) - (min_tips - old_min_tips))

            parent.remove_child(n)
            for c in n.children:
                parent.add_child(c)
