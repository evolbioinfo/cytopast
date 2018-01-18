import logging
import os
from collections import defaultdict, Counter
from functools import reduce
from queue import Queue

import numpy as np
import pandas as pd
from ete3 import Tree, TreeNode

CATEGORIES = 'categories'

TREE_NWK_PASTML_OUTPUT = 'Result_treeIDs.{tips}.taxa.{states}.states.tre'
STATES_TAB_PASTML_OUTPUT = 'Result_states_probs.FULL.{tips}.taxa.{states}.states.txt'

SIZE = 'size'
MIN_NUM_TIPS_INSIDE = 'min_size'
MAX_NUM_TIPS_INSIDE = 'max_size'

MIN_NUM_TIPS_BELOW = 'min_num_tips'
MAX_NUM_TIPS_BELOW = 'max_num_tips'

EDGE_SIZE = 'edge_size'
METACHILD = 'metachild'
FONT_SIZE = 'fontsize'


def name_tree(tree):
    """
    Names unnamed tree nodes.
    :param tree: ete3.Tree
    :return: True if the tree was changed, False otherwise (if all the nodes were already named)
    """
    i = 0
    names = list(n.name for n in tree.traverse())
    need_to_rename = len(set(names)) != len(names)
    for n in tree.traverse():
        if not n.is_leaf() and (not n.name or need_to_rename):
            if n.is_root():
                n.name = 'ROOT'
            else:
                n.name = 'node_{i}'.format(i=i)
            i += 1
    if i > 0:
        logging.info('Named tree internal nodes')
    return i > 0


def apply_pastml(annotation_file, tree_file, pastml, out_dir=None, model='JC'):
    """
    Applies PASTML on the given tree and annotation file.
    :param annotation_file: path to the csv state file: tip_name,state
    :param tree_file: path to the tree nwk file
    :param pastml: path to the PASTML binary
    :return: path to the annotation file produced by PASTML
    """
    tree = read_tree(tree_file)
    n_tips = len(tree.get_leaves())

    df = pd.read_csv(annotation_file, index_col=0, header=None)

    names = df.index.astype(np.str)
    df = df[np.in1d(names, [n.name for n in tree.iter_leaves()])]
    df.to_csv(annotation_file, header=False, index=True)

    states = df[1].unique()
    logging.info('States are {}'.format(states))
    n_states = len([s for s in states if not pd.isnull(s)])
    if out_dir is None:
        out_dir = os.path.dirname(annotation_file)
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    command = 'cd {dir}; {pastml} -a {annotation_file} -t {tree_file} -m {model} -I T > {log_file}'.format(
        dir=out_dir, pastml=pastml, annotation_file=os.path.abspath(annotation_file),
        tree_file=os.path.abspath(tree_file), model=model,
        log_file=os.path.join(out_dir, 'pastml_log.log'))
    logging.info(command)
    os.system(command)

    res_data = os.path.join(out_dir, STATES_TAB_PASTML_OUTPUT).format(tips=n_tips, states=n_states)
    pd.read_table(res_data, sep=', ', header=0, index_col=0).astype(bool).to_csv(res_data, sep=',', index=True)
    return res_data


def pasml_annotations2cytoscape_annotation(cat2file, output, sep='\t'):

    def get_states(name, cat_df):
        row = cat_df.loc[name, :]
        states = cat_df.columns[row]
        return None if len(states) != 1 else states[0]

    cat2df = {cat: pd.read_table(data_path, sep=',', index_col=0, header=0) for (cat, data_path) in
              cat2file.items()}
    df = pd.DataFrame(index=next(iter(cat2df.values())).index, columns=cat2df.keys())
    for cat, cat_df in cat2df.items():
        cat_df.index = cat_df.index.map(str)
        df[cat] = df.index.map(lambda name: get_states(name, cat_df))

    logging.info(df.sample(n=5))
    df.to_csv(output, sep=sep)


def annotate_tree_with_cyto_metadata(tree, data_path, sep='\t', one_state=False):
    df = pd.read_table(data_path, sep=sep, index_col=0, header=0)
    df.index = df.index.map(str)
    df.fillna('', inplace=True)
    tree = read_tree(tree) if not isinstance(tree, TreeNode) else tree

    for n in tree.traverse():
        if not one_state:
            n.add_features(**df.loc[n.name, :].to_dict())
        else:
            data = df.loc[n.name, :]
            data = data[data != False]
            n.add_features(**data.to_dict())
    return tree, sorted(df.columns)


def annotate_tree_with_metadata(tree_path, data_path, sep=','):
    df = pd.read_table(data_path, sep=sep, index_col=0, header=0)
    df.index = df.index.map(str)
    tree = read_tree(tree_path)

    def get_states(name):
        row = df.loc[name, :]
        return df.columns[row > 0]

    for n in tree.traverse():
        states = get_states(n.name)
        n.add_feature('state', str(states[0]) if len(states) == 1 else '')
        n.add_feature(CATEGORIES, tuple(sorted(states)))
    return tree, sorted(df.columns)


def read_tree(tree_path):
    try:
        tree = Tree(tree_path, format=3)
    except:
        try:
            tree = Tree(tree_path, format=1)
        except:
            try:
                tree = Tree(tree_path, format=2)
            except:
                tree = Tree(tree_path, format=0)
    return tree


def get_states(n, categories):
    return ['{}:{}'.format(cat, getattr(n, cat)) for cat in categories if hasattr(n, cat)]


def compress_tree(tree, categories, can_merge_diff_sizes=True, cut=True, name_feature=None):
    for n in tree.traverse():
        if n.is_leaf():
            n.add_feature(MIN_NUM_TIPS_INSIDE, 1)
            n.add_feature(MAX_NUM_TIPS_INSIDE, 1)
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
    parents = [tree]
    while parents:
        collapse_horizontally(lambda _: get_states(_, categories), parents=parents, tips2bin=tips2bin)
        parents = reduce(lambda l1, l2: l1 + l2, (p.children for p in parents))

    if cut:
        tip_sizes = [getattr(_, MAX_NUM_TIPS_INSIDE, 0) * getattr(_, EDGE_SIZE, 1) for _ in tree.get_leaves()]
        if len(tip_sizes) > 15:
            threshold = sorted(tip_sizes)[-15]
            logging.info('Removing tips of size less than {}'.format(threshold))
            remove_small_tips(tree,
                              to_be_removed=lambda _: getattr(_, MAX_NUM_TIPS_INSIDE, 0)
                                                      * getattr(_, EDGE_SIZE, 1) <= threshold)
        remove_mediators(tree, lambda _: get_states(_, categories))
    logging.info('Gonna collapse horizontally, {}merging nodes of different sizes'
                 .format('' if merge_different_sizes else 'not '))
    parents = [tree]
    while parents:
        collapse_horizontally(lambda _: get_states(_, categories), parents=parents, tips2bin=tips2bin)
        parents = reduce(lambda l1, l2: l1 + l2, (p.children for p in parents))

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


def get_scaling_function(y_m, y_M, x_m, x_M):
    # calculate a linear function y = k x + b, where y \in [m, M]
    if x_M <= x_m:
        return lambda _: y_m
    k = (y_M - y_m) / (x_M - x_m)
    b = y_m - k * x_m
    return lambda _: int(k * _ + b)


def collapse_horizontally(get_states, parents, tips2bin=lambda _: _):
    def get_sorted_states(n, add_edge_size=True):
        return tuple(sorted(get_states(n))), tips2bin(getattr(n, MAX_NUM_TIPS_INSIDE, 0)), \
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
        states = get_states(n)
        parent = n.up
        if not parent:
            continue
        parent_states = get_states(parent)
        # merge the child into its parent if their states are the same
        if parent_states == states:
            old_max_tips = getattr(parent, MAX_NUM_TIPS_INSIDE, 0)
            old_min_tips = getattr(parent, MIN_NUM_TIPS_INSIDE, 0)

            merge_features(parent, (parent, n),
                           (MAX_NUM_TIPS_INSIDE, MIN_NUM_TIPS_INSIDE), sum, default_value=0)

            max_tips = getattr(parent, MAX_NUM_TIPS_INSIDE, 0)
            min_tips = getattr(parent, MIN_NUM_TIPS_INSIDE, 0)

            parent.add_feature(MAX_NUM_TIPS_BELOW, getattr(parent, MAX_NUM_TIPS_BELOW) - (max_tips - old_max_tips))
            parent.add_feature(MIN_NUM_TIPS_BELOW, getattr(parent, MIN_NUM_TIPS_BELOW) - (min_tips - old_min_tips))

            parent.remove_child(n)
            for c in n.children:
                parent.add_child(c)


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
        if states == set(parent_states) | set(get_states(child)):
            old_max_tips = getattr(parent, MAX_NUM_TIPS_INSIDE, 0)
            old_min_tips = getattr(parent, MIN_NUM_TIPS_INSIDE, 0)

            merge_features(parent, (parent, n),
                           (MAX_NUM_TIPS_INSIDE, MIN_NUM_TIPS_INSIDE), sum, default_value=0)

            max_tips = getattr(parent, MAX_NUM_TIPS_INSIDE, 0)
            min_tips = getattr(parent, MIN_NUM_TIPS_INSIDE, 0)

            parent.add_feature(MAX_NUM_TIPS_BELOW, getattr(parent, MAX_NUM_TIPS_BELOW) - (max_tips - old_max_tips))
            parent.add_feature(MIN_NUM_TIPS_BELOW, getattr(parent, MIN_NUM_TIPS_BELOW) - (min_tips - old_min_tips))

            parent.remove_child(n)
            for c in n.children:
                parent.add_child(c)
