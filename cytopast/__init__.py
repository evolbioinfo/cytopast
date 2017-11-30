import logging
import os
from collections import defaultdict
from functools import reduce
from queue import Queue

import numpy as np
import pandas as pd
from ete3 import Tree

CATEGORIES = 'categories'

PASTML = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'PASTML')
# PASTML = 'PASTML'
TREE_NWK_PASTML_OUTPUT = 'Result_treeIDs.{tips}.taxa.{states}.states.tre'
STATES_TAB_PASTML_OUTPUT = 'Result_states_probs.FULL.{tips}.taxa.{states}.states.txt'

SIZE = 'size'
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


def apply_pastml(annotation_file, tree_file, pastml=PASTML, out_dir=None, model='JC'):
    """
    Applies PASTML on the given tree and annotation file.
    :param annotation_file: path to the csv state file: tip_name,state
    :param tree_file: path to the tree nwk file
    :param pastml: path to the PASTML binary
    :return: path to the annotation file produced by PASTML
    """
    tree = read_tree(tree_file)

    df = pd.read_csv(annotation_file, index_col=0, header=None)

    names = df.index.astype(np.str)
    nodes = [n for n in tree.iter_leaves() if n.name in names]
    n_tips = len(nodes)
    need_to_prune = len(tree.get_leaves()) > n_tips
    if need_to_prune:
        logging.info('Pruning...')
        tree.prune(nodes, preserve_branch_length=True)
    if need_to_prune or name_tree(tree):
        tree.write(outfile=tree_file, format=3, format_root_node=True)
    df = df[np.in1d(names, [n.name for n in tree.iter_leaves()])]
    df.to_csv(annotation_file, header=False, index=True)

    states = df[1].unique()
    logging.info('States are {}'.format(states))
    n_states = len([s for s in states if not pd.isnull(s)])
    if out_dir is None:
        out_dir = os.path.dirname(annotation_file)
    os.makedirs(out_dir, exist_ok=True)

    command = 'cd {dir}; {pastml} -a {annotation_file} -t {tree_file} -m {model} -I T > {log_file}'.format(
        dir=out_dir, pastml=pastml, annotation_file=annotation_file, tree_file=tree_file, model=model,
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
        return None if len(states) > 1 else states[0]

    cat2df = {cat: pd.read_table(data_path, sep=',', index_col=0, header=0) for (cat, data_path) in
              cat2file.items()}
    df = pd.DataFrame(index=next(iter(cat2df.values())).index, columns=cat2df.keys())
    for cat, cat_df in cat2df.items():
        cat_df.index = cat_df.index.map(str)
        df[cat] = df.index.map(lambda name: get_states(name, cat_df))

    logging.info(df.sample(n=5))
    df.to_csv(output, sep=sep)


def annotate_tree_with_cyto_metadata(tree_path, data_path, sep='\t', one_state=False):
    df = pd.read_table(data_path, sep=sep, index_col=0, header=0)
    df.index = df.index.map(str)
    df.fillna('', inplace=True)
    tree = read_tree(tree_path)

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
        tree = Tree(tree_path, format=1)
    return tree


def compress_tree(tree, categories, can_merge_diff_sizes=True, cut=True, name_feature=None):
    for n in tree.traverse():
        n.add_feature('num_tips', len(n.get_leaves()))
        if n.is_leaf():
            n.add_feature(SIZE, 1)

    def get_states(n):
        return set('{}:{}'.format(cat, getattr(n, cat)) for cat in categories if hasattr(n, cat))

    collapse_vertically(tree, get_states)
    tip_sizes = set(getattr(l, SIZE, 0) for l in tree.iter_leaves())
    merge_different_sizes = max(tip_sizes) / min(tip_sizes) > 10 and can_merge_diff_sizes
    if merge_different_sizes:
        tips2bin = lambda n_tips: int(np.log10(max(1, n_tips)))
    else:
        tips2bin = lambda n_tips: n_tips

    logging.info('Gonna collapse horizontally, {}merging nodes of different sizes'
                 .format('' if merge_different_sizes else 'not '))
    parents = [tree]
    while parents:
        collapse_horizontally(get_states, parents=parents, tips2bin=tips2bin)
        parents = reduce(lambda l1, l2: l1 + l2, (p.children for p in parents))

    def get_size(n):
        return getattr(n, EDGE_SIZE, 1) * (getattr(n, SIZE, 0) + sum(get_size(c) for c in n.children))

    n2size = {n: get_size(n) for n in tree.children}
    sizes = sorted(n2size.values())
    threshold_size = 0 if len(sizes) < 10 else sizes[-10]
    logging.info('Threshold size is set to {}'.format(threshold_size))
    for n, size in n2size.items():
        if size < threshold_size:
            tree.remove_child(n)

    szs = sorted(getattr(n, SIZE, 0) * getattr(n, EDGE_SIZE, 1) for n in tree.iter_leaves())
    if len(szs) > 10 and cut:
        threshold = szs[-10]
        logging.info('Removing tips of size less than {}'.format(threshold))
        remove_small_tips(threshold, tree)

    szs = [getattr(n, SIZE, 0) for n in tree.traverse() if getattr(n, SIZE, 0)]
    max_size = max(szs)
    min_size = min(szs)
    need_log = max_size / min_size > 100
    if need_log:
        max_size = np.log10(max_size)
        min_size = np.log10(min_size)

    e_szs = [getattr(n, EDGE_SIZE, 1) for n in tree.traverse()]
    max_e_size = max(e_szs)
    min_e_size = min(e_szs)
    need_e_log = max_e_size / min_e_size > 100
    if need_e_log:
        max_e_size = np.log10(max_e_size)
        min_e_size = np.log10(min_e_size)

    logging.info('Max cluster size is {}, min is {}'.format(max_size, min_size))

    for n in tree.traverse():
        n_tips = getattr(n, SIZE, 0)
        edge_size = getattr(n, EDGE_SIZE, 1)
        scaled_size = ((np.log10(max(n_tips, 1)) if need_log else n_tips) - min_size) / max(max_size - min_size, 1)
        n.add_feature(SIZE, 20 if n_tips == 0 else int(40 + 360 * scaled_size))
        n.add_feature(FONT_SIZE, 10 if n_tips == 0 else int(10 + 40 * scaled_size))

        n.add_feature('edge_name', str(edge_size) if edge_size > 1 else '')
        scaled_e_size = ((np.log10(edge_size) if need_e_log else edge_size) - min_e_size) \
                        / max(max_e_size - min_e_size, 1)
        n.add_feature(EDGE_SIZE, int(10 * (1 + scaled_e_size)))

        state = getattr(n, name_feature, '') if name_feature is not None else ''
        is_metachild = getattr(n, METACHILD, False)
        real_num_tips = getattr(n, 'num_tips')
        n.state = '{} {}{}'.format(state,
                                   '{}{}'.format('~' if is_metachild and merge_different_sizes else '', n_tips)
                                   if n_tips else '',
                                   ' ({})'.format(real_num_tips - n_tips) if real_num_tips > n_tips else '')
        # if hasattr(n, METACHILD):
        #     n.del_feature(METACHILD)

    return tree


def collapse_horizontally(get_states, parents, tips2bin=lambda n_tips: int(np.log10(max(1, n_tips)))):
    def get_sorted_states(n):
        return tuple(sorted(get_states(n))), tips2bin(getattr(n, SIZE, 0))

    def get_configuration(n):
        queue = Queue()
        queue.put((0, n), block=False)
        config = [(0, get_sorted_states(n))]
        while not queue.empty():
            level, n = queue.get(block=False)
            for (c, states) in sorted(((c, get_sorted_states(c)) for c in n.children), key=lambda kv: kv[1]):
                config.append((level + 1, states))
                queue.put((level + 1, c))
        return tuple(config)

    for p in parents:
        state2children = defaultdict(list)
        for c in p.children:
            state2children[get_configuration(c)].append(c)
        for children in (children for children in state2children.values() if len(children) > 1):
            edge_size = len(children)
            child = children[0]
            for c in children[1:]:
                p.remove_child(c)
            child.add_feature(EDGE_SIZE, edge_size)
            child.add_feature(METACHILD, True)
            child.add_feature(SIZE,
                              int(round(np.mean([getattr(c, SIZE, 0) for c in children]))))
            # queue = Queue()
            # queue.put(child, block=False)
            # while not queue.empty():
            #     child = queue.get(block=False)
            #     for grandchild in child.children:
            #         queue.put(grandchild, block=False)


def remove_small_tips(threshold, tree):
    changed = True
    while changed:
        changed = False
        for l in tree.get_leaves():
            if l.up and getattr(l, SIZE, 0) * getattr(l, EDGE_SIZE, 1) < threshold:
                l.up.remove_child(l)
                changed = True


def collapse_vertically(tree, get_states):
    """
    Collapses a child node into its parent if they are in the same state.
    :param get_states: a function that returns a set of node states
    :param tree: ete3.Tree
    :return: void, modifies the input tree
    """
    queue = Queue()
    for n in tree.traverse('postorder'):
        queue.put(n, block=False)
    while not queue.empty():
        n = queue.get(block=False)
        states = get_states(n)
        parent = n.up
        if not parent:
            continue
        parent_states = get_states(parent)
        # merge the child into its parent if their states are the same
        if parent_states == states:
            parent.add_feature(SIZE, getattr(parent, SIZE, 0) + getattr(n, SIZE, 0))
            parent.remove_child(n)
            for c in n.children:
                parent.add_child(c)
            continue
        # remove intermediate nodes that are just mediators between their parent and child states
        if n.is_leaf() or len(n.children) > 1:
            continue
        child = n.children[0]
        if states == parent_states | get_states(child):
            parent.add_feature(SIZE, getattr(parent, SIZE, 0) + getattr(n, SIZE, 0))
            parent.remove_child(n)
            parent.add_child(child)
