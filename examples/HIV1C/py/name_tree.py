from collections import Counter

from ete3 import Tree
from ete3.parser.newick import write_newick


def read_tree(tree_path):
    for f in (3, 2, 5, 1, 0, 3, 4, 6, 7, 8, 9):
        try:
            return Tree(tree_path, format=f)
        except:
            continue
    raise ValueError('Could not read the tree {}. Is it a valid newick?'.format(tree_path))


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


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_tree', required=True, type=str)
    params = parser.parse_args()

    tr = read_tree(params.input_tree)
    name_tree(tr)

    nwk = write_newick(tr, format_root_node=True)
    with open(params.output_tree, 'w+') as f:
        f.write('%s\n' % nwk)
