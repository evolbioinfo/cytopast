from ete3 import Tree
import pandas as pd


def read_tree(tree_path):
    for f in (3, 2, 5, 1, 0, 3, 4, 6, 7, 8, 9):
        try:
            return Tree(tree_path, format=f)
        except:
            continue
    raise ValueError('Could not read the tree {}. Is it a valid newick?'.format(tree_path))


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--trees', required=True, type=str, nargs='+')
    parser.add_argument('--output', required=True, type=str)
    params = parser.parse_args()

    results = pd.DataFrame(columns=['tree1', 'tree2', 'RF'])
    path2tree = {_: read_tree(_) for _ in params.trees}
    for path1 in params.trees:
        tree1 = read_tree(path1)
        for path2 in params.trees:
            tree2 = read_tree(path2)
            results.append({'tree1': path1, 'tree2': path2,
                            'RF': tree1.robinson_foulds(tree2, unrooted_trees=False, expand_polytomies=False,
                                                        polytomy_size_limit=1000)}, ignore_index=True)
    results.to_csv(params.output, sep='\t', index=False)
