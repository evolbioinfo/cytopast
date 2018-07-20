import pandas as pd

from cytopast import remove_certain_leaves, read_tree

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_tree', required=True, type=str)
    parser.add_argument('--metadata', required=True, type=str)
    parser.add_argument('--out_tree_pattern', required=True, type=str)
    parser.add_argument('--years', required=True, type=int, nargs='+')
    parser.add_argument('--date_column', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_table(params.metadata, header=0, index_col=0)
    df.index = df.index.map(str)
    tree = read_tree(params.in_tree)

    for year in sorted(params.years, key=lambda _: -_):
        tree = remove_certain_leaves(tree, to_remove=lambda node: pd.isnull(df.loc[node.name, params.date_column])
                                                                  or df.loc[node.name, params.date_column] > year)
        tree.write(outfile=params.out_tree_pattern.format(year), format=3)
