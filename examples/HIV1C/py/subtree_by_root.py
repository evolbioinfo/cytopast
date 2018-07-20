from cytopast import read_tree

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_tree', required=True, type=str)
    parser.add_argument('--root', required=True, type=str)
    parser.add_argument('--out_tree', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.in_tree)
    for _ in tree.traverse():
        if _.name == params.root:
            _.write(outfile=params.out_tree, format=3)
