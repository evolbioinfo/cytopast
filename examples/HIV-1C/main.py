import os

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE1 = os.path.join(DATA_DIR, 'collen_tree.nwk')
TREE2 = os.path.join(DATA_DIR, 'tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data.tab')

if '__main__' == __name__:
    columns = [' Region ', 'Location', 'RT:M184V', 'RT:K103N']
    for model in ('F81', 'JC'):
        for tree, label in ((TREE2, 'zero'), (TREE1, 'poly'),):
            print('Processing model {}, tree {}'.format(model, label))
            pastml_pipeline(data=STATES_INPUT, tree=tree,
                            columns=columns,
                            html_compressed=os.path.join(DATA_DIR, 'map_{}_{}_{}.html'.
                                                         format(model, label, '_'.join(columns))),
                            html=os.path.join(DATA_DIR, 'tree_{}_{}_{}.html'.
                                              format(model, label, '_'.join(columns))), model=model,
                            name_column=' Region ', verbose=True)
