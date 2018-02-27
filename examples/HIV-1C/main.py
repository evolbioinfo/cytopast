import os

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE = os.path.join(DATA_DIR, 'collen_tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data.tab')

if '__main__' == __name__:
    columns = [' Region ', 'Location', 'RT:M184V', 'RT:K103N']
    for model in ('F81', 'JC'):
            print('Processing model {}'.format(model))
            pastml_pipeline(data=STATES_INPUT, tree=TREE,
                            columns=columns,
                            html_compressed=os.path.join(DATA_DIR, 'map_{}_{}.html'.
                                                         format(model, '_'.join(columns))),
                            html=os.path.join(DATA_DIR, 'tree_{}_{}.html'.
                                              format(model, '_'.join(columns))), model=model,
                            name_column=' Region ', verbose=True)
