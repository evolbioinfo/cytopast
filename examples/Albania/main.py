import logging
import os

from cytopast.pastml_analyser import pastml

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.txt')
STATES_INPUT_PERCENT = os.path.join(DATA_DIR, 'Annotation.Albanian.5chars.randomized_{}_percent.txt')

if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    for model in ('F81', 'JC'):
        res_tree = TREE_NWK
        pastml(data=STATES_INPUT, tree=res_tree,
               html_compressed=os.path.join(DATA_DIR, 'albania_compressed_{}.html'.format(model)),
               data_sep=',', work_dir=DATA_DIR,
               pastml_exe=os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', '..', 'PASTML'), model=model)

    # for percent in (5, 10, 25):
    #     df = pd.read_csv(STATES_INPUT, index_col=0, header=0)
    #     random_index = df.sample(frac=percent / 100).index
    #     new_values = np.random.choice(df['Country'].unique(), size=len(random_index))
    #     df.loc[random_index, 'Country'] = new_values
    #     randomized_state_file = STATES_INPUT_PERCENT.format(percent)
    #     df.to_csv(randomized_state_file, header=True, index=True)
    #
    #     pastml(data=randomized_state_file, tree=res_tree,
    #            html_compressed=os.path.join(DATA_DIR, 'albania_compressed_{}_{}.html'.format(percent, model)),
    #            data_sep=',', work_dir=os.path.join(DATA_DIR, '{}'.format(percent)),
    #            pastml_exe=os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', '..', 'PASTML'), model=model)
