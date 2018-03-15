import glob
import logging
import os

from cytopast import pasml_annotations2cytoscape_annotation
from cytopast.pastml_analyser import _past_vis

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'tree.nwk')


def visualise(tree, res_data, html_compressed, html=None, name=None):
    columns = ['Location']
    col2annotation_files = {columns[0]: res_data}
    res_annotations = os.path.join(DATA_DIR, 'combined_annotations_{}.tab'.format(name))
    pasml_annotations2cytoscape_annotation(col2annotation_files, res_annotations)
    _past_vis(tree, res_annotations, html_compressed, html, columns=columns)


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    for state_file in glob.glob(os.path.join(DATA_DIR, 'state_Location.csv.pastml.out*')):
        name = state_file.replace(os.path.join(DATA_DIR, 'state_Location.csv.pastml.out.'), '')
        visualise(TREE_NWK, state_file, html_compressed=os.path.join(DATA_DIR, 'map_{}.html'.format(name)), 
                  html=os.path.join(DATA_DIR, 'tree_{}.html'.format(name)), name=name)
