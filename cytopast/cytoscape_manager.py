import json

from colour import Color
from jinja2 import Environment, PackageLoader

FONT_SIZE = 'fontsize'

DEFAULT_SIZE = 50

STATE = 'state'

SIZE = 'size'
EDGE_SIZE = 'edge_size'

DATA = 'data'
ID = 'id'
NAME = 'name'
EDGES = 'edges'
NODES = 'nodes'
ELEMENTS = 'elements'

INTERACTION = 'interaction'

TARGET = 'target'
SOURCE = 'source'


def tree2json(tree, add_fake_nodes=True):
    nodes, edges = [], []
    if add_fake_nodes:
        max_dist = max(n.dist for n in tree.traverse())
        dist_step = max_dist / 10 if max_dist < 10 or max_dist > 100 else 1
    for n in tree.traverse('preorder'):
        features = {feature: getattr(n, feature) for feature in n.features}
        features[ID] = n.name
        features[NAME] = n.state
        del features['state']
        if SIZE not in n.features:
            features[SIZE] = DEFAULT_SIZE
        if FONT_SIZE not in n.features:
            features[FONT_SIZE] = 10
        if 'shape' not in n.features:
            features['shape'] = 'ellipse'
        nodes.append(get_node(**features))
        for child in n.children:
            source_name = n.name
            edge_size = getattr(child, EDGE_SIZE, None)
            if hasattr(child, 'edge_name'):
                edge_name = getattr(child, 'edge_name')
            else:
                edge_name = '%g' % child.dist
            if not edge_size:
                edge_size = 6
            if add_fake_nodes and int(child.dist / dist_step) > 0:
                target_name = 'fake_node_{}'.format(child.name)
                nodes.append(get_node(**{ID: target_name, NAME: '', 'is_fake': 1, SIZE: 1, 'shape': 'rectangle',
                                         FONT_SIZE: 1}))
                edges.append(get_edge(**{SOURCE: source_name, TARGET: target_name, INTERACTION: 'none',
                                         SIZE: edge_size, NAME: ''}))
                source_name = target_name
            edge_data = {SOURCE: source_name, TARGET: child.name, INTERACTION: 'triangle',
                         SIZE: edge_size, NAME: edge_name}
            if add_fake_nodes:
                edge_data['minLen'] = int(child.dist / dist_step)
            edges.append(get_edge(**edge_data))

    json_dict = {NODES: nodes, EDGES: edges}
    return json_dict


def json2cyjs(json_dict, out_cyjs, graph_name='Tree'):
    json_dict = {DATA: {NAME: graph_name}, ELEMENTS: json_dict}
    with open(out_cyjs, 'w+') as fp:
        json.dump(json_dict, fp)


def save_as_cytoscape_html(tree, out_html, categories, graph_name='Untitled', dist_step=None, layout='dagre'):
    """
    Converts a tree to an html representation using Cytoscape.js.

    If categories are specified they are visualised as pie-charts inside the nodes,
    given that each node contains features corresponding to these categories with values being the percentage.
    For instance, given categories ['A', 'B', 'C'], a node with features {'A': 50, 'B': 50}
    will have a half-half pie-chart (half-colored in a colour of A, and half B).

    If dist_step is specified, the edges are rescaled accordingly to their dist (node.dist / dist_step),
    otherwise all edges are drawn of the same length.

    :param dist_step: if specified, the edges are rescaled accordingly to their dist (node.dist / dist_step),
    otherwise all edges are drawn of the same length.
    :param graph_name: str, the name of the web-page
    :param categories: a list of categories for the pie-charts inside the nodes
    :param tree: ete3.Tree
    :param out_html: path where to save the resulting html file.
    """
    json_dict = tree2json(tree, add_fake_nodes=layout == 'dagre')

    env = Environment(loader=PackageLoader('cytopast', 'templates'))
    template = env.get_template('pie_tree.js')
    graph = template.render(name2colour=[(name, Color(hue=i / len(categories), saturation=.8, luminance=.5).get_hex())
                                         for (i, name) in enumerate(categories, start=1)], elements=json_dict,
                            layout=layout)
    template = env.get_template('index.html')
    page = template.render(graph=graph, title=graph_name)

    with open(out_html, 'w+') as fp:
        fp.write(page)


def get_node(**data):
    return {DATA: data}


def get_edge(**data):
    return {DATA: data}
