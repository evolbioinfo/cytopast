import json
from queue import Queue

from colour import Color
from jinja2 import Environment, PackageLoader

DEFAULT_EDGE_SIZE = 6

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


def tree2json(tree, add_fake_nodes=True, categories=None, name_feature=STATE):
    nodes, edges = [], []
    features_to_keep = [SIZE, 'shape', FONT_SIZE]
    if categories:
        features_to_keep += categories

    sorted_categories = sorted(categories)

    categories = set(categories) if categories else set()
    if add_fake_nodes:
        max_dist = max(n.dist for n in tree.traverse())
        dist_step = max_dist / 10 if max_dist < 10 or max_dist > 100 else 1

    node2tooltip = {n: ' / '.join(cat for cat in sorted_categories if cat in n.features) for n in tree.traverse()}
    queue = Queue()
    queue.put(tree, block=False)
    node2id = {}
    i = 0
    while not queue.empty():
        n = queue.get(block=False)
        node2id[n] = i
        i += 1    
        for c in sorted(n.children, key=lambda c: (node2tooltip[c], str(getattr(c, name_feature, c.name)))):
            queue.put(c, block=False)
    
    for n, n_id in sorted(node2id.items(), key=lambda ni: ni[1]):
        if n == tree and add_fake_nodes and int(n.dist / dist_step) > 0:
            edge_name = getattr(n, 'edge_name', '%g' % n.dist)
            source_name = 'fake_node_{}'.format(n_id)
            nodes.append(get_node({ID: source_name, NAME: '', 'is_fake': 1, SIZE: 1, 'shape': 'rectangle',
                                   FONT_SIZE: 1}))
            edges.append(get_edge(**{SOURCE: source_name, TARGET: n_id, INTERACTION: 'triangle',
                                     SIZE: getattr(n, EDGE_SIZE, DEFAULT_EDGE_SIZE), NAME: edge_name}))
        features = {feature: getattr(n, feature) for feature in n.features if feature in features_to_keep}
        features[ID] = n_id
        features[NAME] = str(getattr(n, name_feature, n.name))
        if not (set(features.keys()) & categories):
            features[str(n.state)] = 100
        if SIZE not in n.features:
            features[SIZE] = DEFAULT_SIZE
        if FONT_SIZE not in n.features:
            features[FONT_SIZE] = 10
        if 'shape' not in n.features:
            features['shape'] = 'ellipse'
        features['tooltip'] = node2tooltip[n]
        nodes.append(get_node(features))
        for child in sorted(n.children, key=lambda c: node2id[c]):
            source_name = n_id
            edge_size = getattr(child, EDGE_SIZE, DEFAULT_EDGE_SIZE)
            edge_name = getattr(child, 'edge_name', '%g' % child.dist)
            if add_fake_nodes and int(child.dist / dist_step) > 0:
                target_name = 'fake_node_{}'.format(node2id[child])
                nodes.append(get_node({ID: target_name, NAME: '', 'is_fake': 1, SIZE: 1, 'shape': 'rectangle',
                                       FONT_SIZE: 1}))
                edges.append(get_edge(**{SOURCE: source_name, TARGET: target_name, INTERACTION: 'none',
                                         SIZE: edge_size, NAME: ''}))
                source_name = target_name
            edge_data = {SOURCE: source_name, TARGET: node2id[child], INTERACTION: 'triangle',
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


def save_as_cytoscape_html(tree, out_html, categories, graph_name='Untitled', layout='dagre', name_feature=STATE,
                           name2colour=None, add_fake_nodes=True):
    """
    Converts a tree to an html representation using Cytoscape.js.

    If categories are specified they are visualised as pie-charts inside the nodes,
    given that each node contains features corresponding to these categories with values being the percentage.
    For instance, given categories ['A', 'B', 'C'], a node with features {'A': 50, 'B': 50}
    will have a half-half pie-chart (half-colored in a colour of A, and half B).

    If dist_step is specified, the edges are rescaled accordingly to their dist (node.dist / dist_step),
    otherwise all edges are drawn of the same length.

    otherwise all edges are drawn of the same length.
    :param graph_name: str, the name of the web-page
    :param categories: a list of categories for the pie-charts inside the nodes
    :param tree: ete3.Tree
    :param out_html: path where to save the resulting html file.
    """
    json_dict = tree2json(tree, add_fake_nodes=add_fake_nodes, categories=categories, name_feature=name_feature)

    env = Environment(loader=PackageLoader('cytopast', 'templates'))
    template = env.get_template('pie_tree.js')
    if name2colour is None:
        name2colour = [(name, Color(hue=i / len(categories), saturation=.8, luminance=.5).get_hex())
                       for (i, name) in enumerate(categories, start=1)]
    name2colour = sorted(name2colour, key=lambda nk: nk[0])
    graph = template.render(name2colour=name2colour, elements=json_dict, layout=layout)
    template = env.get_template('index.html')
    page = template.render(graph=graph, title=graph_name)

    with open(out_html, 'w+') as fp:
        fp.write(page)


def get_node(data, position=None):
    res = {DATA: data}
    if position:
        res['position'] = position
    return res


def get_edge(**data):
    return {DATA: data}
