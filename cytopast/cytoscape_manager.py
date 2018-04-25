import os
from queue import Queue
import numpy as np

from jinja2 import Environment, PackageLoader

from cytopast import METACHILD, MAX_NUM_TIPS_INSIDE, get_states, DATE, NODE_SIZE, FONT_SIZE, EDGE_NAME, NODE_NAME, \
    EDGE_SIZE, TIPS_INSIDE, TIPS_BELOW, sum_len_values

TOOLTIP = 'tooltip'

DEFAULT_FONT_SIZE = 10

COLOR = 'color'

SHAPE = 'shape'

FAKE_NODE_SHAPE = 'rectangle'

MERGED_NODE_SHAPE = 'ellipse'

INNER_NODE_SHAPE = 'rectangle'

TIP_SHAPE = 'ellipse'

DEFAULT_EDGE_SIZE = 10
DEFAULT_EDGE_COLOR = '#909090'
SPECIAL_EDGE_COLOR = '#383838'

DEFAULT_SIZE = 50

STATE = 'state'

DATA = 'data'
ID = 'id'
EDGES = 'edges'
NODES = 'nodes'
ELEMENTS = 'elements'

INTERACTION = 'interaction'

TARGET = 'target'
SOURCE = 'source'


def get_fake_node(n, n_id):
    attributes = {ID: n_id, 'is_fake': 1, SHAPE: FAKE_NODE_SHAPE}
    for feature in n.features:
        if feature.startswith(NODE_NAME):
            attributes[feature] = ''
        if feature.startswith(NODE_SIZE):
            attributes[feature] = 1
        if feature.startswith(FONT_SIZE):
            attributes[feature] = 1
    return _get_node(attributes)


def get_node(n, n_id, tooltip='', clazz=None):
    features = {feature: getattr(n, feature) for feature in n.features if feature in [SHAPE, DATE]
                or feature.startswith('node_')}
    features[ID] = n_id
    if SHAPE not in n.features:
        features[SHAPE] = TIP_SHAPE if n.is_leaf() and getattr(n, MAX_NUM_TIPS_INSIDE, 1) == 1 \
            else INNER_NODE_SHAPE if getattr(n, MAX_NUM_TIPS_INSIDE, 0) == 0 else MERGED_NODE_SHAPE
    features[TOOLTIP] = tooltip
    return _get_node(features, clazz=_clazz_list2css_class(clazz))


def get_edge(n, source_name, target_name, **kwargs):
    features = {SOURCE: source_name, TARGET: target_name, INTERACTION: 'triangle', COLOR: DEFAULT_EDGE_COLOR}
    features.update({feature: getattr(n, feature) for feature in n.features if feature.startswith('edge_')
                     or feature == DATE})
    features.update(kwargs)
    return _get_edge(**features)


def get_scaling_function(y_m, y_M, x_m, x_M):
    """
    Returns a linear function y = k x + b, where y \in [m, M]
    :param y_m:
    :param y_M:
    :param x_m:
    :param x_M:
    :return:
    """
    if x_M <= x_m:
        return lambda _: y_m
    k = (y_M - y_m) / (x_M - x_m)
    b = y_m - k * x_m
    return lambda _: int(k * _ + b)


def set_cyto_features(n, tips_inside, tips_below, size_scaling, e_size_scaling, font_scaling, transform_size,
                      transform_e_size, state, suffix=''):
    min_n_tips = min(sum_len_values(_) for _ in tips_inside) if tips_inside else 0
    max_n_tips = max(sum_len_values(_) for _ in tips_inside) if tips_inside else 0

    min_n_tips_below = min(sum_len_values(_) for _ in tips_below) if tips_below else 0
    max_n_tips_below = max(sum_len_values(_) for _ in tips_below) if tips_below else 0

    n.add_feature('{}{}'.format(NODE_NAME, suffix),
                  '{}{}{}'.format(state,
                                  (' {}'.format('{}-{}'.format(min_n_tips, max_n_tips)
                                                if min_n_tips != max_n_tips else min_n_tips))
                                  if max_n_tips > 0 else '',
                                  ' ({})'.format('{}-{}'.format(min_n_tips_below, max_n_tips_below)
                                                 if min_n_tips_below != max_n_tips_below else min_n_tips_below)
                                  if max_n_tips_below > 0 else ''))

    n.add_feature('{}{}'.format(NODE_SIZE, suffix), 20 if max_n_tips == 0 else size_scaling(transform_size(max_n_tips)))
    n.add_feature('{}{}'.format(FONT_SIZE, suffix), 10 if max_n_tips == 0 else font_scaling(transform_size(max_n_tips)))

    # If it is a compressed tree, set horizontal size as label
    if max_n_tips > 0 or max_n_tips_below > 0:
        edge_size = max(len(tips_inside), 1)
        n.add_feature('{}{}'.format(EDGE_NAME, suffix), str(edge_size) if edge_size != 1 else '')
        n.add_feature('{}{}'.format(EDGE_SIZE, suffix), e_size_scaling(transform_e_size(edge_size)))
    # else set dist as label
    else:
        n.add_feature('{}{}'.format(EDGE_NAME, suffix), '%g' % round(n.dist, 2))
        n.add_feature('{}{}'.format(EDGE_SIZE, suffix), e_size_scaling(transform_e_size(1)))


def _tree2json(tree, categories, add_fake_nodes, name_feature, n2tooltip, sort_key, min_date=0, max_date=0):
    e_size_scaling, font_scaling, size_scaling, transform_e_size, transform_size = get_size_transformations(tree)

    for n in tree.traverse():
        state = getattr(n, name_feature, '') if name_feature is not None else ''

        tips_inside, tips_below = getattr(n, TIPS_INSIDE, []), getattr(n, TIPS_BELOW, [])
        if isinstance(tips_inside, dict):
            tips_inside = [tips_inside]
        if isinstance(tips_below, dict):
            tips_below = [tips_below]
        set_cyto_features(n, tips_inside, tips_below, size_scaling, e_size_scaling, font_scaling, transform_size,
                          transform_e_size, state, suffix='')

        if min_date != max_date:
            for year in range(int(max_date), int(min_date) - 1, -1):
                tips_inside = [e for e in
                               ({k: _ for (k, _) in year2tips.items() if k <= year} for year2tips in tips_inside)
                               if e]
                tips_below = [e for e in
                              ({k: _ for (k, _) in year2tips.items() if k <= year} for year2tips in tips_below)
                              if e]
                set_cyto_features(n, tips_inside, tips_below, size_scaling, e_size_scaling, font_scaling,
                                  transform_size,
                                  transform_e_size, state, suffix='_{}'.format(year))

    clazzes = set()
    nodes, edges = [], []

    if add_fake_nodes:
        max_dist = max(n.dist for n in tree.traverse())
        dist_step = max_dist / 10 if max_dist < 10 or max_dist > 100 else 1

    node2tooltip = {n: n2tooltip(n, categories) for n in tree.traverse()}
    queue = Queue()
    queue.put(tree, block=False)
    node2id = {}
    i = 0
    while not queue.empty():
        n = queue.get(block=False)
        node2id[n] = i
        i += 1
        for c in sorted(n.children, key=lambda _: sort_key(_, name_feature, node2tooltip)):
            queue.put(c, block=False)

    for n, n_id in sorted(node2id.items(), key=lambda ni: ni[1]):
        if n == tree and add_fake_nodes and int(n.dist / dist_step) > 0:
            fake_id = 'fake_node_{}'.format(n_id)
            nodes.append(get_fake_node(n, fake_id))
            edges.append(get_edge(n, fake_id, n_id))

        clazz = tuple('{}_{}'.format(cat, getattr(n, cat)) for cat in categories if hasattr(n, cat))
        if clazz:
            clazzes.add(clazz)
        nodes.append(get_node(n, n_id, tooltip=node2tooltip[n], clazz=clazz))

        for child in sorted(n.children, key=lambda _: node2id[_]):
            edge_attributes = {COLOR: (SPECIAL_EDGE_COLOR if getattr(child, METACHILD, False) else DEFAULT_EDGE_COLOR)}
            source_name = n_id
            if add_fake_nodes and int(child.dist / dist_step) > 0:
                target_name = 'fake_node_{}'.format(node2id[child])
                nodes.append(get_fake_node(child, target_name))
                fake_edge_attributes = {k: '' for k in child.features if k.startswith(EDGE_NAME)}
                fake_edge_attributes.update(edge_attributes)
                fake_edge_attributes[INTERACTION] = 'none'
                edges.append(get_edge(child, source_name, target_name, **fake_edge_attributes))
                source_name = target_name
            if add_fake_nodes:
                edge_attributes['minLen'] = int(child.dist / dist_step)
            edges.append(get_edge(child, source_name, node2id[child], **edge_attributes))

    json_dict = {NODES: nodes, EDGES: edges}
    return json_dict, sorted(clazzes)


def get_size_transformations(tree):
    n_sizes = [getattr(n, MAX_NUM_TIPS_INSIDE) for n in tree.traverse() if getattr(n, MAX_NUM_TIPS_INSIDE, False)]
    max_size = max(n_sizes) if n_sizes else 1
    min_size = min(n_sizes) if n_sizes else 1
    need_log = max_size / min_size > 100
    transform_size = lambda _: np.power(np.log10(_ + 9) if need_log else _, 1 / 2)

    e_szs = [len(getattr(n, TIPS_INSIDE)) for n in tree.traverse() if getattr(n, TIPS_INSIDE, False)]
    max_e_size = max(e_szs) if e_szs else 1
    min_e_size = min(e_szs) if e_szs else 1
    need_e_log = max_e_size / min_e_size > 100
    transform_e_size = lambda _: np.log10(_) if need_e_log else _

    size_scaling = get_scaling_function(y_m=30, y_M=30 * min(8, int(max_size / min_size)),
                                        x_m=transform_size(min_size), x_M=transform_size(max_size))
    font_scaling = get_scaling_function(y_m=10, y_M=10 * min(3, int(max_size / min_size)),
                                        x_m=transform_size(min_size), x_M=transform_size(max_size))
    e_size_scaling = get_scaling_function(y_m=10, y_M=10 * min(3, int(max_e_size / min_e_size)),
                                          x_m=transform_e_size(min_e_size), x_M=transform_e_size(max_e_size))

    return e_size_scaling, font_scaling, size_scaling, transform_e_size, transform_size


def save_as_cytoscape_html(tree, out_html, categories, layout='dagre', name_feature=STATE,
                           name2colour=None, add_fake_nodes=True,
                           n2tooltip=lambda n, categories: ', '.join(get_states(n, categories)),
                           sort_key=lambda n, name_feature, node2tooltip:
                           (str(getattr(n, name_feature, '')) if name_feature else '', node2tooltip[n],
                            -getattr(n, NODE_SIZE, 0), n.name), min_date=0, max_date=0):
    """
    Converts a tree to an html representation using Cytoscape.js.

    If categories are specified they are visualised as pie-charts inside the nodes,
    given that each node contains features corresponding to these categories with values being the percentage.
    For instance, given categories ['A', 'B', 'C'], a node with features {'A': 50, 'B': 50}
    will have a half-half pie-chart (half-colored in a colour of A, and half B).

    If dist_step is specified, the edges are rescaled accordingly to their dist (node.dist / dist_step),
    otherwise all edges are drawn of the same length.

    otherwise all edges are drawn of the same length.
    :param name_feature: str, a node feature whose value will be used as a label
    :param sort_key: a function, that given a tree node, the name_feature and a node to tooltip dict, 
    returns a key to be used for sorting nodes on the same level in the tree.
    :param n2tooltip: dict, TreeNode to str mapping tree nodes to tooltips.
    :param add_fake_nodes: bool, if to add fake nodes, needed for showing branch lengths. 
    :param layout: str, name of the layout for Cytoscape.js 
    :param name2colour: dict, str to str, category name to HEX colour mapping 
    :param categories: a list of categories for the pie-charts inside the nodes
    :param tree: ete3.Tree
    :param out_html: path where to save the resulting html file.
    """
    graph_name = os.path.splitext(os.path.basename(out_html))[0]

    json_dict, clazzes \
        = _tree2json(tree, categories=categories, add_fake_nodes=add_fake_nodes, name_feature=name_feature,
                     n2tooltip=n2tooltip, sort_key=sort_key, min_date=min_date, max_date=max_date)
    env = Environment(loader=PackageLoader('cytopast'))
    template = env.get_template('pie_tree.js')

    clazz2css = {}
    for clazz_list in clazzes:
        n = len(clazz_list)
        css = ''
        for i, cat in enumerate(clazz_list, start=1):
            css += """
                'pie-{i}-background-color': "{colour}",
                'pie-{i}-background-size': '{percent}\%',
            """.format(i=i, percent=round(100 / n, 2), colour=name2colour[cat])
        clazz2css[_clazz_list2css_class(clazz_list)] = css
    graph = template.render(clazz2css=clazz2css.items(), elements=json_dict, layout=layout, title=graph_name,
                            min_date=min_date, max_date=max_date)
    slider = env.get_template('time_slider.html').render(min_date=min_date, max_date=max_date) \
        if min_date != max_date else ''

    template = env.get_template('index.html')
    page = template.render(graph=graph, title=graph_name, slider=slider)

    os.makedirs(os.path.abspath(os.path.dirname(out_html)), exist_ok=True)
    with open(out_html, 'w+') as fp:
        fp.write(page)


def _clazz_list2css_class(clazz_list):
    if not clazz_list:
        return None
    return ''.join(c for c in '-'.join(clazz_list) if c.isalnum() or '-' == c)


def _get_node(data, position=None, clazz=None):
    res = {DATA: data}
    if clazz:
        res['classes'] = clazz
    if position:
        res['position'] = position
    return res


def _get_edge(**data):
    return {DATA: data}
