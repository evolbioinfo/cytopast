import logging

from cytopast.cytoscape_manager import save_as_cytoscape_html

from cytopast.colour_generator import get_enough_colours, WHITE

from cytopast import REASONABLE_NUMBER_OF_TIPS, DATE, compress_tree


def visualize(tree, columns, name_column=None, html=None, html_compressed=None,
              tip_size_threshold=REASONABLE_NUMBER_OF_TIPS, min_date=0, max_date=0):
    one_column = len(columns) == 1

    column2values = {}
    for feature in columns:
        column2values[feature] = annotate(tree, feature, unique=one_column)

    if not name_column and one_column:
        name_column = columns[0]

    name2colour = {}
    for cat in columns:
        unique_values = column2values[cat]
        num_unique_values = len(unique_values)
        colours = get_enough_colours(num_unique_values)
        for value, col in zip(unique_values, colours):
            name2colour['{}_{}'.format(value, True) if one_column else '{}_{}'.format(cat, value)] = col
        logging.info('Mapped states to colours for {} as following: {} -> {}'.format(cat, unique_values, colours))
        # let ambiguous values be white
        if not one_column:
            name2colour['{}_'.format(cat)] = WHITE

    # set internal node dates to min of its tips' dates
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            if not hasattr(n, DATE):
                n.add_feature(DATE, 0)
        else:
            n.add_feature(DATE, min(getattr(_, DATE) for _ in n))

    def get_category_str(n):
        if one_column:
            return '{}: {}'.format(columns[0], ' or '.join('{}'.format(_)
                                                           for _ in column2values[columns[0]]
                                                           if hasattr(n, _) and getattr(n, _, '') != ''))
        return '<br>'.join('{}: {}'.format(_, getattr(n, _))
                           for _ in columns if hasattr(n, _) and getattr(n, _, '') != '')

    if html:
        save_as_cytoscape_html(tree, html, categories=column2values[columns[0]] if one_column else columns,
                               name2colour=name2colour,
                               n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               name_feature='name', min_date=min_date, max_date=max_date, is_compressed=False)

    if html_compressed:
        tree = compress_tree(tree, categories=column2values[columns[0]] if one_column else columns,
                             tip_size_threshold=tip_size_threshold)
        save_as_cytoscape_html(tree, html_compressed, categories=column2values[columns[0]] if one_column else columns,
                               name2colour=name2colour, n2tooltip={n: get_category_str(n) for n in tree.traverse()},
                               min_date=min_date, max_date=max_date, name_feature=name_column, is_compressed=True)
    return tree


def annotate(tree, feature, unique=True):
    all_states = set()
    for node in tree.traverse():
        possible_states = getattr(node, feature)
        if isinstance(possible_states, list):
            node.add_feature(feature, '')
        else:
            all_states.add(possible_states)
        if unique:
            for state in (possible_states if isinstance(possible_states, list) else [possible_states]):
                node.add_feature(state, True)
                all_states.add(state)
    return sorted(all_states)
