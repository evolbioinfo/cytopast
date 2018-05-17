import colorsys

NUM2COLOURS = {
    1: ['#1b9e77'],
    2: ['#d95f02', '#1b9e77'],
    3: ['#a6cee3', '#1f78b4', '#b2df8a'],
    4: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3'],
    5: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00'],
    6: ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f'],
    7: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f'],
    8: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00'],
    9: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6'],
    10: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'],
    11: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
         '#ffff99'],
    12: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
         '#ffff99', '#b15928']
}

WHITE = '#ffffff'


def get_enough_colours(num_unique_values):
    """
    Generates and returns an array of `num_unique_values` HEX colours.
    :param num_unique_values: int, number of colours to be generated.
    :return: array of str, containing colours in HEX format.
    """
    if num_unique_values in NUM2COLOURS:
        return NUM2COLOURS[num_unique_values]
    return ['#%02x%02x%02x' % tuple(rgb) for rgb in
            (map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*hsv))
             for hsv in ((_ / num_unique_values, 0.5, 0.5) for _ in range(1, num_unique_values + 1)))]
