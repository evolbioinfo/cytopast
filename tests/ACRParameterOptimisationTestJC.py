import os
import unittest

import numpy as np
import pandas as pd

from cytopast import read_tree
from pypastml import acr, MPPA, JC, get_personalized_feature_name, LH, LH_SF

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')


class ACRParameterOptimisationTestJC(unittest.TestCase):

    def setUp(self):
        self.feature = 'Country'
        df = pd.read_csv(STATES_INPUT, index_col=0, header=0)[[self.feature]]
        self.tree = read_tree(TREE_NWK)
        self.acr_result = acr(self.tree, df, prediction_method=MPPA, model=JC)[0]

    def test_likelihood(self):
        self.assertAlmostEqual(self.acr_result.likelihood, -121.873, places=3,
                               msg='Likelihood was supposed to be the {:.3f}, got {:3f}'
                               .format(-121.873, self.acr_result.likelihood))

    def test_sf(self):
        self.assertAlmostEqual(self.acr_result.sf, 4.951, places=3,
                               msg='SF was supposed to be the {:.3f}, got {:3f}'
                               .format(4.951, self.acr_result.sf))

    def test_frequencies(self):
        value = self.acr_result.frequencies
        expected_value = np.ones(len(value), np.float64) / len(value)
        self.assertListEqual(value.tolist(), expected_value.tolist(),
                             msg='Frequencies were supposed to be the {}, got {}'.format(expected_value, value))

    def test_frequencies_sum_to_1(self):
        value = self.acr_result.frequencies.sum()
        self.assertAlmostEqual(value, 1, places=3,
                               msg='Frequencies were supposed to sum to 1, not to {:3f}'.format(value))

    def test_likelihood_same_for_all_nodes(self):
        """
        Tests if marginal likelihoods were correctly calculated
        by comparing the likelihoods of all the nodes (should be all the same).
        """
        lh_feature = get_personalized_feature_name(self.feature, LH)
        lh_sf_feature = get_personalized_feature_name(self.feature, LH_SF)

        for node in self.tree.traverse():
            if not node.is_root() and not (node.is_leaf() and node.dist == 0):
                node_loglh = np.log10(getattr(node, lh_feature).sum()) - getattr(node, lh_sf_feature)
                parent_loglh = np.log10(getattr(node.up, lh_feature).sum()) - getattr(node.up, lh_sf_feature)
                self.assertAlmostEqual(node_loglh, parent_loglh, places=2,
                                       msg='Likelihoods of {} and {} were supposed to be the same.'
                                       .format(node.name, node.up.name))
