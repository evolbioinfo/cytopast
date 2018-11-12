import os
import unittest

import numpy as np
import pandas as pd

from cytopast import read_tree
from pypastml import get_personalized_feature_name
from pypastml.acr import acr
from pypastml.ml import LH, LH_SF, MPPA, EFT, MARGINAL_PROBS

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')


class ACRParameterOptimisationTestEFT(unittest.TestCase):

    def setUp(self):
        self.feature = 'Country'
        df = pd.read_csv(STATES_INPUT, index_col=0, header=0)[[self.feature]]
        self.tree = read_tree(TREE_NWK)
        self.acr_result = acr(self.tree, df, prediction_method=MPPA, model=EFT)[0]

    def test_likelihood(self):
        self.assertAlmostEqual(-123.173, self.acr_result.likelihood, places=3,
                               msg='Likelihood was supposed to be the {:.3f}, got {:3f}'
                               .format(-123.173, self.acr_result.likelihood))

    def test_sf(self):
        self.assertAlmostEqual(5.38, self.acr_result.sf, places=3,
                               msg='SF was supposed to be the {:.3f}, got {:3f}'
                               .format(5.38, self.acr_result.sf))

    def test_frequencies(self):
        for loc, expected_value in {'Africa': 0.325, 'Albania': 0.208, 'EastEurope': 0.065,
                                    'Greece': 0.26, 'WestEurope': 0.143}.items():
            value = self.acr_result.frequencies[np.where(self.acr_result.states == loc)][0]
            self.assertAlmostEqual(value, expected_value, places=3,
                                   msg='Frequency of {} was supposed to be the {:.3f}, got {:3f}'
                                   .format(loc, expected_value, value))

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

    def test_marginal_probs_root(self):
        mps = getattr(self.tree, get_personalized_feature_name(self.feature, MARGINAL_PROBS))
        expected_values = {'Africa': 0.70653324, 'Albania': 0.03733449, 'EastEurope': 0.04927248,
                           'Greece': 0.04667803, 'WestEurope': 0.16018175}
        for loc, expected_value in expected_values.items():
            value = mps[np.where(self.acr_result.states == loc)][0]
            self.assertAlmostEqual(value, expected_value, places=3,
                                   msg='{}: Marginal probability of {} was supposed to be the {:.3f}, got {:3f}'
                                   .format(self.tree.name, loc, expected_value, value))

    def test_marginal_probs_internal_node(self):
        for node in self.tree.traverse():
            if 'node_4' == node.name:
                mps = getattr(node, get_personalized_feature_name(self.feature, MARGINAL_PROBS))
                expected_values = {'Africa': 0.43054907, 'Albania': 0.00376214, 'EastEurope': 0.00185369,
                                   'Greece': 0.00470286, 'WestEurope': 0.55913223}
                for loc, expected_value in expected_values.items():
                    value = mps[np.where(self.acr_result.states == loc)][0]
                    self.assertAlmostEqual(value, expected_value, places=3,
                                           msg='{}: Marginal probability of {} was supposed to be the {:.3f}, got {:3f}'
                                           .format(node.name, loc, expected_value, value))
                break

    def test_marginal_probs_tip(self):
        for node in self.tree:
            if '02ALAY1660' == node.name:
                mps = getattr(node, get_personalized_feature_name(self.feature, MARGINAL_PROBS))
                expected_values = {'Africa': 0, 'Albania': 1, 'EastEurope': 0, 'Greece': 0, 'WestEurope': 0}
                for loc, expected_value in expected_values.items():
                    value = mps[np.where(self.acr_result.states == loc)][0]
                    self.assertAlmostEqual(value, expected_value, places=3,
                                           msg='{}: Marginal probability of {} was supposed to be the {:.3f}, got {:3f}'
                                           .format(node.name, loc, expected_value, value))
                break

