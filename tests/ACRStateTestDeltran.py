import os
import unittest

import pandas as pd

from cytopast import read_tree
from pypastml.acr import acr
from pypastml.parsimony import DELTRAN

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')


class ACRStateTestDeltran(unittest.TestCase):

    def setUp(self):
        self.feature = 'Country'
        df = pd.read_csv(STATES_INPUT, index_col=0, header=0)[[self.feature]]
        self.tree = read_tree(TREE_NWK)
        self.acr_result = acr(self.tree, df, prediction_method=DELTRAN)[0]

    def test_num_steps(self):
        self.assertEqual(32, self.acr_result.steps,
                         msg='Was supposed to have {} parsimonious steps, got {}.'.format(32, self.acr_result.steps))

    def test_state_root(self):
        expected_state = 'Africa'
        state = getattr(self.tree, self.feature)
        self.assertEqual(expected_state, state,
                         msg='Root state was supposed to be {}, got {}.'.format(expected_state, state))

    def test_state_resolved_node_129(self):
        expected_state = 'Greece'
        for node in self.tree.traverse():
            if 'node_129' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break

    def test_state_resolved_node_25(self):
        expected_state = 'Africa'
        for node in self.tree.traverse():
            if 'node_25' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break

    def test_state_resolved_node_21(self):
        expected_state = 'Africa'
        for node in self.tree.traverse():
            if 'node_21' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                    .format(node.name, expected_state, state))
                break

    def test_state_unresolved_node_32(self):
        expected_state = {'WestEurope', 'Greece'}
        for node in self.tree.traverse():
            if 'node_32' == node.name:
                state = set(getattr(node, self.feature))
                self.assertSetEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                    .format(node.name, expected_state, state))
                break

    def test_state_zero_tip(self):
        expected_state = 'Albania'
        for node in self.tree.traverse():
            if '01ALAY1715' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break

    def test_state_tip(self):
        expected_state = 'WestEurope'
        for node in self.tree:
            if '94SEAF9671' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break
