import os
import unittest

import pandas as pd

from cytopast import read_tree
from pypastml import acr, MPPA, F81

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')


class ACRStateTestMPPAF81(unittest.TestCase):

    def setUp(self):
        self.feature = 'Country'
        df = pd.read_csv(STATES_INPUT, index_col=0, header=0)[[self.feature]]
        self.tree = read_tree(TREE_NWK)
        acr(self.tree, df, prediction_method=MPPA, model=F81)

    def test_state_root(self):
        expected_state = 'Africa'
        state = getattr(self.tree, self.feature)
        self.assertEqual(expected_state, state,
                         msg='Root state was supposed to be {}, got {}.'.format(expected_state, state))

    def test_num_unresolved_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if isinstance(getattr(node, self.feature), list):
                num += 1
        self.assertEqual(6, num, msg='Was supposed to have {} unresolved nodes, got {}.'.format(6, num))

    def test_num_African_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if not isinstance(getattr(node, self.feature), list) and 'Africa' == getattr(node, self.feature):
                num += 1
        self.assertEqual(113, num, msg='Was supposed to have {} African nodes, got {}.'.format(113, num))

    def test_num_Albanian_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if not isinstance(getattr(node, self.feature), list) and 'Albania' == getattr(node, self.feature):
                num += 1
        self.assertEqual(50, num, msg='Was supposed to have {} Albanian nodes, got {}.'.format(50, num))

    def test_num_Greek_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if not isinstance(getattr(node, self.feature), list) and 'Greece' == getattr(node, self.feature):
                num += 1
        self.assertEqual(65, num, msg='Was supposed to have {} Greek nodes, got {}.'.format(65, num))

    def test_num_WE_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if not isinstance(getattr(node, self.feature), list) and 'WestEurope' == getattr(node, self.feature):
                num += 1
        self.assertEqual(27, num, msg='Was supposed to have {} West European nodes, got {}.'.format(27, num))

    def test_state_unresolved_internal_node(self):
        expected_state = {'Africa', 'Greece'}
        for node in self.tree.traverse():
            if 'node_79' == node.name:
                state = set(getattr(node, self.feature))
                self.assertSetEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                    .format(node.name, expected_state, state))
                break

    def test_state_node_32(self):
        expected_state = {'WestEurope', 'Greece'}
        for node in self.tree.traverse():
            if 'node_32' == node.name:
                state = set(getattr(node, self.feature))
                self.assertSetEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                    .format(node.name, expected_state, state))
                break

    def test_state_resolved_internal_node(self):
        expected_state = 'Greece'
        for node in self.tree.traverse():
            if 'node_80' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
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
