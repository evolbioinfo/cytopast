import os
import unittest

import pandas as pd

from cytopast import read_tree
from pypastml import acr, EFT, JOINT

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')


class ACRStateTestJointEFT(unittest.TestCase):

    def setUp(self):
        self.feature = 'Country'
        df = pd.read_csv(STATES_INPUT, index_col=0, header=0)[[self.feature]]
        self.tree = read_tree(TREE_NWK)
        acr(self.tree, df, prediction_method=JOINT, model=EFT)

    def test_num_African_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if 'Africa' == getattr(node, self.feature):
                num += 1
        self.assertEqual(115, num, msg='Was supposed to have {} African nodes, got {}.'.format(115, num))

    def test_num_Albanian_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if 'Albania' == getattr(node, self.feature):
                num += 1
        self.assertEqual(50, num, msg='Was supposed to have {} Albanian nodes, got {}.'.format(50, num))

    def test_num_Greek_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if 'Greece' == getattr(node, self.feature):
                num += 1
        self.assertEqual(66, num, msg='Was supposed to have {} Greek nodes, got {}.'.format(66, num))

    def test_num_WE_nodes(self):
        num = 0
        for node in self.tree.traverse():
            if 'WestEurope' == getattr(node, self.feature):
                num += 1
        self.assertEqual(30, num, msg='Was supposed to have {} West European nodes, got {}.'.format(30, num))

    def test_state_root(self):
        expected_state = 'Africa'
        state = getattr(self.tree, self.feature)
        self.assertEqual(expected_state, state,
                         msg='Root state was supposed to be {}, got {}.'.format(expected_state, state))

    def test_state_resolved_internal_node_1(self):
        expected_state = 'Africa'
        for node in self.tree.traverse():
            if 'node_79' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break

    def test_state_resolved_internal_node_2(self):
        expected_state = 'Greece'
        for node in self.tree.traverse():
            if 'node_80' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break

    def test_state_resolved_internal_node_3(self):
        expected_state = 'WestEurope'
        for node in self.tree.traverse():
            if 'node_25' == node.name:
                state = getattr(node, self.feature)
                self.assertEqual(expected_state, state, msg='{} state was supposed to be {}, got {}.'
                                 .format(node.name, expected_state, state))
                break

    def test_state_resolved_internal_node_4(self):
        expected_state = 'Africa'
        for node in self.tree.traverse():
            if 'node_48' == node.name:
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
