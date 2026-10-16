import unittest
import numpy as np
from helper_functions import neighbor_joining, Node
from test_helper_functions import is_isomorhipic

class TestNeighborJoining(unittest.TestCase):

    def test_two_taxa(self):
        """Two taxa should join directly to the root."""
        D = np.array([
            [0, 4],
            [4, 0]
        ])
        labels = ["A", "B"]

        root = neighbor_joining(D, labels)

        self.assertEqual(root.name, "root")
        self.assertEqual(root.left.name, "A")
        self.assertEqual(root.right.name, "B")
        self.assertAlmostEqual(root.left_distance, 2.0)
        self.assertAlmostEqual(root.right_distance, 2.0)

    def test_three_taxa(self):
        """Check the tree and branch lengths for three taxa."""
        distances = np.array(
            [
                [0, 5, 9],
                [5, 0, 10],
                [9, 10, 0],
            ],
            dtype=float,
        )

        result = neighbor_joining(distances, list("ABC"))

        # Expected topology: A and B join first, then C
        expected = Node(
            "root",
            left=Node("C", None, 0, None, 0),
            left_distance=3.5,
            right=Node(
                "(A,B)",
                left=Node("A", None, 0, None, 0),
                left_distance=2.0,
                right=Node("B", None, 0, None, 0),
                right_distance=3.0,
            ),
            right_distance=3.5,
        )

        self.assertTrue(is_isomorhipic(result, expected))

    def test_four_taxa(self):
        """Check a four-taxon tree"""
        distances = np.array(
            [
                [0, 5, 9, 9],
                [5, 0, 10, 10],
                [9, 10, 0, 8],
                [9, 10, 8, 0],
            ],
            dtype=float,
        )

        result = neighbor_joining(distances, list("ABCD"))

        expected = Node(
            "root",
            left=Node(
                "(A,B)",
                left=Node("A", None, 0, None, 0),
                left_distance=2.0,
                right=Node("B", None, 0, None, 0),
                right_distance=3.0,
            ),
            left_distance=1.5,
            right=Node(
                "(C,D)",
                left=Node("C", None, 0, None, 0),
                left_distance=4.0,
                right=Node("D", None, 0, None, 0),
                right_distance=4.0,
            ),
            right_distance=1.5,
        )

        self.assertTrue(is_isomorhipic(result, expected))

    def test_tree_is_binary(self):
        """Every internal node should have exactly two children."""
        D = np.array([
            [0, 5, 9, 9],
            [5, 0, 10, 10],
            [9, 10, 0, 8],
            [9, 10, 8, 0]
        ])
        labels = ["A", "B", "C", "D"]

        root = neighbor_joining(D, labels)

        def check(node):
            if node.left is None and node.right is None:
                return

            self.assertIsNotNone(node.left)
            self.assertIsNotNone(node.right)

            check(node.left)
            check(node.right)

        check(root)

    def test_invalid_input(self):
        """A distance matrix with fewer than two taxa should fail."""
        with self.assertRaises(ValueError):
            neighbor_joining(np.array([[0]]), ["A"])

    def test_two_taxa_unequal_distance(self):
        """Branch lengths should each be half the pairwise distance."""
        D = np.array([
            [0, 3],
            [3, 0]
        ])
        root = neighbor_joining(D, ["X", "Y"])

        self.assertAlmostEqual(root.left_distance, 1.5)
        self.assertAlmostEqual(root.right_distance, 1.5)

    def test_leaf_labels_preserved(self):
        """All input labels should appear exactly once."""
        D = np.array([
            [0, 5, 9],
            [5, 0, 10],
            [9, 10, 0]
        ])
        labels = ["A", "B", "C"]

        root = neighbor_joining(D, labels)

        found = []

        def collect(node):
            if node.left is None and node.right is None:
                found.append(node.name)
                return

            collect(node.left)
            collect(node.right)

        collect(root)

        self.assertCountEqual(found, labels)
        self.assertEqual(len(found), len(set(found)))


if __name__ == "__main__":
    unittest.main()