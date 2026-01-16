"""Tests for classification task."""

import pytest
from haplogrep3.models import (
    Sample,
    Polymorphism,
    Phylotree,
    PhyloTreeNode,
    Haplogroup,
)
from haplogrep3.tasks import ClassificationTask
from haplogrep3.distance import Distance


def create_test_tree() -> Phylotree:
    """Create a simple test phylotree."""
    # Build a simple tree:
    # root
    #   └── H (263G, 750G)
    #       └── H1 (263G, 750G, 3010A)
    #           └── H1a (263G, 750G, 3010A, 3915A)

    h1a = PhyloTreeNode(
        haplogroup=Haplogroup(name="H1a"),
        polymorphisms=[Polymorphism.from_string("3915A")],
        children=[],
    )

    h1 = PhyloTreeNode(
        haplogroup=Haplogroup(name="H1"),
        polymorphisms=[Polymorphism.from_string("3010A")],
        children=[h1a],
    )
    h1a.parent = h1

    h = PhyloTreeNode(
        haplogroup=Haplogroup(name="H"),
        polymorphisms=[
            Polymorphism.from_string("263G"),
            Polymorphism.from_string("750G"),
        ],
        children=[h1],
    )
    h1.parent = h

    root = PhyloTreeNode(
        haplogroup=Haplogroup(name="root"),
        polymorphisms=[],
        children=[h],
    )
    h.parent = root

    return Phylotree(
        id="test-tree",
        name="Test Tree",
        root=root,
    )


class TestClassificationTask:
    def test_classify_exact_match(self):
        tree = create_test_tree()
        task = ClassificationTask(tree, Distance.KULCZYNSKI, hits=1)

        # Sample with exact H1a polymorphisms
        sample = Sample(
            id="test1",
            polymorphisms=[
                Polymorphism.from_string("263G"),
                Polymorphism.from_string("750G"),
                Polymorphism.from_string("3010A"),
                Polymorphism.from_string("3915A"),
            ],
        )

        results = task.classify([sample])

        assert len(results) == 1
        assert results[0].classification is not None
        assert results[0].classification.haplogroup.name == "H1a"
        assert results[0].classification.quality == 1.0

    def test_classify_partial_match(self):
        tree = create_test_tree()
        task = ClassificationTask(tree, Distance.KULCZYNSKI, hits=1)

        # Sample with H polymorphisms only
        sample = Sample(
            id="test2",
            polymorphisms=[
                Polymorphism.from_string("263G"),
                Polymorphism.from_string("750G"),
            ],
        )

        results = task.classify([sample])

        assert len(results) == 1
        assert results[0].classification is not None
        # Should match H best (perfect match for H)
        assert results[0].classification.haplogroup.name == "H"
        assert results[0].classification.quality == 1.0

    def test_classify_multiple_hits(self):
        tree = create_test_tree()
        task = ClassificationTask(tree, Distance.KULCZYNSKI, hits=3)

        sample = Sample(
            id="test3",
            polymorphisms=[
                Polymorphism.from_string("263G"),
                Polymorphism.from_string("750G"),
                Polymorphism.from_string("3010A"),
            ],
        )

        results = task.classify([sample])

        assert len(results) == 1
        result = results[0].classification
        assert result is not None
        # Top hit should be H1
        assert result.haplogroup.name == "H1"
        # Should have other hits
        assert len(result.other_results) >= 1

    def test_classify_no_match(self):
        tree = create_test_tree()
        task = ClassificationTask(tree, Distance.KULCZYNSKI, hits=1)

        # Sample with no matching polymorphisms
        sample = Sample(
            id="test4",
            polymorphisms=[
                Polymorphism.from_string("16519C"),
            ],
        )

        results = task.classify([sample])

        assert len(results) == 1
        result = results[0].classification
        assert result is not None
        # Quality should be low
        assert result.quality < 0.5
