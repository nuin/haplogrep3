"""Tests for distance metrics."""

import pytest
from mtclassify.models import Polymorphism
from mtclassify.distance import (
    KulczynskiDistance,
    HammingDistance,
    JaccardDistance,
    Distance,
    get_distance_metric,
)


def make_poly_set(*positions: int) -> set[Polymorphism]:
    """Helper to create a set of polymorphisms."""
    return {Polymorphism(position=p, mutation=Polymorphism.from_string(f"{p}G").mutation) for p in positions}


class TestKulczynskiDistance:
    def test_perfect_match(self):
        metric = KulczynskiDistance()
        sample = make_poly_set(263, 315, 750)
        expected = make_poly_set(263, 315, 750)

        score = metric.calculate(sample, expected)
        assert score == 1.0

    def test_no_overlap(self):
        metric = KulczynskiDistance()
        sample = make_poly_set(263, 315)
        expected = make_poly_set(750, 1438)

        score = metric.calculate(sample, expected)
        assert score == 0.0

    def test_partial_match(self):
        metric = KulczynskiDistance()
        sample = make_poly_set(263, 315, 750)
        expected = make_poly_set(263, 315, 1438, 4769)

        score = metric.calculate(sample, expected)
        # found=2, sample=3, expected=4
        # ratio1 = 2/4 = 0.5, ratio2 = 2/3 = 0.667
        # score = 0.5 * (0.5 + 0.667) = 0.5833
        assert 0.58 < score < 0.59

    def test_empty_sets(self):
        metric = KulczynskiDistance()
        assert metric.calculate(set(), set()) == 0.0
        assert metric.calculate(make_poly_set(263), set()) == 0.0
        assert metric.calculate(set(), make_poly_set(263)) == 0.0


class TestHammingDistance:
    def test_perfect_match(self):
        metric = HammingDistance()
        sample = make_poly_set(263, 315, 750)
        expected = make_poly_set(263, 315, 750)

        score = metric.calculate(sample, expected)
        assert score == 1.0

    def test_no_overlap(self):
        metric = HammingDistance()
        sample = make_poly_set(263, 315)
        expected = make_poly_set(750, 1438)

        score = metric.calculate(sample, expected)
        # missing=2, extra=2, max_size=2
        # score = 1 - 4/2 = negative, clamped to 0
        assert score == 0.0


class TestJaccardDistance:
    def test_perfect_match(self):
        metric = JaccardDistance()
        sample = make_poly_set(263, 315, 750)
        expected = make_poly_set(263, 315, 750)

        score = metric.calculate(sample, expected)
        assert score == 1.0

    def test_no_overlap(self):
        metric = JaccardDistance()
        sample = make_poly_set(263, 315)
        expected = make_poly_set(750, 1438)

        score = metric.calculate(sample, expected)
        assert score == 0.0

    def test_partial_overlap(self):
        metric = JaccardDistance()
        sample = make_poly_set(263, 315, 750)
        expected = make_poly_set(263, 315, 1438, 4769)

        score = metric.calculate(sample, expected)
        # intersection=2, union=5
        # score = 2/5 = 0.4
        assert score == 0.4


class TestGetDistanceMetric:
    def test_get_kulczynski(self):
        metric = get_distance_metric(Distance.KULCZYNSKI)
        assert isinstance(metric, KulczynskiDistance)

    def test_get_hamming(self):
        metric = get_distance_metric(Distance.HAMMING)
        assert isinstance(metric, HammingDistance)

    def test_get_jaccard(self):
        metric = get_distance_metric(Distance.JACCARD)
        assert isinstance(metric, JaccardDistance)
