"""Jaccard distance metric."""

from haplogrep3.models.polymorphism import Polymorphism
from .base import DistanceMetric


class JaccardDistance(DistanceMetric):
    """Jaccard similarity coefficient.

    Calculates:

        Q = |A ∩ B| / |A ∪ B|

    Where:
        A = sample polymorphisms
        B = expected polymorphisms for haplogroup

    Returns 1.0 for identical sets, 0.0 for no overlap.
    """

    @property
    def name(self) -> str:
        return "Jaccard"

    def calculate(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
        weights: dict[int, float] | None = None,
    ) -> float:
        """Calculate Jaccard similarity coefficient.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup
            weights: Optional position weights (not used)

        Returns:
            Quality score between 0 and 1
        """
        if not expected_polymorphisms and not sample_polymorphisms:
            return 1.0

        if not expected_polymorphisms or not sample_polymorphisms:
            return 0.0

        intersection = sample_polymorphisms & expected_polymorphisms
        union = sample_polymorphisms | expected_polymorphisms

        if not union:
            return 1.0

        return len(intersection) / len(union)
