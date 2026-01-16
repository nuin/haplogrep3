"""Kulczynski distance metric (default for haplogrep)."""

from haplogrep3.models.polymorphism import Polymorphism
from .base import DistanceMetric


class KulczynskiDistance(DistanceMetric):
    """Kulczynski similarity coefficient.

    The Kulczynski measure is the default distance metric for haplogroup
    classification. It calculates:

        Q = (1/2) * (|A ∩ B| / |A| + |A ∩ B| / |B|)

    Where:
        A = sample polymorphisms
        B = expected polymorphisms for haplogroup
        |A ∩ B| = number of matching polymorphisms

    This can be simplified to:
        Q = |A ∩ B| * (|A| + |B|) / (2 * |A| * |B|)

    Or equivalently:
        Q = 0.5 * (found/expected + found/sample)

    Returns 1.0 for perfect match, 0.0 for no overlap.
    """

    @property
    def name(self) -> str:
        return "Kulczynski"

    def calculate(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
        weights: dict[int, float] | None = None,
    ) -> float:
        """Calculate Kulczynski similarity coefficient.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup
            weights: Optional position weights (not used in basic Kulczynski)

        Returns:
            Quality score between 0 and 1
        """
        if not expected_polymorphisms or not sample_polymorphisms:
            return 0.0

        found = self.get_found_polymorphisms(sample_polymorphisms, expected_polymorphisms)
        found_count = len(found)

        if found_count == 0:
            return 0.0

        sample_count = len(sample_polymorphisms)
        expected_count = len(expected_polymorphisms)

        # Kulczynski formula: average of two ratios
        # ratio1 = found / expected (what fraction of expected did we find)
        # ratio2 = found / sample (what fraction of sample matches expected)
        ratio1 = found_count / expected_count
        ratio2 = found_count / sample_count

        return 0.5 * (ratio1 + ratio2)
