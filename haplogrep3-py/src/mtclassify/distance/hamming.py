"""Hamming distance metric."""

from mtclassify.models.polymorphism import Polymorphism
from .base import DistanceMetric


class HammingDistance(DistanceMetric):
    """Hamming distance-based similarity.

    Calculates similarity based on the number of differences:

        Q = 1 - (|missing| + |extra|) / max(|sample|, |expected|)

    Where:
        missing = expected polymorphisms not in sample
        extra = sample polymorphisms not in expected

    Returns 1.0 for perfect match, lower for more differences.
    """

    @property
    def name(self) -> str:
        return "Hamming"

    def calculate(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
        weights: dict[str, float] | None = None,
        hotspots: set[str] | None = None,
    ) -> float:
        """Calculate Hamming-based similarity.

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

        missing = self.get_missing_polymorphisms(sample_polymorphisms, expected_polymorphisms)
        extra = self.get_remaining_polymorphisms(sample_polymorphisms, expected_polymorphisms)

        total_differences = len(missing) + len(extra)
        max_size = max(len(sample_polymorphisms), len(expected_polymorphisms))

        if max_size == 0:
            return 1.0

        # Convert distance to similarity
        return max(0.0, 1.0 - (total_differences / max_size))
