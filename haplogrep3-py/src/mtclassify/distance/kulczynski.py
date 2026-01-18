"""Kulczynski distance metric (default for haplogrep)."""

from mtclassify.models.polymorphism import Polymorphism
from .base import DistanceMetric


class KulczynskiDistance(DistanceMetric):
    """Kulczynski similarity coefficient.

    The Kulczynski measure is the default distance metric for haplogroup
    classification. When weights are provided, it calculates:

        Q = 0.5 * (weightFound / weightSample + weightFound / weightExpected)

    Where:
        weightFound = sum of weights for matching polymorphisms
        weightSample = sum of weights for all sample polymorphisms
        weightExpected = sum of weights for all expected polymorphisms

    Without weights, it falls back to unweighted counts:
        Q = 0.5 * (found/sample + found/expected)

    Returns 1.0 for perfect match, 0.0 for no overlap.
    """

    @property
    def name(self) -> str:
        return "Kulczynski"

    def calculate(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
        weights: dict[str, float] | None = None,
        hotspots: set[str] | None = None,
    ) -> float:
        """Calculate Kulczynski similarity coefficient.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup
            weights: Optional polymorphism weights (poly string -> weight)
            hotspots: Optional set of hotspot polymorphism strings to exclude from sample weight

        Returns:
            Quality score between 0 and 1
        """
        if not expected_polymorphisms or not sample_polymorphisms:
            return 0.0

        found = self.get_found_polymorphisms(sample_polymorphisms, expected_polymorphisms)

        if not found:
            return 0.0

        if weights:
            # Weighted Kulczynski (matching Java haplogrep behavior)
            hotspots = hotspots or set()
            weight_found = sum(weights.get(str(p), 1.0) for p in found)
            # Exclude hotspots from sample weight calculation (Java behavior)
            weight_sample = sum(
                weights.get(str(p), 1.0)
                for p in sample_polymorphisms
                if str(p) not in hotspots
            )
            weight_expected = sum(weights.get(str(p), 1.0) for p in expected_polymorphisms)

            if weight_sample == 0 or weight_expected == 0:
                return 0.0

            # ratio1 = found weight / sample weight (what fraction of sample matches)
            # ratio2 = found weight / expected weight (what fraction of expected did we find)
            ratio1 = weight_found / weight_sample
            ratio2 = weight_found / weight_expected

            return 0.5 * ratio1 + 0.5 * ratio2
        else:
            # Unweighted fallback
            found_count = len(found)
            sample_count = len(sample_polymorphisms)
            expected_count = len(expected_polymorphisms)

            ratio1 = found_count / sample_count
            ratio2 = found_count / expected_count

            return 0.5 * (ratio1 + ratio2)
