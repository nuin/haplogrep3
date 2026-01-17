"""Kimura 2-parameter distance metric."""

from haplogrep3.models.polymorphism import Polymorphism, Mutation
from .base import DistanceMetric


# Transition pairs (purine-purine or pyrimidine-pyrimidine)
TRANSITIONS = {
    (Mutation.A, Mutation.G),
    (Mutation.G, Mutation.A),
    (Mutation.C, Mutation.T),
    (Mutation.T, Mutation.C),
}

# Transversion pairs (purine-pyrimidine)
# All other combinations are transversions


def is_transition(ref: str, alt: Mutation) -> bool:
    """Check if a mutation is a transition.

    Transitions: A<->G, C<->T (purine-purine or pyrimidine-pyrimidine)
    Transversions: all other changes

    Args:
        ref: Reference nucleotide
        alt: Alternate nucleotide (Mutation)

    Returns:
        True if transition, False if transversion
    """
    try:
        ref_mut = Mutation(ref.upper())
        return (ref_mut, alt) in TRANSITIONS
    except ValueError:
        return False


class KimuraDistance(DistanceMetric):
    """Kimura 2-parameter distance-based similarity.

    The Kimura 2-parameter model accounts for different rates of
    transitions (A<->G, C<->T) and transversions (all others).

    Transitions are more common than transversions, so they are
    weighted less heavily in the distance calculation.

    For haplogroup classification, we adapt this to:
        - Missing transitions are less penalized than missing transversions
        - The quality score favors matches and penalizes mismatches
    """

    # Weight for transversions (higher = more penalty)
    TRANSVERSION_WEIGHT = 2.0
    # Weight for transitions
    TRANSITION_WEIGHT = 1.0

    @property
    def name(self) -> str:
        return "Kimura"

    def calculate(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
        weights: dict[str, float] | None = None,
        hotspots: set[str] | None = None,
    ) -> float:
        """Calculate Kimura-weighted similarity.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup
            weights: Optional position weights

        Returns:
            Quality score between 0 and 1
        """
        if not expected_polymorphisms and not sample_polymorphisms:
            return 1.0

        if not expected_polymorphisms or not sample_polymorphisms:
            return 0.0

        found = self.get_found_polymorphisms(sample_polymorphisms, expected_polymorphisms)
        missing = self.get_missing_polymorphisms(sample_polymorphisms, expected_polymorphisms)

        # Calculate weighted scores
        found_score = 0.0
        total_expected_weight = 0.0

        for poly in expected_polymorphisms:
            weight = self._get_weight(poly)
            total_expected_weight += weight
            if poly in found:
                found_score += weight

        if total_expected_weight == 0:
            return 0.0

        return found_score / total_expected_weight

    def _get_weight(self, poly: Polymorphism) -> float:
        """Get the weight for a polymorphism based on mutation type.

        Args:
            poly: Polymorphism to weight

        Returns:
            Weight value
        """
        if poly.reference and is_transition(poly.reference, poly.mutation):
            return self.TRANSITION_WEIGHT
        return self.TRANSVERSION_WEIGHT
