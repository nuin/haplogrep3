"""Base class for distance metrics."""

from abc import ABC, abstractmethod
from typing import Protocol

from mtclassify.models.polymorphism import Polymorphism


class DistanceMetric(ABC):
    """Abstract base class for distance/similarity metrics.

    Distance metrics calculate how well a sample's polymorphisms
    match a haplogroup's expected polymorphisms.
    """

    @abstractmethod
    def calculate(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
        weights: dict[str, float] | None = None,
        hotspots: set[str] | None = None,
    ) -> float:
        """Calculate the quality/similarity score.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup
            weights: Optional polymorphism weights (poly string -> weight)
            hotspots: Optional set of hotspot polymorphism strings

        Returns:
            Quality score between 0 and 1 (1 = perfect match)
        """
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name of this metric."""
        pass

    def get_found_polymorphisms(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
    ) -> set[Polymorphism]:
        """Get expected polymorphisms that were found in sample.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup

        Returns:
            Intersection of sample and expected polymorphisms
        """
        return sample_polymorphisms & expected_polymorphisms

    def get_missing_polymorphisms(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
    ) -> set[Polymorphism]:
        """Get expected polymorphisms that were NOT found in sample.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup

        Returns:
            Expected polymorphisms not in sample
        """
        return expected_polymorphisms - sample_polymorphisms

    def get_remaining_polymorphisms(
        self,
        sample_polymorphisms: set[Polymorphism],
        expected_polymorphisms: set[Polymorphism],
    ) -> set[Polymorphism]:
        """Get sample polymorphisms that are not in expected set.

        These are "private" mutations or mutations from other haplogroups.

        Args:
            sample_polymorphisms: Set of polymorphisms observed in sample
            expected_polymorphisms: Set of expected polymorphisms for haplogroup

        Returns:
            Sample polymorphisms not in expected
        """
        return sample_polymorphisms - expected_polymorphisms
