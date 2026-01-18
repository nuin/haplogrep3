"""Distance metrics for haplogroup classification."""

from enum import Enum
from .base import DistanceMetric
from .kulczynski import KulczynskiDistance
from .hamming import HammingDistance
from .jaccard import JaccardDistance
from .kimura import KimuraDistance


class Distance(str, Enum):
    """Available distance metrics."""
    KULCZYNSKI = "kulczynski"
    HAMMING = "hamming"
    JACCARD = "jaccard"
    KIMURA = "kimura"


def get_distance_metric(distance: Distance) -> DistanceMetric:
    """Get distance metric instance by enum value.

    Args:
        distance: Distance enum value

    Returns:
        DistanceMetric instance
    """
    metrics = {
        Distance.KULCZYNSKI: KulczynskiDistance(),
        Distance.HAMMING: HammingDistance(),
        Distance.JACCARD: JaccardDistance(),
        Distance.KIMURA: KimuraDistance(),
    }
    return metrics[distance]


__all__ = [
    "Distance",
    "DistanceMetric",
    "KulczynskiDistance",
    "HammingDistance",
    "JaccardDistance",
    "KimuraDistance",
    "get_distance_metric",
]
