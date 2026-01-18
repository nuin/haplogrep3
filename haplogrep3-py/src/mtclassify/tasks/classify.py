"""Classification task for haplogroup assignment."""

from typing import Optional
import heapq

from mtclassify.models import (
    Sample,
    Phylotree,
    PhyloTreeNode,
    Haplogroup,
    Polymorphism,
    RankedResult,
    ClassificationResult,
)
from mtclassify.distance import Distance, DistanceMetric, get_distance_metric


class ClassificationTask:
    """Classifies samples against a phylogenetic tree to determine haplogroups."""

    def __init__(
        self,
        phylotree: Phylotree,
        distance: Distance = Distance.KULCZYNSKI,
        hits: int = 1,
    ):
        """Initialize the classification task.

        Args:
            phylotree: Phylogenetic tree to classify against
            distance: Distance metric to use
            hits: Number of top hits to return per sample
        """
        self.phylotree = phylotree
        self.distance = distance
        self.hits = max(1, hits)
        self.metric = get_distance_metric(distance)

    def classify(self, samples: list[Sample]) -> list[Sample]:
        """Classify a list of samples.

        Args:
            samples: List of samples to classify

        Returns:
            List of samples with classification results populated
        """
        for sample in samples:
            sample.classification = self._classify_sample(sample)
        return samples

    def _classify_sample(self, sample: Sample) -> ClassificationResult:
        """Classify a single sample.

        Args:
            sample: Sample to classify

        Returns:
            ClassificationResult with ranked haplogroup matches
        """
        sample_polys = sample.polymorphism_set

        # Get weights and hotspots from phylotree for weighted distance calculation
        weights = self.phylotree.get_weights() if self.phylotree.weights_file else None
        hotspots = self.phylotree.hotspots if self.phylotree.hotspots else None

        # Score all haplogroups
        scores: list[tuple[float, PhyloTreeNode]] = []

        for node in self.phylotree.root.get_all_nodes():
            if node.haplogroup.name == "root":
                continue

            expected_polys = set(node.get_all_polymorphisms())
            if not expected_polys:
                continue

            quality = self.metric.calculate(sample_polys, expected_polys, weights, hotspots)
            scores.append((quality, node))

        # Get top N hits
        top_scores = heapq.nlargest(self.hits, scores, key=lambda x: x[0])

        if not top_scores:
            # No matches found - return unknown haplogroup
            return ClassificationResult(
                top_result=RankedResult(
                    haplogroup=Haplogroup(name="?"),
                    quality=0.0,
                ),
            )

        # Build results
        results = []
        for quality, node in top_scores:
            expected_polys = set(node.get_all_polymorphisms())
            found = self.metric.get_found_polymorphisms(sample_polys, expected_polys)
            missing = self.metric.get_missing_polymorphisms(sample_polys, expected_polys)
            remaining = self.metric.get_remaining_polymorphisms(sample_polys, expected_polys)

            result = RankedResult(
                haplogroup=node.haplogroup,
                quality=quality,
                expected_polymorphisms=sorted(expected_polys, key=lambda p: p.sort_key),
                found_polymorphisms=sorted(found, key=lambda p: p.sort_key),
                missing_polymorphisms=sorted(missing, key=lambda p: p.sort_key),
                remaining_polymorphisms=sorted(remaining, key=lambda p: p.sort_key),
            )
            results.append(result)

        return ClassificationResult(
            top_result=results[0],
            other_results=results[1:] if len(results) > 1 else [],
        )


def classify_samples(
    samples: list[Sample],
    phylotree: Phylotree,
    distance: Distance = Distance.KULCZYNSKI,
    hits: int = 1,
) -> list[Sample]:
    """Convenience function to classify samples.

    Args:
        samples: List of samples to classify
        phylotree: Phylogenetic tree
        distance: Distance metric
        hits: Number of top hits

    Returns:
        Classified samples
    """
    task = ClassificationTask(phylotree, distance, hits)
    return task.classify(samples)
