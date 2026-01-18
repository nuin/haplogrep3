"""Sample model representing classified mtDNA samples."""

from typing import Optional
from pydantic import BaseModel, Field

from .polymorphism import Polymorphism
from .haplogroup import Haplogroup


class RankedResult(BaseModel):
    """A classification result with haplogroup and quality score."""
    haplogroup: Haplogroup
    quality: float = Field(ge=0.0, le=1.0, description="Quality score (0-1)")
    expected_polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Polymorphisms expected for this haplogroup"
    )
    found_polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Expected polymorphisms found in sample"
    )
    missing_polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Expected polymorphisms not found in sample"
    )
    remaining_polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Sample polymorphisms not in haplogroup definition"
    )

    def __str__(self) -> str:
        return f"{self.haplogroup.name} ({self.quality:.2%})"


class ClassificationResult(BaseModel):
    """Complete classification result for a sample."""
    top_result: RankedResult
    other_results: list[RankedResult] = Field(default_factory=list)

    @property
    def haplogroup(self) -> Haplogroup:
        """Get the top-hit haplogroup."""
        return self.top_result.haplogroup

    @property
    def quality(self) -> float:
        """Get the top-hit quality score."""
        return self.top_result.quality


class Sample(BaseModel):
    """Represents an mtDNA sample with polymorphisms and classification results."""
    id: str = Field(description="Sample identifier")
    polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Observed polymorphisms in the sample"
    )
    range_start: int = Field(default=1, ge=1, le=16569)
    range_end: int = Field(default=16569, ge=1, le=16569)
    classification: Optional[ClassificationResult] = Field(
        default=None,
        description="Classification result after running haplogroup assignment"
    )

    @property
    def range_str(self) -> str:
        """Get range as string format."""
        return f"{self.range_start}-{self.range_end}"

    @property
    def polymorphism_set(self) -> set[Polymorphism]:
        """Get polymorphisms as a set for comparison."""
        return set(self.polymorphisms)

    def get_polymorphism_strings(self) -> list[str]:
        """Get polymorphisms as sorted list of strings."""
        sorted_polys = sorted(self.polymorphisms, key=lambda p: p.sort_key)
        return [str(p) for p in sorted_polys]
