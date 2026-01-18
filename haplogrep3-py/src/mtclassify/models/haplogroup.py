"""Haplogroup model representing mtDNA haplogroups."""

from typing import Optional
from pydantic import BaseModel, Field

from .polymorphism import Polymorphism


class Haplogroup(BaseModel):
    """Represents an mtDNA haplogroup.

    A haplogroup is a genetic population group defined by specific
    polymorphisms in the mitochondrial DNA.
    """
    name: str = Field(description="Haplogroup name (e.g., 'H', 'H1a', 'L3e2b')")
    expected_polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Polymorphisms that define this haplogroup"
    )

    def __str__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Haplogroup):
            return False
        return self.name.lower() == other.name.lower()

    def is_super_haplogroup(self, other: "Haplogroup") -> bool:
        """Check if this haplogroup is a superhaplogroup of another.

        For example, 'H' is a superhaplogroup of 'H1a'.

        Args:
            other: The potential subhaplogroup

        Returns:
            True if this is a superhaplogroup of other
        """
        return other.name.lower().startswith(self.name.lower())
