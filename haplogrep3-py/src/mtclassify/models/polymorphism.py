"""Polymorphism model representing mtDNA variants."""

from enum import Enum
from typing import Optional
from pydantic import BaseModel, Field


class Mutation(str, Enum):
    """Nucleotide mutation types."""
    A = "A"
    T = "T"
    G = "G"
    C = "C"
    # IUPAC ambiguity codes
    R = "R"  # A or G
    Y = "Y"  # C or T
    K = "K"  # G or T
    M = "M"  # A or C
    S = "S"  # G or C
    W = "W"  # A or T
    H = "H"  # A or C or T
    B = "B"  # C or G or T
    V = "V"  # A or C or G
    D = "D"  # A or G or T
    N = "N"  # any
    DEL = "d"  # deletion


class Polymorphism(BaseModel):
    """Represents a single mtDNA polymorphism (variant).

    Format examples: 263G, 315.1C (insertion), 522d (deletion)
    Position is 1-based (mtDNA positions 1-16569).
    """
    position: int = Field(ge=1, le=16569, description="mtDNA position (1-based)")
    insertion_position: Optional[int] = Field(
        default=None,
        ge=1,
        description="Insertion position for indels (e.g., 315.1 has insertion_position=1)"
    )
    mutation: Mutation = Field(description="The variant nucleotide")
    reference: Optional[str] = Field(default=None, description="Reference nucleotide")

    @classmethod
    def from_string(cls, poly_str: str) -> "Polymorphism":
        """Parse polymorphism from string format like '263G', '315.1C', '522d'.

        Args:
            poly_str: String representation of polymorphism

        Returns:
            Polymorphism instance
        """
        poly_str = poly_str.strip()

        # Handle deletion
        if poly_str.endswith('d'):
            position = int(poly_str[:-1])
            return cls(position=position, mutation=Mutation.DEL)

        # Find where the position ends and mutation begins
        i = 0
        while i < len(poly_str) and (poly_str[i].isdigit() or poly_str[i] == '.'):
            i += 1

        position_str = poly_str[:i]
        mutation_str = poly_str[i:].upper()

        # Handle insertion position (e.g., 315.1)
        if '.' in position_str:
            parts = position_str.split('.')
            position = int(parts[0])
            insertion_position = int(parts[1])
        else:
            position = int(position_str)
            insertion_position = None

        return cls(
            position=position,
            insertion_position=insertion_position,
            mutation=Mutation(mutation_str[0]) if mutation_str else Mutation.N
        )

    def to_string(self) -> str:
        """Convert to standard string format."""
        if self.insertion_position is not None:
            return f"{self.position}.{self.insertion_position}{self.mutation.value}"
        return f"{self.position}{self.mutation.value}"

    def __str__(self) -> str:
        return self.to_string()

    def __hash__(self) -> int:
        return hash((self.position, self.insertion_position, self.mutation))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Polymorphism):
            return False
        return (
            self.position == other.position
            and self.insertion_position == other.insertion_position
            and self.mutation == other.mutation
        )

    @property
    def sort_key(self) -> tuple[int, int]:
        """Key for sorting polymorphisms by position."""
        return (self.position, self.insertion_position or 0)
