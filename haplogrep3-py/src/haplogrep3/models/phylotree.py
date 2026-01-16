"""Phylotree model for haplogroup classification."""

from typing import Optional
from pydantic import BaseModel, Field

from .polymorphism import Polymorphism
from .haplogroup import Haplogroup


class PhyloTreeNode(BaseModel):
    """A node in the phylogenetic tree representing a haplogroup."""
    haplogroup: Haplogroup
    polymorphisms: list[Polymorphism] = Field(
        default_factory=list,
        description="Polymorphisms defining this node (branch-specific)"
    )
    children: list["PhyloTreeNode"] = Field(
        default_factory=list,
        description="Child nodes (subhaplogroups)"
    )
    parent: Optional["PhyloTreeNode"] = Field(
        default=None,
        exclude=True,  # Avoid circular serialization
        description="Parent node"
    )

    model_config = {"arbitrary_types_allowed": True}

    def get_all_polymorphisms(self) -> list[Polymorphism]:
        """Get all polymorphisms from root to this node.

        Returns:
            List of all polymorphisms accumulated from root to this node
        """
        polys = list(self.polymorphisms)
        if self.parent:
            polys = self.parent.get_all_polymorphisms() + polys
        return polys

    def find_node(self, haplogroup_name: str) -> Optional["PhyloTreeNode"]:
        """Find a node by haplogroup name.

        Args:
            haplogroup_name: Name of haplogroup to find

        Returns:
            PhyloTreeNode if found, None otherwise
        """
        if self.haplogroup.name.lower() == haplogroup_name.lower():
            return self
        for child in self.children:
            result = child.find_node(haplogroup_name)
            if result:
                return result
        return None

    def get_all_nodes(self) -> list["PhyloTreeNode"]:
        """Get all nodes in the subtree rooted at this node.

        Returns:
            List of all nodes including this one
        """
        nodes = [self]
        for child in self.children:
            nodes.extend(child.get_all_nodes())
        return nodes


class Phylotree(BaseModel):
    """Represents a complete phylogenetic tree for haplogroup classification."""
    id: str = Field(description="Tree identifier (e.g., 'phylotree-fu-rcrs@1.2')")
    name: str = Field(description="Human-readable name")
    version: str = Field(default="")
    root: PhyloTreeNode = Field(description="Root node of the tree")
    reference_fasta: Optional[str] = Field(default=None, description="Path to reference FASTA")
    weights_file: Optional[str] = Field(default=None, description="Path to weights file")
    hotspots: set[int] = Field(default_factory=set, description="Known hotspot positions")

    model_config = {"arbitrary_types_allowed": True}

    def get_haplogroup(self, name: str) -> Optional[Haplogroup]:
        """Find a haplogroup by name.

        Args:
            name: Haplogroup name to find

        Returns:
            Haplogroup if found, None otherwise
        """
        node = self.root.find_node(name)
        return node.haplogroup if node else None

    def get_node(self, haplogroup_name: str) -> Optional[PhyloTreeNode]:
        """Find a tree node by haplogroup name.

        Args:
            haplogroup_name: Name of haplogroup to find

        Returns:
            PhyloTreeNode if found, None otherwise
        """
        return self.root.find_node(haplogroup_name)

    def get_all_haplogroups(self) -> list[Haplogroup]:
        """Get all haplogroups in the tree.

        Returns:
            List of all haplogroups
        """
        return [node.haplogroup for node in self.root.get_all_nodes()]

    def get_polymorphisms_for_haplogroup(self, name: str) -> list[Polymorphism]:
        """Get all polymorphisms that define a haplogroup.

        Args:
            name: Haplogroup name

        Returns:
            List of polymorphisms from root to haplogroup
        """
        node = self.root.find_node(name)
        if node:
            return node.get_all_polymorphisms()
        return []
