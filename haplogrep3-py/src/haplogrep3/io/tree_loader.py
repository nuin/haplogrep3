"""Phylotree loader for loading haplogroup trees from various sources."""

import os
import re
from pathlib import Path
from typing import Optional

import yaml

from haplogrep3.models import Phylotree, PhyloTreeNode, Haplogroup, Polymorphism


class PhylotreeLoader:
    """Loads phylogenetic trees from files or remote repositories."""

    # Default tree repository URL
    DEFAULT_REPOSITORY = "https://raw.githubusercontent.com/genepi/haplogrep-trees/main"

    def __init__(self, trees_dir: Optional[Path] = None):
        """Initialize the tree loader.

        Args:
            trees_dir: Directory containing local tree files
        """
        self.trees_dir = trees_dir or Path.home() / ".haplogrep3" / "trees"
        self.trees_dir.mkdir(parents=True, exist_ok=True)

    def load(self, tree_id: str) -> Phylotree:
        """Load a phylotree by ID.

        Args:
            tree_id: Tree identifier (e.g., 'phylotree-fu-rcrs@1.2')

        Returns:
            Loaded Phylotree instance
        """
        # Check if tree_id is a file path
        if os.path.exists(tree_id):
            return self._load_from_file(Path(tree_id))

        # Check if tree exists locally
        tree_path = self._find_local_tree(tree_id)
        if tree_path:
            return self._load_from_file(tree_path)

        raise FileNotFoundError(
            f"Tree '{tree_id}' not found. Install it with: haplogrep3 install-tree {tree_id}"
        )

    def _find_local_tree(self, tree_id: str) -> Optional[Path]:
        """Find a tree in the local trees directory.

        Args:
            tree_id: Tree identifier

        Returns:
            Path to tree YAML file if found
        """
        # Look for tree.yaml in tree directory
        tree_dir = self.trees_dir / tree_id
        yaml_file = tree_dir / "tree.yaml"
        if yaml_file.exists():
            return yaml_file

        # Also check for direct YAML files
        yaml_file = self.trees_dir / f"{tree_id}.yaml"
        if yaml_file.exists():
            return yaml_file

        return None

    def _load_from_file(self, path: Path) -> Phylotree:
        """Load a phylotree from a YAML file.

        Args:
            path: Path to the tree YAML file

        Returns:
            Phylotree instance
        """
        with open(path) as f:
            data = yaml.safe_load(f)

        tree_dir = path.parent

        # Load the mtphyl tree file
        tree_file = tree_dir / data.get("tree", "tree.txt")
        root = self._parse_mtphyl_tree(tree_file)

        return Phylotree(
            id=data.get("id", path.stem),
            name=data.get("name", data.get("id", path.stem)),
            version=data.get("version", ""),
            root=root,
            reference_fasta=str(tree_dir / data["fasta"]) if "fasta" in data else None,
            weights_file=str(tree_dir / data["weights"]) if "weights" in data else None,
            hotspots=set(data.get("hotspots", [])),
        )

    def _parse_mtphyl_tree(self, tree_file: Path) -> PhyloTreeNode:
        """Parse an mtphyl format tree file.

        The mtphyl format uses indentation to represent the tree hierarchy:
        - Each line is a haplogroup with its defining mutations
        - Indentation (tabs) indicates parent-child relationships
        - Format: [tabs]HaplogroupName[tab]mutation1[tab]mutation2...

        Args:
            tree_file: Path to the mtphyl tree file

        Returns:
            Root PhyloTreeNode
        """
        root = PhyloTreeNode(
            haplogroup=Haplogroup(name="root"),
            polymorphisms=[],
            children=[],
        )

        if not tree_file.exists():
            return root

        # Stack to track parent nodes at each indentation level
        stack: list[tuple[int, PhyloTreeNode]] = [(-1, root)]

        with open(tree_file) as f:
            for line in f:
                if not line.strip():
                    continue

                # Count leading tabs for indentation level
                indent = 0
                for char in line:
                    if char == '\t':
                        indent += 1
                    else:
                        break

                # Parse the line: haplogroup name and mutations
                parts = line.strip().split('\t')
                if not parts:
                    continue

                haplogroup_name = parts[0]
                mutations = []
                for mut_str in parts[1:]:
                    mut_str = mut_str.strip()
                    if mut_str:
                        try:
                            mutations.append(Polymorphism.from_string(mut_str))
                        except (ValueError, IndexError):
                            pass  # Skip invalid mutations

                # Create new node
                node = PhyloTreeNode(
                    haplogroup=Haplogroup(name=haplogroup_name),
                    polymorphisms=mutations,
                    children=[],
                )

                # Find parent by popping stack until we find lower indent
                while stack and stack[-1][0] >= indent:
                    stack.pop()

                # Add as child of parent
                if stack:
                    parent = stack[-1][1]
                    parent.children.append(node)
                    node.parent = parent

                # Push this node onto stack
                stack.append((indent, node))

        return root

    def list_available(self) -> list[str]:
        """List available local trees.

        Returns:
            List of tree IDs
        """
        trees = []
        if self.trees_dir.exists():
            for item in self.trees_dir.iterdir():
                if item.is_dir() and (item / "tree.yaml").exists():
                    trees.append(item.name)
                elif item.suffix == ".yaml":
                    trees.append(item.stem)
        return sorted(trees)


def load_phylotree(tree_id: str, trees_dir: Optional[Path] = None) -> Phylotree:
    """Convenience function to load a phylotree.

    Args:
        tree_id: Tree identifier or path
        trees_dir: Optional directory containing trees

    Returns:
        Loaded Phylotree instance
    """
    loader = PhylotreeLoader(trees_dir)
    return loader.load(tree_id)
