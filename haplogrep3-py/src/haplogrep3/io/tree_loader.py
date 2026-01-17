"""Phylotree loader for loading haplogroup trees from various sources."""

import os
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path
from typing import Optional

import yaml

from haplogrep3.models import Phylotree, PhyloTreeNode, Haplogroup, Polymorphism

# Module-level cache for loaded trees
_tree_cache: dict[str, Phylotree] = {}


def get_cached_tree(tree_path: str) -> Optional[Phylotree]:
    """Get a tree from the cache if available."""
    return _tree_cache.get(tree_path)


def cache_tree(tree_path: str, tree: Phylotree) -> None:
    """Cache a loaded tree."""
    _tree_cache[tree_path] = tree


def clear_tree_cache() -> None:
    """Clear the tree cache."""
    _tree_cache.clear()


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

    def _load_from_file(self, path: Path, use_cache: bool = True) -> Phylotree:
        """Load a phylotree from a YAML file.

        Args:
            path: Path to the tree YAML file
            use_cache: Whether to use/update the cache

        Returns:
            Phylotree instance
        """
        cache_key = str(path.resolve())

        # Check cache first
        if use_cache:
            cached = get_cached_tree(cache_key)
            if cached is not None:
                return cached

        with open(path) as f:
            data = yaml.safe_load(f)

        tree_dir = path.parent

        # Load the tree file (XML or mtphyl format)
        tree_filename = data.get("tree", "tree.txt")
        tree_file = tree_dir / tree_filename

        if tree_filename.endswith(".xml"):
            root = self._parse_xml_tree(tree_file)
        else:
            root = self._parse_mtphyl_tree(tree_file)

        tree = Phylotree(
            id=data.get("id", path.stem),
            name=data.get("name", data.get("id", path.stem)),
            version=str(data.get("version", "")),
            root=root,
            reference_fasta=str(tree_dir / data["fasta"]) if "fasta" in data else None,
            weights_file=str(tree_dir / data["weights"]) if "weights" in data else None,
            hotspots=set(data.get("hotspots", [])),
        )

        # Cache the loaded tree
        if use_cache:
            cache_tree(cache_key, tree)

        return tree

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

    def _parse_xml_tree(self, tree_file: Path) -> PhyloTreeNode:
        """Parse an XML format tree file.

        The XML format uses nested haplogroup elements with poly children:
        <phylotree>
          <haplogroup name="H">
            <details>
              <poly>263G</poly>
            </details>
            <haplogroup name="H1">...</haplogroup>
          </haplogroup>
        </phylotree>

        Args:
            tree_file: Path to the XML tree file

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

        tree = ET.parse(tree_file)
        xml_root = tree.getroot()

        def parse_haplogroup(element: ET.Element, parent: PhyloTreeNode) -> PhyloTreeNode:
            """Recursively parse haplogroup elements."""
            name = element.get("name", "unknown")

            # Get polymorphisms from details/poly elements
            mutations = []
            details = element.find("details")
            if details is not None:
                for poly in details.findall("poly"):
                    if poly.text:
                        try:
                            mutations.append(Polymorphism.from_string(poly.text.strip()))
                        except (ValueError, IndexError):
                            pass

            node = PhyloTreeNode(
                haplogroup=Haplogroup(name=name),
                polymorphisms=mutations,
                children=[],
                parent=parent,
            )

            # Recursively parse child haplogroups
            for child_elem in element.findall("haplogroup"):
                child_node = parse_haplogroup(child_elem, node)
                node.children.append(child_node)

            return node

        # Parse all top-level haplogroups
        for hg_elem in xml_root.findall(".//haplogroup"):
            # Only process top-level haplogroups (those directly under phylotree or without a haplogroup parent)
            parent_tag = None
            for parent in xml_root.iter():
                if hg_elem in parent:
                    parent_tag = parent.tag
                    break

            if parent_tag == "phylotree":
                node = parse_haplogroup(hg_elem, root)
                root.children.append(node)

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
