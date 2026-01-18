"""Tests for Polymorphism model."""

import pytest
from mtclassify.models import Polymorphism, Mutation


class TestPolymorphism:
    def test_from_string_simple(self):
        poly = Polymorphism.from_string("263G")
        assert poly.position == 263
        assert poly.mutation == Mutation.G
        assert poly.insertion_position is None

    def test_from_string_insertion(self):
        poly = Polymorphism.from_string("315.1C")
        assert poly.position == 315
        assert poly.insertion_position == 1
        assert poly.mutation == Mutation.C

    def test_from_string_deletion(self):
        poly = Polymorphism.from_string("522d")
        assert poly.position == 522
        assert poly.mutation == Mutation.DEL

    def test_to_string(self):
        poly = Polymorphism(position=263, mutation=Mutation.G)
        assert str(poly) == "263G"

    def test_to_string_insertion(self):
        poly = Polymorphism(position=315, insertion_position=1, mutation=Mutation.C)
        assert str(poly) == "315.1C"

    def test_equality(self):
        poly1 = Polymorphism.from_string("263G")
        poly2 = Polymorphism.from_string("263G")
        assert poly1 == poly2

    def test_inequality(self):
        poly1 = Polymorphism.from_string("263G")
        poly2 = Polymorphism.from_string("263A")
        assert poly1 != poly2

    def test_hash(self):
        poly1 = Polymorphism.from_string("263G")
        poly2 = Polymorphism.from_string("263G")
        assert hash(poly1) == hash(poly2)

        # Can be used in sets
        poly_set = {poly1, poly2}
        assert len(poly_set) == 1

    def test_sort_key(self):
        polys = [
            Polymorphism.from_string("315.2C"),
            Polymorphism.from_string("263G"),
            Polymorphism.from_string("315.1C"),
        ]
        sorted_polys = sorted(polys, key=lambda p: p.sort_key)
        assert sorted_polys[0].position == 263
        assert sorted_polys[1].insertion_position == 1
        assert sorted_polys[2].insertion_position == 2
