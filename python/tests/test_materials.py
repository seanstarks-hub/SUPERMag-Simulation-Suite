"""Tests for the materials database."""

import pytest
from supermag.materials import get_material, list_materials, MATERIALS


class TestGetMaterial:
    def test_nb(self):
        nb = get_material("Nb")
        assert nb["Tc"] == 9.2
        assert nb["xi_S"] == 38.0
        assert nb["type"] == "superconductor"

    def test_fe(self):
        fe = get_material("Fe")
        assert fe["E_ex"] == 256.0
        assert fe["type"] == "ferromagnet"

    def test_unknown(self):
        with pytest.raises(KeyError, match="not found"):
            get_material("Unobtanium")

    def test_returns_copy(self):
        """Modifying returned dict should not affect database."""
        nb = get_material("Nb")
        nb["Tc"] = 0.0
        assert get_material("Nb")["Tc"] == 9.2


class TestListMaterials:
    def test_structure(self):
        mats = list_materials()
        assert "superconductor" in mats
        assert "ferromagnet" in mats

    def test_known_materials(self):
        mats = list_materials()
        assert "Nb" in mats["superconductor"]
        assert "Fe" in mats["ferromagnet"]
