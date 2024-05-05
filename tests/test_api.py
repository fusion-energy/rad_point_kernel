"""Python API tests for rad_point_kernel."""

import math

import pytest

import rad_point_kernel as pkc


# --- Source tests ---


class TestSource:
    def test_photon_source(self):
        s = pkc.Source("photon", 662e3)
        assert s.particle == "photon"

    def test_neutron_source(self):
        s = pkc.Source("neutron", 14.1e6)
        assert s.particle == "neutron"

    def test_spectrum_source(self):
        s = pkc.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
        assert s.particle == "photon"

    def test_invalid_particle(self):
        with pytest.raises(ValueError):
            pkc.Source("electron", 1e6)

    def test_repr(self):
        s = pkc.Source("neutron", 14.1e6)
        assert "neutron" in repr(s)


# --- Material tests ---


class TestMaterial:
    def test_element_material(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        assert iron.density == 7.874
        assert iron.fraction == "mass"

    def test_nuclide_material(self):
        mat = pkc.Material(composition={"Li6": 0.6, "Li7": 0.4}, density=0.534)
        assert mat.density == 0.534

    def test_formula_material(self):
        water = pkc.Material(composition={"H2O": 1.0}, density=1.0)
        assert water.density == 1.0

    def test_multi_element(self):
        concrete = pkc.Material(
            composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
            density=2.3, fraction="mass",
        )
        assert concrete.density == 2.3

    def test_atom_fraction(self):
        water = pkc.Material(composition={"H": 2.0, "O": 1.0}, density=1.0, fraction="atom")
        assert water.fraction == "atom"

    def test_volume_mix(self):
        concrete = pkc.Material(composition={"O": 0.53, "Si": 0.34, "Ca": 0.13}, density=2.3)
        steel = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        mixed = pkc.Material.volume_mix(concrete, 0.97, steel, 0.03)
        expected = 0.97 * 2.3 + 0.03 * 7.874
        assert mixed.density == pytest.approx(expected, rel=1e-10)

    def test_invalid_density(self):
        with pytest.raises(ValueError):
            pkc.Material(composition={"Fe": 1.0}, density=-1.0)

    def test_invalid_fraction(self):
        with pytest.raises(ValueError):
            pkc.Material(composition={"Fe": -0.5}, density=7.874)


