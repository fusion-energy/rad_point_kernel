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


# --- Layer tests ---


class TestLayer:
    def test_void_layer(self):
        layer = pkc.Layer(thickness=50)
        assert layer.thickness == 50.0
        assert layer.has_material is False

    def test_material_layer(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layer = pkc.Layer(thickness=10, material=iron)
        assert layer.thickness == 10.0
        assert layer.has_material is True

    def test_invalid_thickness(self):
        with pytest.raises(ValueError):
            pkc.Layer(thickness=-5)

    def test_zero_thickness(self):
        layer = pkc.Layer(thickness=0)
        assert layer.thickness == 0.0


# --- Transmission tests ---


class TestTransmission:
    def test_void_neutron(self):
        s = pkc.Source("neutron", 14.06e6)
        frac = pkc.calculate_transmission([pkc.Layer(thickness=100)], s)
        assert frac == pytest.approx(1.0, abs=1e-15)

    def test_void_photon(self):
        s = pkc.Source("photon", 1e6)
        frac = pkc.calculate_transmission([pkc.Layer(thickness=100)], s)
        assert frac == pytest.approx(1.0, abs=1e-15)

    def test_material_less_than_one(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        s = pkc.Source("neutron", 14.06e6)
        frac = pkc.calculate_transmission([pkc.Layer(thickness=10, material=iron)], s)
        assert 0.0 < frac < 1.0

    def test_thicker_less_transmission(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        s = pkc.Source("neutron", 14.06e6)
        thin = pkc.calculate_transmission([pkc.Layer(thickness=5, material=iron)], s)
        thick = pkc.calculate_transmission([pkc.Layer(thickness=50, material=iron)], s)
        assert thin > thick

    def test_multi_layer_equals_single(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        s = pkc.Source("neutron", 14.06e6)
        single = pkc.calculate_transmission([pkc.Layer(thickness=10, material=iron)], s)
        double = pkc.calculate_transmission(
            [pkc.Layer(thickness=5, material=iron), pkc.Layer(thickness=5, material=iron)], s
        )
        assert single == pytest.approx(double, rel=1e-12)

    def test_spectrum(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        s = pkc.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
        frac = pkc.calculate_transmission([pkc.Layer(thickness=10, material=iron)], s)
        assert 0.0 < frac < 1.0

    def test_different_energies(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=10, material=iron)]
        s1 = pkc.Source("neutron", 2e6)
        s2 = pkc.Source("neutron", 14.06e6)
        assert pkc.calculate_transmission(layers, s1) != pytest.approx(
            pkc.calculate_transmission(layers, s2), rel=0.01
        )


# --- Flux tests ---
