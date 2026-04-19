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


class TestFlux:
    def test_inverse_square_law(self):
        s = pkc.Source("neutron", 14.06e6)
        r10 = pkc.calculate_flux(1e12, [pkc.Layer(thickness=10)], s)
        r100 = pkc.calculate_flux(1e12, [pkc.Layer(thickness=100)], s)
        assert r10.uncollided_flux / r100.uncollided_flux == pytest.approx(100.0, rel=1e-12)

    def test_void_flux_formula(self):
        s = pkc.Source("neutron", 14.06e6)
        result = pkc.calculate_flux(1e12, [pkc.Layer(thickness=200)], s)
        expected = 1e12 / (4 * math.pi * 200**2)
        assert result.uncollided_flux == pytest.approx(expected, rel=1e-12)

    def test_transmission_in_result(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("neutron", 14.06e6)
        flux = pkc.calculate_flux(1e12, layers, s)
        trans = pkc.calculate_transmission(layers, s)
        assert flux.transmission_fraction == pytest.approx(trans, rel=1e-12)


# --- Dose tests ---


class TestDose:
    def test_neutron_dose(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=50), pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("neutron", 14.06e6)
        result = pkc.calculate_dose(1e12, layers, s, "AP")
        assert result.dose_rate > 0

    def test_photon_dose(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=50), pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("photon", 1e6)
        result = pkc.calculate_dose(1e12, layers, s, "AP")
        assert result.dose_rate > 0

    def test_all_geometries(self):
        s = pkc.Source("neutron", 14.06e6)
        for geo in ["AP", "PA", "RLAT", "LLAT", "ROT", "ISO"]:
            result = pkc.calculate_dose(1e12, [pkc.Layer(thickness=100)], s, geo)
            assert result.dose_rate > 0, f"Failed for geometry {geo}"

    def test_invalid_geometry(self):
        s = pkc.Source("neutron", 14.06e6)
        with pytest.raises(ValueError):
            pkc.calculate_dose(1e12, [pkc.Layer(thickness=100)], s, "INVALID")

    def test_spectrum_dose(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
        result = pkc.calculate_dose(1e12, layers, s, "AP")
        assert result.dose_rate > 0


# --- Buildup model tests ---


class TestBuildup:
    def test_no_buildup(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        s = pkc.Source("neutron", 14.06e6)
        result = pkc.calculate_flux(1e12, [pkc.Layer(thickness=10, material=iron)], s)
        assert result.buildup_factor == pytest.approx(1.0)

    def test_linear_buildup(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("neutron", 14.06e6)
        buildup = pkc.BuildupModel.linear(a=0.5)
        r_no = pkc.calculate_flux(1e12, layers, s)
        r_bu = pkc.calculate_flux(1e12, layers, s, buildup=buildup)
        assert r_bu.uncollided_flux > r_no.uncollided_flux

    def test_buildup_result_auto_resolve(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("neutron", 14.06e6)

        br = pkc.BuildupResult()
        br.buildup = {"flux": 2.5}
        result = pkc.calculate_flux(1e12, layers, s, buildup=br)
        assert result.buildup_factor == pytest.approx(2.5)

    def test_interpolation_result_as_buildup(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=10, material=iron)]
        s = pkc.Source("photon", 662e3)

        ir = pkc.InterpolationResult(value=3.0, sigma=0.1)
        result = pkc.calculate_flux(1e12, layers, s, buildup=ir)
        assert result.buildup_factor == pytest.approx(3.0)


# --- Material fraction tests ---


class TestMaterialFractions:
    def test_atom_vs_mass_fraction_water(self):
        water_atom = pkc.Material(composition={"H": 2.0, "O": 1.0}, density=1.0, fraction="atom")
        water_mass = pkc.Material(composition={"H": 0.111898, "O": 0.888102}, density=1.0, fraction="mass")
        s = pkc.Source("neutron", 14.06e6)
        frac_atom = pkc.calculate_transmission([pkc.Layer(thickness=20, material=water_atom)], s)
        frac_mass = pkc.calculate_transmission([pkc.Layer(thickness=20, material=water_mass)], s)
        assert frac_atom == pytest.approx(frac_mass, rel=1e-3)


# --- Secondary photon dose tests ---


class TestSecondaryPhotonDose:
    def test_void_no_gammas(self):
        result = pkc.calculate_secondary_photon_dose_rate(1e12, [pkc.Layer(thickness=100)], pkc.Source("neutron", 14.06e6), "AP")
        assert result.secondary_photon_dose_rate == 0.0
        assert result.neutron_dose_rate > 0.0

    def test_material_produces_gammas(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=50), pkc.Layer(thickness=10, material=iron)]
        result = pkc.calculate_secondary_photon_dose_rate(1e12, layers, pkc.Source("neutron", 14.06e6), "AP")
        assert result.secondary_photon_dose_rate > 0.0

    def test_total_is_sum(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=50), pkc.Layer(thickness=10, material=iron)]
        result = pkc.calculate_secondary_photon_dose_rate(1e12, layers, pkc.Source("neutron", 14.06e6), "AP")
        assert result.total_dose_rate == pytest.approx(
            result.neutron_dose_rate + result.secondary_photon_dose_rate, rel=1e-12
        )

    def test_neutron_dose_matches_standalone(self):
        iron = pkc.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [pkc.Layer(thickness=50), pkc.Layer(thickness=10, material=iron)]
        coupled = pkc.calculate_secondary_photon_dose_rate(1e12, layers, pkc.Source("neutron", 14.06e6), "AP")
        s = pkc.Source("neutron", 14.06e6)
        standalone = pkc.calculate_dose(1e12, layers, s, "AP")
        assert coupled.neutron_dose_rate == pytest.approx(standalone.dose_rate, rel=1e-10)


# --- BuildupResult tests ---


class TestBuildupResult:
    def test_to_dict_round_trip(self):
        r = pkc.BuildupResult()
        r.mc = {"dose-AP": 1.5e-11}
        r.mc_std_dev = {"dose-AP": 1e-13}
        r.pk = {"dose-AP": 1.2e-11}
        r.buildup = {"dose-AP": 1.25}
        r.optical_thickness = 3.5
        d = r.to_dict()
        r2 = pkc.BuildupResult.from_dict(d)
        assert r2.mc == r.mc
        assert r2.buildup == r.buildup
        assert r2.optical_thickness == r.optical_thickness

    def test_save_load_round_trip(self, tmp_path):
        r1 = pkc.BuildupResult()
        r1.mc = {"flux": 1e-6}
        r1.buildup = {"flux": 1.25}
        r2 = pkc.BuildupResult()
        r2.mc = {"flux": 5e-7}
        r2.buildup = {"flux": 1.3}
        path = tmp_path / "cache.json"
        pkc.BuildupResult.save([r1, r2], path)
        loaded = pkc.BuildupResult.load(path)
        assert len(loaded) == 2
        assert loaded[0].mc == r1.mc

    def test_empty_result(self):
        r = pkc.BuildupResult()
        assert r.mc == {}
        assert r.optical_thickness == 0.0


# --- BuildupTable tests ---


def _make_results(b_values, quantity="dose-AP"):
    results = []
    for b in b_values:
        r = pkc.BuildupResult()
        r.mc = {quantity: b * 1e-11}
        r.mc_std_dev = {quantity: 0.001 * 1e-11}
        r.pk = {quantity: 1e-11}
        r.buildup = {quantity: b}
        results.append(r)
    return results


class TestBuildupTable:
    def test_1d_interpolation(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        r = table.interpolate(t=20)
        assert r.value == pytest.approx(1.3, abs=0.01)
        assert not r.is_extrapolated

    def test_1d_extrapolation_detected(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        r = table.interpolate(t=50, warn=False)
        assert r.is_extrapolated
        assert "t" in r.extrapolated_axes

    def test_extrapolation_larger_sigma(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        assert table.interpolate(t=100, warn=False).sigma > table.interpolate(t=20).sigma

    def test_2d(self):
        points = [{"a": 10, "b": 10}, {"a": 10, "b": 20}, {"a": 20, "b": 10}, {"a": 20, "b": 20}]
        table = pkc.BuildupTable(points, _make_results([1.1, 1.2, 1.3, 1.4]))
        r = table.interpolate(a=15, b=15)
        assert 1.0 < r.value < 1.5

    def test_available_quantities(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3], "flux"))
        assert "flux" in table.available_quantities

    def test_default_quantity(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        assert table.interpolate(t=15).value > 0

    def test_invalid_quantity_raises(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        with pytest.raises(ValueError, match="not available"):
            table.interpolate(t=15, quantity="flux")

    def test_wrong_axes_raises(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        with pytest.raises(ValueError, match="Expected axes"):
            table.interpolate(wrong=15)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            pkc.BuildupTable([{"t": 10}], _make_results([1.2]))

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            pkc.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2]))

    def test_inconsistent_axes_raises(self):
        with pytest.raises(ValueError, match="same axes"):
            pkc.BuildupTable([{"t": 10}, {"x": 20}], _make_results([1.2, 1.3]))

    def test_axis_ranges(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 30}], _make_results([1.2, 1.3]))
        assert table.axis_ranges == {"t": (10.0, 30.0)}

    def test_repr(self):
        table = pkc.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        assert "BuildupTable" in repr(table)
