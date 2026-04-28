"""Python API tests for rad_point_kernel."""

import math

import pytest

import rad_point_kernel as rpk


# Source tests


class TestSource:
    def test_photon_source(self):
        s = rpk.Source("photon", 662e3)
        assert s.particle == "photon"

    def test_neutron_source(self):
        s = rpk.Source("neutron", 14.1e6)
        assert s.particle == "neutron"

    def test_spectrum_source(self):
        s = rpk.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
        assert s.particle == "photon"

    def test_invalid_particle(self):
        with pytest.raises(ValueError):
            rpk.Source("electron", 1e6)

    def test_repr(self):
        s = rpk.Source("neutron", 14.1e6)
        assert "neutron" in repr(s)


# Material tests


class TestMaterial:
    def test_element_material(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        assert iron.density == 7.874
        assert iron.fraction == "mass"

    def test_nuclide_material(self):
        mat = rpk.Material(composition={"Li6": 0.6, "Li7": 0.4}, density=0.534)
        assert mat.density == 0.534

    def test_formula_material(self):
        water = rpk.Material(composition={"H2O": 1.0}, density=1.0)
        assert water.density == 1.0

    def test_multi_element(self):
        concrete = rpk.Material(
            composition={"H": 0.01, "O": 0.53, "Si": 0.34, "Ca": 0.04, "Al": 0.03, "Fe": 0.01},
            density=2.3, fraction="mass",
        )
        assert concrete.density == 2.3

    def test_atom_fraction(self):
        water = rpk.Material(composition={"H": 2.0, "O": 1.0}, density=1.0, fraction="atom")
        assert water.fraction == "atom"

    def test_volume_mix(self):
        concrete = rpk.Material(composition={"O": 0.53, "Si": 0.34, "Ca": 0.13}, density=2.3)
        steel = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        mixed = rpk.Material.volume_mix(concrete, 0.97, steel, 0.03)
        expected = 0.97 * 2.3 + 0.03 * 7.874
        assert mixed.density == pytest.approx(expected, rel=1e-10)

    def test_invalid_density(self):
        with pytest.raises(ValueError):
            rpk.Material(composition={"Fe": 1.0}, density=-1.0)

    def test_invalid_fraction(self):
        with pytest.raises(ValueError):
            rpk.Material(composition={"Fe": -0.5}, density=7.874)


# Layer tests


class TestLayer:
    def test_void_layer(self):
        layer = rpk.Layer(thickness=50)
        assert layer.thickness == 50.0
        assert layer.has_material is False

    def test_material_layer(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layer = rpk.Layer(thickness=10, material=iron)
        assert layer.thickness == 10.0
        assert layer.has_material is True

    def test_invalid_thickness(self):
        with pytest.raises(ValueError):
            rpk.Layer(thickness=-5)

    def test_zero_thickness(self):
        layer = rpk.Layer(thickness=0)
        assert layer.thickness == 0.0


# Transmission tests


class TestTransmission:
    def test_void_neutron(self):
        s = rpk.Source("neutron", 14.06e6)
        frac = rpk.calculate_transmission([rpk.Layer(thickness=100)], s)
        assert frac == pytest.approx(1.0, abs=1e-15)

    def test_void_photon(self):
        s = rpk.Source("photon", 1e6)
        frac = rpk.calculate_transmission([rpk.Layer(thickness=100)], s)
        assert frac == pytest.approx(1.0, abs=1e-15)

    def test_material_less_than_one(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        s = rpk.Source("neutron", 14.06e6)
        frac = rpk.calculate_transmission([rpk.Layer(thickness=10, material=iron)], s)
        assert 0.0 < frac < 1.0

    def test_thicker_less_transmission(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        s = rpk.Source("neutron", 14.06e6)
        thin = rpk.calculate_transmission([rpk.Layer(thickness=5, material=iron)], s)
        thick = rpk.calculate_transmission([rpk.Layer(thickness=50, material=iron)], s)
        assert thin > thick

    def test_multi_layer_equals_single(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        s = rpk.Source("neutron", 14.06e6)
        single = rpk.calculate_transmission([rpk.Layer(thickness=10, material=iron)], s)
        double = rpk.calculate_transmission(
            [rpk.Layer(thickness=5, material=iron), rpk.Layer(thickness=5, material=iron)], s
        )
        assert single == pytest.approx(double, rel=1e-12)

    def test_spectrum(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        s = rpk.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
        frac = rpk.calculate_transmission([rpk.Layer(thickness=10, material=iron)], s)
        assert 0.0 < frac < 1.0

    def test_different_energies(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=10, material=iron)]
        s1 = rpk.Source("neutron", 2e6)
        s2 = rpk.Source("neutron", 14.06e6)
        assert rpk.calculate_transmission(layers, s1) != pytest.approx(
            rpk.calculate_transmission(layers, s2), rel=0.01
        )


# Flux tests


class TestFlux:
    def test_inverse_square_law(self):
        s = rpk.Source("neutron", 14.06e6)
        r10 = rpk.calculate_flux([rpk.Layer(thickness=10)], s).scale(1e12)
        r100 = rpk.calculate_flux([rpk.Layer(thickness=100)], s).scale(1e12)
        assert r10.flux / r100.flux == pytest.approx(100.0, rel=1e-12)

    def test_void_flux_formula(self):
        s = rpk.Source("neutron", 14.06e6)
        result = rpk.calculate_flux([rpk.Layer(thickness=200)], s).scale(1e12)
        expected = 1e12 / (4 * math.pi * 200**2)
        assert result.flux == pytest.approx(expected, rel=1e-12)

    def test_transmission_in_result(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("neutron", 14.06e6)
        flux = rpk.calculate_flux(layers, s)
        trans = rpk.calculate_transmission(layers, s)
        assert flux.transmission_fraction == pytest.approx(trans, rel=1e-12)

    def test_scale_returns_new(self):
        s = rpk.Source("neutron", 14.06e6)
        layers = [rpk.Layer(thickness=200)]
        per_particle = rpk.calculate_flux(layers, s)
        scaled = per_particle.scale(strength=1e12)
        assert per_particle.source_strength == 1.0
        assert scaled.source_strength == 1e12
        assert scaled.flux == pytest.approx(
            per_particle.flux * 1e12, rel=1e-12
        )

    def test_scale_replaces_not_compounds(self):
        s = rpk.Source("photon", 662e3)
        layers = [rpk.Layer(thickness=10)]
        first = rpk.calculate_flux(layers, s).scale(1e12)
        second = first.scale(1e10)
        assert second.source_strength == 1e10
        per_particle = rpk.calculate_flux(layers, s)
        assert second.flux == pytest.approx(
            per_particle.flux * 1e10, rel=1e-12
        )

    def test_scale_invalid(self):
        s = rpk.Source("neutron", 14.06e6)
        result = rpk.calculate_flux([rpk.Layer(thickness=100)], s)
        with pytest.raises(ValueError):
            result.scale(0.0)
        with pytest.raises(ValueError):
            result.scale(-1.0)


# Dose tests


class TestDose:
    def test_neutron_dose(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=50), rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("neutron", 14.06e6)
        result = rpk.calculate_dose(layers, s, "AP")
        assert result.dose > 0

    def test_photon_dose(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=50), rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("photon", 1e6)
        result = rpk.calculate_dose(layers, s, "AP")
        assert result.dose > 0

    def test_all_geometries(self):
        s = rpk.Source("neutron", 14.06e6)
        for geo in ["AP", "PA", "RLAT", "LLAT", "ROT", "ISO"]:
            result = rpk.calculate_dose([rpk.Layer(thickness=100)], s, geo)
            assert result.dose > 0, f"Failed for geometry {geo}"

    def test_invalid_geometry(self):
        s = rpk.Source("neutron", 14.06e6)
        with pytest.raises(ValueError):
            rpk.calculate_dose([rpk.Layer(thickness=100)], s, "INVALID")

    def test_spectrum_dose(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("photon", [(1173e3, 1.0), (1333e3, 1.0)])
        result = rpk.calculate_dose(layers, s, "AP")
        assert result.dose > 0

    def test_scale_dose(self):
        s = rpk.Source("photon", 662e3)
        layers = [rpk.Layer(thickness=100)]
        per_particle = rpk.calculate_dose(layers, s, "AP")
        scaled = per_particle.scale(strength=1e10)
        assert scaled.dose == pytest.approx(per_particle.dose * 1e10, rel=1e-12)


# Buildup model tests


class TestBuildup:
    def test_no_buildup(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        s = rpk.Source("neutron", 14.06e6)
        result = rpk.calculate_flux([rpk.Layer(thickness=10, material=iron)], s)
        assert result.buildup_factor == pytest.approx(1.0)

    def test_linear_buildup(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("neutron", 14.06e6)
        buildup = rpk.BuildupModel.linear(a=0.5)
        r_no = rpk.calculate_flux(layers, s)
        r_bu = rpk.calculate_flux(layers, s, buildup=buildup)
        assert r_bu.flux > r_no.flux

    def test_buildup_result_auto_resolve(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("neutron", 14.06e6)

        br = rpk.BuildupResult()
        br.buildup = {"flux": 2.5}
        result = rpk.calculate_flux(layers, s, buildup=br)
        assert result.buildup_factor == pytest.approx(2.5)

    def test_interpolation_result_as_buildup(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=10, material=iron)]
        s = rpk.Source("photon", 662e3)

        ir = rpk.InterpolationResult(value=3.0, sigma=0.1)
        result = rpk.calculate_flux(layers, s, buildup=ir)
        assert result.buildup_factor == pytest.approx(3.0)


# Material fraction tests


class TestMaterialFractions:
    def test_atom_vs_mass_fraction_water(self):
        water_atom = rpk.Material(composition={"H": 2.0, "O": 1.0}, density=1.0, fraction="atom")
        water_mass = rpk.Material(composition={"H": 0.111898, "O": 0.888102}, density=1.0, fraction="mass")
        s = rpk.Source("neutron", 14.06e6)
        frac_atom = rpk.calculate_transmission([rpk.Layer(thickness=20, material=water_atom)], s)
        frac_mass = rpk.calculate_transmission([rpk.Layer(thickness=20, material=water_mass)], s)
        assert frac_atom == pytest.approx(frac_mass, rel=1e-3)


# Secondary photon dose tests


class TestSecondaryPhotonDose:
    def test_void_no_gammas(self):
        result = rpk.calculate_secondary_photon_dose([rpk.Layer(thickness=100)], rpk.Source("neutron", 14.06e6), "AP")
        assert result.secondary_photon_dose == 0.0
        assert result.neutron_dose > 0.0

    def test_material_produces_gammas(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=50), rpk.Layer(thickness=10, material=iron)]
        result = rpk.calculate_secondary_photon_dose(layers, rpk.Source("neutron", 14.06e6), "AP")
        assert result.secondary_photon_dose > 0.0

    def test_total_is_sum(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=50), rpk.Layer(thickness=10, material=iron)]
        result = rpk.calculate_secondary_photon_dose(layers, rpk.Source("neutron", 14.06e6), "AP")
        assert result.total_dose == pytest.approx(
            result.neutron_dose + result.secondary_photon_dose, rel=1e-12
        )

    def test_neutron_dose_matches_standalone(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=50), rpk.Layer(thickness=10, material=iron)]
        coupled = rpk.calculate_secondary_photon_dose(layers, rpk.Source("neutron", 14.06e6), "AP")
        s = rpk.Source("neutron", 14.06e6)
        standalone = rpk.calculate_dose(layers, s, "AP")
        assert coupled.neutron_dose == pytest.approx(standalone.dose, rel=1e-10)

    def test_scale_secondary_photon(self):
        iron = rpk.Material(composition={"Fe": 1.0}, density=7.874)
        layers = [rpk.Layer(thickness=50), rpk.Layer(thickness=10, material=iron)]
        per_particle = rpk.calculate_secondary_photon_dose(layers, rpk.Source("neutron", 14.06e6), "AP")
        scaled = per_particle.scale(strength=1e12)
        assert scaled.source_strength == 1e12
        assert scaled.neutron_dose == pytest.approx(per_particle.neutron_dose * 1e12, rel=1e-12)
        assert scaled.secondary_photon_dose == pytest.approx(per_particle.secondary_photon_dose * 1e12, rel=1e-12)
        assert scaled.total_dose == pytest.approx(per_particle.total_dose * 1e12, rel=1e-12)


# BuildupResult tests


class TestBuildupResult:
    def test_to_dict_round_trip(self):
        r = rpk.BuildupResult()
        r.mc = {"dose-AP": 1.5e-11}
        r.mc_std_dev = {"dose-AP": 1e-13}
        r.pk = {"dose-AP": 1.2e-11}
        r.buildup = {"dose-AP": 1.25}
        r.optical_thickness = 3.5
        d = r.to_dict()
        r2 = rpk.BuildupResult.from_dict(d)
        assert r2.mc == r.mc
        assert r2.buildup == r.buildup
        assert r2.optical_thickness == r.optical_thickness

    def test_save_load_round_trip(self, tmp_path):
        r1 = rpk.BuildupResult()
        r1.mc = {"flux": 1e-6}
        r1.buildup = {"flux": 1.25}
        r2 = rpk.BuildupResult()
        r2.mc = {"flux": 5e-7}
        r2.buildup = {"flux": 1.3}
        path = tmp_path / "cache.json"
        rpk.BuildupResult.save([r1, r2], path)
        loaded = rpk.BuildupResult.load(path)
        assert len(loaded) == 2
        assert loaded[0].mc == r1.mc

    def test_empty_result(self):
        r = rpk.BuildupResult()
        assert r.mc == {}
        assert r.optical_thickness == 0.0
        assert r.source_strength == 1.0

    def test_scale(self):
        r = rpk.BuildupResult()
        r.mc = {"dose-AP": 1.5e-11}
        r.mc_std_dev = {"dose-AP": 1e-13}
        r.pk = {"dose-AP": 1.2e-11}
        r.buildup = {"dose-AP": 1.25}
        r.optical_thickness = 3.5

        scaled = r.scale(strength=1e10)
        assert scaled.source_strength == 1e10
        assert scaled.mc["dose-AP"] == pytest.approx(1.5e-11 * 1e10)
        assert scaled.pk["dose-AP"] == pytest.approx(1.2e-11 * 1e10)
        # Buildup is a ratio, must be unchanged
        assert scaled.buildup["dose-AP"] == pytest.approx(1.25)
        # Original is unchanged
        assert r.source_strength == 1.0
        assert r.mc["dose-AP"] == pytest.approx(1.5e-11)

    def test_scale_replaces_not_compounds(self):
        r = rpk.BuildupResult()
        r.mc = {"flux": 1e-6}
        r.pk = {"flux": 1e-6}
        r.buildup = {"flux": 1.0}
        once = r.scale(strength=1e12)
        twice = once.scale(strength=1e10)
        assert twice.source_strength == 1e10
        assert twice.mc["flux"] == pytest.approx(1e-6 * 1e10)

    def test_scale_round_trips_in_json(self, tmp_path):
        r = rpk.BuildupResult()
        r.mc = {"dose-AP": 1.5e-11}
        r.buildup = {"dose-AP": 1.25}
        scaled = r.scale(strength=1e10)
        path = tmp_path / "cache.json"
        rpk.BuildupResult.save([scaled], path)
        loaded = rpk.BuildupResult.load(path)
        assert loaded[0].source_strength == 1e10

    def test_scale_invalid(self):
        r = rpk.BuildupResult()
        r.mc = {"flux": 1e-6}
        with pytest.raises(ValueError):
            r.scale(strength=0.0)
        with pytest.raises(ValueError):
            r.scale(strength=-1.0)


# BuildupTable tests


def _make_results(b_values, quantity="dose-AP"):
    results = []
    for b in b_values:
        r = rpk.BuildupResult()
        r.mc = {quantity: b * 1e-11}
        r.mc_std_dev = {quantity: 0.001 * 1e-11}
        r.pk = {quantity: 1e-11}
        r.buildup = {quantity: b}
        results.append(r)
    return results


class TestBuildupTable:
    def test_1d_interpolation(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        r = table.interpolate(t=20)
        assert r.value == pytest.approx(1.3, abs=0.01)
        assert not r.is_extrapolated

    def test_1d_extrapolation_detected(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        r = table.interpolate(t=50, warn=False)
        assert r.is_extrapolated
        assert "t" in r.extrapolated_axes

    def test_extrapolation_larger_sigma(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        assert table.interpolate(t=100, warn=False).sigma > table.interpolate(t=20).sigma

    def test_2d(self):
        points = [{"a": 10, "b": 10}, {"a": 10, "b": 20}, {"a": 20, "b": 10}, {"a": 20, "b": 20}]
        table = rpk.BuildupTable(points, _make_results([1.1, 1.2, 1.3, 1.4]))
        r = table.interpolate(a=15, b=15)
        assert 1.0 < r.value < 1.5

    def test_available_quantities(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3], "flux"))
        assert "flux" in table.available_quantities

    def test_default_quantity(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        assert table.interpolate(t=15).value > 0

    def test_invalid_quantity_raises(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        with pytest.raises(ValueError, match="not available"):
            table.interpolate(t=15, quantity="flux")

    def test_wrong_axes_raises(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        with pytest.raises(ValueError, match="Expected axes"):
            table.interpolate(wrong=15)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            rpk.BuildupTable([{"t": 10}], _make_results([1.2]))

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            rpk.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2]))

    def test_inconsistent_axes_raises(self):
        with pytest.raises(ValueError, match="same axes"):
            rpk.BuildupTable([{"t": 10}, {"x": 20}], _make_results([1.2, 1.3]))

    def test_axis_ranges(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 30}], _make_results([1.2, 1.3]))
        assert table.axis_ranges == {"t": (10.0, 30.0)}

    def test_repr(self):
        table = rpk.BuildupTable([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        assert "BuildupTable" in repr(table)


def _build_paired_result(geo, n_mc, p_mc, n_std, p_std, pk_n):
    r = rpk.BuildupResult()
    n_name = f"dose-{geo}"
    p_name = f"dose-{geo}-coupled-photon"
    r.mc = {n_name: n_mc, p_name: p_mc}
    r.mc_std_dev = {n_name: n_std, p_name: p_std}
    r.pk = {n_name: pk_n}
    if pk_n > 0:
        r.buildup = {n_name: n_mc / pk_n}
    return r, n_name, p_name


class TestSynthesizeDoseTotals:
    def test_synthesizes_total_for_paired_dose(self):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r, n_name, p_name = _build_paired_result("AP", 1.0e-12, 4.0e-13, 5e-14, 3e-14, 8e-13)
        parsed = [(q, _parse_quantity(q)) for q in [n_name, p_name]]

        _synthesize_dose_totals(r, parsed)

        total = "dose-AP-total"
        assert r.mc[total] == pytest.approx(1.4e-12)
        assert r.mc_std_dev[total] == pytest.approx(math.sqrt(5e-14**2 + 3e-14**2))
        assert r.pk[total] == pytest.approx(8e-13)
        assert r.buildup[total] == pytest.approx(1.4e-12 / 8e-13)

    def test_no_synthesis_when_only_neutron_dose(self):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r = rpk.BuildupResult()
        r.mc = {"dose-AP": 1e-12}
        r.pk = {"dose-AP": 8e-13}
        parsed = [("dose-AP", _parse_quantity("dose-AP"))]

        _synthesize_dose_totals(r, parsed)
        assert "dose-AP-total" not in r.mc

    def test_no_synthesis_when_only_coupled_photon(self):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r = rpk.BuildupResult()
        r.mc = {"dose-AP-coupled-photon": 4e-13}
        parsed = [("dose-AP-coupled-photon", _parse_quantity("dose-AP-coupled-photon"))]

        _synthesize_dose_totals(r, parsed)
        assert "dose-AP-total" not in r.mc

    def test_flux_is_not_totalled(self):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r = rpk.BuildupResult()
        r.mc = {"flux": 1e-6, "flux-coupled-photon": 2e-7}
        parsed = [(q, _parse_quantity(q)) for q in ["flux", "flux-coupled-photon"]]

        _synthesize_dose_totals(r, parsed)
        assert "flux-total" not in r.mc

    @pytest.mark.parametrize("geo", ["AP", "PA", "RLAT", "LLAT", "ROT", "ISO"])
    def test_works_for_all_icrp_geometries(self, geo):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r, n_name, p_name = _build_paired_result(geo, 1e-12, 5e-13, 0, 0, 8e-13)
        parsed = [(q, _parse_quantity(q)) for q in [n_name, p_name]]

        _synthesize_dose_totals(r, parsed)
        assert f"dose-{geo}-total" in r.mc
        assert r.mc[f"dose-{geo}-total"] == pytest.approx(1.5e-12)

    def test_two_geometries_get_separate_totals(self):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r = rpk.BuildupResult()
        r.mc = {
            "dose-AP": 1e-12, "dose-AP-coupled-photon": 4e-13,
            "dose-PA": 8e-13, "dose-PA-coupled-photon": 3e-13,
        }
        r.mc_std_dev = {k: 0.0 for k in r.mc}
        r.pk = {"dose-AP": 8e-13, "dose-PA": 7e-13}
        r.buildup = {"dose-AP": 1.25, "dose-PA": 1.14}
        parsed = [
            (q, _parse_quantity(q)) for q in
            ["dose-AP", "dose-AP-coupled-photon", "dose-PA", "dose-PA-coupled-photon"]
        ]

        _synthesize_dose_totals(r, parsed)
        assert r.mc["dose-AP-total"] == pytest.approx(1.4e-12)
        assert r.mc["dose-PA-total"] == pytest.approx(1.1e-12)
        assert r.buildup["dose-AP-total"] == pytest.approx(1.4e-12 / 8e-13)
        assert r.buildup["dose-PA-total"] == pytest.approx(1.1e-12 / 7e-13)

    def test_no_buildup_when_pk_is_zero(self):
        from rad_point_kernel.buildup import _parse_quantity, _synthesize_dose_totals

        r, n_name, p_name = _build_paired_result("AP", 1e-12, 4e-13, 0, 0, 0.0)
        parsed = [(q, _parse_quantity(q)) for q in [n_name, p_name]]

        _synthesize_dose_totals(r, parsed)
        assert "dose-AP-total" in r.mc
        assert "dose-AP-total" not in r.buildup
