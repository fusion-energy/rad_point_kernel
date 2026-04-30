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


# BuildupFit tests


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


class TestBuildupFit:
    def test_1d_interpolation(self):
        # Shin-Ishii is a smooth least-squares fit, not an interpolator;
        # don't expect to pass through training points exactly. Sub-5%
        # match on a smooth 3-point dataset is fine.
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        r = table.interpolate(t=20)
        assert r.value == pytest.approx(1.3, abs=0.05)
        assert not r.is_extrapolated

    def test_1d_extrapolation_detected(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        r = table.interpolate(t=50, warn=False)
        assert r.is_extrapolated
        assert "t" in r.extrapolated_axes

    def test_sigma_not_exposed(self):
        # BuildupFit doesn't expose predictive sigma (NaN sentinel); see TODO.md.
        import math
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}, {"t": 30}], _make_results([1.2, 1.3, 1.25]))
        assert math.isnan(table.interpolate(t=20).sigma)
        assert math.isnan(table.interpolate(t=100, warn=False).sigma)

    def test_2d(self):
        points = [{"a": 10, "b": 10}, {"a": 10, "b": 20}, {"a": 20, "b": 10}, {"a": 20, "b": 20}]
        table = rpk.BuildupFit(points, _make_results([1.1, 1.2, 1.3, 1.4]))
        r = table.interpolate(a=15, b=15)
        assert 1.0 < r.value < 1.5

    def test_available_quantities(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3], "flux"))
        assert "flux" in table.available_quantities

    def test_default_quantity(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        assert table.interpolate(t=15).value > 0

    def test_multi_quantity_default_raises(self):
        results = []
        for b_neutron, b_total in [(1.2, 1.5), (1.3, 1.7)]:
            r = rpk.BuildupResult()
            r.mc = {"dose-AP": b_neutron * 1e-11, "dose-AP-total": b_total * 1e-11}
            r.mc_std_dev = {"dose-AP": 1e-13, "dose-AP-total": 1e-13}
            r.pk = {"dose-AP": 1e-11, "dose-AP-total": 1e-11}
            r.buildup = {"dose-AP": b_neutron, "dose-AP-total": b_total}
            results.append(r)
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}], results)

        with pytest.raises(ValueError, match="multiple quantities"):
            table.interpolate(t=15)

        # Explicit quantity= still works.
        assert table.interpolate(t=15, quantity="dose-AP-total").value > 0

    def test_invalid_quantity_raises(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        with pytest.raises(ValueError, match="not available"):
            table.interpolate(t=15, quantity="flux")

    def test_wrong_axes_raises(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        with pytest.raises(ValueError, match="Expected axes"):
            table.interpolate(wrong=15)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            rpk.BuildupFit([{"t": 10}], _make_results([1.2]))

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            rpk.BuildupFit([{"t": 10}, {"t": 20}], _make_results([1.2]))

    def test_inconsistent_axes_raises(self):
        with pytest.raises(ValueError, match="same axes"):
            rpk.BuildupFit([{"t": 10}, {"x": 20}], _make_results([1.2, 1.3]))

    def test_axis_ranges(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 30}], _make_results([1.2, 1.3]))
        assert table.axis_ranges == {"t": (10.0, 30.0)}

    def test_repr(self):
        table = rpk.BuildupFit([{"t": 10}, {"t": 20}], _make_results([1.2, 1.3]))
        assert "BuildupFit" in repr(table)

    def test_coupled_photon_dispatch(self):
        # Build BuildupResults with both `dose-AP` (B(0)=1) and
        # `dose-AP-coupled-photon` (B(0)=0). Verify BuildupFit fits both
        # quantities and dispatches to the right form: Shin-Ishii for the
        # primary dose, Power x saturator for the coupled photons.
        thicks = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        # Construct synthetic curves matching the expected forms
        results = []
        for t in thicks:
            B_n = 1.0 + 0.5 * t                         # rises from 1, B(0)=1
            B_c = 0.04 * t ** 1.1 * (1 - 2.71828 ** (-0.6 * t))  # rises from 0
            r = rpk.BuildupResult()
            r.optical_thickness = t
            r.mc = {"dose-AP": B_n * 1e-11, "dose-AP-coupled-photon": B_c * 1e-11}
            r.mc_std_dev = {"dose-AP": 1e-13, "dose-AP-coupled-photon": 1e-13}
            r.pk = {"dose-AP": 1e-11}
            r.buildup = {"dose-AP": B_n}
            r.synthesize_dose_totals()
            # synthesize_dose_totals adds dose-AP-coupled-photon and
            # dose-AP-total to pk and buildup using pk_neutron as reference.
            results.append(r)

        fit = rpk.BuildupFit([{"t": t} for t in thicks], results)
        avail = sorted(fit.available_quantities)
        assert "dose-AP" in avail
        assert "dose-AP-coupled-photon" in avail
        assert "dose-AP-total" in avail

        # B(0) ~= 1 for primary (Shin-Ishii)
        b_n = fit.interpolate(t=1e-6, quantity="dose-AP", warn=False).value
        assert abs(b_n - 1.0) < 0.05

        # B(0) ~= 0 for coupled photon (Power x saturator)
        b_c = fit.interpolate(t=1e-6, quantity="dose-AP-coupled-photon", warn=False).value
        assert abs(b_c) < 1e-3, f"expected B(0)~0, got {b_c}"

    def test_total_is_sum_of_components_on_u_shape_data(self):
        # Regression for: 1D dose-AP-total used to be fit independently
        # of dose-AP and dose-AP-coupled-photon. On poly-style data the
        # total is U-shaped (neutron-dominated near the surface, secondary
        # photon-dominated at depth), and Shin-Ishii (sum of two decaying
        # exponentials) cannot represent a U, so the direct B_total fit
        # used to drift >30% off the MC anchors and explode by ~60x in
        # extrapolation. Fix: when both component fits are available,
        # dose-{geo}-total is now derived as their sum, preserving the
        # B_t = B_n + B_p identity that holds in the MC inputs.
        t_anchors = [25, 50, 75, 125, 175, 250]
        # U-shaped MC: B_n decays, B_p grows as a power, sum dips near
        # 125 cm and climbs again at 250 cm.
        b_n = [math.exp(-0.040 * (t - 25)) * 0.97 for t in t_anchors]
        b_p = [2.52e-5 * t**1.976 for t in t_anchors]
        b_t = [n + p for n, p in zip(b_n, b_p)]
        # Sanity: sum is U-shaped over the anchor range.
        assert b_t[0] > b_t[3] and b_t[-1] > b_t[3]

        results = []
        for n, p, t in zip(b_n, b_p, b_t):
            r = rpk.BuildupResult()
            r.mc = {
                "dose-AP": n * 1e-11,
                "dose-AP-coupled-photon": p * 1e-11,
                "dose-AP-total": t * 1e-11,
            }
            r.mc_std_dev = {
                "dose-AP": 1e-13,
                "dose-AP-coupled-photon": 1e-13,
                "dose-AP-total": 1e-13,
            }
            r.pk = {
                "dose-AP": 1e-11,
                "dose-AP-coupled-photon": 1e-11,
                "dose-AP-total": 1e-11,
            }
            r.buildup = {"dose-AP": n, "dose-AP-coupled-photon": p, "dose-AP-total": t}
            results.append(r)

        fit = rpk.BuildupFit([{"t": t} for t in t_anchors], results)

        # 1) Identity holds at every training anchor (machine precision).
        #    Under the old code dose-AP-total came from an independent
        #    Shin-Ishii fit and drifted up to ~34% off the sum; this is
        #    the assertion that would have failed.
        for t, n_mc, p_mc, t_mc in zip(t_anchors, b_n, b_p, b_t):
            bn = fit.interpolate(t=t, quantity="dose-AP").value
            bp = fit.interpolate(t=t, quantity="dose-AP-coupled-photon").value
            bt = fit.interpolate(t=t, quantity="dose-AP-total").value
            assert bt == pytest.approx(bn + bp, rel=1e-12), (
                f"At t={t}, total {bt} != sum {bn + bp}"
            )
            # And the sum (and therefore total) stays close to MC; under
            # the old code the direct total fit was off by 30%+ at the
            # tails of this U-shape.
            assert bt == pytest.approx(t_mc, rel=0.10), (
                f"At t={t}, total {bt} drifted from MC {t_mc}"
            )

        # 2) Extrapolation stays bounded. The old direct fit on this
        #    U-shape exploded (~60x) at deep thickness; the composite
        #    stays equal to the sum of components (each of which
        #    extrapolates within its own form).
        for t in [10, 400, 500]:
            bn = fit.interpolate(t=t, quantity="dose-AP", warn=False).value
            bp = fit.interpolate(
                t=t, quantity="dose-AP-coupled-photon", warn=False
            ).value
            bt = fit.interpolate(
                t=t, quantity="dose-AP-total", warn=False
            ).value
            assert bt == pytest.approx(bn + bp, rel=1e-12)


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
        r, _, _ = _build_paired_result("AP", 1.0e-12, 4.0e-13, 5e-14, 3e-14, 8e-13)
        r.synthesize_dose_totals()

        total = "dose-AP-total"
        assert r.mc[total] == pytest.approx(1.4e-12)
        assert r.mc_std_dev[total] == pytest.approx(math.sqrt(5e-14**2 + 3e-14**2))
        assert r.pk[total] == pytest.approx(8e-13)
        assert r.buildup[total] == pytest.approx(1.4e-12 / 8e-13)

    def test_no_synthesis_when_only_neutron_dose(self):
        r = rpk.BuildupResult()
        r.mc = {"dose-AP": 1e-12}
        r.pk = {"dose-AP": 8e-13}
        r.synthesize_dose_totals()
        assert "dose-AP-total" not in r.mc

    def test_no_synthesis_when_only_coupled_photon(self):
        r = rpk.BuildupResult()
        r.mc = {"dose-AP-coupled-photon": 4e-13}
        r.synthesize_dose_totals()
        assert "dose-AP-total" not in r.mc

    def test_flux_is_not_totalled(self):
        r = rpk.BuildupResult()
        r.mc = {"flux": 1e-6, "flux-coupled-photon": 2e-7}
        r.synthesize_dose_totals()
        assert "flux-total" not in r.mc

    @pytest.mark.parametrize("geo", ["AP", "PA", "RLAT", "LLAT", "ROT", "ISO"])
    def test_works_for_all_icrp_geometries(self, geo):
        r, _, _ = _build_paired_result(geo, 1e-12, 5e-13, 0, 0, 8e-13)
        r.synthesize_dose_totals()
        assert f"dose-{geo}-total" in r.mc
        assert r.mc[f"dose-{geo}-total"] == pytest.approx(1.5e-12)

    def test_two_geometries_get_separate_totals(self):
        r = rpk.BuildupResult()
        r.mc = {
            "dose-AP": 1e-12, "dose-AP-coupled-photon": 4e-13,
            "dose-PA": 8e-13, "dose-PA-coupled-photon": 3e-13,
        }
        r.mc_std_dev = {k: 0.0 for k in r.mc}
        r.pk = {"dose-AP": 8e-13, "dose-PA": 7e-13}
        r.buildup = {"dose-AP": 1.25, "dose-PA": 1.14}
        r.synthesize_dose_totals()
        assert r.mc["dose-AP-total"] == pytest.approx(1.4e-12)
        assert r.mc["dose-PA-total"] == pytest.approx(1.1e-12)
        assert r.buildup["dose-AP-total"] == pytest.approx(1.4e-12 / 8e-13)
        assert r.buildup["dose-PA-total"] == pytest.approx(1.1e-12 / 7e-13)

    def test_no_buildup_when_pk_is_zero(self):
        r, _, _ = _build_paired_result("AP", 1e-12, 4e-13, 0, 0, 0.0)
        r.synthesize_dose_totals()
        assert "dose-AP-total" in r.mc
        assert "dose-AP-total" not in r.buildup

    def test_idempotent(self):
        r, _, _ = _build_paired_result("AP", 1e-12, 4e-13, 0, 0, 8e-13)
        r.synthesize_dose_totals()
        r.synthesize_dose_totals()  # second call must not break anything
        assert r.mc["dose-AP-total"] == pytest.approx(1.4e-12)


class TestBuildupResultValidation:
    def test_setter_rejects_non_float_value(self):
        r = rpk.BuildupResult()
        with pytest.raises(TypeError, match="must be float"):
            r.mc = {"dose-AP": "not a number"}

    def test_setter_rejects_non_float_in_buildup(self):
        r = rpk.BuildupResult()
        with pytest.raises(TypeError, match="BuildupResult.buildup"):
            r.buildup = {"dose-AP": [1.0]}

    def test_from_dict_rejects_empty_dict(self):
        with pytest.raises(ValueError, match="malformed"):
            rpk.BuildupResult.from_dict({})

    def test_from_dict_rejects_dict_with_no_known_keys(self):
        with pytest.raises(ValueError, match="malformed"):
            rpk.BuildupResult.from_dict({"junk": 1, "other": 2})

    def test_from_dict_accepts_partial_payload(self):
        # Building from buildup-only is valid (e.g. for resolve_buildup_py).
        r = rpk.BuildupResult.from_dict({"buildup": {"dose-AP": 1.5}})
        assert r.buildup["dose-AP"] == 1.5

    def test_from_dict_rejects_non_float_inside(self):
        with pytest.raises(TypeError, match="must be float"):
            rpk.BuildupResult.from_dict({"mc": {"dose-AP": "string"}})


class TestTotalDoseRequest:
    def test_expands_total_to_both_halves(self):
        from rad_point_kernel.buildup import _expand_total_requests

        s = rpk.Source("neutron", 14e6)
        out = _expand_total_requests(["dose-AP-total"], s)
        assert out == ["dose-AP", "dose-AP-coupled-photon"]

    def test_dedups_when_halves_already_listed(self):
        from rad_point_kernel.buildup import _expand_total_requests

        s = rpk.Source("neutron", 14e6)
        out = _expand_total_requests(
            ["dose-AP", "dose-AP-total", "dose-AP-coupled-photon"], s
        )
        assert out == ["dose-AP", "dose-AP-coupled-photon"]

    def test_preserves_other_quantities(self):
        from rad_point_kernel.buildup import _expand_total_requests

        s = rpk.Source("neutron", 14e6)
        out = _expand_total_requests(["flux", "dose-AP-total"], s)
        assert out == ["flux", "dose-AP", "dose-AP-coupled-photon"]

    def test_rejects_total_for_photon_source(self):
        from rad_point_kernel.buildup import _expand_total_requests

        s = rpk.Source("photon", 1e6)
        with pytest.raises(ValueError, match="requires a neutron source"):
            _expand_total_requests(["dose-AP-total"], s)

    def test_rejects_invalid_geometry(self):
        from rad_point_kernel.buildup import _expand_total_requests

        s = rpk.Source("neutron", 14e6)
        with pytest.raises(ValueError, match="must be one of"):
            _expand_total_requests(["dose-OBLIQUE-total"], s)


class TestDoseGeometryValidation:
    def test_calculate_dose_rejects_lowercase_geometry(self):
        layers = [rpk.Layer(thickness=10)]
        s = rpk.Source("photon", 662e3)
        with pytest.raises(ValueError, match="invalid geometry"):
            rpk.calculate_dose(layers=layers, source=s, geometry="ap")

    def test_calculate_dose_rejects_unknown_geometry(self):
        layers = [rpk.Layer(thickness=10)]
        s = rpk.Source("photon", 662e3)
        with pytest.raises(ValueError, match="invalid geometry"):
            rpk.calculate_dose(layers=layers, source=s, geometry="OBLIQUE")

    def test_calculate_dose_rejects_geometry_with_whitespace(self):
        layers = [rpk.Layer(thickness=10)]
        s = rpk.Source("photon", 662e3)
        with pytest.raises(ValueError, match="invalid geometry"):
            rpk.calculate_dose(layers=layers, source=s, geometry="AP ")

    def test_calculate_dose_accepts_all_six_geometries(self):
        layers = [rpk.Layer(thickness=10)]
        s = rpk.Source("photon", 662e3)
        for geo in ("AP", "PA", "RLAT", "LLAT", "ROT", "ISO"):
            result = rpk.calculate_dose(layers=layers, source=s, geometry=geo)
            assert result.dose > 0
