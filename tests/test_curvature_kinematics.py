#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Deterministic verification of the implemented SA-RC and SST-RC algebra."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.curvature_validation.kinematics import (  # noqa: E402
    ZERO,
    add,
    correction_state,
    from_openfoam_gradient,
    hodge_dual,
    rotate_tensor,
)


SIMPLE_SHEAR = ((0, 2, 0), (0, 0, 0), (0, 0, 0))
PURE_STRAIN = ((2, 0, 0), (0, -2, 0), (0, 0, 0))
SOLID_ROTATION = ((0, -2, 0), (2, 0, 0), (0, 0, 0))
Q_Z_90 = ((0, -1, 0), (1, 0, 0), (0, 0, 1))


class KinematicStates(unittest.TestCase):
    def test_uniform_flow_has_zero_rates_and_zero_production_response(self) -> None:
        for family, expected_factor in (("sst", 0.0), ("sa", -1.0)):
            state = correction_state(ZERO, family=family, base_production=0.0)
            self.assertEqual(state.strain_rate, 0.0)
            self.assertEqual(state.rotation_rate, 0.0)
            self.assertEqual(state.indicator, 0.0)
            self.assertEqual(state.factor, expected_factor)
            self.assertEqual(state.corrected_production, 0.0)

    def test_simple_shear_is_neutral_without_strain_material_derivative(self) -> None:
        for family in ("sst", "sa"):
            state = correction_state(SIMPLE_SHEAR, family=family, base_production=7.0)
            self.assertAlmostEqual(state.strain_rate, 2.0)
            self.assertAlmostEqual(state.rotation_rate, 2.0)
            self.assertAlmostEqual(state.r_star, 1.0)
            self.assertAlmostEqual(state.indicator, 0.0)
            self.assertAlmostEqual(state.factor, 1.0)
            self.assertAlmostEqual(state.corrected_production, 7.0)
            self.assertFalse(state.clipped)

    def test_pure_strain_exposes_sst_upper_clip_and_sa_asymptote(self) -> None:
        sst = correction_state(PURE_STRAIN, family="sst")
        sa = correction_state(PURE_STRAIN, family="sa")
        self.assertAlmostEqual(sst.raw_factor, 3.0)
        self.assertEqual(sst.factor, 1.25)
        self.assertTrue(sst.clipped)
        self.assertAlmostEqual(sa.factor, 3.0)

    def test_solid_body_rotation_has_no_spurious_production(self) -> None:
        for family in ("sst", "sa"):
            state = correction_state(SOLID_ROTATION, family=family, base_production=0.0)
            self.assertEqual(state.strain_rate, 0.0)
            self.assertGreater(state.rotation_rate, 0.0)
            self.assertEqual(state.corrected_production, 0.0)

    def test_zero_rotation_and_zero_strain_are_regularised(self) -> None:
        zero_rotation = correction_state(PURE_STRAIN, family="sst")
        zero_strain = correction_state(SOLID_ROTATION, family="sst")
        self.assertTrue(math.isfinite(zero_rotation.r_star))
        self.assertTrue(math.isfinite(zero_strain.r_star))
        self.assertTrue(math.isfinite(zero_rotation.factor))
        self.assertTrue(math.isfinite(zero_strain.factor))

    def test_small_and_large_rates_remain_finite(self) -> None:
        for scale in (1.0e-14, 1.0e14):
            gradient = tuple(tuple(scale * item for item in row) for row in SIMPLE_SHEAR)
            for family in ("sst", "sa"):
                state = correction_state(gradient, family=family, omega=max(scale, 1.0e-30))
                self.assertTrue(math.isfinite(state.factor))
                self.assertAlmostEqual(state.factor, 1.0)

    def test_near_zero_denominator_is_finite(self) -> None:
        gradient = ((0, 1e-40, 0), (0, 0, 0), (0, 0, 0))
        dsdt = ((1e-40, 0, 0), (0, -1e-40, 0), (0, 0, 0))
        for family in ("sst", "sa"):
            state = correction_state(gradient, dsdt, family=family, omega=1e-40)
            self.assertTrue(math.isfinite(state.indicator))
            self.assertTrue(math.isfinite(state.factor))

    def test_convex_concave_idealisation_has_opposite_response(self) -> None:
        # Opposite DS/Dt signs model opposite streamline-curvature senses for
        # identical local S and W.  The invariant must be odd in DS/Dt.
        dsdt = ((1.0, 0.5, 0), (0.5, -1.0, 0), (0, 0, 0))
        opposite = tuple(tuple(-item for item in row) for row in dsdt)
        for family in ("sst", "sa"):
            first = correction_state(SIMPLE_SHEAR, dsdt, family=family)
            second = correction_state(SIMPLE_SHEAR, opposite, family=family)
            self.assertAlmostEqual(first.indicator, -second.indicator)
            self.assertNotAlmostEqual(first.factor, second.factor)
            self.assertLess((first.factor - 1.0) * (second.factor - 1.0), 0.0)

    def test_rotation_sign_reversal_changes_indicator_sign(self) -> None:
        dsdt = ((1.0, 0.5, 0), (0.5, -1.0, 0), (0, 0, 0))
        reversed_rotation = ((0, 0, 0), (2, 0, 0), (0, 0, 0))
        first = correction_state(SIMPLE_SHEAR, dsdt, family="sa")
        second = correction_state(reversed_rotation, dsdt, family="sa")
        self.assertAlmostEqual(first.strain_rate, second.strain_rate)
        self.assertAlmostEqual(first.rotation_rate, second.rotation_rate)
        self.assertAlmostEqual(first.indicator, -second.indicator)

    def test_clipping_boundaries(self) -> None:
        high = correction_state(PURE_STRAIN, family="sst", fr_max=3.0)
        at_limit = correction_state(PURE_STRAIN, family="sst", fr_max=high.raw_factor)
        below = correction_state(PURE_STRAIN, family="sst", fr_max=high.raw_factor - 1e-12)
        self.assertFalse(at_limit.clipped)
        self.assertTrue(below.clipped)
        self.assertAlmostEqual(at_limit.factor - below.factor, 1e-12, places=13)

    def test_continuity_around_regular_state(self) -> None:
        factors = []
        for epsilon in (-1e-8, 0.0, 1e-8):
            gradient = ((0, 2 + epsilon, 0), (0, 0, 0), (0, 0, 0))
            factors.append(correction_state(gradient, family="sa").factor)
        self.assertLess(max(factors) - min(factors), 1e-12)

    def test_reference_frame_basis_invariance(self) -> None:
        dsdt = ((1.0, 0.5, 0), (0.5, -1.0, 0), (0, 0, 0))
        for family in ("sst", "sa"):
            original = correction_state(SIMPLE_SHEAR, dsdt, family=family)
            transformed = correction_state(
                rotate_tensor(SIMPLE_SHEAR, Q_Z_90),
                rotate_tensor(dsdt, Q_Z_90),
                family=family,
            )
            self.assertAlmostEqual(original.strain_rate, transformed.strain_rate)
            self.assertAlmostEqual(original.rotation_rate, transformed.rotation_rate)
            self.assertAlmostEqual(original.indicator, transformed.indicator)
            self.assertAlmostEqual(original.factor, transformed.factor)

    def test_zero_frame_rotation_exactly_recovers_inertial_state(self) -> None:
        dsdt = ((1.0, 0.5, 0), (0.5, -1.0, 0), (0, 0, 0))
        for family in ("sst", "sa"):
            inertial = correction_state(SIMPLE_SHEAR, dsdt, family=family)
            framed = correction_state(
                SIMPLE_SHEAR, dsdt, family=family, frame_omega=(0, 0, 0)
            )
            self.assertEqual(inertial, framed)

    def test_openfoam_absolute_and_published_relative_mrf_agree(self) -> None:
        frame_omega = (0.0, 0.0, 0.4)
        absolute_gradient = add(SIMPLE_SHEAR, hodge_dual(frame_omega))
        dsdt = ((1.0, 0.5, 0), (0.5, -1.0, 0), (0, 0, 0))

        for family in ("sst", "sa"):
            published = correction_state(
                SIMPLE_SHEAR,
                dsdt,
                family=family,
                frame_omega=frame_omega,
                gradient_is_relative=True,
            )
            openfoam = correction_state(
                absolute_gradient,
                dsdt,
                family=family,
                frame_omega=frame_omega,
            )
            for first, second in zip(published.strain, openfoam.strain):
                for first_value, second_value in zip(first, second):
                    self.assertAlmostEqual(first_value, second_value)
            for first, second in zip(published.rotation, openfoam.rotation):
                for first_value, second_value in zip(first, second):
                    self.assertAlmostEqual(first_value, second_value)
            for first, second in zip(
                published.material_derivative_strain,
                openfoam.material_derivative_strain,
            ):
                for first_value, second_value in zip(first, second):
                    self.assertAlmostEqual(first_value, second_value)
            self.assertAlmostEqual(published.indicator, openfoam.indicator)
            self.assertAlmostEqual(published.factor, openfoam.factor)

    def test_rotating_basis_term_has_published_sign(self) -> None:
        state = correction_state(
            PURE_STRAIN, family="sa", frame_omega=(0, 0, 0.5)
        )
        self.assertEqual(
            state.material_derivative_strain,
            ((0.0, 2.0, 0.0), (2.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
        )

    def test_openfoam_gradient_storage_conversion(self) -> None:
        stored = tuple(zip(*SIMPLE_SHEAR))
        self.assertEqual(from_openfoam_gradient(stored), SIMPLE_SHEAR)

    def test_disabled_correction_exactly_recovers_base_state(self) -> None:
        for family in ("sst", "sa"):
            state = correction_state(
                SIMPLE_SHEAR,
                family=family,
                base_production=7.25,
                turbulent_viscosity=3.5,
                working_variable=2.25,
                enabled=False,
            )
            self.assertEqual(state.factor, 1.0)
            self.assertEqual(state.corrected_production, state.base_production)
            self.assertEqual(state.turbulent_viscosity, 3.5)
            self.assertEqual(state.working_variable, 2.25)

    def test_nan_and_infinity_fail_fast_in_oracle(self) -> None:
        for bad in (math.nan, math.inf, -math.inf):
            with self.assertRaises(ValueError):
                correction_state(((bad, 0, 0), (0, 0, 0), (0, 0, 0)), family="sst")
            with self.assertRaises(ValueError):
                correction_state(ZERO, family="sa", cr1=bad)
            with self.assertRaises(ValueError):
                correction_state(ZERO, family="sa", frame_omega=(0, 0, bad))


class SourceContract(unittest.TestCase):
    """Guard the oracle against unnoticed divergence from production code."""

    @classmethod
    def setUpClass(cls) -> None:
        source_root = ROOT / "src/curvatureCorrectionModels"
        cls.source = (source_root / "common/curvatureFunctions.C").read_text()
        cls.model_sources = "\n".join(
            (source_root / path).read_text()
            for path in (
                "kOmegaSSTCC/kOmegaSSTCC.C",
                "SpalartAllmarasCC/SpalartAllmarasCC.C",
            )
        )
        cls.model_headers = "\n".join(
            (source_root / path).read_text()
            for path in (
                "kOmegaSSTCC/kOmegaSSTCC.H",
                "SpalartAllmarasCC/SpalartAllmarasCC.H",
            )
        )

    def test_sst_denominator_and_clipping_contract(self) -> None:
        self.assertIn("max(sqr(strainRate), 0.09*sqr(omega))", self.source)
        self.assertIn("rotationRate*D2*sqrt(D2)", self.source)
        self.assertIn("min(max(fRotation, frMin), frMax)", self.source)

    def test_sa_denominator_and_unbounded_contract(self) -> None:
        self.assertIn("0.5*(sqr(strainRate) + sqr(rotationRate))", self.source)
        self.assertIn("(2.0*WS && DSDt)/max(sqr(D2)", self.source)
        self.assertIn("SA-RC is deliberately not bounded", self.source)

    def test_material_derivative_contains_time_and_convection(self) -> None:
        self.assertIn("symm(fvc::grad(fvc::ddt(U)))", self.source)
        self.assertIn(
            'fvc::div(phi, strain, "div(phi,curvatureStrain)")',
            self.source,
        )
        self.assertIn("- strain*fvc::div(phi)", self.source)

    def test_safe_defaults_and_diagnostics_contract(self) -> None:
        self.assertIn('type == word("kOmegaSST2003RC")', self.model_sources)
        self.assertIn(
            'type == word("SpalartAllmarasNoFt2RC")', self.model_sources
        )
        self.assertIn('TypeName("kOmegaSST2003RC")', self.model_headers)
        self.assertIn('TypeName("SpalartAllmarasNoFt2RC")', self.model_headers)
        self.assertIn("makeRASModel(kOmegaSST2003RC)", (ROOT / "src/curvatureCorrectionModels/kOmegaSSTCC/makekOmegaSSTCC.C").read_text())
        self.assertIn("makeRASModel(SpalartAllmarasNoFt2RC)", (ROOT / "src/curvatureCorrectionModels/SpalartAllmarasCC/makeSpalartAllmarasCC.C").read_text())
        self.assertIn("ft2 = scalar(0)", self.model_sources)
        self.assertIn("requires ft2 false", self.model_sources)
        self.assertEqual(self.model_sources.count("Switch::getOrAddToDict"), 4)
        self.assertEqual(self.model_sources.count("\"writeCurvatureDiagnostics\""), 4)
        for field in (
            "curvatureIndicator_",
            "productionBase_",
            "productionCorrected_",
            "productionRatio_",
            "fr1Clipped_",
        ):
            self.assertGreaterEqual(self.model_sources.count(field), 2)
        self.assertGreaterEqual(self.model_sources.count("std::isfinite"), 2)
        self.assertIn("requireFinite(rHat", self.source)
        self.assertIn("requireFinite(rTilde", self.source)

    def test_mrf_source_contract(self) -> None:
        self.assertIn("mesh.csorted<IOMRFZoneList>()", self.source)
        self.assertIn("const volTensorField frameRotation(*OmegaRot)", self.source)
        self.assertIn("(frameRotation & strain) - (strain & frameRotation)", self.source)
        self.assertEqual(
            self.model_sources.count("updateMRFAngularVelocity(OmegaRot_)"), 2
        )


if __name__ == "__main__":
    unittest.main()
