#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Independent scalar/tensor oracle for the implemented RC algebra.

The module intentionally uses only the Python standard library.  Tensor
components follow the published convention ``grad_u[i][j] = dU_i/dx_j``;
``from_openfoam_gradient`` converts OpenFOAM's transposed storage convention.
It does not discretise the material derivative: deterministic tests supply
``DS/Dt`` directly so the algebra can be isolated from a mesh and a scheme.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence

Tensor = tuple[tuple[float, float, float], ...]
Vector = tuple[float, float, float]


def tensor(rows: Sequence[Sequence[float]]) -> Tensor:
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        raise ValueError("a 3x3 tensor is required")
    result = tuple(tuple(float(value) for value in row) for row in rows)
    if not all(math.isfinite(value) for row in result for value in row):
        raise ValueError("tensor input contains NaN or infinity")
    return result


def vector(values: Sequence[float]) -> Vector:
    if len(values) != 3:
        raise ValueError("a 3-component vector is required")
    result = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in result):
        raise ValueError("vector input contains NaN or infinity")
    return result


ZERO: Tensor = tensor(((0, 0, 0), (0, 0, 0), (0, 0, 0)))


def transpose(value: Tensor) -> Tensor:
    return tensor(tuple(tuple(value[j][i] for j in range(3)) for i in range(3)))


def from_openfoam_gradient(value: Sequence[Sequence[float]]) -> Tensor:
    """Convert OpenFOAM ``grad(U)_ij=dU_j/dx_i`` to published indexing."""
    return transpose(tensor(value))


def symmetric(value: Tensor) -> Tensor:
    return tensor(
        tuple(tuple(0.5 * (value[i][j] + value[j][i]) for j in range(3)) for i in range(3))
    )


def antisymmetric(value: Tensor) -> Tensor:
    return tensor(
        tuple(tuple(0.5 * (value[i][j] - value[j][i]) for j in range(3)) for i in range(3))
    )


def add(first: Tensor, second: Tensor) -> Tensor:
    return tensor(
        tuple(tuple(first[i][j] + second[i][j] for j in range(3)) for i in range(3))
    )


def subtract(first: Tensor, second: Tensor) -> Tensor:
    return tensor(
        tuple(tuple(first[i][j] - second[i][j] for j in range(3)) for i in range(3))
    )


def hodge_dual(value: Sequence[float]) -> Tensor:
    x, y, z = vector(value)
    return tensor(((0, -z, y), (z, 0, -x), (-y, x, 0)))


def matmul(first: Tensor, second: Tensor) -> Tensor:
    return tensor(
        tuple(
            tuple(sum(first[i][k] * second[k][j] for k in range(3)) for j in range(3))
            for i in range(3)
        )
    )


def double_contract(first: Tensor, second: Tensor) -> float:
    return sum(first[i][j] * second[i][j] for i in range(3) for j in range(3))


def rate(value: Tensor) -> float:
    return math.sqrt(2.0 * double_contract(value, value))


@dataclass(frozen=True)
class CorrectionState:
    strain: Tensor
    rotation: Tensor
    material_derivative_strain: Tensor
    strain_rate: float
    rotation_rate: float
    r_star: float
    indicator: float
    raw_factor: float
    factor: float
    clipped: bool
    base_production: float
    corrected_production: float
    turbulent_viscosity: float
    working_variable: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def correction_state(
    grad_u: Sequence[Sequence[float]],
    material_derivative_strain: Sequence[Sequence[float]] = ZERO,
    *,
    family: str,
    frame_omega: Sequence[float] = (0.0, 0.0, 0.0),
    gradient_is_relative: bool = False,
    omega: float = 1.0,
    cr1: float = 1.0,
    cr2: float | None = None,
    cr3: float = 1.0,
    fr_min: float = 0.0,
    fr_max: float = 1.25,
    base_production: float = 1.0,
    turbulent_viscosity: float = 1.0,
    working_variable: float = 1.0,
    enabled: bool = True,
    small_rate: float = 1.0e-30,
    small_fourth_rate: float = 1.0e-120,
) -> CorrectionState:
    """Evaluate the RC algebra in inertial or rotating-frame form.

    OpenFOAM MRF stores an absolute velocity gradient, while the published
    equations use a relative gradient plus the frame rotation.  Set
    ``gradient_is_relative`` only for the latter representation.
    """
    gradient = tensor(grad_u)
    dsdt = tensor(material_derivative_strain)
    frame_rotation = hodge_dual(frame_omega)
    scalars: Iterable[float] = (
        omega,
        cr1,
        cr3,
        fr_min,
        fr_max,
        base_production,
        turbulent_viscosity,
        working_variable,
        small_rate,
        small_fourth_rate,
    )
    if not all(math.isfinite(value) for value in scalars):
        raise ValueError("scalar input contains NaN or infinity")
    if small_rate <= 0 or small_fourth_rate <= 0:
        raise ValueError("regularisation values must be positive")

    family_key = family.lower()
    if family_key not in {"sst", "sa"}:
        raise ValueError("family must be 'sst' or 'sa'")
    if cr2 is None:
        cr2 = 2.0 if family_key == "sst" else 12.0
    if not math.isfinite(cr2):
        raise ValueError("Cr2 contains NaN or infinity")
    if family_key == "sst" and (fr_min < 0 or fr_min > fr_max):
        raise ValueError("invalid SST clipping interval")

    strain = symmetric(gradient)
    rotation = antisymmetric(gradient)
    if gradient_is_relative:
        rotation = add(rotation, frame_rotation)
    dsdt = add(
        dsdt,
        subtract(matmul(frame_rotation, strain), matmul(strain, frame_rotation)),
    )
    strain_rate = rate(strain)
    rotation_rate = rate(rotation)
    r_star = strain_rate / max(rotation_rate, small_rate)
    ws = matmul(rotation, strain)
    numerator = 2.0 * double_contract(ws, dsdt)

    if family_key == "sst":
        d2 = max(strain_rate * strain_rate, 0.09 * omega * omega)
        denominator = max(rotation_rate * d2 * math.sqrt(d2), small_fourth_rate)
    else:
        d2 = 0.5 * (strain_rate * strain_rate + rotation_rate * rotation_rate)
        denominator = max(d2 * d2, small_fourth_rate)
    indicator = numerator / denominator
    raw_factor = (
        (1.0 + cr1)
        * (2.0 * r_star / (1.0 + r_star))
        * (1.0 - cr3 * math.atan(cr2 * indicator))
        - cr1
    )
    factor = min(max(raw_factor, fr_min), fr_max) if family_key == "sst" else raw_factor
    clipped = family_key == "sst" and factor != raw_factor
    if not enabled:
        factor = 1.0
        clipped = False

    return CorrectionState(
        strain=strain,
        rotation=rotation,
        material_derivative_strain=dsdt,
        strain_rate=strain_rate,
        rotation_rate=rotation_rate,
        r_star=r_star,
        indicator=indicator,
        raw_factor=raw_factor,
        factor=factor,
        clipped=clipped,
        base_production=base_production,
        corrected_production=base_production * factor,
        turbulent_viscosity=turbulent_viscosity,
        working_variable=working_variable,
    )


def rotate_tensor(value: Sequence[Sequence[float]], transform: Sequence[Sequence[float]]) -> Tensor:
    """Apply ``Q A Q^T`` for an orthogonal change of Cartesian basis."""
    q = tensor(transform)
    return matmul(matmul(q, tensor(value)), transpose(q))
