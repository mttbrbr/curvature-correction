#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate a periodic full-cylinder O-grid for the Xiao rotating pipe."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

RADIUS = 0.5
LENGTH = 15.0 * RADIUS
LEVELS = {
    "smoke": (64, 3, 9),
    "quick": (128, 6, 18),
    "medium": (256, 12, 36),
}


def vertex(plane: int, local: int) -> int:
    return 8 * plane + local


def expected_cells(level: str) -> int:
    axial, core, radial = LEVELS[level]
    return axial * (core * core + 4 * core * radial)


def mesh_text(level: str) -> str:
    axial, core, radial = LEVELS[level]
    inner = 0.4 * RADIUS / math.sqrt(2.0)
    outer = RADIUS / math.sqrt(2.0)
    cross_section = (
        (inner, inner), (-inner, inner), (-inner, -inner), (inner, -inner),
        (outer, outer), (-outer, outer), (-outer, -outer), (outer, -outer),
    )
    quads = (
        (2, 3, 0, 1),
        (1, 0, 4, 5),
        (2, 1, 5, 6),
        (3, 2, 6, 7),
        (0, 3, 7, 4),
    )
    lines = [
        "/* SPDX-License-Identifier: GPL-3.0-or-later */",
        "FoamFile { format ascii; class dictionary; object blockMeshDict; }",
        "scale 1;", "", "vertices", "(",
    ]
    for x in (0.0, LENGTH):
        for y, z in cross_section:
            lines.append(f"    ({x:.12g} {y:.12g} {z:.12g})")
    lines.extend((");", "", "blocks", "("))
    for block, (q0, q1, q2, q3) in enumerate(quads):
        vertices = (
            vertex(0, q0), vertex(1, q0), vertex(1, q1), vertex(0, q1),
            vertex(0, q3), vertex(1, q3), vertex(1, q2), vertex(0, q2),
        )
        transverse = core if block == 0 else radial
        grading = "1" if block == 0 else "0.2"
        lines.append(
            "    hex (" + " ".join(map(str, vertices)) + ") rotating "
            f"({axial} {core} {transverse}) simpleGrading (1 1 {grading})"
        )
    lines.extend((");", "", "edges", "("))
    for plane, x in enumerate((0.0, LENGTH)):
        for first, second in ((4, 5), (5, 6), (6, 7), (7, 4)):
            lines.append(
                f"    arc {vertex(plane, first)} {vertex(plane, second)} "
                f"origin ({x:.12g} 0 0)"
            )
    lines.extend((");", "", "boundary", "("))
    for name, plane, order, neighbour in (
        ("periodicIn", 0, (0, 3, 2, 1), "periodicOut"),
        ("periodicOut", 1, (0, 1, 2, 3), "periodicIn"),
    ):
        lines.extend((
            f"    {name}", "    {", "        type cyclic;",
            f"        neighbourPatch {neighbour};", "        faces", "        (",
        ))
        for quad in quads:
            lines.append("            (" + " ".join(
                str(vertex(plane, quad[index])) for index in order
            ) + ")")
        lines.extend(("        );", "    }"))
    lines.extend(("    wall", "    {", "        type wall;", "        faces", "        ("))
    for q0, q1, q2, q3 in quads[1:]:
        lines.append(
            f"            ({vertex(0, q3)} {vertex(1, q3)} "
            f"{vertex(1, q2)} {vertex(0, q2)})"
        )
    lines.extend(("        );", "    }", ");", "", "mergePatchPairs ();", ""))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--level", choices=tuple(LEVELS), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(mesh_text(args.level), encoding="utf-8")
    print(f"Wrote {args.output} ({expected_cells(args.level)} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
