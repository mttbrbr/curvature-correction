#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate the non-singular full-cylinder O-grid for the rotating pipe."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


X_PLANES = (0.0, 9.0, 10.5, 10.809, 12.5710162184)
AXIAL_CELLS = {
    "smoke": (600, 100, 21, 118),
    "quick": (800, 150, 32, 163),
    "medium": (1275, 240, 51, 260),
}
TRANSVERSE_CELLS = {
    "smoke": (3, 9),
    "quick": (4, 12),
    "medium": (6, 18),
}


def vertex(plane: int, local: int) -> int:
    return 8 * plane + local


def mesh_text(level: str) -> str:
    axial = AXIAL_CELLS[level]
    n_core, n_radial = TRANSVERSE_CELLS[level]
    radius = 0.03
    core_corner_radius = 0.012
    a = core_corner_radius / math.sqrt(2.0)
    outer = radius / math.sqrt(2.0)

    # Four inner-square and four circular-boundary vertices, counter-clockwise.
    cross_section = (
        (a, a),
        (-a, a),
        (-a, -a),
        (a, -a),
        (outer, outer),
        (-outer, outer),
        (-outer, -outer),
        (outer, -outer),
    )
    # Counter-clockwise quadrilaterals: square core, then four annular sectors.
    quads = (
        (2, 3, 0, 1),
        (1, 0, 4, 5),
        (2, 1, 5, 6),
        (3, 2, 6, 7),
        (0, 3, 7, 4),
    )

    lines = [
        "/* SPDX-License-Identifier: GPL-3.0-or-later */",
        "FoamFile",
        "{",
        "    format ascii;",
        "    class dictionary;",
        "    object blockMeshDict;",
        "}",
        "",
        "scale 1;",
        "",
        "vertices",
        "(",
    ]
    for x in X_PLANES:
        for y, z in cross_section:
            lines.append(f"    ({x:.12g} {y:.12g} {z:.12g})")
    lines.extend((");", "", "blocks", "("))

    for slab, n_axial in enumerate(axial):
        zone = " rotating" if slab > 0 else ""
        for block, quad in enumerate(quads):
            q0, q1, q2, q3 = quad
            vertices = (
                vertex(slab, q0),
                vertex(slab + 1, q0),
                vertex(slab + 1, q1),
                vertex(slab, q1),
                vertex(slab, q3),
                vertex(slab + 1, q3),
                vertex(slab + 1, q2),
                vertex(slab, q2),
            )
            radial = n_core if block == 0 else n_radial
            radial_grading = "1" if block == 0 else "0.2"
            vertex_list = " ".join(str(value) for value in vertices)
            lines.append(
                f"    hex ({vertex_list}){zone} "
                f"({n_axial} {n_core} {radial}) "
                f"simpleGrading (1 1 {radial_grading})"
            )
    lines.extend((");", "", "edges", "("))
    for plane, x in enumerate(X_PLANES):
        for first, second in ((4, 5), (5, 6), (6, 7), (7, 4)):
            lines.append(
                f"    arc {vertex(plane, first)} {vertex(plane, second)} "
                f"origin ({x:.12g} 0 0)"
            )
    lines.extend((");", "", "boundary", "("))

    lines.extend(("    inlet", "    {", "        type patch;", "        faces", "        ("))
    for q0, q1, q2, q3 in quads:
        lines.append(
            f"            ({vertex(0, q0)} {vertex(0, q3)} "
            f"{vertex(0, q2)} {vertex(0, q1)})"
        )
    lines.extend(("        );", "    }"))

    last = len(X_PLANES) - 1
    lines.extend(("    outlet", "    {", "        type patch;", "        faces", "        ("))
    for q0, q1, q2, q3 in quads:
        lines.append(
            f"            ({vertex(last, q0)} {vertex(last, q1)} "
            f"{vertex(last, q2)} {vertex(last, q3)})"
        )
    lines.extend(("        );", "    }"))

    for patch_name, slabs in (("stationaryWall", range(0, 1)), ("rotatingWall", range(1, 4))):
        lines.extend(
            (
                f"    {patch_name}",
                "    {",
                "        type wall;",
                "        faces",
                "        (",
            )
        )
        for slab in slabs:
            for q0, q1, q2, q3 in quads[1:]:
                lines.append(
                    f"            ({vertex(slab, q3)} {vertex(slab + 1, q3)} "
                    f"{vertex(slab + 1, q2)} {vertex(slab, q2)})"
                )
        lines.extend(("        );", "    }"))

    lines.extend((");", "", "mergePatchPairs ();", ""))
    return "\n".join(lines)


def expected_cells(level: str) -> int:
    n_core, n_radial = TRANSVERSE_CELLS[level]
    cells_per_plane = n_core * n_core + 4 * n_core * n_radial
    return sum(AXIAL_CELLS[level]) * cells_per_plane


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--level", choices=sorted(AXIAL_CELLS), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(mesh_text(args.level), encoding="utf-8")
    print(f"Wrote {args.output} ({expected_cells(args.level)} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
