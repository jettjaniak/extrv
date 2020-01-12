#pragma once

#include "types.h"

namespace helpers {
    /**
     * Compute (r_cir, alpha) parametrization from (x, y) coordinates of point laying on a sphere.
     * @param lig_x x ligand coordinate
     * @param lig_y y ligand coordinate
     * @return (r_cir, alpha), where alpha is in [-π, π]
     */
    pair<double, double> parametrize_ligand(double lig_x, double lig_y);
}