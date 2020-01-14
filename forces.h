#pragma once

#include "types.h"
#include "Parameters.h"


namespace forces {
    /**
     * Computes forces acting on the sphere from all sources other than bonds.
     *
     * @param shear shear rate of fluid flow
     * @param h distance between sphere and surface
     * @param p model's parameters
     */
    forces_t non_bond_forces(double shear, double h, const Parameters *p);
}