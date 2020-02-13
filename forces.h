#pragma once

#include "types.h"
#include "SimulationSettings.h"


namespace forces {
    /**
     * Computes forces acting on the sphere from all sources other than bonds.
     *
     * @param shear_rate shear rate of fluid flow
     * @param h distance between sphere and surface
     * @param p model's parameters
     */
    forces_t non_bond_forces(double shear_rate, double h, const Parameters *p);

    /**
     * Computes forces that shear flow exerts on the sphere.
     *
     * @param shear_rate shear rate of fluid flow
     * @param h distance between sphere and surface
     * @param p model's parameters
     */
    forces_t shear_forces(double shear_rate, double h, const Parameters *p);

    /**
     * Computes repulsive force between sphere and surface.
     *
     * @param h distance between sphere and surface
     * @param p model's parameters
     */
    double repulsive_force(double h, const Parameters *p);

    /**
     * Compute gravitational force acting in y direction, with correct sign.
     * Archimedes rule.
     *
     * @param p model's parameters
     */
    double grav_force(const Parameters *p);
}