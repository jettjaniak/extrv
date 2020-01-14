#pragma once

#include "types.h"
#include "Parameters.h"


namespace velocities {
    /**
     * Computes velocities from forces.
     *
     * In this model, for a given distance between surface and sphere,
     * sphere velocities (translational and angular) can be determined
     * from forces acting on the sphere.
     *
     * @param h sphere's height above the surface (distance)
     * @param f forces acting on the sphere
     * @param p model parameters
     *
     * @return computed velocities
     *
     */
    velocities_t compute_velocities(double h, const forces_t & f, const Parameters* p);
}