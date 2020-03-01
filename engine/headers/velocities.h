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
    velocities_t compute_velocities(double h, const forces_t &f, const Parameters *p);

    /**
     * Function related to hydrodynamics used in velocity computation.
     */
    double lambda_fun(double h_over_r);

    /**
     * Function related to hydrodynamics used in velocity computation.
     */
    double d_fun(double t_t_val, double f_r_val, double f_t_val, double t_r_val);
}