#pragma once

#include "types.h"
#include "Parameters.h"
#include "helpers.h"

/**
 * Each ligand has xyz coordinates in the model.
 * Coordinate z is defined only up to a sign.
 * We take a circle that is a slice of a sphere at z = [lig. z coord.].
 * r_cir is a radius of that circle.
 * alpha is an angle from the bottommost circle point to the ligand.
 */
class Ligand {
private:
    // radius of the circle
    double r_cir;
    // alpha - alpha_0, where alpha is the angle
    // from the bottommost circle point to the ligand
    // alpha_inc is in [-π, π]
    double alpha_inc;
    // value 0 indicates no bonding, higher values indicate bonding to different receptors
    int bond_state = 0;

    // ligand's parameters
    LigandParameters* lig_p;

public:
    // x coordinate of receptor currently bonded to the ligand,
    // valid only when ligand is bonded
    double bd_rec_x = 0.0;

    /**
     * Ligand constructor.
     *
     * @param lig_x x coordinate of point on the sphere
     * @param lig_y y coordinate of point on the sphere
     * @param lig_p_ ligand's parameters
     */
    Ligand(double lig_x, double lig_y, LigandParameters* lig_p_)
    {
        auto r_alpha_pair = helpers::parametrize_ligand(lig_x, lig_y);
        r_cir = r_alpha_pair.first;
        alpha_inc = r_alpha_pair.second;

        lig_p = lig_p_;
    }

    /**
     * Computes bonding probability and draws if bonding will happen.
     * If there are more than one receptor to bond to, it chooses one accordingly to probability distribution.
     * If bonding happens, it saves all bonding information except for bd_rec_x.
     * It has to be updated after the sphere is moved.
     *
     * @param h distance from sphere to surface
     * @param alpha_0 sphere's rotation
     * @param dt simulation time step
     * @param p model's parameters
     * @param generator random number generator
     * @return true if bonding will happen, false otherwise
     */
    bool prepare_bonding(double h, double alpha_0, double dt, Parameters* p, generator_t generator)
    {
        // TODO
        return false;
    }

    /**
     * Computes bond rupture probability and draws if rupture will happen.
     *
     * @param h distance from sphere to surface
     * @param alpha_0 sphere's rotation
     * @param dt simulation time step
     * @param p model's parameters
     * @param generator random number generator
     * @return true if bond rupture will happen, false otherwise
     */
    bool prepare_rupture(double h, double alpha_0, double dt, Parameters* p, generator_t generator)
    {
        // TODO
        return false;
    }

    /**
     * Computes forces that bond exerts on the sphere.
     *
     * @param h distance from sphere to surface
     * @param alpha_0 sphere's rotation
     * @param p model's parameters
     * @return forces and torques exerted on the sphere
     */
    forces_t bond_force(double h, double alpha_0, Parameters* p)
    {
        // TODO
        forces_t f;
        return f;
    }
};



