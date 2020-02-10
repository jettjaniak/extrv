#pragma once

#include "types.h"
#include "Parameters.h"

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
    int prepared_bond_state = 0;
    // x coordinate of receptor currently bonded to the ligand,
    // valid only when ligand is bonded
    double bd_rec_x = 0.0;

    // ligand's parameters
    LigandParameters* lig_p;

    /**
     * Computes x coordinate of ligand.
     * @param alpha_0 sphere's rotation
     */
    double x_pos(double alpha_0);

    /**
     * Computes y coordinate of ligand.
     * @param alpha_0 sphere's rotation
     */
    double y_pos(double alpha_0);

    /**
     * Compute distance from ligand to surface.
     *
     * @param h distance from sphere to surface
     * @param alpha_0 sphere's rotation
     * @param p model's parameters
     * @return distance from ligand to surface
     */
    double surface_dist(double h, double alpha_0, Parameters* p);

public:
    /**
     * Ligand constructor.
     *
     * @param lig_x x coordinate of point on the sphere
     * @param lig_y y coordinate of point on the sphere
     * @param lig_p_ ligand's parameters
     */
    Ligand(double lig_x, double lig_y, LigandParameters* lig_p_);

    /**
     * Computes bonding probability and draws if bonding will happen.
     * If ligand is bonded (bond_state != 0), always returns false.
     * If there are more than one receptor to bond to, it chooses one accordingly to probability distribution.
     * If bonding happens, it saves all bonding information except for bd_rec_x - it has to be updated after
     * the sphere is moved, because we assume bonding will happen at the end of time step.
     *
     * @param h distance from sphere to surface
     * @param alpha_0 sphere's rotation
     * @param dt simulation time step
     * @param p model's parameters
     * @param generator random number generator
     * @return true if bonding will happen, false otherwise
     */
    bool prepare_binding(double h, double alpha_0, double dt, Parameters *p, generator_t generator);

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
    bool prepare_rupture(double h, double alpha_0, double dt, Parameters* p, generator_t generator);

    /**
     * Creates previously prepared bond.
     * It will set `bd_rec_x` and `bond_state`.
     *
     * @param alpha_0 sphere's rotation
     */
    void bond(double alpha_0);

    /**
     * Ruptures current bond.
     * It will set `bond_state`.
     */
    void rupture();


    void move_bd_rec(double x_dist);

    /**
     * Computes forces that bond exerts on the sphere.
     *
     * @param h distance from sphere to surface
     * @param alpha_0 sphere's rotation
     * @param p model's parameters
     * @return forces and torques exerted on the sphere
     */
    forces_t bond_forces(double h, double alpha_0, Parameters* p);
};



