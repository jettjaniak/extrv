#pragma once

#include "types.h"
#include "Parameters.h"
#include "AbstractBondType.h"

/**
 * Each ligand has xyz coordinates in the model.
 * Coordinate z is defined only up to a sign.
 * We take a circle that is a slice of a sphere at z = [lig. z coord.].
 * r_cir is a radius of that circle.
 * alpha is an angle from the bottommost circle point to the ligand.
 */
struct Ligand {
    /// radius of the circle in μm
    double r_cir;
    /**
     * alpha - rot, where alpha is the angle
     * from the bottommost circle point to the ligand
     * rot_inc is in [0, 2π]
     */
    double rot_inc;
    /// value 0 indicates no bonding, higher values indicate bonding to different receptors
    int bond_state = 0;

    vector<double> binding_rates;

    /**
     * x coordinate of receptor currently bonded to the ligand,
     * valid only when ligand is bonded
     */
    double bd_rec_x = INFTY;

    /// type of ligand
    Parameters::LigandType* lig_type;

    /**
     * Ligand constructor.
     *
     * @param lig_xy x coordinate of point on the sphere
     * @param lig_y y coordinate of point on the sphere
     * @param lig_type ligand's parameters
     * @param p_ model's parameters
     */
    Ligand(xy_t lig_xy, Parameters::LigandType *lig_type);

    /**
     * Computes x coordinate of ligand.
     * @param rot sphere's rotation
     */
    double x_pos(double rot) const;

    /**
     * Computes y coordinate of ligand.
     * @param rot sphere's rotation
     */
    double y_pos(double rot) const;

    /**
     * Compute distance from ligand to surface.
     *
     * @param h distance from sphere to surface
     * @param rot sphere's rotation
     * @return distance from ligand to surface
     */
    double surface_dist(double h, double rot);

    /**
     * Compute bond length.
     * Should abort if ligand is not bonded (bond_state == 0).
     *
     * @param h distance from sphere to surface
     * @param rot sphere's rotation
     * @return bond length in μm
     */
    double bond_length(double h, double rot);

    /**
     * Returns parameters of current bond.
     */
    AbstractBondType* get_curr_bond_type();

    double update_binding_rates(double h, double rot);

    double rupture_rate(double h, double rot);

    /**
     * Computes bonding probability and draws if bonding will happen.
     * If ligand is bonded (bond_state != 0), always returns false.
     * If there are more than one receptor to bond to, it chooses one accordingly to probability distribution.
     * If bonding happens, it saves all bonding information except for bd_rec_x - it has to be updated after
     * the sphere is moved, because we assume bonding will happen at the end of time step.
     *
     * @param h distance from sphere to surface
     * @param rot sphere's rotation
     * @param dt simulation time step
     * @param generator random number generator
     * @return true if bonding will happen, false otherwise
     */
    bool prepare_binding(double h, double rot, double dt, generator_t &generator);

    /**
     * Computes bond rupture probability and draws if rupture will happen.
     * Should abort if ligand is not bonded (bond_state == 0).
     *
     * @param h distance from sphere to surface
     * @param rot sphere's rotation
     * @param dt simulation time step
     * @param generator random number generator
     * @return true if bond rupture will happen, false otherwise
     */
    bool prepare_rupture(double h, double rot, double dt, generator_t &generator);

    /**
     * Creates previously prepared bond.
     * It will set `bd_rec_x` and `bond_state` to appropriate values
     * and `prepared_bond_state` to -1.
     *
     * @param rot sphere's rotation
     */
    void bond(double rot);

    /**
     * Ruptures current bond.
     * It will set `bond_state` to 0 and `bd_rec_x` to infinity.
     */
    void rupture();

    /**
     * Change position of receptor bonded to this ligand.
     *
     * In our frame of reference cell center is at (0, 0),
     * so when we move in x direction by x_dist we have to move receptor by - x_dist.
     *
     * @param x_dist how much cell moved in x direction in μm
     */
    void move_bd_rec(double x_dist);

    /**
     * Computes forces that bond exerts on the sphere.
     *
     * @param h distance from sphere to surface
     * @param rot sphere's rotation
     * @return forces and torques exerted on the sphere
     */
    forces_t bond_forces(double h, double rot);
};



