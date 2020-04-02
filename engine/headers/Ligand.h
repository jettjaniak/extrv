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
     * TODO
     */
    double surface_dist(const array<double, 3> & pos);

    /**
     * Compute bond length.
     * Should abort if ligand is not bonded (bond_state == 0).
     *
     * TODO
     */
    double bond_length(const array<double, 3> & pos);

    /**
     * Returns parameters of current bond.
     */
    AbstractBondType* get_curr_bond_type();

    double update_binding_rates(const array<double, 3> & pos);

    double rupture_rate(const array<double, 3> & pos);

    /**
     * TODO
     */
    void bond(double rot, double dist, generator_t &generator);

    /**
     * Ruptures current bond.
     * It will set `bond_state` to 0 and `bd_rec_x` to infinity.
     */
    void rupture();

    /**
     * Computes forces that bond exerts on the sphere.
     *
     * TODO
     */
    forces_t bond_forces(const array<double, 3> & pos);
};



