#pragma once

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"
#include "helpers.h"
#include "History.h"


struct SimulationState {
    /// distance from sphere to surface in μm
    double h;
    /// sphere's rotation in radians
    double rot = 0.0;

    /// ligands on sphere
    vector<Ligand> ligands;

    // indices of leftmost and rightmost ligands that are close enough to surface to bond
    size_t left_lig_ind = 0;
    size_t right_lig_ind = 0;

    /// indices of bonded ligands
    set<size_t> bd_lig_ind;

    /// random number generator
    generator_t generator;

    /// model parameters
    Parameters* p;

    /**
     * Initialize simulation.
     *
     * @param h_0 starting height above surface
     * @param p model parameters
     * @param seed number that random number generator will be seeded with
     */
    SimulationState(double h_0, Parameters* p, unsigned int seed);

    /**
     * Do one step of simulation.
     *
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     */
    void simulate_one_step(double dt, double shear_rate);

    /**
     * Do n_steps of simulation.
     *
     * @param n_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     */
    void simulate(size_t n_steps, double dt, double shear_rate, bool stop_if_no_bonds = false);

    /**
     * Do n_steps of simulation and return history updated every save_every steps.
     *
     * @param n_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     * @param save_every number of steps between history updates
     * @return simulation history
     */
    History simulate_with_history(size_t n_steps, double dt, double shear, bool stop_if_no_bonds = false,
            size_t save_every = 1000);

    /// reseed random number generator
    void reseed(unsigned int seed);

    /// Check if given ligand is close enough to surface to consider its binding.
    bool surf_dist_small(size_t lig_ind);

    /**
     * Update value of one side of range of ligands that are close to surface.
     *
     * @param curr_ind one side of range
     * @param step +1 for right side of range, -1 for left
     * @return updated index
     */
    void update_one_side_of_range(size_t & curr_ind, int step);
};