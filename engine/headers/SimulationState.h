#pragma once

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"
#include "helpers.h"
#include "History.h"


struct SimulationState {
    /// local height, distance and rotation
    Position pos;

    /// sphere's rotation in radians
    double global_rot = 0.0;
    /// distance traveled in direction of flow (x)
    double global_dist = 0.0;

    double shear_rate = 0.0;

    /// ligands on sphere
    vector<Ligand> ligands;

    /// range of ligands' rot_incs for which ligands can bond at given step
    double left_rot_inc = 0.0;
    double right_rot_inc = 0.0;

    /// indices of leftmost and rightmost ligands that are close enough to surface to bond
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
    void simulate_one_step();

    /**
     * Do n_steps of simulation.
     *
     * @param n_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     */
    void simulate(size_t n_steps);

    /**
     * Do n_steps of simulation and return history updated every save_every steps.
     *
     * @param n_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     * @param save_every number of steps between history updates
     * @return simulation history
     */
    History simulate_with_history(size_t n_steps, size_t save_every);

    /// reseed random number generator
    void reseed(unsigned int seed);

    void update_rot_inc_range();

    void update_rot_inc_ind();

    void rhs(const Position & x, Position & dxdt, double /*t*/);
};