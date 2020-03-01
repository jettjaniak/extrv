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
    void simulate(size_t n_steps, double dt, double shear_rate);

    /**
     * Do n_steps of simulation and return history updated every save_every steps.
     *
     * @param n_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     * @param save_every number of steps between history updates
     * @return simulation history
     */
    History simulate_with_history(size_t n_steps, double dt, double shear, size_t save_every = 1000);

    /// reseed random number generator
    void reseed(unsigned int seed);
};