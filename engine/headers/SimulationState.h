#pragma once

#include "types.h"
#include "Ligand.h"
#include "Settings.h"
#include "helpers.h"
#include "History.h"


struct SimulationState {
    // distance from sphere to surface
    double h;
    // sphere's rotation
    double rot = 0.0;

    // TODO: documentation
    vector<Ligand> ligands;
    // indices of bonded ligands
    set<size_t> bd_lig_ind;

    // random number generator
    generator_t generator;

    // TODO: documentation
    Settings* settings;

    // TODO: documentation
    SimulationState(double h_0, Settings* settings_, unsigned int seed);

    // TODO: documentation
    void simulate_one_step(double dt, double shear);

    // TODO: documentation
    void simulate(size_t n_steps, double dt, double shear);

    // TODO: documentation
    History simulate_with_history(size_t n_steps, double dt, double shear, size_t save_every = 1000);

    // TODO: generate seed if seed = 0 (timestamp won't work in multiprocessing!)
    void reseed(unsigned int seed);
};