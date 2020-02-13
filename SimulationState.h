#pragma once

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"

class SimulationState {
public:
    // distance from sphere to surface
    double h;
    // sphere's rotation
    double alpha_0 = 0.0;

    vector<Ligand> ligands;
    // indices of bonded ligands
    set<size_t> bd_lig_ind;

    // random number generator
    generator_t generator;

    SimulationSettings* settings;

    // TODO: check seed type for default_random_engine
    SimulationState(double h_0, SimulationSettings* settings_, size_t seed);

    void simulate_one_step(double dt, double shear);

    void reseed(size_t seed) {
        generator = std::default_random_engine(seed);
    }
};



