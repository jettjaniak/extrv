#pragma once

#include "types.h"
#include "Ligand.h"
#include "Settings.h"

struct Stats {
    vector<size_t> n_bd_lig_vec;

    Stats() = default;

    explicit Stats(int n_lig_types) {
        n_bd_lig_vec.resize(n_lig_types);
    }

};

class SimulationState {
public:
    // distance from sphere to surface
    double h;
    // sphere's rotation
    double rot = 0.0;

    vector<Ligand> ligands;
    // indices of bonded ligands
    set<size_t> bd_lig_ind;

    // random number generator
    generator_t generator;

    Settings* settings;
    Stats stats;

    SimulationState() = default;
    SimulationState(double h_0, Settings* settings_, unsigned int seed);

    void simulate_one_step(double dt, double shear);

    void reseed(unsigned int seed) {
        generator = generator_t{seed};
    }

    void update_stats();
};



