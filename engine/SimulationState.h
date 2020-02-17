#pragma once

#include "types.h"
#include "Ligand.h"
#include "SimulationSettings.h"

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

    SimulationState(double h_0, SimulationSettings* settings_, unsigned int seed);

    void simulate_one_step(double dt, double shear);

    void reseed(unsigned int seed) {
        generator = generator_t{seed};
    }
};


//struct SimulationStats {
//    size_t n_all_bonds = 0;
//
//    void update(const SimulationState& ss) {
//
//    }
//};

