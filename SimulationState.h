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
    // model's parameters
    Parameters* p;

    // TODO: check seed type for default_random_engine
    SimulationState(double h_0, Parameters* p_, size_t seed) {
        h = h_0;
        p = p_;
        generator = generator_t(seed);
    }
};



