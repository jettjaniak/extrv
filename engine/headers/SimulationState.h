#pragma once

#include "types.h"
#include "Ligand.h"
#include "Settings.h"
#include "helpers.h"


struct SimulationState {
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

    SimulationState(double h_0, Settings* settings_, unsigned int seed);

    void simulate_one_step(double dt, double shear);

    // TODO: simulate, simulate_with_history

    // TODO: generate seed if seed = 0 (timestamp won't work in multiprocessing!)
    void reseed(unsigned int seed);
};


struct History {
    // TODO: documentation, documentation, documentation
    struct active_traj_t {
        size_t start_i;
        vector<xy_t> positions;

        active_traj_t(size_t start_i_) {
            start_i = start_i_;
        }
    };

    struct final_traj_t {
        size_t start_i;
        size_t n_of_pos;
        double* pos_ptr;

        final_traj_t(const active_traj_t& active_traj) {
            start_i = active_traj.start_i;
            n_of_pos = active_traj.positions.size();
            pos_ptr = (double*) active_traj.positions.data();
        }
    };

    const SimulationState* s;
    size_t hist_i = 0;
    // bonded ligand indices from previous update
    set<size_t> prev_blis;
    // just bonded ligand indices
    vector<size_t> bon_lis;
    // just ruptured ligand indices
    vector<size_t> rup_lis;

    // here we store ligands that are continuously bonded, key: ligand index
    unordered_map<size_t, active_traj_t> active_trajs_map;
    // here we move values from active_trajs_map after rupture
    vector<final_traj_t> final_trajs_vec;

    explicit History(const SimulationState* s_);

    void update();

    void finish();
};