#pragma once

#include "types.h"
#include "Ligand.h"
#include "Settings.h"
#include "helpers.h"


struct SimulationState;

class History {  // TODO: documentation, documentation, documentation
public:
    struct BondTrajectory {
        size_t start_i;
        vector<xy_t> positions;

        explicit BondTrajectory(size_t start_i);
    };

    const SimulationState* s;
    vector<double> h;
    vector<double> rot;
    // here we move values from active_trajs_map after rupture or finish
    vector<BondTrajectory> bond_trajectories;

    explicit History(const SimulationState *s);

    void update();

    // before calling this function simulation results are not complete
    void finish();

private:
    size_t hist_i = 0;
    // bonded ligand indices from previous update
    set<size_t> prev_blis;
    // here we store ligands that are continuously bonded, key: ligand index
    unordered_map<size_t, BondTrajectory> active_trajs_map;

    void update_bond_trajectories();
};

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