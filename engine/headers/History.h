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

    // TODO: documentation
    explicit History(const SimulationState *s);

    // TODO: documentation
    void update();

    // before calling this function simulation results are not complete
    void finish();

private:
    size_t hist_i = 0;
    // bonded ligand indices from previous update
    set<size_t> prev_blis;
    // here we store ligands that are continuously bonded, key: ligand index
    unordered_map<size_t, BondTrajectory> active_trajs_map;

    // TODO: documentation
    void update_bond_trajectories();
};