#pragma once

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"
#include "helpers.h"


struct AbstractSimulationState;

/**
 * Builds and represents simulation history.
 */
class History {
public:
    /**
     * Stores index for which given bond was first observed and
     * subsequent positions of bonded ligand.
     */
    struct BondTrajectory {
        // TODO: add bond type
        size_t start_i;
        vector<xy_t> positions;

        explicit BondTrajectory(size_t start_i);
    };

    vector<double> time;
    /// vector of heights above surface in μm
    vector<double> h;
    /// vector of sphere's rotation in radians
    vector<double> rot;
    /// distances traveled in direction of fluid flow (x)
    vector<double> dist;
    /// here we move values from active_trajs_map after rupture or finish
    vector<BondTrajectory> bond_trajectories;

    History() = default;

    /// update history using current data from simulation state
    void update(const AbstractSimulationState* s);

    /// before calling this function simulation results are not complete
    void finish();

private:
    /// number of calls to update
    size_t hist_i = 0;
    /// bonded ligand indices from previous update
    set<size_t> prev_blis;
    /// here we store ligands that are continuously bonded, key: ligand index
    unordered_map<size_t, BondTrajectory> active_trajs_map;

    /**
     * Find which bonds are new, which ruptured, update positions etc.
     */
    void update_bond_trajectories(const AbstractSimulationState *s);
};