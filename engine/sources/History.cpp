#include "History.h"

#include <utility>
#include <algorithm>

#include "velocities.h"
#include "SimulationState.h"


History::BondTrajectory::BondTrajectory(size_t start_i) : start_i(start_i) {}

void History::update(const SimulationState *s) {
    update_bond_trajectories(s);
    time.push_back(s->time);
    h.push_back(s->h());
    rot.push_back(s->global_rot());
    dist.push_back(s->global_dist());
    hist_i++;
}

void History::finish() {
    for (auto &value_p : active_trajs_map)
        bond_trajectories.emplace_back(value_p.second);
    active_trajs_map.clear();
}

void History::update_bond_trajectories(const SimulationState *s) {

    const set<size_t> & curr_blis = s->bd_lig_ind;
    vector<size_t>::iterator it;

    // just bonded ligand indices
    vector<size_t> bon_lis(curr_blis.size());
    // bon_lis = curr_blis \ prev_blis
    it = std::set_difference(
            curr_blis.cbegin(), curr_blis.cend(),
            prev_blis.cbegin(), prev_blis.cend(),
            bon_lis.begin()
    );
    bon_lis.resize(it - bon_lis.begin());

    // just ruptured ligand indices
    vector<size_t> rup_lis(prev_blis.size());
    // rup_lis = prev_blis \ curr_blis
    it = std::set_difference(
            prev_blis.cbegin(), prev_blis.cend(),
            curr_blis.cbegin(), curr_blis.cend(),
            rup_lis.begin()
    );
    rup_lis.resize(it - rup_lis.begin());

    for (size_t bli : bon_lis)
        active_trajs_map.emplace(bli, BondTrajectory(hist_i));

    for (size_t rli : rup_lis) {
        auto active_traj_map_it = active_trajs_map.find(rli);
        // not found, shouldn't happen
        if (active_traj_map_it == active_trajs_map.end())
            abort();
        auto &active_traj = (*active_traj_map_it).second;
        bond_trajectories.emplace_back(active_traj);
        active_trajs_map.erase(active_traj_map_it);
    }

    for (size_t cli : curr_blis) {
        auto active_traj_it = active_trajs_map.find(cli);
        // not found, shouldn't happen
        if (active_traj_it == active_trajs_map.end())
            abort();
        vector<xy_t> &active_traj_pos = (*active_traj_it).second.positions;
        const Ligand &ligand = s->ligands[cli];
        active_traj_pos.emplace_back(ligand.x_pos(s->ode_x[POS_ROT]), ligand.y_pos(s->ode_x[POS_ROT]));
    }

    prev_blis = curr_blis;
}