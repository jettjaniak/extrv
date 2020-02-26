#include "../headers/SimulationState.h"

#include <algorithm>
#include <utility>
#include <iostream>

// TODO: add includes to CMake and Cython, so we can just type "[name].h"
#include "../headers/forces.h"
#include "../headers/velocities.h"
#include "../headers/helpers.h"


void SimulationState::simulate_one_step(double dt, double shear) {

    ///////////////////////////////////
    // Compute forces and velocities //
    ///////////////////////////////////

    forces_t f = forces::non_bond_forces(shear, h, settings->p);
    // add forces of each bond
    for (auto & bd_i : bd_lig_ind)
        f += ligands[bd_i].bond_forces(h, rot);

    velocities_t v = velocities::compute_velocities(h, f, settings->p);


    //////////////////////////////
    // Prepare changes in bonds //
    //////////////////////////////

    vector<size_t> new_bondings_lig_ind;
    // We try with all ligands, including those already bonded!
    // But `prepare_binding` will check it.
    for (size_t i = 0; i < ligands.size(); i++)
        if (ligands[i].prepare_binding(h, rot, dt, generator))
            new_bondings_lig_ind.push_back(i);

    // Here we try only with bonded ligands.
    vector<size_t> new_ruptures_lig_ind;
    for (auto & bd_i : bd_lig_ind)
        if(ligands[bd_i].prepare_rupture(h, rot, dt, generator))
            new_ruptures_lig_ind.push_back(bd_i);


    //////////
    // Move //
    //////////

    // update height
    h += dt * v.v_y;
    if (h < 0) abort();  // TODO: use throw / configure CLion to catch aborts

    // move surface in x direction (sphere's center is always at origin),
    // i.e. move bonded receptors on surface in opposite direction
    for (auto & bd_i : bd_lig_ind)
        ligands[bd_i].move_bd_rec(- dt * v.v_x);

    // rotate sphere
    rot += dt * v.o_z;


    ////////////////////////
    // Apply bond changes //
    ////////////////////////

    for (auto & new_bon_i : new_bondings_lig_ind) {
        bd_lig_ind.insert(new_bon_i);
        ligands[new_bon_i].bond(rot);
    }

    for (auto & new_rup_i : new_ruptures_lig_ind) {
        bd_lig_ind.erase(new_rup_i);
        ligands[new_rup_i].rupture();
    }
}

SimulationState::SimulationState(double h_0, Settings* settings_, unsigned int seed) {
    h = h_0;
    settings = settings_;

    reseed(seed);

    LigandType* lig_p;
    size_t n_of_lig;
    for (auto& lig_type : settings->lig_types_and_nrs) {
        lig_p = lig_type.first;
        n_of_lig = lig_type.second;
        for (size_t i = 0; i < n_of_lig; i++) {
            xyz_t lig_xyz = helpers::draw_from_uniform_dist_on_sphere(settings->p->r_c, generator);
            // TODO: do it wisely, use EPS_PROB
            if (std::abs(lig_xyz.z) < 0.1) {
                xy_t lig_xy{lig_xyz};
                ligands.emplace_back(lig_xy, lig_p, settings->p);
            }
        }
        // TODO: sort ligands by r_cir
        // TODO @Kajetan: indicate which ligands have chance of bonding,
        //   by specifying range of indices (or iterators)
    }
}


History::History(const SimulationState *s_) {
    s = s_;
}

void History::update() {
    const set<size_t> & curr_blis = s->bd_lig_ind;

    vector<size_t>::iterator it;

    // bon_lis = curr_blis \ prev_blis
    bon_lis.resize(curr_blis.size());
    it = std::set_difference(
            curr_blis.cbegin(), curr_blis.cend(),
            prev_blis.cbegin(), prev_blis.cend(),
            bon_lis.begin()
    );
    bon_lis.resize(it - bon_lis.begin());

    // rup_lis = prev_blis \ curr_blis
    rup_lis.resize(prev_blis.size());
    it = std::set_difference(
            prev_blis.cbegin(), prev_blis.cend(),
            curr_blis.cbegin(), curr_blis.cend(),
            rup_lis.begin()
    );
    rup_lis.resize(it - rup_lis.begin());

    for (size_t bli : bon_lis)
        active_trajs_map.emplace(bli, active_traj_t(hist_i));

    for (size_t rli : rup_lis) {
        auto active_traj_map_it = active_trajs_map.find(rli);
        // not found, shouldn't happen
        if (active_traj_map_it == active_trajs_map.end())
            abort();
        auto &active_traj = (*active_traj_map_it).second;
        final_trajs_vec.emplace_back(active_traj);
        active_trajs_map.erase(active_traj_map_it);
    }

    for (size_t cli : curr_blis) {
        auto active_traj_it = active_trajs_map.find(cli);
        // not found, shouldn't happen
        if (active_traj_it == active_trajs_map.end())
            abort();
        vector<xy_t> &active_traj_pos = (*active_traj_it).second.positions;
        const Ligand &ligand = s->ligands[cli];
        active_traj_pos.emplace_back(ligand.x_pos(s->rot), ligand.y_pos(s->rot));
    }

    prev_blis = curr_blis;
    hist_i++;
}

void History::finish() {
    for (auto &value_p : active_trajs_map)
        final_trajs_vec.emplace_back(value_p.second);
    active_trajs_map.clear();
}

