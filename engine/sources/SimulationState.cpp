#include "SimulationState.h"

#include <algorithm>
#include <utility>
#include <iostream>

#include "forces.h"
#include "velocities.h"
#include "helpers.h"


SimulationState::SimulationState(double h_0, Parameters* p, unsigned int seed) : h(h_0), p(p) {
    reseed(seed);

    Parameters::LigandType* lig_type;
    size_t n_of_lig;
    for (auto& lig_type_and_nr : p->lig_types_and_nrs) {
        lig_type = lig_type_and_nr.first;
        n_of_lig = lig_type_and_nr.second;
        std::cout << "n_of_lig = " << n_of_lig << std::endl;
        for (size_t i = 0; i < n_of_lig; i++) {
            xy_t lig_xy{helpers::draw_from_uniform_dist_on_sphere(p->r_cell, generator)};
            Ligand new_ligand {lig_xy, lig_type};
            // If ligand will be always far from surface it will never bind, so we ignore it.
            double min_surf_dist = p->r_cell - new_ligand.r_cir;
            if (min_surf_dist < lig_type->max_surf_dist)
                ligands.push_back(new_ligand);
        }
    }
    std::cout << "before sort" << std::endl;
    std::sort(ligands.begin(), ligands.end(),
         [](const Ligand & a, const Ligand & b) -> bool {
            // TODO: make rot_inc from 0 to 2 pi
            return std::fmod(a.rot_inc + 2 * PI, 2 * PI) < std::fmod(b.rot_inc + 2 * PI, 2 * PI);
         });

}

void SimulationState::simulate_one_step(double dt, double shear_rate) {

    update_one_side_of_range(right_lig_ind, 1);
    update_one_side_of_range(left_lig_ind, -1);

    ///////////////////////////////////
    // Compute forces and velocities //
    ///////////////////////////////////

    forces_t f = forces::non_bond_forces(shear_rate, h, p);
    // add forces of each bond
    for (auto & bd_i : bd_lig_ind)
        f += ligands[bd_i].bond_forces(h, rot);

    velocities_t v = velocities::compute_velocities(h, f, p);


    //////////////////////////////
    // Prepare changes in bonds //
    //////////////////////////////

    vector<size_t> new_bondings_lig_ind;
    // We try with all ligands, including those already bonded!
    // But `prepare_binding` will check it.
    for (size_t i = left_lig_ind; i != (right_lig_ind + 1) % ligands.size(); i = (i + 1) % ligands.size())
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

    double x_dist = dt * v.v_x;
    dist += x_dist;
    // move surface in x direction (sphere's center is always at origin),
    // i.e. move bonded receptors on surface in opposite direction
    for (auto & bd_i : bd_lig_ind)
        ligands[bd_i].move_bd_rec(x_dist);

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

void SimulationState::simulate(size_t n_steps, double dt, double shear_rate) {
    for (int i = 0; i < n_steps; ++i)
        simulate_one_step(dt, shear_rate);
}

History SimulationState::simulate_with_history(size_t n_steps, double dt, double shear, size_t save_every) {
    History hist;
    hist.update(this);
    for (int i = 0; i < n_steps; ++i) {
        simulate_one_step(dt, shear);
        if (i % save_every == 0)
            hist.update(this);
    }
    hist.finish();
    return hist;
}

void SimulationState::reseed(unsigned int seed) {
    generator = generator_t{seed};
}

bool SimulationState::surf_dist_small(size_t lig_ind) {
    return ligands[lig_ind].surface_dist(h, rot) < p->max_surf_dist;
}

void SimulationState::update_one_side_of_range(size_t & curr_ind, int step) {
    size_t next_ind = (ligands.size() + curr_ind + step) % ligands.size();
    bool curr_ok = surf_dist_small(curr_ind);
    bool next_ok = surf_dist_small(next_ind);

    // As long as current and next ligands are close to surface we can widen the range.
    while (curr_ok && next_ok) {
        curr_ind = next_ind;
        next_ind = (ligands.size() + curr_ind + step) % ligands.size();
        curr_ok = next_ok;
        next_ok = surf_dist_small(next_ind);
    }

    // If current ligand is not close to surface we should shrink the range.
    while (!curr_ok ) {
        curr_ind = (ligands.size() + curr_ind - step) % ligands.size();
        curr_ok = surf_dist_small(curr_ind);
    }
}


