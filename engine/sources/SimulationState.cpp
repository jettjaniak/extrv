#include "SimulationState.h"

#include <algorithm>
#include <utility>
#include <iostream>

#include "forces.h"
#include "velocities.h"
#include "helpers.h"


SimulationState::SimulationState(double h_0, Parameters* p, unsigned int seed) : h(h_0), p(p) {
    reseed(seed);

    Parameters::LigandType* lig_p;
    size_t n_of_lig;
    for (auto& lig_type : p->lig_types_and_nrs) {
        lig_p = lig_type.first;
        n_of_lig = lig_type.second;
        for (size_t i = 0; i < n_of_lig; i++) {
            xyz_t lig_xyz = helpers::draw_from_uniform_dist_on_sphere(p->r_cell, generator);
            // TODO: do it wisely, use EPS_PROB
            if (std::abs(lig_xyz.z) < 0.1) {
                xy_t lig_xy{lig_xyz};
                ligands.emplace_back(lig_xy, lig_p);
            }
        }
        // TODO: sort ligands by alpha_inc
        // TODO @Kajetan: indicate which ligands have chance of bonding,
        //   by specifying range of indices (or iterators)
    }
}

void SimulationState::simulate_one_step(double dt, double shear_rate) {

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