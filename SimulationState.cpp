#include "SimulationState.h"

#include "forces.h"
#include "velocities.h"
#include "helpers.h"


void SimulationState::simulate_one_step(double dt, double shear) {

    ///////////////////////////////////
    // Compute forces and velocities //
    ///////////////////////////////////

    forces_t f = forces::non_bond_forces(shear, h, p);
    // add forces of each bond
    for (auto & bd_i : bd_lig_ind)
        f += ligands[bd_i].bond_forces(h, alpha_0);

    velocities_t v = velocities::compute_velocities(h, f, p);


    //////////////////////////////
    // Prepare changes in bonds //
    //////////////////////////////

    vector<size_t> new_bondings_lig_ind;
    // We try with all ligands, including those already bonded!
    // But `prepare_binding` will check it.
    for (size_t i = 0; i < ligands.size(); i++)
        if (ligands[i].prepare_binding(h, alpha_0, dt, generator))
            new_bondings_lig_ind.push_back(i);

    // Here we try only with bonded ligands.
    vector<size_t> new_ruptures_lig_ind;
    for (auto & bd_i : bd_lig_ind)
        if(ligands[bd_i].prepare_rupture(h, alpha_0, dt, generator))
            new_ruptures_lig_ind.push_back(bd_i);


    //////////
    // Move //
    //////////

    // update height
    h += dt * v.v_y;

    // move surface in x direction (sphere's center is always at origin),
    // i.e. move bonded receptors on surface in opposite direction
    for (auto & bd_i : bd_lig_ind)
        ligands[bd_i].move_bd_rec(- dt * v.v_x);

    // rotate sphere
    alpha_0 += dt * v.o_z;


    ////////////////////////
    // Apply bond changes //
    ////////////////////////

    for (auto & new_bon_i : new_bondings_lig_ind) {
        bd_lig_ind.insert(new_bon_i);
        ligands[new_bon_i].bond(alpha_0);
    }

    for (auto & new_rup_i : new_ruptures_lig_ind) {
        bd_lig_ind.erase(new_rup_i);
        ligands[new_rup_i].rupture();
    }

}

SimulationState::SimulationState(double h_0, Parameters *p_, size_t seed) {
    h = h_0;
    p = p_;
    reseed(seed);

    LigandParameters* lig_p;
    size_t n_of_lig;
    for (auto& lig_type : p->lig_types) {
        lig_p = lig_type.first;
        n_of_lig = lig_type.second;
        for (size_t i; i < n_of_lig; i++) {
            xyz_t lig_xyz = helpers::draw_from_uniform_dist_on_sphere(0, p->r_c);
            // TODO: check if |z| is small
            // TODO: keep model's parameters in separate class
            ligands.emplace_back(lig_xyz.x, lig_xyz.y, lig_p, p);
        }
        // TODO: sort ligands by r_cir
    }

}
