#include "../headers/SimulationState.h"

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
            // TODO: do it wisely
            if (std::abs(lig_xyz.z) < 0.3 && lig_xyz.y < -4.45) {
                xy_t lig_xy{lig_xyz};
                ligands.emplace_back(lig_xy, lig_p, settings->p);
            }
        }
        // TODO: sort ligands by r_cir
    }
}

void Stats::update() {
    bd_ligs_ind_and_xy.clear();
    // Reset number of bonds of each ligand type to zero.
    std::fill(n_bd_lig_vec.begin(), n_bd_lig_vec.end(), 0);
    for (auto lig_i : s->bd_lig_ind) {
        auto & ligand = s->ligands[lig_i];
        xy_t lig_xy_pos = xy_t{ligand.x_pos(s->rot), ligand.y_pos(s->rot)};
        bd_ligs_ind_and_xy.emplace_back(lig_i, lig_xy_pos);
        int this_lig_type_ind = ligand.lig_type->index_in_settings;
        n_bd_lig_vec[this_lig_type_ind]++;
    }
}

Stats::Stats(const SimulationState* s_) {
    s = s_;
    n_bd_lig_vec.resize(s->settings->lig_types_and_nrs.size());
    update();
}
