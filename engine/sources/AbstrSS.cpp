#include "AbstrSS.h"

#include <algorithm>
#include <utility>
#include <iostream>


#include "forces.h"
#include "velocities.h"
#include "helpers.h"


AbstrSS::AbstrSS(double h_0, Parameters* p, unsigned int seed) : p(p) {

    pos[POS_LOG_H] = log(h_0);
    reseed(seed);

    Parameters::LigandType* lig_type;
    size_t n_of_lig;
    for (auto& lig_type_and_nr : p->lig_types_and_nrs) {
        lig_type = lig_type_and_nr.first;
        n_of_lig = lig_type_and_nr.second;
        for (size_t i = 0; i < n_of_lig; i++) {
            xy_t lig_xy{helpers::draw_from_uniform_dist_on_sphere(p->r_cell, generator)};
            Ligand new_ligand {lig_xy, lig_type};
            // If ligand will be always far from surface it will never bind, so we ignore it.
            double min_surf_dist = p->r_cell - new_ligand.r_cir;
            if (min_surf_dist < lig_type->max_surf_dist)
                ligands.push_back(new_ligand);
        }
    }
    std::sort(ligands.begin(), ligands.end(),
         [](const Ligand & a, const Ligand & b) -> bool {
            return a.rot_inc < b.rot_inc;
         });
}

double AbstrSS::h() const {
    return exp(pos[POS_LOG_H]);
}

double AbstrSS::rot() const {
    return global_rot + pos[POS_ROT];
}

double AbstrSS::dist() const {
    return global_dist + pos[POS_DIST];
}

void AbstrSS::simulate_one_step() {

    // Optimization - track only ligands that are close to surface

    update_rot_inc_range();
    update_rot_inc_ind();

    after_right_lig_ind = (right_lig_ind + 1) % ligands.size();

    if (left_lig_ind <= right_lig_ind)
        n_active_lig = right_lig_ind - left_lig_ind + 1;
    else
        n_active_lig = (right_lig_ind + ligands.size() - left_lig_ind + 1) % ligands.size();

    // just for debug
    // check_rot_ind();


    ///////////////////
    //  Event rates  //
    ///////////////////

    size_t rates_i = 0;
    rates.resize(n_active_lig + bd_lig_ind.size());

    // binding events
    for (size_t i = left_lig_ind; i != after_right_lig_ind; i = (i + 1) % ligands.size(), rates_i++) {
        // ignore bonded ligands
        if (ligands[i].bond_state != 0) {
            rates[rates_i] = 0.0;
            continue;
        }
        rates[rates_i] = ligands[i].update_binding_rates(pos);
    }

    // rupture events
    for (const auto &i : bd_lig_ind) {
        rates[rates_i] = ligands[i].rupture_rate(pos);
        rates_i++;
    }

    double dt_bonds = compute_dt_bonds();
    if (dt_bonds != -INFTY)
        try_dt = std::min(try_dt, dt_bonds);

    // ODE step
    double step_done_with_dt = do_ode_step();
    diag.add_dt(step_done_with_dt);  // diagnostics

    if (step_done_with_dt >= dt_bonds) {
        vector<size_t> event_nrs = compute_event_nrs(step_done_with_dt);

        vector<size_t> ligs_to_erase;
        vector<size_t> ligs_to_insert;
        for (const auto &event_nr : event_nrs) {
            // Apply event

            // binding event
            if (event_nr < n_active_lig) {
                size_t lig_nr = (left_lig_ind + event_nr) % ligands.size();
                // bd_lig_ind.insert(lig_nr);
                ligs_to_insert.push_back(lig_nr);
                ligands[lig_nr].bond(pos[POS_ROT], pos[POS_DIST], generator);
            }
            // rupture event
            else {
                auto it = bd_lig_ind.begin();
                std::advance(it, event_nr - n_active_lig);
                ligands[*it].rupture();
                // bd_lig_ind.erase(*it);
                ligs_to_erase.push_back(*it);
            }
        }

        for (const auto &lig_to_erase : ligs_to_erase)
            bd_lig_ind.erase(lig_to_erase);

        diag.n_bonds_created += ligs_to_insert.size();
        for (const auto &lig_to_insert : ligs_to_insert) {
            bd_lig_ind.insert(lig_to_insert);
        }

        // Update positions after events

        double rot_reminder = std::fmod(pos[POS_ROT], 2 * PI);
        global_rot += pos[POS_ROT] - rot_reminder;
        pos[POS_ROT] = rot_reminder;

        for (const auto &i : bd_lig_ind)
            ligands[i].bd_rec_x -= pos[POS_DIST];
        global_dist += pos[POS_DIST];
        pos[POS_DIST] = 0.0;

        // ODE has changed
        reset_stepper();
    }
}

void AbstrSS::simulate(double max_time, size_t max_steps) {
    for (size_t i = 0; i < max_steps && time < max_time; ++i) {
        simulate_one_step();
    }
}

History AbstrSS::simulate_with_history(double max_time, size_t max_steps, double save_every) {
    History hist;
    hist.update(this);

    double last_time = time;
    for (size_t i = 1; i <= max_steps && time < max_time; ++i) {
        simulate_one_step();
        if (time - last_time >= save_every) {
            last_time = time;
            hist.update(this);
        }
    }
    hist.finish();
    return hist;
}

void AbstrSS::reseed(unsigned int seed) {
    generator = generator_t{seed};
}

void AbstrSS::update_rot_inc_range() {
    if (p->max_surf_dist <= exp(pos[POS_LOG_H])) {
        left_rot_inc = right_rot_inc = - pos[POS_ROT];
        return;
    }
    // rot + rot_inc, in [0, π]
    double beta = acos(1 - (p->max_surf_dist - exp(pos[POS_LOG_H])) / p->r_cell);
    right_rot_inc = beta - pos[POS_ROT];
    left_rot_inc = 2 * PI - beta - pos[POS_ROT];
    // projecting on [0, 2π]
    right_rot_inc = std::fmod(right_rot_inc + 2 * PI, 2 * PI);
    left_rot_inc = std::fmod(left_rot_inc + 2 * PI, 2 * PI);
}

// TODO: move it somewhere else
void AbstrSS::update_rot_inc_ind() {
    // |----[------]----|
    // 0    L      R   2 PI
    if (left_rot_inc <= right_rot_inc) {
        // R
        Ligand * right_lig = &ligands[right_lig_ind];
        if (right_lig->rot_inc <= right_rot_inc) {
            size_t new_right_lig_ind = helpers::cyclic_add(right_lig_ind, 1, ligands.size());
            Ligand * new_right_lig = &ligands[new_right_lig_ind];
            while (new_right_lig->rot_inc >= right_lig->rot_inc &&  // We didn't jump through 2 PI
                   new_right_lig->rot_inc <= right_rot_inc)         // We didn't jump through R
            {
                right_lig_ind = new_right_lig_ind;

                right_lig = &ligands[right_lig_ind];
                new_right_lig_ind = helpers::cyclic_add(right_lig_ind, 1, ligands.size());
                new_right_lig = &ligands[new_right_lig_ind];
            }
            if (right_lig->rot_inc < left_rot_inc) {
                // no ligands in [L, R]
                left_lig_ind = helpers::cyclic_add(right_lig_ind, 1, ligands.size());
                return;
            }
        } else {
            size_t new_right_lig_ind = helpers::cyclic_add(right_lig_ind, -1, ligands.size());
            Ligand * new_right_lig = &ligands[new_right_lig_ind];
            while (right_lig->rot_inc > right_rot_inc &&            // We are too far right
                   new_right_lig->rot_inc <= right_lig->rot_inc &&  // We didn't jump through 0
                   new_right_lig->rot_inc >= left_rot_inc)          // We didn't jump through L
            {
                right_lig_ind = new_right_lig_ind;

                right_lig = &ligands[right_lig_ind];
                new_right_lig_ind = helpers::cyclic_add(right_lig_ind, -1, ligands.size());
                new_right_lig = &ligands[new_right_lig_ind];
            }
            if (right_lig->rot_inc > right_rot_inc) {
                // no ligands in [L, R]
                left_lig_ind = helpers::cyclic_add(right_lig_ind, 1, ligands.size());
                return;
            }
        }

        // L
        Ligand * left_lig = &ligands[left_lig_ind];
        if (left_lig->rot_inc >= left_rot_inc) {
            size_t new_left_lig_ind = helpers::cyclic_add(left_lig_ind, -1, ligands.size());
            Ligand * new_left_lig = &ligands[new_left_lig_ind];
            while (new_left_lig->rot_inc <= left_lig->rot_inc &&  // We didn't jump through 0
                   new_left_lig->rot_inc >= left_rot_inc)         // We didn't jump through L
            {
                left_lig_ind = new_left_lig_ind;

                left_lig = &ligands[left_lig_ind];
                new_left_lig_ind = helpers::cyclic_add(left_lig_ind, -1, ligands.size());
                new_left_lig = &ligands[new_left_lig_ind];
            }
            if (left_lig->rot_inc > right_rot_inc) {
                // no ligands in [L, R]
                right_lig_ind = helpers::cyclic_add(left_lig_ind, -1, ligands.size());
                return;
            }
        } else {
            size_t new_left_lig_ind = helpers::cyclic_add(left_lig_ind, 1, ligands.size());
            Ligand * new_left_lig = &ligands[new_left_lig_ind];
            while (left_lig->rot_inc < left_rot_inc &&            // We are too far left
                   new_left_lig->rot_inc >= left_lig->rot_inc &&  // We didn't jump through 2 PI
                   new_left_lig->rot_inc <= right_rot_inc)        // We didn't jump through R
            {
                left_lig_ind = new_left_lig_ind;

                left_lig = &ligands[left_lig_ind];
                new_left_lig_ind = helpers::cyclic_add(left_lig_ind, 1, ligands.size());
                new_left_lig = &ligands[new_left_lig_ind];
            }
            if (left_lig->rot_inc < left_rot_inc) {
                // no ligands in [L, R]
                right_lig_ind = helpers::cyclic_add(left_lig_ind, -1, ligands.size());
                return;
            }
        }
    // |----]------[----|
    // 0    R      L   2 PI
    } else {
        // R
        Ligand * right_lig = &ligands[right_lig_ind];
        // Try to find largest index for which rot_inc <= R

        // Case 1.
        // |--x-]------[----|
        // 0    R      L   2 PI
        if (right_lig->rot_inc <= right_rot_inc) {
            for (size_t i = right_lig_ind + 1; i < ligands.size(); i++) {
                if (ligands[i].rot_inc <= right_rot_inc)
                    right_lig_ind = i;
                else
                    break;
            }
        // Case 2.
        // |----]--x---[----|
        // 0    R      L   2 PI
        } else {
            // size_t is unsigned, hence workaround with i_plus_1
            size_t i;
            for (size_t i_plus_1 = right_lig_ind; i_plus_1 != 0; i_plus_1--) {
                i = i_plus_1 - 1;
                right_lig_ind = i;
                // We are satisfied with first ligand with rot_inc <= R
                if (ligands[i].rot_inc <= right_rot_inc)
                    break;
            }
            // There are no ligands with rot_inc <= R
            if (ligands[right_lig_ind].rot_inc > right_rot_inc) {
                // We should look for rightmost ligand with rot_inc >= L
                if (ligands.back().rot_inc >= left_rot_inc) {
                    right_lig_ind = ligands.size() - 1;
                } else {
                    // no ligands in [0, R] or [L, 2 PI]
                    right_lig_ind = 0;
                    left_lig_ind = 1;
                    return;
                }
            }
        }

        // L
        Ligand * left_lig = &ligands[left_lig_ind];
        // Try to find the smallest index for which rot_inc >= L

        // Case 1.
        // |----]------[-x--|
        // 0    R      L   2 PI
        if (left_lig->rot_inc >= left_rot_inc) {
            size_t i;
            for (size_t i_plus_1 = left_lig_ind; i_plus_1 != 0; i_plus_1--) {
                i = i_plus_1 - 1;
                if (ligands[i].rot_inc >= left_rot_inc)
                    left_lig_ind = i;
                else
                    break;
            }
        // Case 2.
        // |----]----x-[----|
        // 0    R      L   2 PI
        } else {
            for (size_t i = left_lig_ind + 1; i < ligands.size(); i++) {
                left_lig_ind = i;
                if (ligands[i].rot_inc >= left_rot_inc)
                    break;
            }
            // There are no ligands with rot_ind >= L
            if (ligands[left_lig_ind].rot_inc < left_rot_inc) {
                // We should look for leftmost ligand with rot_inc <= R
                if (ligands.front().rot_inc <= right_rot_inc) {
                    left_lig_ind = 0;
                } else {
                    // no ligands in [0, R] or [L, 2 PI]
                    right_lig_ind = 0;
                    left_lig_ind = 1;
                    return;
                }
            }
        }
    }
}

void AbstrSS::rhs(const array<double, 3> &x, array<double, 3> &dxdt, double /*t*/) {
    forces_t f = forces::non_bond_forces(shear_rate, exp(x[POS_LOG_H]), p);
    // add forces of each bond
    for (auto & bd_i : bd_lig_ind)
        f += ligands[bd_i].bond_forces(pos);

    dxdt = velocities::compute_velocities(x[POS_LOG_H], f, p);
}

void AbstrSS::check_rot_ind() {
    if (helpers::cyclic_add(left_lig_ind, -1, ligands.size()) != right_lig_ind) {
        Ligand & left_lig = ligands[left_lig_ind];
        Ligand & right_lig = ligands[right_lig_ind];
        if (left_rot_inc <= right_rot_inc) {
            if (left_lig.rot_inc < left_rot_inc || left_lig.rot_inc > right_rot_inc ||
                right_lig.rot_inc < left_rot_inc || right_lig.rot_inc > right_rot_inc) {
                std::cout << "wrong rot_inc indices";
            }
        } else {
            if ((left_lig.rot_inc > right_rot_inc && left_lig.rot_inc < left_rot_inc) ||
                (right_lig.rot_inc > right_rot_inc && right_lig.rot_inc < left_rot_inc)) {
                std::cout << "wrong rot_inc indices";
            }
        }
    }
}

void AbstrSS::Diagnostic::add_dt(double dt) {
    auto i = size_t(-log10(dt));
    if (dt_freq.size() < i + 1) {
        if (i + 1 > dt_freq.max_size())
            std::cout << i + 1 << " > max_size, dt = " << dt << std::endl;
        dt_freq.resize(i + 1);
    }
    dt_freq[i]++;
}
