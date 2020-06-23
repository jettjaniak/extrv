#include "SimulationState.h"

#include <algorithm>
#include <iostream>

#include "forces.h"
#include "velocities.h"
#include "helpers.h"

#include <boost/numeric/odeint/stepper/generation.hpp>


SimulationState::SimulationState(
        double h_0, Parameters* p, unsigned int seed,
        double max_dt, double ode_abs_err, double ode_rel_err,
        double rate_integral_tol_) :

        p(p),
        max_dt(max_dt),
        ode_abs_err(ode_abs_err),
        ode_rel_err(ode_rel_err),
        rate_integral_tol(rate_integral_tol_)
{
    reset_stepper();
    try_dt = max_dt;
    reseed(seed);
    update_randomness();
    if (rate_integral_tol == -1.0)
        rate_integral_tol = ode_abs_err;

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

    std::fill(ode_x.begin(), ode_x.end(), 0.0);
    ode_x[POS_LOG_H] = log(h_0);
}

double SimulationState::h() const {
    return exp(ode_x[POS_LOG_H]);
}

double SimulationState::rot() const {
    return ode_x[POS_ROT];
}

double SimulationState::dist() const {
    return ode_x[POS_DIST];
}

double SimulationState::global_rot() const {
    return cumulated_rot + rot();
}

double SimulationState::global_dist() const {
    return cumulated_dist + dist();
}

void SimulationState::simulate_one_step() {

    // ODE step
    double rate_integral = do_ode_step();

    if (rate_integral > rate_integral_value - rate_integral_tol) {
        vector<double> weights(ligands.size(), 0.0);

        update_rot_inc_range();
        update_rot_inc_ind();

        after_right_lig_ind = (right_lig_ind + 1) % ligands.size();

        for (size_t i = left_lig_ind; i != after_right_lig_ind; i = (i + 1) % ligands.size()) {
            // ignore bonded ligands
            if (ligands[i].bond_state != 0) {
                continue;
            }
            // TODO: now it only works with one type of bond per ligand
            weights[i] = ligands[i].update_binding_rates(h(), rot());
        }
        for (const auto &i : bd_lig_ind) {
            weights[i] = ligands[i].rupture_rate(h(), rot(), dist());
        }

        std::discrete_distribution<size_t> which_event_distribution(weights.begin(), weights.end());
        size_t lig_nr = which_event_distribution(generator);

        ode_x[POS_SUM_RATE_INT] = 0.0;
        // binding event
        if (ligands[lig_nr].bond_state == 0) {
            bd_lig_ind.insert(lig_nr);
            ligands[lig_nr].bond(rot(), dist(), generator);
            diag.n_bonds_created++;
        }
        // rupture event
        else {
            ligands[lig_nr].rupture();
            bd_lig_ind.erase(lig_nr);
        }

        // Update positions after events

        double rot_reminder = std::fmod(rot(), 2 * PI);
        cumulated_rot += rot() - rot_reminder;
        ode_x[POS_ROT] = rot_reminder;

        for (const auto &i : bd_lig_ind)
            ligands[i].bd_rec_x -= dist();
        cumulated_dist += dist();
        ode_x[POS_DIST] = 0.0;

        update_randomness();
        // ODE has changed
        reset_stepper();
    }
}

void SimulationState::simulate(double max_time, size_t max_steps) {
    for (size_t i = 0; i < max_steps && time < max_time; ++i) {
        simulate_one_step();
    }
}

History SimulationState::simulate_with_history(double max_time, size_t max_steps, double save_every) {
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

void SimulationState::reseed(unsigned int seed) {
    generator = generator_t{seed};
}

void SimulationState::update_rot_inc_range() {
    if (p->max_surf_dist <= h()) {
        left_rot_inc = right_rot_inc = - rot();
        return;
    }
    // rot + rot_inc, in [0, π]
    double beta = acos(1 - (p->max_surf_dist - h()) / p->r_cell);
    right_rot_inc = beta - rot();
    left_rot_inc = 2 * PI - beta - rot();
    // projecting on [0, 2π]
    right_rot_inc = std::fmod(right_rot_inc + 2 * PI, 2 * PI);
    left_rot_inc = std::fmod(left_rot_inc + 2 * PI, 2 * PI);
}

// TODO: move it somewhere else
void SimulationState::update_rot_inc_ind() {
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

void SimulationState::rhs(const array<double, 4> &x, array<double, 4> &dxdt, double /*t*/) {
    std::fill(dxdt.begin(), dxdt.end(), 0.0);

    update_rot_inc_range();
    update_rot_inc_ind();

    after_right_lig_ind = (right_lig_ind + 1) % ligands.size();

    if (x[POS_LOG_H] > MAX_LOG_H)
        return;
    double h = exp(x[POS_LOG_H]);
    // Don't bother computing anything, it will be caught later.
    if (std::isnan(h) or std::isinf(h))
        return;

    dxdt[POS_SUM_RATE_INT] = 0.0;
    for (size_t i = left_lig_ind; i != after_right_lig_ind; i = (i + 1) % ligands.size()) {
        // ignore bonded ligands
        if (ligands[i].bond_state != 0) {
            continue;
        }
        // TODO: now it only works with one type of bond per ligand
        dxdt[POS_SUM_RATE_INT] += ligands[i].update_binding_rates(h, x[POS_ROT]);
    }

    // rupture events
    for (const auto &i : bd_lig_ind) {
        dxdt[POS_SUM_RATE_INT] += ligands[i].rupture_rate(h, x[POS_ROT], x[POS_DIST]);
    }

    forces_t f = forces::non_bond_forces(shear_rate, h, p);
    // add forces of each bond
    for (auto & bd_i : bd_lig_ind)
        f += ligands[bd_i].bond_forces(h, x[POS_ROT], x[POS_DIST]);

    array<double, 3> vel = velocities::compute_velocities(h, f, p);
    dxdt[POS_LOG_H] = vel[POS_LOG_H];
    dxdt[POS_ROT] = vel[POS_ROT];
    dxdt[POS_DIST] = vel[POS_DIST];

}

void SimulationState::check_rot_ind() {
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

void SimulationState::Diagnostic::add_dt(double dt) {
    auto i = size_t(-log10(dt));
    if (dt_freq.size() < i + 1) {
        if (i + 1 > dt_freq.max_size())
            std::cout << i + 1 << " > max_size, dt = " << dt << std::endl;
        dt_freq.resize(i + 1);
    }
    dt_freq[i]++;
}

double SimulationState::do_ode_step() {
    try_dt = std::min(try_dt, max_dt);
    // TODO: define as class parameter or define operator()
    namespace pl = std::placeholders;
    static auto rhs_system = std::bind(&SimulationState::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);

    controlled_step_result result;
    double dt_inout, time_inout;
    static array<double, 4> ode_x_inout;
    double old_rate_integral, rate_integral;

    while (true) {
        dt_inout = try_dt;
        ode_x_inout = ode_x;
        time_inout = time;
        result = stepper.try_step(rhs_system, ode_x_inout, time_inout, dt_inout);
        if (result == fail) {
            // solver made dt_inout smaller
            try_dt = dt_inout;
            continue;
        }
        if (ode_x[POS_LOG_H] > MAX_LOG_H || helpers::values_not_ok(ode_x_inout.begin(), ode_x_inout.end())) {
            reset_stepper();
            try_dt /= 2;
            diag.n_pos_not_ok++;
            continue;
        }
        rate_integral = ode_x_inout[POS_SUM_RATE_INT];
        if (rate_integral > rate_integral_value + rate_integral_tol) {
            diag.n_rate_int_too_big++;
            reset_stepper();
            old_rate_integral = ode_x[POS_SUM_RATE_INT];
            double diff_curr = rate_integral - old_rate_integral;
            double diff_opt = rate_integral_value - old_rate_integral;
            try_dt *= diff_opt / diff_curr;
            continue;
        }
        break;
    }
    diag.add_dt(try_dt);

//    static long n_eq = 0;
//    for (int i = 0; i < 3; i++) {
//        if (ode_x[i] == ode_x_inout[i]) {
//            std::cout << ++n_eq << ", dt = " << try_dt << std::endl;
//            std::cout << "pos before: [";
//            for (const auto &item : ode_x) {
//                std::cout << item << ", ";
//            }
//            std::cout << "]" << std::endl;
//
//            std::cout << "pos after:  [";
//            for (const auto &item : ode_x_inout) {
//                std::cout << item << ", ";
//            }
//            std::cout << "]" << std::endl;
//
////            std::cout << "dx / dt:    [";
////            for (const auto &item : dxdt) {
////                std::cout << item << ", ";
////            }
////            std::cout << "]" << std::endl << std::endl;
//            std::cout << std::endl;
//            break;
//        }
//    }

    try_dt = dt_inout;  // possibly larger after successful step
    ode_x = ode_x_inout;
    time = time_inout;

    return rate_integral;
}

void SimulationState::reset_stepper() {
    stepper = make_controlled<error_stepper_type>(ode_abs_err, ode_rel_err);
}

void SimulationState::update_randomness() {
    const double u = helpers::draw_from_uniform_dist(generator);
    rate_integral_value = -log(u);
}