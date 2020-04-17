#include "SimulationState.h"

#include <algorithm>
#include <utility>
#include <iostream>

#include "forces.h"
#include "velocities.h"
#include "helpers.h"

#include <boost/numeric/odeint/stepper/generation.hpp>


SimulationState::SimulationState(
        double h_0, Parameters* p, unsigned int seed,
        double max_dt, double max_dt_with_bonds,
        double abs_err, double rel_err) :

        p(p),
        max_dt(max_dt),
        max_dt_with_bonds(max_dt_with_bonds),
        abs_err(abs_err),
        rel_err(rel_err)
{
    reset_stepper();
    try_dt = max_dt;
    reseed(seed);
    u = helpers::draw_from_uniform_dist(generator);

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

    ode_x.resize(ligands.size() + 1 + 3, 0.0);
    h_ode_i = ligands.size() + 1 + POS_H;
    rot_ode_i = ligands.size() + 1 + POS_ROT;
    dist_ode_i = ligands.size() + 1 + POS_DIST;
//    ode_x.resize(3);
//    h_ode_i = POS_H;
//    rot_ode_i = POS_ROT;
//    dist_ode_i = POS_DIST;
    ode_x[h_ode_i] = h_0;
}

double SimulationState::h() const {
    return ode_x[h_ode_i];
}

double SimulationState::rot() const {
    return ode_x[rot_ode_i];
}

double SimulationState::dist() const {
    return ode_x[dist_ode_i];
}

double SimulationState::global_rot() const {
    return cumulated_rot + rot();
}

double SimulationState::global_dist() const {
    return cumulated_dist + dist();
}

void SimulationState::simulate_one_step() {

    // ODE step
    double step_done_with_dt = do_ode_step();
    diag.add_dt(step_done_with_dt);  // diagnostics


    if (std::abs(1 - exp(- ode_x[ligands.size()]) - u) < 1e-3) {
        std::discrete_distribution<size_t> which_event_distribution(ode_x.begin(), ode_x.end() - 4);
        size_t lig_nr = which_event_distribution(generator);
        std::fill(ode_x.begin(), ode_x.end() - 3, 0.0);
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

        double rot_reminder = std::fmod(ode_x[rot_ode_i], 2 * PI);
        cumulated_rot += ode_x[rot_ode_i] - rot_reminder;
        ode_x[rot_ode_i] = rot_reminder;

        for (const auto &i : bd_lig_ind)
            ligands[i].bd_rec_x -= ode_x[dist_ode_i];
        cumulated_dist += ode_x[dist_ode_i];
        ode_x[dist_ode_i] = 0.0;


        // ODE has changed
        reset_stepper();

        u = helpers::draw_from_uniform_dist(generator);
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

void SimulationState::rhs(const vector<double> &x, vector<double> &dxdt, double /*t*/) {
    std::fill(dxdt.begin(), dxdt.end(), 0.0);

    update_rot_inc_range();
    update_rot_inc_ind();

    after_right_lig_ind = (right_lig_ind + 1) % ligands.size();

    if (left_lig_ind <= right_lig_ind)
        n_active_lig = right_lig_ind - left_lig_ind + 1;
    else
        n_active_lig = (right_lig_ind + ligands.size() - left_lig_ind + 1) % ligands.size();

    double any_event_rate = 0.0;
    for (size_t i = left_lig_ind; i != after_right_lig_ind; i = (i + 1) % ligands.size()) {
        // ignore bonded ligands
        if (ligands[i].bond_state != 0) {
            continue;
        }
        dxdt[i] = ligands[i].update_binding_rates(x[h_ode_i], x[rot_ode_i]);
        any_event_rate += dxdt[i];
    }

    // rupture events
    for (const auto &i : bd_lig_ind) {
        dxdt[i] = ligands[i].rupture_rate(x[h_ode_i], x[rot_ode_i], x[dist_ode_i]);
        any_event_rate += dxdt[i];
    }
    dxdt[ligands.size()] = any_event_rate;

    forces_t f = forces::non_bond_forces(shear_rate, x[h_ode_i], p);
    // add forces of each bond
    for (auto & bd_i : bd_lig_ind)
        f += ligands[bd_i].bond_forces(x[h_ode_i], x[rot_ode_i], x[dist_ode_i]);

    array<double, 3> vel = velocities::compute_velocities(x[h_ode_i], f, p);
    dxdt[h_ode_i] = vel[POS_H];
    dxdt[rot_ode_i] = vel[POS_ROT];
    dxdt[dist_ode_i] = vel[POS_DIST];
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
    auto rhs_system = std::bind(&SimulationState::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);

    controlled_step_result result;
    double dt_inout = 0.0, time_inout = 0.0;
    vector<double> pos_inout(3);

    bool step_done = false;
    while (!step_done) {
        dt_inout = try_dt;
        pos_inout = ode_x;
        time_inout = time;
        result = stepper.try_step(rhs_system, pos_inout, time_inout, dt_inout);
        if (result == fail) {
            // solver made dt_inout smaller
            try_dt = dt_inout;
            continue;
        }
        double cdf = 1 - exp(- pos_inout[ligands.size()]);
        if (cdf > u + 1e-3) {
            reset_stepper();
            try_dt /= 2;
            continue;
        }
        step_done = true;
    }

    double step_done_with_dt = try_dt;
    try_dt = dt_inout;  // possibly larger after successful step
    ode_x = pos_inout;
    time = time_inout;

    return step_done_with_dt;
}

void SimulationState::reset_stepper() {
    stepper = make_controlled<error_stepper_type>(abs_err, rel_err);
}

double SimulationState::compute_dt_bonds() {
    //////////////////////////////////////////////
    //  Gillespie algorithm - first event time  //
    //////////////////////////////////////////////

    double any_event_rate = 0.0;
    for (const auto &rate : rates)
        any_event_rate += rate;

    // TODO: just for debug, remove
    static constexpr double max_rate = 1e80;
    if (any_event_rate > max_rate)
        std::cout << "EVENT RATE > " << max_rate << std::endl;

    if (any_event_rate > 0.0) {
        std::exponential_distribution<double> dt_bonds_distribution(any_event_rate);
        return dt_bonds_distribution(generator);
    }
    return INFTY;
}

vector<size_t> SimulationState::compute_event_nrs(double) {
    std::discrete_distribution<size_t> which_event_distribution(rates.begin(), rates.end());
    return {which_event_distribution(generator)};
}