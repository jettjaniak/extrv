#pragma once

//#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"
#include "helpers.h"
#include "History.h"

using namespace boost::numeric::odeint;

struct SimulationState {
    struct Diagnostic {
        /// stores frequencies of -int(log10(try_dt))
        vector<size_t> dt_freq;
        /// how many times position was NaN or outside allowed range
        size_t n_pos_not_ok = 0;

        void add_dt(double dt);
    };
    Diagnostic diag;

    /// log height, rotation and distance
    array<double, 3> pos {};
    double try_dt = MAX_DT;

    double time = 0.0;
    /// sphere's rotation in radians
    double global_rot = 0.0;
    /// distance traveled in direction of flow (x)
    double global_dist = 0.0;

    double shear_rate = 0.0;

    /// ligands on sphere
    vector<Ligand> ligands;

    /// range of ligands' rot_incs for which ligands can bond at given step
    double left_rot_inc = 0.0;
    double right_rot_inc = 0.0;

    /// indices of leftmost and rightmost ligands that are close enough to surface to bond
    size_t left_lig_ind = 1;
    size_t right_lig_ind = 0;
    size_t after_right_lig_ind = 1;
    size_t n_active_lig = 0;

    vector<double> rates;

    /// indices of bonded ligands
    set<size_t> bd_lig_ind;

    /// random number generator
    generator_t generator;

    /// model parameters
    Parameters* p;

    typedef runge_kutta_dopri5<array<double, 3>> error_stepper_type;
    typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    controlled_stepper_type stepper;

    /**
     * Initialize simulation.
     *
     * @param h_0 starting height above surface
     * @param p model parameters
     * @param seed number that random number generator will be seeded with
     */
    SimulationState(double h_0, Parameters* p, unsigned int seed);

    /**
     * Do one step of simulation.
     *
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     */
    void simulate_one_step();

    /**
     * Do n_steps of simulation.
     *
     * @param max_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     */
    void simulate(double max_time, size_t max_steps);

    /**
     * Do n_steps of simulation and return history updated every save_every steps.
     *
     * @param max_steps number of simulation steps
     * @param dt time step in seconds
     * @param shear_rate fluid shear rate in 1/s
     * @param save_every number of steps between history updates
     * @return simulation history
     */
    History simulate_with_history(double max_time, size_t max_steps, double save_every);

    /// reseed random number generator
    void reseed(unsigned int seed);

    void update_rot_inc_range();

    void update_rot_inc_ind();

    void rhs(const array<double, 3> & x, array<double, 3> & dxdt, double /*t*/);
};