#pragma once

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"
#include "helpers.h"
#include "History.h"

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>

using namespace boost::numeric::odeint;

typedef runge_kutta_dopri5<array<double, 4>> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> adaptive_stepper_type;


struct SimulationState {
    struct Diagnostic {
        /// stores frequencies of -int(log10(try_dt))
        vector<size_t> dt_freq;
        size_t n_bonds_created = 0;
        size_t n_pos_not_ok = 0;
        size_t n_rate_int_too_big = 0;

        void add_dt(double dt);
    };
    Diagnostic diag;

    /// reaction rate and cell position
    array<double, 4> ode_x;

    /// time step that we try to use
    double try_dt;
    /// maximal time step
    double max_dt;
    /// limit on absolute error of ODE solver
    double ode_abs_err;
    /// limit on relative error of ODE solver
    double ode_rel_err;

    /// absolute tolerance value of rate integral
    double rate_integral_tol;
    /// if integral of first event rate is close to it, the event fires
    double rate_integral_value;

    double time = 0.0;
    /// sphere's rotation in radians
    double cumulative_rot = 0.0;
    /// distance traveled in direction of flow (x)
    double cumulative_dist = 0.0;

    /// fluid flow shear rate in 1/s
    double shear_rate = 0.0;

    /// ligands on sphere
    vector<Ligand> ligands;
    /// indices of bonded ligands
    set<size_t> bd_lig_ind;
    /// model parameters
    Parameters* p;

    /// random number generator
    generator_t generator;
    /// ODE stepper
    adaptive_stepper_type stepper;

    /**
     * Initialize simulation.
     *
     * @param h_0 starting height above surface
     * @param p model parameters
     * @param seed number that random number generator will be seeded with
     */
    SimulationState(double h_0, Parameters* p, unsigned int seed,
                    double max_dt = 0.1, double ode_abs_err = 1e-10, double ode_rel_err = 0.0,
                    double rate_integral_tol = -1.0);


    double h() const;
    double rot() const;
    double dist() const;
    double global_rot() const;
    double global_dist() const;

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


    /// range of ligands' rot_incs for which ligands can bond at given step
    double left_rot_inc = 0.0;
    double right_rot_inc = 0.0;
    /// update left_ and right_rot_inc
    void update_rot_inc_range();

    /// indices of leftmost and rightmost ligands that are close enough to surface to bond
    size_t left_lig_ind = 1;
    size_t right_lig_ind = 0;
    size_t after_right_lig_ind = 1;
    /// update left_ and right_lig_ind
    void update_rot_inc_ind();

    // sanity check for the above
    void check_rot_ind();

    /// ODE's RHS
    void rhs(const array<double, 4> &x, array<double, 4> &dxdt, double /*t*/);

    // make new ODE stepper in case the previous got in trouble
    void reset_stepper();

    //
    void do_ode_step();

    /// sample rate integral value for which random event will trigger
    void update_randomness();
};