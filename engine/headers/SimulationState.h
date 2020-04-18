#pragma once

#include "types.h"
#include "Ligand.h"
#include "Parameters.h"
#include "helpers.h"
#include "History.h"

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>

using namespace boost::numeric::odeint;

typedef runge_kutta_dopri5<vector<double>> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> adaptive_stepper_type;


struct SimulationState {
    struct Diagnostic {
        /// stores frequencies of -int(log10(try_dt))
        vector<size_t> dt_freq;
        size_t n_bonds_created = 0;

        void add_dt(double dt);
    };
    Diagnostic diag;

    /// height, rotation and distance
    vector<double> ode_x;
    size_t log_h_ode_i;
    size_t rot_ode_i;
    size_t dist_ode_i;

    double try_dt;
    double u;

    double time = 0.0;
    /// sphere's rotation in radians
    double cumulated_rot = 0.0;
    /// distance traveled in direction of flow (x)
    double cumulated_dist = 0.0;

    /// fluid flow shear rate in 1/s
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

    /// rates of binding and rupture reactions
    vector<double> rates;

    /// indices of bonded ligands
    set<size_t> bd_lig_ind;

    /// random number generator
    generator_t generator;

    /// model parameters
    Parameters* p;

    adaptive_stepper_type stepper;

    /// maximal time step
    double max_dt;
    /// maximal time step
    double max_dt_with_bonds;
    /// limit on absolute error of ODE solver
    double abs_err;
    /// limit on relative error of ODE solver
    double rel_err;

    /**
     * Initialize simulation.
     *
     * @param h_0 starting height above surface
     * @param p model parameters
     * @param seed number that random number generator will be seeded with
     */
    SimulationState(double h_0, Parameters* p, unsigned int seed,
                    double max_dt = 0.1, double max_dt_with_bonds = 1e-4,
                    double abs_err = 1e-10, double rel_err = 1e-6);


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

    void update_rot_inc_range();
    void update_rot_inc_ind();
    void check_rot_ind();

    /// ODE's RHS
    void rhs(const vector<double> & x, vector<double> & dxdt, double /*t*/);

    void reset_stepper();
    double do_ode_step();
    double compute_dt_bonds();
    vector<size_t> compute_event_nrs(double dt);
};