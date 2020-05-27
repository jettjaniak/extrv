import time

import numpy as np
from scipy import integrate
import multiprocessing
import pickle
import os
from datetime import datetime
from collections import defaultdict, namedtuple


from extrv_engine import EulerGillSS, EulerProbSS, RKGillSS, RKProbSS, AdapGillSS, AdapProbSS, SlipBondType, Parameters

CONST_DT = 1e-6
FALLING_TIME = 1
ROLLING_TIME = 10
N_TRIALS = 50
# SOLVERS = ['euler_gill', 'euler_prob', 'rk_gill', 'rk_prob', 'adap_gill', 'adap_prob']
SOLVERS = [
    'euler_gill, dt <= 1e-6', 'euler_prob, dt = 1e-6',
    'rk_gill, dt <= 1e-5', 'rk_prob, dt = 1e-6',
    'adap_gill, dt <= 0.1 (with bonds 1e-4)', 'adap_prob, dt <= 1e-6'
]
EULER_GILL_DT = 1e-6
EULER_PROB_DT = 1e-6
RK_GILL_DT = 1e-5
RK_PROB_DT = 1e-6
ADAP_GILL_MAX_DT = 0.1
ADAP_GILL_MAX_DT_WITH_BONDS = 1e-4
ADAP_PROB_MAX_DT = 1e-6
ADAP_PROB_MAX_DT_WITH_BONDS = 1e-6

# SHEAR_RATES = [0.1, 0.9]
SHEAR_RATES = [0.1, 0.3, 0.5, 0.7, 0.9]

ADAP_ABS_ERR = 1e-10
ADAP_REL_ERR = 1e-6
INITIAL_HEIGHT = 0.03
SAVE_EVERY = 1e-4

MAX_STEPS_FALLING = int(2 * FALLING_TIME * 1e6)
MAX_STEPS_ROLLING = int(2 * ROLLING_TIME * 1e6)

SimulationStats = namedtuple("SimulationStats", ['mean_h', 'rot', 'dist', 'n_bonds', 'mean_bond_ls', 'comp_time'])


def setup_parameters():
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        rep_0=5075,
        rep_scale=1145
    )

    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=27,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=750,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6
    )
    psgl = Parameters.LigandType()
    psgl.add_bond_type(psgl_plus_esel_bond)
    p.add_ligands(lig_type=psgl, n_of_lig=20000)
    # Don't let them be garbage collected!
    return p, psgl, psgl_plus_esel_bond


def many_tests(solver_name, pool):
    for shear_rate in SHEAR_RATES:
        one_test(solver_name, shear_rate, pool)


def one_test(solver_name, shear_rate, pool):
    seeds = np.random.randint(uint_info.max, size=N_TRIALS, dtype='uint32')
    # seeds = np.arange(n_trials) + 100
    for seed in seeds:
        pool.apply_async(iteration_wrapper, (solver_name, shear_rate, seed),
                         callback=iteration_callback, error_callback=iteration_error_callback)


def iteration_wrapper(*args, **kwargs):
    try:
        return iteration(*args, **kwargs)
    except Exception as error:
        print("ERROR RAISED!")
        raise error


def iteration(solver_name, shear_rate=0, seed=0, save_every=SAVE_EVERY):
    p, psgl, psgl_plus_esel_bond = setup_parameters()

    if solver_name.startswith('euler_gill'):
        ss = EulerGillSS(INITIAL_HEIGHT, p, seed, dt=EULER_GILL_DT)
    elif solver_name.startswith('euler_prob'):
        ss = EulerProbSS(INITIAL_HEIGHT, p, seed, dt=EULER_PROB_DT)
    elif solver_name.startswith('rk_gill'):
        ss = RKGillSS(INITIAL_HEIGHT, p, seed, dt=RK_GILL_DT)
    elif solver_name.startswith('rk_prob'):
        ss = RKProbSS(INITIAL_HEIGHT, p, seed, dt=RK_PROB_DT)
    elif solver_name.startswith('adap_gill'):
        ss = AdapGillSS(INITIAL_HEIGHT, p, seed,
                        max_dt=ADAP_GILL_MAX_DT, max_dt_with_bonds=ADAP_GILL_MAX_DT_WITH_BONDS,
                        abs_err=ADAP_ABS_ERR, rel_err=ADAP_REL_ERR)
    elif solver_name.startswith('adap_prob'):
        ss = AdapProbSS(INITIAL_HEIGHT, p, seed,
                        max_dt=ADAP_PROB_MAX_DT, max_dt_with_bonds=ADAP_PROB_MAX_DT_WITH_BONDS,
                        abs_err=ADAP_ABS_ERR, rel_err=ADAP_REL_ERR)
    else:
        raise ValueError

    # falling
    ss.simulate(FALLING_TIME, MAX_STEPS_FALLING)
    n_bonds_before_falling = ss.diag.n_bonds_created - len(ss.bd_lig_ind)

    # rolling
    ss.shear_rate = shear_rate

    comp_start_time = time.time()
    sim_hist = ss.simulate_with_history(FALLING_TIME + ROLLING_TIME, MAX_STEPS_ROLLING, save_every=save_every)
    comp_end_time = time.time()

    print("simulation done for solver", solver_name, "seed", seed, "shear", shear_rate)

    n_bonds_since_falling = ss.diag.n_bonds_created - n_bonds_before_falling
    mean_bonds_ls = 0.0
    if n_bonds_since_falling != 0:
        n_short_bonds = n_bonds_since_falling - len(sim_hist.bond_trajectories)
        mean_short_bond_lifespan = save_every / 2

        all_bonds_lifespan = n_short_bonds * mean_short_bond_lifespan
        for bond_traj in sim_hist.bond_trajectories:
            start_i = bond_traj.start_i
            end_i = min(bond_traj.start_i + len(bond_traj.positions), len(sim_hist.time) - 1)
            bond_lifespan = sim_hist.time[end_i] - sim_hist.time[start_i]
            all_bonds_lifespan += bond_lifespan
        mean_bonds_ls = all_bonds_lifespan / n_bonds_since_falling


    sim_stats = SimulationStats(
        mean_h=integrate.simps(sim_hist.h, sim_hist.time) / (sim_hist.time[-1] - sim_hist.time[0]),
        rot=sim_hist.rot[-1] - sim_hist.rot[0],
        dist=sim_hist.dist[-1] - sim_hist.dist[0],
        n_bonds=n_bonds_since_falling,
        mean_bond_ls=mean_bonds_ls,
        comp_time=comp_end_time - comp_start_time
    )

    return solver_name, shear_rate, sim_stats
    # return solver_name, shear_rate, sim_stats, ss, sim_hist


def iteration_callback(result):
    solver_name, shear_rate, sim_stats = result
    test_results[solver_name][shear_rate].append(sim_stats)


def iteration_error_callback(error):
    print("ERROR CALLBACK CALLED!")
    raise error


if __name__ == '__main__':
    uint_info = np.iinfo(np.uint32)
    pool = multiprocessing.Pool()

    test_results = dict()
    for solver_name in SOLVERS:
        test_results[solver_name] = defaultdict(list)
        many_tests(solver_name, pool)

    pool.close()
    pool.join()

    directory = '../results/diff_shear/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    date_str = datetime.now().strftime("%d-%m-%Y_%Hh%Mm%Ss")
    with open(f'{directory}/{date_str}.pickle', 'wb') as file:
        pickle.dump(test_results, file)