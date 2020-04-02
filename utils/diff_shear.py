import numpy as np
import multiprocessing
import pickle
import os
from datetime import datetime
from collections import defaultdict


from extrv_engine import SimulationState, SlipBondType, Parameters

DT = 1e-6
FALLING_TIME = 1
ROLLING_TIME = 10

N_TRIALS = 100
SHEAR_RATES = [0.1, 0.3, 0.5, 0.7, 0.9]

# DT = 1e-5
# FALLING_TIME = 0.1
# ROLLING_TIME = 1
#
# N_TRIALS = 10
# SHEAR_RATES = [0.3, 0.7]

INITIAL_HEIGHT = 0.03


def setup_parameters():
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        rep_0=5074.616349947093,
        rep_scale=1146.409200818891
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


def many_tests(pool):
    for shear_rate in SHEAR_RATES:
        one_test(shear_rate, pool)


def one_test(shear_rate, pool):
    seeds = np.random.randint(uint_info.max, size=N_TRIALS, dtype='uint32')
    # seeds = np.arange(n_trials) + 100
    for seed in seeds:
        pool.apply_async(iteration_wrapper, (shear_rate, seed),
                         callback=iteration_callback, error_callback=iteration_error_callback)


def iteration_wrapper(*args, **kwargs):
    try:
        return iteration(*args, **kwargs)
    except Exception as error:
        print("ERROR RAISED!")
        raise error


def iteration(shear_rate=0, seed=0, save_every=1000):
    p, psgl, psgl_plus_esel_bond = setup_parameters()
    n_steps = int(ROLLING_TIME / DT)
    n_steps_falling = int(FALLING_TIME / DT)

    ss = SimulationState(INITIAL_HEIGHT, p, seed)

    ss.simulate(n_steps_falling, DT, shear=0.0)  # falling
    sim_hist = ss.simulate_with_history(n_steps, DT, shear=shear_rate, save_every=save_every)
    print("simulation done for seed", seed, "shear", shear_rate)

    # TODO: pickle History object
    vel = np.diff(sim_hist.dist) / (DT * save_every)
    return shear_rate, np.mean(vel)


def iteration_callback(result):
    shear_rate, mean_v = result
    test_results[shear_rate].append(mean_v)


def iteration_error_callback(error):
    print("ERROR CALLBACK CALLED!")
    raise error


if __name__ == '__main__':
    uint_info = np.iinfo(np.uint32)
    test_results = defaultdict(list)

    pool = multiprocessing.Pool()

    many_tests(pool)

    pool.close()
    pool.join()

    directory = '../results/diff_shear/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    date_str = datetime.now().strftime("%d-%m-%Y_%Hh%Mm%Ss")
    with open(f'{directory}/{date_str}.pickle', 'wb') as file:
        pickle.dump(test_results, file)
