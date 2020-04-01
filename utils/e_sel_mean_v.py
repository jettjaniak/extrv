import numpy as np
import multiprocessing
import pickle
import os
from datetime import datetime
from collections import defaultdict

from extrv_engine import SlipBondType, Parameters, SimulationState

REC_DENS_0 = 750
BINDING_RATE_0 = 0.06


def setup_parameters(fold_change):
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        rep_0=5074.616349947093,
        rep_scale=1146.409200818891,
        abs_err=1e-10,
        rel_err=1e-6
    )
    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=27,
        spring_const=100,
        binding_rate_0=BINDING_RATE_0,
        rec_dens=REC_DENS_0 * fold_change,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6
    )
    psgl = Parameters.LigandType()
    psgl.add_bond_type(psgl_plus_esel_bond)
    p.add_ligands(lig_type=psgl, n_of_lig=20000)
    # Don't let them be garbage collected!
    return p, psgl, psgl_plus_esel_bond


def iteration_wrapper(fold_change=1., seed=0, max_time=5, falling_time=1, max_dt=10**-6,
              shear=1, initial_height=0.03, save_every=1000):
    try:
        return iteration(fold_change, seed, max_time, falling_time, max_dt,
              shear, initial_height, save_every)
    except Exception as e:
        print("ERROR")
        raise e


def iteration(fold_change=1., seed=0, shear=1, initial_height=0.03, save_every=1e-4, max_time=5.0):
    p, psgl, psgl_plus_esel_bond = setup_parameters(fold_change)
    ss = SimulationState(initial_height, p, seed)

    # ss.simulate(max_time=1.0, max_steps=int(1e6))  # falling
    # ss.shear_rate = shear
    sim_hist = ss.simulate_with_history(max_time=1.0, max_steps=int(1e6), save_every=save_every)
    ss.shear_rate = shear
    sim_hist = ss.simulate_with_history(max_time=max_time, max_steps=int(max_time * 1e6), save_every=save_every)
    print("simulation done for seed", seed)

    # vel = np.diff(sim_hist.dist) / np.diff(sim_hist.time)
    # # discard first 20% of velocity data
    # vel = vel[len(vel) // 5:]
    return fold_change, 0, sim_hist
    # return fold_change, np.mean(vel), sim_hist


def iteration_callback(result):
    fold_change, mean_v = result
    test_results[fold_change].append(mean_v)


def iteration_error_callback(error):
    print("error callback called!".upper())
    raise error


def one_test(fold_change, n_trials, pool, max_time):
    seeds = np.random.randint(uint_info.max, size=n_trials, dtype='uint')
    # seeds = np.arange(n_trials) + 100
    for seed in seeds:
        pool.apply_async(iteration_wrapper, (fold_change, seed, max_time),
                         callback=iteration_callback, error_callback=iteration_error_callback)


def many_tests(n_trials, pool, max_time):
    for fold_change in FOLD_CHANGES:
        one_test(fold_change, n_trials, pool, max_time)


if __name__ == '__main__':
    N_TRIALS = 100
    MAX_TIME = 5
    FOLD_CHANGES = [
        0.714285714285713,
        1.07142857142856,
        1.78571428571428,
        2.14285714285713,
        3.57142857142856,
        4, 5, 6, 7, 8, 9, 10,
    ]

    uint_info = np.iinfo(np.uint32)

    test_results = defaultdict(list)

    pool = multiprocessing.Pool()

    many_tests(N_TRIALS, pool, MAX_TIME)

    pool.close()
    pool.join()

    directory = f'../results/e_sel_mean_v/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    date_str = datetime.now().strftime("%d-%m-%Y_%Hh%Mm%Ss")
    with open(f'{directory}/{date_str}.pickle', 'wb') as file:
        pickle.dump(test_results, file)
