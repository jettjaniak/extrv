import numpy as np
import multiprocessing
import pickle
import os
from datetime import datetime
from collections import defaultdict

from extrv_engine import SlipBondType, Parameters, SimulationState


def setup_parameters(fold_change):
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        # TODO: remove, use nonspecific
        f_rep_0=1e3,
        tau=5
    )
    psgl = Parameters.LigandType()
    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=27,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=750 * fold_change,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6
    )
    psgl.add_bond_type(psgl_plus_esel_bond)
    p.add_ligands(lig_type=psgl, n_of_lig=20000)
    return p


def iteration(fold_change=1., seed=0, max_time=5, falling_time=1, max_dt=10**-6,
              shear=1, initial_height=0.03, save_every=1000):
    p = setup_parameters(fold_change)
    k_on_0 = p.bond_types[0].rec_dens * p.bond_types[0].binding_rate_0
    dt = min(0.1 / k_on_0, max_dt)
    n_steps = int(max_time / dt)
    n_steps_falling = int(falling_time / dt)

    ss = SimulationState(initial_height, p, seed)
    ss.simulate(n_steps_falling, dt, shear=0.0)  # falling
    result = ss.simulate_with_history(n_steps, dt, shear, save_every=save_every)

    # discard first 20% of velocity data
    # TODO: velocity in History
    return fold_change, np.mean(result['v'][(n_steps // save_every) // 5:, 0])


def iteration_callback(result):
    fold_change, mean_v = result
    test_results[fold_change].append(mean_v)


def iteration_error_callback(error):
    print("error callback called!".upper())


def one_test(fold_change, n_trials, pool, max_time):
    seeds = np.random.randint(uint_info.max, size=n_trials, dtype='uint')
    for seed in seeds:
        pool.apply_async(iteration, (fold_change, seed, max_time),
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
