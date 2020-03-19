import multiprocessing
import pickle
import os
import numpy as np
from collections import Counter
from datetime import datetime

from extrv_engine import Parameters, SimulationState, SlipBondType

N_STEPS_FALLING = int(1e5)
N_STEPS_TEST = int(5e6)

p = Parameters(
    r_c=4.5,
    visc=0.01,
    temp=310,
    dens_diff=0.05,
    # not really using it
    f_rep_0=1e3,
    tau=5
)

# For this parameters simulation makes sense.
psgl = Parameters.LigandType()
psgl_plus_esel_bond = SlipBondType(
    eq_bond_len=20,
    spring_const=2,
    binding_rate_0=1.2e5,
    rec_dens=1,
    react_compl_slip=0.75,
    rup_rate_0_slip=0.01
)
psgl.add_bond_type(psgl_plus_esel_bond)
p.add_ligands(psgl, 10000)

FORCES = [2,]# 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8]
COEFF = 1e-5 * psgl_plus_esel_bond.binding_rate_0
DT = 0.1 / psgl_plus_esel_bond.binding_rate_0


def iteration(test_nr, force, seed):
    ss = SimulationState(h_0=0.0225, p=p, seed=seed)
    ss.simulate(n_steps=N_STEPS_FALLING, dt=DT, shear=0.0)  # falling
    ss.simulate(n_steps=N_STEPS_TEST, dt=DT, shear=force * COEFF, stop_if_no_bonds=True)

    detached = False if ss.bd_lig_ind else True
    print("iteration returning", )
    return test_nr, force, detached, seed


def iteration_wrapper(test_nr, force, seed):
    try:
        iteration(test_nr, force, seed)
    except:
        print("error", test_nr)


def iteration_callback(result):
    print("iteration callback called")
    test_nr, force, detached, seed = result
    # print(f"name: {name}, test nr {test_nr}, force {force}, seed {seed}: ", end="")
    if detached:
        # print("detached.")
        test_results[0][test_nr][force] += 1
    # else:
    #     print("not detached.")


def iteration_error_callback(error):
    print("error callback called!".upper())
    raise error


def one_test(test_nr, n_trials, pool):
    # TODO: seed type
    seeds = np.random.randint(uint_info.max, size=n_trials, dtype='uint')
    # seeds = np.arange(n_trials) + 100
    for force in FORCES:
        # each seed represents one trial
        for seed in seeds:
            pool.apply_async(iteration, (test_nr, force, seed),
                             callback=iteration_callback, error_callback=iteration_error_callback)


def many_tests(n_tests, n_trials, pool):
    for i in range(n_tests):
        one_test(i, n_trials, pool)


if __name__ == '__main__':
    N_TRIALS = 4
    N_TESTS = 1
    uint_info = np.iinfo(np.uint32)
    test_results = ([Counter() for _ in range(N_TESTS)], N_TRIALS)

    pool = multiprocessing.Pool()

    many_tests(N_TESTS, N_TRIALS, pool)

    pool.close()
    pool.join()

    directory = f'results/test96_same_cell/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    date_str = datetime.now().strftime("%d-%m-%Y_%Hh%Mm%Ss")
    with open(f'{directory}/{date_str}.pickle', 'wb') as file:
        pickle.dump(test_results, file)