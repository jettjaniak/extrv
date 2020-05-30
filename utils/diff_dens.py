import time
import numpy as np
from scipy import integrate
import multiprocessing
import pickle
import os
from datetime import datetime


from extrv_engine import SimulationState
from utils.testing_utils import SimulationStats, setup_parameters

FALLING_TIME = 1
ROLLING_TIME = 10
N_TRIALS = 30
MAX_DT = 0.1
SHEAR_RATE = 100
REC_DENS_0 = 1000
BINDING_RATE_0 = 0.09

WALL_ADHESINS_DENSITY_FOLD_CHANGES = [1, 2, 3, 4, 5]

ABS_ERR = 1e-12
REL_ERR = 0
INITIAL_HEIGHT = 0.03
SAVE_EVERY = 1e-4

MAX_STEPS_FALLING = int(2 * FALLING_TIME * 1e6)
MAX_STEPS_ROLLING = int(2 * ROLLING_TIME * 1e6)


def one_test(fold_change, pool_, seeds_, abs_err, rel_err):
    for seed in seeds_:
        pool_.apply_async(iteration_wrapper, (fold_change, seed, abs_err, rel_err),
                          callback=iteration_callback, error_callback=iteration_error_callback)


def iteration_wrapper(*args, **kwargs):
    try:
        fold_change, sim_stats, ss, sim_hist = iteration(*args, **kwargs)
        return fold_change, sim_stats
    except Exception as error:
        print("ERROR RAISED!")
        raise error


def iteration(fold_change=1, seed=0, abs_err=ABS_ERR, rel_err=REL_ERR, save_every=SAVE_EVERY):
    print("simulation starting for fold change", fold_change, "seed", seed)
    p, psgl, psgl_plus_esel_bond = setup_parameters(rec_dens=REC_DENS_0 * fold_change)

    # falling
    ss = SimulationState(INITIAL_HEIGHT, p, seed, max_dt=MAX_DT,
                         ode_abs_err=abs_err, ode_rel_err=rel_err)
    ss.simulate(FALLING_TIME, MAX_STEPS_FALLING)

    n_bonds_before_falling = ss.diag.n_bonds_created - len(ss.bd_lig_ind)

    # rolling
    ss.shear_rate = SHEAR_RATE

    comp_start_time = time.time()
    sim_hist = ss.simulate_with_history(FALLING_TIME + ROLLING_TIME, MAX_STEPS_ROLLING, save_every=save_every)
    comp_end_time = time.time()

    print("simulation is done  for fold change", fold_change, "seed", seed)

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

    return fold_change, sim_stats, ss, sim_hist


def iteration_callback(result):
    fold_change, sim_stats = result
    test_results[fold_change].append(sim_stats)


def iteration_error_callback(error):
    print("ERROR CALLBACK CALLED!")
    raise error


def generate_good_seeds(n_seeds):
    def seed_test(seed_):
        p, psgl, psgl_plus_esel_bond = setup_parameters(rec_dens=REC_DENS_0)
        ss = SimulationState(INITIAL_HEIGHT, p, seed_, max_dt=MAX_DT,
                             ode_abs_err=ABS_ERR, ode_rel_err=REL_ERR)
        ss.simulate(FALLING_TIME, MAX_STEPS_FALLING)
        return ss.diag.n_bonds_created > 0

    uint_info = np.iinfo(np.uint32)
    n_good = 0
    while n_good < n_seeds:
        seed = np.random.randint(uint_info.max, dtype='uint32')
        if seed_test(seed):
            n_good += 1
            yield seed


if __name__ == '__main__':
    directory = '../results/diff_dens/diff_err/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    seeds = tuple(generate_good_seeds(N_TRIALS))

    for abs_err_exp in [-12, -3, -10, -8, -4]:
        abs_err = 10 ** abs_err_exp
        # We use dict comprehension instead of default dict to have ordered keys.
        test_results = {fold_change: [] for fold_change in WALL_ADHESINS_DENSITY_FOLD_CHANGES}

        pool = multiprocessing.Pool()
        for fold_change in WALL_ADHESINS_DENSITY_FOLD_CHANGES:
            one_test(fold_change, pool, seeds, abs_err, 0)
        pool.close()
        pool.join()

        date_str = datetime.now().strftime("%d.%m.%Y_%H-%M-%S")
        with open(f'{directory}/abs1e{abs_err_exp}_rel0_{date_str}.pickle', 'wb') as file:
            pickle.dump(test_results, file)
