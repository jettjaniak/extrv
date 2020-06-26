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
MAX_DT = 0.1
SHEAR_RATE = 10
REC_DENS_0 = 750
BINDING_RATE_0 = 0.06

WALL_ADHESINS_DENSITY_FOLD_CHANGES = [1, 2, 3, 4, 5]

ABS_ERR = 1e-12
REL_ERR = 0
INITIAL_HEIGHT = 0.03
SAVE_EVERY = 1e-4

MAX_STEPS_FALLING = int(2 * FALLING_TIME * 1e6)
MAX_STEPS_ROLLING = int(2 * ROLLING_TIME * 1e6)


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


def iteration_wrapper(*args, **kwargs):
    try:
        iteration(*args, **kwargs)
    except Exception as error:
        print("ERROR RAISED!")
        raise error


def iteration_error_callback(error):
    print("ERROR CALLBACK CALLED!")
    raise error


def iteration(stats_list, fold_change=1, seed=0, abs_err=ABS_ERR, save_every=SAVE_EVERY):
    print("simulation starting for fold change", fold_change, "seed", seed)
    p, psgl, psgl_plus_esel_bond = setup_parameters(rec_dens=REC_DENS_0 * fold_change)

    # falling
    ss = SimulationState(INITIAL_HEIGHT, p, seed, max_dt=MAX_DT,
                         ode_abs_err=abs_err, ode_rel_err=REL_ERR)
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

    stats_list.append(sim_stats)
    return fold_change, sim_stats, ss, sim_hist


def for_err_test(err_base, err_exp, nice_inc=0):
    os.nice(nice_inc)
    err = err_base ** err_exp
    with multiprocessing.Manager() as manager:
        # We use dict comprehension instead of default dict to have ordered keys.
        test_results = {fold_change: [] for fold_change in WALL_ADHESINS_DENSITY_FOLD_CHANGES}
        test_results = manager.dict(test_results)
        
        pool = multiprocessing.Pool()
        for fold_change in WALL_ADHESINS_DENSITY_FOLD_CHANGES:
            stats_list = manager.list()
            test_results[fold_change] = stats_list
            for seed in seeds:
                pool.apply_async(iteration_wrapper, (stats_list, fold_change, seed, err),
                                 error_callback=iteration_error_callback)
        pool.close()
        pool.join()

        tr_cast = {}
        for fold_change in test_results.keys():
            tr_cast[fold_change] = list(test_results[fold_change])

    date_str = datetime.now().strftime("%d.%m.%Y_%H-%M-%S")
    with open(f'{directory}/err{err_base}^{err_exp}_{date_str}.pickle', 'wb') as file:
        pickle.dump(tr_cast, file)


if __name__ == '__main__':
    directory = '../results/diff_dens/diff_err/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # 80 seeds
    GOOD_SEEDS = (1731712957, 3190037783, 3759517182, 4177028560, 1830405957,
                  3762189420, 1546919209, 752713666, 2277760481, 2225307367, 771879583, 2686016408, 3091179815,
                  54641801, 822526847, 3485949959, 1808968316, 2740946722, 2128909877, 1735928740, 709544941,
                  256163407, 2950738251, 67318337, 3705137511, 23055023, 3076940865, 2962143592, 3516257955,
                  3987088566, 997047223, 3296440446, 3213720315, 501563541, 3609636494, 3947717566, 1205502233,
                  3141042261, 2198693213, 1654518378, 8356930, 682264254, 1174718341, 556715953, 3943414791,
                  3560480401, 1483454036, 1038552105, 325198726, 712250723, 783014756, 2310304458, 1895812771,
                  2484051420, 622190755, 1641609963, 3466489807, 1441286071, 2419453959, 1533929822, 2464008272,
                  4191568628, 1751682214, 2288410354, 3596124070, 1820690884, 2414023370, 2264240200, 3133142593,
                  3102573477, 2513463839, 1460214160, 629244395, 2998528036, 1651710757, 156828413, 2688925137,
                  1635904256, 1031466595, 2716957902)

    seeds = GOOD_SEEDS[0:30]

    processes = []
    nice_inc_ = 5
    for err_exp_ in [-13]:
        proc_kwargs = dict(err_base=2, err_exp=err_exp_, nice_inc=nice_inc_)
        proc = multiprocessing.Process(target=for_err_test, kwargs=proc_kwargs)
        nice_inc_ += 1
        processes.append(proc)
        proc.start()

    for proc in processes:
        proc.join()
