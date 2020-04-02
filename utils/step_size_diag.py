import numpy as np
import matplotlib.pyplot as plt

from extrv_engine import SimulationState
from params import setup_parameters

INITIAL_HEIGHT = 0.03
MAX_SEED = np.iinfo(np.uint32).max


def plot_distribution(max_abs_err_exp=-6, min_abs_err_exp=-10,
                      max_rel_err_exp=-3, min_rel_err_exp=-6):
    seed = np.random.randint(MAX_SEED, dtype='uint32')
    print("seed", seed)

    n_abs_errs = int(max_abs_err_exp - min_abs_err_exp) + 1
    abs_errs = np.logspace(min_abs_err_exp, max_abs_err_exp, n_abs_errs)
    n_rel_errs = int(max_rel_err_exp - min_rel_err_exp) + 1
    rel_errs = np.logspace(min_rel_err_exp, max_rel_err_exp, n_rel_errs)
    dt_freqs = []
    max_len = 0
    for rel_err in rel_errs:
        for abs_err in abs_errs:
            print(f"trying for abs_err={abs_err}, rel_err={rel_err}...", end=' ')
            p, lig, bond = setup_parameters(abs_err=abs_err, rel_err=rel_err)
            ss = SimulationState(INITIAL_HEIGHT, p, seed)
            ss.simulate(max_time=1.0, max_steps=int(1e6))
            ss.shear_rate = 0.5
            ss.simulate(max_time=3.0, max_steps=int(1e7))
            dt_freqs.append(list(ss.diag.dt_freq))
            max_len = max(max_len, len(ss.diag.dt_freq))
            cumsum = np.cumsum(ss.diag.dt_freq)
            cumsum_normalized = cumsum / cumsum[-1]
            print(f"done")
            plt.plot(cumsum_normalized, label=f"abs_err={abs_err}, rel_err={rel_err}")

    max_ticks = min(10, max_len)
    plt.xticks(range(0, max_len, max_len // max_ticks),
               ["1e-" + str(exp) for exp in range(0, max_len, max_len // max_ticks)])
    plt.title("CDFs of $dt$ for different ODE error limits")
    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_distribution(max_abs_err_exp=-9, max_rel_err_exp=-5)
    # 1528777241