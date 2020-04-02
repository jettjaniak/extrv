import numpy as np
import matplotlib.pyplot as plt

from extrv_engine import SimulationState
from params import setup_parameters

INITIAL_HEIGHT = 0.03
MAX_SEED = np.iinfo(np.uint32).max


def plot_distribution(max_abs_err_exp=-6, min_abs_err_exp=-10):
    seed = np.random.randint(MAX_SEED, dtype='uint32')
    print("seed", seed)

    n_abs_errs = int(max_abs_err_exp - min_abs_err_exp) + 1
    abs_errs = np.logspace(min_abs_err_exp, max_abs_err_exp, n_abs_errs)
    dt_freqs = []
    max_len = 0
    for abs_err in abs_errs:
        p, lig, bond = setup_parameters(abs_err=abs_err)
        ss = SimulationState(INITIAL_HEIGHT, p, seed)
        ss.simulate(max_time=1.0, max_steps=int(1e6))
        ss.shear_rate = 0.5
        ss.simulate(max_time=3.0, max_steps=int(1e7))
        dt_freqs.append(list(ss.diag.dt_freq))
        max_len = max(max_len, len(ss.diag.dt_freq))
        print("done for abs_err", abs_err)

    plt.figure()
    plt.title("CDFs of $dt$ for different ODE error limits")

    for dt_freq, abs_err in zip(dt_freqs, abs_errs):
        cumsum = np.cumsum(dt_freq)
        cumsum_normalized = cumsum / cumsum[-1]
        plt.plot(cumsum_normalized, label="abs_err=" + str(abs_err))
        plt.legend()
        plt.show(block=False)

    max_ticks = min(10, max_len)
    plt.xticks(range(0, max_len, max_len // max_ticks),
               ["1e-" + str(exp) for exp in range(0, max_len, max_len // max_ticks)])

    plt.show()


if __name__ == '__main__':
    plot_distribution(-8)
    # 1528777241