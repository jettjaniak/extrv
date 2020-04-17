import numpy as np
import pickle
import matplotlib.pyplot as plt
from utils.diff_shear import SimulationStats

import matplotlib
matplotlib.rcParams.update({
    'errorbar.capsize': 4,
    'font.size': 15
})


def standard_error_of_mean(values):
    return np.sqrt(np.var(values) / len(values))


def plot_one_field(ax, field, solver, test_results):
    res = test_results[solver]
    x = []
    y = []
    err = []
    for shear, stats_list in res.items():
        values = [getattr(stat, field) for stat in stats_list]
        x.append(shear)
        y.append(np.mean(values))
        err.append(standard_error_of_mean(values))

    x = np.array(x)
    y = np.array(y)
    err = np.array(err)

    sort_ind = np.argsort(x)
    x = x[sort_ind]
    y = y[sort_ind]
    err = err[sort_ind]

    ax.errorbar(x, y, err, label=solver)


def plot(all_tests_results):
    fig, axes = plt.subplots(1, 1, sharex=True)
    # fig, axes = plt.subplots(2, 3, sharex=True)
    axes = [axes]
    # axes = axes.flatten()

    for ax, field in zip(axes[-1:], SimulationStats._fields[-1:]):
        for solver_name in all_tests_results.keys():
            plot_one_field(ax, field, solver_name, all_tests_results)
        ax.title.set_text(field)

    # plt.xlabel("shear rate $(1/s)$")
    # plt.ylabel("mean velocity $(\\mu m/s)$")
    # plt.tight_layout()
    plt.yscale('log')
    # plt.legend()