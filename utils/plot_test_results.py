import numpy as np
import matplotlib.pyplot as plt
from utils.diff_shear import SimulationStats

import matplotlib
matplotlib.rcParams.update({
    'errorbar.capsize': 10,
    'font.size': 15
})


def standard_error_of_mean(values):
    return np.sqrt(np.var(values) / len(values))


def plot_one_field(ax, field, test_results):
    variables = []
    means = []
    errors = []
    for variable, stats_list in test_results.items():
        values = [getattr(stat, field) for stat in stats_list]
        variables.append(variable)
        ax.scatter(
            np.repeat(variable, len(values)), values,
            s=20, c=[(1, 0, 0, 0.3)]
        )
        means.append(np.mean(values))
        errors.append(standard_error_of_mean(values))

    variables = np.array(variables)
    means = np.array(means)
    errors = np.array(errors)

    ax.errorbar(variables, means, errors, linewidth=2, zorder=10, color='black')


def plot(test_results, title=''):
    # fig, axes = plt.subplots(1, 1, sharex=True)
    fig, axes = plt.subplots(2, 3, sharex=True)
    # axes = [axes]
    axes = axes.flatten()
    if title:
        fig.suptitle(title)

    for ax, field in zip(axes, SimulationStats._fields):
        plot_one_field(ax, field, test_results)
        ax.title.set_text(field)

    # plt.xlabel("shear rate $(1/s)$")
    # plt.ylabel("mean velocity $(\\mu m/s)$")
    plt.tight_layout()
    # plt.yscale('log')
    # plt.legend()