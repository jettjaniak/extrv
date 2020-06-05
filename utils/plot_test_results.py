import math
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from utils.diff_shear import SimulationStats

import matplotlib
matplotlib.rcParams.update({
    'errorbar.capsize': 10,
    'font.size': 18
})


def standard_error_of_mean(values):
    return np.sqrt(np.var(values) / len(values))


def plot_one_field(ax, field, test_results, color='black', label=''):
    variables = []
    means = []
    errors = []
    for variable, stats_list in test_results.items():
        values = [getattr(stat, field) for stat in stats_list]
        variables.append(variable)
        # ax.scatter(
        #     np.repeat(variable, len(values)), values,
        #     s=20, c=[(1, 0, 0, 0.3)]
        # )
        means.append(np.mean(values))
        errors.append(standard_error_of_mean(values))

    variables = np.array(variables)
    means = np.array(means)
    errors = np.array(errors)

    # ax.errorbar(variables, means, errors, linewidth=2, zorder=10, color='black')
    if label:
        ax.plot(variables, means, linewidth=2, zorder=10, color=color, label=label)
    else:
        ax.plot(variables, means, linewidth=2, zorder=10, color=color)



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


def plot_multi(test_results_dict, title=''):
    # fig, axes = plt.subplots(1, 1, sharex=True)
    fig, axes = plt.subplots(2, 3, sharex=True)
    # axes = [axes]
    axes = axes.flatten()
    if title:
        fig.suptitle(title)

    for ax, field in zip(axes, SimulationStats._fields):
        ax.title.set_text(field)
        for tr_item, color in zip(test_results_dict.items(), plt.get_cmap('tab20').colors):
            label, test_result = tr_item
            plot_one_field(ax, field, test_result, color=color, label=label)

    plt.tight_layout()
    plt.legend()


FIELDS_MAP = dict(
    dist='distance traveled in flow direction',
    rot='cell rotation',
    n_bonds='number of all formed bonds',
    mean_h='mean distance from wall'
)

FIELDS_Y_LABEL_MAP = dict(
    dist='distance in $\\mu m$',
    rot='absolute rotation in radians',
    n_bonds='number of bonds',
    mean_h='mean distance in $\\mu m$'
)


def plot_diff_err_or_dt(results_dir: str, err_or_dt: str = 'err'):
    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(6, 9)
    axes = [
        fig.add_subplot(gs[0:3, 0:4]),
        fig.add_subplot(gs[0:3, 4:8]),
        fig.add_subplot(gs[3:6, 0:4]),
        fig.add_subplot(gs[3:6, 4:8])
    ]
    cb_ax = fig.add_subplot(gs[1:-1, -1:])

    err_or_dt_dict = {}
    min_exp = float('inf')
    max_exp = -float('inf')
    for tr_fn in os.listdir(results_dir):
        tr_path = os.path.join(results_dir, tr_fn)
        with open(tr_path, 'rb') as tr_file:
            tr = pickle.load(tr_file)
        tr_arr = []
        for fold_change, stats_l in tr.items():
            for stat in stats_l:
                tr_row = [fold_change] + [abs(getattr(stat, field)) for field in FIELDS_MAP.keys()]
                tr_arr.append(tr_row)
        tr_df = pd.DataFrame(tr_arr, columns=["fold change", *FIELDS_MAP.values()])
        slice_pos = 5 if err_or_dt == 'err' else 4
        tr_exp = int(tr_fn.split('_', maxsplit=1)[0][slice_pos:])
        err_or_dt_dict[tr_exp] = tr_df
        min_exp = min(min_exp, tr_exp)
        max_exp = max(max_exp, tr_exp)
    err_or_dt_dict = {tr_err_exp: tr for tr_err_exp, tr in sorted(err_or_dt_dict.items())}
    diff_exp = max_exp - min_exp

    cm_name = 'inferno' if err_or_dt == 'err' else 'viridis'
    cm_ = plt.get_cmap(cm_name)
    # clip 10% of highest values
    vmax = max_exp + diff_exp * 0.1
    cm_norm = matplotlib.colors.Normalize(vmin=min_exp, vmax=vmax)

    def cm(exp):
        return cm_(cm_norm(exp))

    for field, y_label, ax in zip(FIELDS_MAP.values(), FIELDS_Y_LABEL_MAP.values(), axes):
        ax.title.set_text(field)
        ax.grid(True)
        ax.set_xlabel("wall adehsins density (1000/$\\mu m^2$)")
        ax.set_ylabel(y_label)

    for exp, df in err_or_dt_dict.items():
        means = df.groupby('fold change').mean()
        for field, ax in zip(FIELDS_MAP.values(), axes):
            ax.plot(means[field], color=cm(exp), linewidth=2)

    cb_ax.imshow(
        np.linspace(min_exp, max_exp, 100).reshape(-1, 1),
        cmap=cm_name, norm=cm_norm, aspect=12/diff_exp, origin='lower',
        extent=(-0.5, 0.5, min_exp, max_exp)
    )
    cb_ax.set_xticks([])
    cb_ax.set_yticks(range(min_exp, max_exp + 1),)
    cb_ax.set_yticklabels(["$10^{" + str(exp) + "}$" for exp in range(min_exp, max_exp + 1)])
    # dummy_ax.set_visible(False)
    # plt.colorbar(img, orientation="vertical", cax=cb_ax)
    # cb_ax.set_ylim((0, 1))ibib

    plt.tight_layout()
