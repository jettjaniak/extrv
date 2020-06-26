import os
import pickle
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import seaborn as sns
from utils.diff_shear import SimulationStats
from utils.testing_utils import FIELDS_MAP, FIELDS_Y_LABEL_MAP, DF_COLUMNS, get_err_or_dt_dict

import matplotlib
matplotlib.rcParams.update({
    'errorbar.capsize': 10,
    'font.size': 17
})


def standard_error_of_mean(values):
    return np.sqrt(np.var(values) / len(values))


def standard_error(values):
    return np.sqrt(np.var(values))


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


def plot_diff_err_or_dt(results_dir: str, err_or_dt: str = 'err', base=2, trim=1, min_exp_bound=-float('inf'), ignore_exp=()):
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(6, 10)
    axes = [
        fig.add_subplot(gs[0:3, 0:4]),
        fig.add_subplot(gs[0:3, 4:8]),
        fig.add_subplot(gs[3:6, 0:4]),
        fig.add_subplot(gs[3:6, 4:8])
    ]
    cb_ax = fig.add_subplot(gs[1:-1, -2:])
    cb_ax.set_ylabel("error tolerance $\\epsilon$")

    err_or_dt_dict, min_exp, max_exp = get_err_or_dt_dict(results_dir, err_or_dt, base, min_exp_bound, ignore_exp)
    diff_exp = max_exp - min_exp

    cm_name = 'inferno' if err_or_dt == 'err' else 'viridis'
    cm_ = plt.get_cmap(cm_name)
    # clip 10% of highest values
    vmax = max_exp + diff_exp * 0.1
    # noinspection PyUnresolvedReferences
    cm_norm = matplotlib.colors.Normalize(vmin=min_exp, vmax=vmax)

    def cm(exp_):
        return cm_(cm_norm(exp_))

    for i, (field_name, y_label, ax) in enumerate(zip(FIELDS_MAP.values(), FIELDS_Y_LABEL_MAP.values(), axes)):
        ax.title.set_text(field_name)
        ax.grid(True)
        if i > 1:
            ax.set_xlabel("wall adehsins density (750/$\\mu m^2$)")
        ax.set_ylabel(y_label)

    ylims = []
    i = 0
    for exp, df in err_or_dt_dict.items():
        i -= 1
        means = df.groupby('fold change').mean()
        errors = df.groupby('fold change').apply(standard_error_of_mean)
        field_0 = tuple(FIELDS_MAP.keys())[0]
        xs = means[field_0].index
        xs_ls = np.linspace(xs.min(), xs.max(), 300)
        ylims.append([])

        for field, ax in zip(FIELDS_MAP.keys(), axes):
            color = cm(exp)
            transparent_color = list(color[:3]) + [0.15]
            means_spline = interpolate.pchip(xs, means[field])
            means_smooth = means_spline(xs_ls)
            ax.plot(xs_ls, means_smooth, color=color, linewidth=2, zorder=50 - i)
            ax.scatter(xs, means[field], color=color, zorder=100 - i, s=8**2)

            errors_spline = interpolate.pchip(xs, errors[field])
            errors_smooth = errors_spline(xs_ls)

            errors_high_smooth = means_smooth + errors_smooth
            errors_low_smooth = means_smooth - errors_smooth

            ax.fill_between(xs_ls, errors_high_smooth, errors_low_smooth, color=transparent_color, zorder=30 - i)
            # ax.plot(xs_ls, errors_high_smooth, color=color, linestyle='--', zorder=5)
            # ax.plot(xs_ls, errors_low_smooth, color=cm(exp), linestyle='--', zorder=5)

            # ax.errorbar(means.index, means[field], errors[field], color=cm(exp), linewidth=2)
            ylims[-1].append(ax.get_ylim())

    try:
        for ax, ylim in zip(axes, ylims[- trim - 1]):
            ax.set_ylim(ylim)
    except IndexError:
        print("Can't trim that much!")

    cb_ax.imshow(
        np.linspace(min_exp, max_exp, 100).reshape(-1, 1),
        cmap=cm_name, norm=cm_norm, aspect=12/diff_exp, origin='lower',
        extent=(-0.5, 0.5, min_exp, max_exp)
    )
    cb_ax.set_xticks([])
    # yticks = range(int(np.ceil(min_exp)), np.floor(max_exp + 1)
    # cb_ax.set_yticks(range(min_exp, max_exp + 1),)
    cb_ax.set_yticklabels(["${" + str(base) + "}^{" + str(exp) + "}$" for exp in cb_ax.get_yticks()])
    # dummy_ax.set_visible(False)
    # plt.colorbar(img, orientation="vertical", cax=cb_ax)
    # cb_ax.set_ylim((0, 1))

    plt.tight_layout()


def plot_dt_and_err(dt_exp=-19, err_exp=-13, our_cmap_val=0.9, their_cmap_val=0.35):

    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(6, 8)
    axes = [
        fig.add_subplot(gs[0:3, 0:4]),
        fig.add_subplot(gs[0:3, 4:8]),
        fig.add_subplot(gs[3:6, 0:4]),
        fig.add_subplot(gs[3:6, 4:8])
    ]

    dt_dict, _, _ = get_err_or_dt_dict('results/diff_dens/diff_dt', 'dt', 2)
    dt_df = dt_dict[dt_exp]
    err_dict, _, _ = get_err_or_dt_dict('results/diff_dens/diff_err', 'err', 2)
    err_df = err_dict[err_exp]
    
    alg_field = "algorithm"
    dt_label = "original algorithm, step size $dt = 2^" + "{" + str(dt_exp) + "}$"
    err_label = "our algorithm, error tolerance $\\epsilon = 2^" + "{" + str(err_exp) + "}$"
    
    dt_df[alg_field] = dt_label
    err_df[alg_field] = err_label
    full_df = dt_df.append(err_df)

    their_cmap = plt.get_cmap('viridis')
    our_cmap = plt.get_cmap('inferno')
    palette = [
        their_cmap(their_cmap_val),
        our_cmap(our_cmap_val)
    ]
    for i, (field, ax) in enumerate(zip(FIELDS_MAP.keys(), axes)):
        sns.boxplot(x="fold change", y=field, hue=alg_field,
                    data=full_df, palette=palette, ax=ax)
        if i != 0:
            ax._remove_legend(7)
        else:
            ax.legend(fontsize=15)

    for i, (field_name, y_label, ax) in enumerate(zip(FIELDS_MAP.values(), FIELDS_Y_LABEL_MAP.values(), axes)):
        ax.title.set_text(field_name)
        if i > 1:
            ax.set_xlabel("wall adehsins density (750/$\\mu m^2$)")
        else:
            ax.set_xlabel(None)
        ax.set_ylabel(y_label)

    plt.tight_layout()