#!/usr/bin/env python3
from textwrap import dedent
import numpy as np
import numexpr as ne
import scipy.stats as stats
import uncertainties as u
from uncertainties import unumpy as unp
from uncertainties import ufloat
import matplotlib as mpl
mpl.use('Agg')
import atlas_mpl_style as ampl
import matplotlib.pyplot as plt  # noqa
import pandas as pd  # noqa
import seaborn  # noqa

pd.options.mode.chained_assignment = 'raise'
samples = ['2b2j', '4b', '4j']
channels = ['res']  # , 'inter', 'boost']
n_tags = [2, 3]
n_bins = 30

sigSel = 'region == "SIG"'
sdba_sel = 'region == "SDBA"'
sdb_all_sel = 'region != "SIG"'
variables = {'m_HH': '$m_{hh}$', 'pt_H0': '$p_T(h_0)$', 'pt_H1': '$p_T(h_1)$'}

print("Reading samples...")
ntuple2t = {
    sample: pd.read_csv(f"baseline_noPU_atlas_qcd/ntuple/{sample}-fullNTuple-2-tag.dat")
    for sample in samples
}
print("2 tag read")
ntuple3t = {
    sample: pd.read_csv(f"baseline_noPU_atlas_qcd/ntuple/{sample}-fullNTuple-3-tag.dat")
    for sample in samples
}
print("3 tag read")
ntuple4t = {
    sample: pd.read_csv(f"baseline_noPU_atlas_qcd/ntuple/{sample}-fullNTuple-4-tag.dat")
    for sample in samples
}
print("4 tag read")

plt.style.use('paper')
seaborn.set_palette('deep')
plt.ioff()


def weighted_chisquare(f_obs, f_exp):
    "Calculate weighted chi-square using method in arXiv:physics/0605123"
    w1 = unp.nominal_values(f_obs)
    w2 = unp.nominal_values(f_exp)
    s1 = unp.std_devs(f_obs)  # noqa
    s2 = unp.std_devs(f_exp)  # noqa
    W1 = np.sum(w1)  # noqa
    W2 = np.sum(w2)  # noqa
    X2 = ne.evaluate(
        "sum((W1*w2 - W2*w1)**2 / (W1**2 * s2**2 + W2**2 * s1**2))")
    p = stats.chi2.sf(X2, np.size(w1) - 1)
    return (X2, p)


def consistent(a, b):
    total_err = np.sqrt(a.s**2 + b.s**2)
    if np.abs(a.n - b.n) <= total_err:
        return True
    else:
        return False


def make_hist(a, weights, bins=n_bins, range_=None):
    if isinstance(weights, pd.Series):
        hasUnc = isinstance(weights.values[0], u.UFloat)
        wgts = weights.values
        # Hopefully working with a numpy array speeds things up a bit
    else:
        hasUnc = isinstance(weights[0], u.UFloat)
        wgts = weights
    if hasUnc:
        raise Exception(
            'make_hist does not support weights with uncertainties (too slow)')
    else:
        y, x = np.histogram(a, bins=bins, weights=wgts, range=range_)
        w2, _ = np.histogram(a, bins=bins, weights=(wgts**2), range=range_)
        y = unp.uarray(y, std_devs=np.sqrt(w2))
    return x, y


def plot_hist(bins, values, label, color, ax=None):
    if ax is None:
        ax = plt.gca()
    vals = unp.nominal_values(values)
    err_high = vals + unp.std_devs(values)
    err_low = vals - unp.std_devs(values)

    # double up
    vals = np.stack((vals, vals)).ravel(1)
    err_high = np.stack((err_high, err_high)).ravel(1)
    err_low = np.stack((err_low, err_low)).ravel(1)

    x = np.stack((bins[:-1], bins[1:])).ravel(1)

    ax.fill_between(x, err_high, err_low, color=color, alpha=0.3)
    ax.plot(x, vals, label=label, color=color, lw=1)
    return ax


def plot_ratio(bins, real, est, ax=None):
    if ax is None:
        ax = plt.gca()

    real_vals = unp.nominal_values(real)
    real_vals[real_vals == 0.] = np.nan
    ratio = (est - real_vals) / real_vals

    # double up real
    real = np.stack((real, real)).ravel(1)
    real_vals = np.stack((real_vals, real_vals)).ravel(1)
    real_err = np.nan_to_num(unp.std_devs(real) / real_vals)
    x = np.stack((bins[:-1], bins[1:])).ravel(1)
    binCenters = (bins[:-1] + bins[1:]) / 2
    ax.set_ylim(-2 * np.max(real_err), 2 * np.max(real_err))
    ax.fill_between(x, real_err, -real_err, color="grey", alpha=0.5)
    ax.axhline(color="grey")
    ax.errorbar(
        binCenters,
        unp.nominal_values(ratio),
        fmt='o',
        yerr=unp.std_devs(ratio),
        color="C0")
    return ax


def plot_matrix(analysis_channel, n_tag, sdbSel):
    if n_tag == 2:
        ntuple = ntuple2t
    elif n_tag == 3:
        ntuple = ntuple3t
    else:
        raise ValueError(f"n_tag invalid : {n_tag} /= 2 or 3")

    # Sanity filters
    sfilt = (f'analysis_channel == "{analysis_channel}" and abs(m_HH) < 2000'
             f' and abs(pt_H0) < 900 and abs(pt_H1) < 900')
    all_bkg = pd.concat(ntuple.values()).query(sfilt)
    all_bkg.m_HH = np.abs(all_bkg.m_HH)
    sig = all_bkg.query(sigSel)
    sdba = all_bkg.query(sdbSel)

    all_bkg4t = pd.concat(ntuple4t.values()).query(sfilt)
    all_bkg4t.m_HH = np.abs(all_bkg4t.m_HH)
    sig4t = all_bkg4t.query(sigSel)
    sdba4t = all_bkg4t.query(sdbSel).copy()

    plot_count = 1
    fig = plt.figure(figsize=(16, 10), dpi=600, facecolor='paper:bg')
    for bin_by in variables:
        rescale_bin_range = (min(sig[bin_by].min(), sdba[bin_by].min(),
                                 sig4t[bin_by].min(), sdba4t[bin_by].min()),
                             max(sig[bin_by].max(), sdba[bin_by].max(),
                                 sig4t[bin_by].max(), sdba4t[bin_by].max()))
        bins = np.linspace(rescale_bin_range[0] - 0.1,
                           rescale_bin_range[1] + 0.1, n_bins)
        _, sig_hist = make_hist(sig[bin_by], sig.weight, bins)
        _, sdba_hist = make_hist(sdba[bin_by], sdba.weight, bins)

        mask = unp.nominal_values(sdba_hist) == 0.0
        sdba_hist[mask] = ufloat(np.nan, 0)

        scale_factor = sig_hist / sdba_hist
        index = np.digitize(sdba4t[bin_by], bins) - 1
        sdba4t.loc[:, f'{bin_by}_sf'] = np.choose(
            index, unp.nominal_values(scale_factor))
        sdba4t.loc[:, f'{bin_by}_sf_err'] = np.choose(
            index, unp.std_devs(scale_factor))
        for hist_of in variables:
            bin_range = (min(sig4t[hist_of].min(), sdba4t[hist_of].min()),
                         max(sig4t[hist_of].max(), sdba4t[hist_of].max()))
            if hist_of != bin_by:
                hist_bins = np.linspace(bin_range[0], bin_range[1], n_bins)
            else:
                hist_bins = bins

            x, real_4t = make_hist(sig4t[hist_of], sig4t.weight, hist_bins)
            _, est_4t = make_hist(sdba4t[hist_of],
                                  sdba4t.weight * sdba4t[f'{bin_by}_sf'], x)

            ax = plt.subplot(len(variables), len(variables), plot_count)
            total_est_yield = np.sum(est_4t)
            total_real_yield = np.sum(real_4t)
            norm = total_real_yield.n / total_est_yield.n
            plot_hist(x, real_4t, "Real", 'C1', ax)
            plot_hist(x, est_4t, "Estimated", 'C0', ax=ax)
            if bin_by == hist_of:
                trans = ax.get_xaxis_transform()
                for blank, (start, end) in zip(mask, zip(bins[:-1], bins[1:])):
                    if blank:
                        rect = mpl.patches.Rectangle(
                            (start, 0),
                            width=(end - start),
                            height=1,
                            transform=trans,
                            color='grey',
                            ec='grey',
                            alpha=0.8)
                        ax.add_patch(rect)
            ax.set_xlabel(f'{variables[hist_of]} / GeV')
            ax.set_title(f"{variables[hist_of]} scaled in {variables[bin_by]}")
            ax.legend(loc=1)
            x2, x2_pval = weighted_chisquare(
                f_obs=real_4t, f_exp=est_4t * norm)
            ax.text(
                0.45,
                0.95,
                dedent(f'''\
                    $\\chi^{{2}} = {x2:6.3g}$   $P = {x2_pval:6.3f}$
                    {"Real Yield: ":17}${total_real_yield:L}$
                    {"Estimated Yield: ":17}${total_est_yield:L}$'''),
                horizontalalignment='center',
                verticalalignment='top',
                multialignment='center',
                transform=ax.transAxes,
                fontsize=8,
                color='black',
                bbox={
                    "fc": ("on:red" if x2_pval <= 0.25 else
                           ("on:orange" if x2_pval <= 0.50 else
                            ("on:yellow" if x2_pval <= 0.75 else
                              "on:green"))),
                    "ec":
                        "none",
                    "alpha":
                        0.5
                })
            plot_count += 1
    plt.tight_layout(pad=1.0)

    fig.savefig(
        f'bkg_estimates/crosstab-{analysis_channel}-channel-{n_tag}-tag'
        f'{"(all-sideband)" if sdbSel == sdb_all_sel else ""}.svg',
        facecolor='paper:bg')
    print(
        dedent(f'''\
                  Completed:
                      Channel: {analysis_channel}
                      Tags: {n_tag}
                      Sideband: {"A only" if sdbSel==sdba_sel else "All"}
        '''))


for analysis_channel in channels:
    for n_tag in n_tags:
        for sdbSel in [sdba_sel, sdb_all_sel]:
            plot_matrix(analysis_channel, n_tag, sdbSel)
