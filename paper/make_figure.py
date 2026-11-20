"""Make paper/figure.png for the JOSS paper: four of astrodatapy's 29 quantities, every dataset converted to
one Hubble parameter (and, for stellar masses, one IMF) by the package.

Usage: python make_figure.py        (needs astrodatapy and matplotlib; writes figure.png next to this script)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astrodatapy.number_density import number_density
from astrodatapy.evolution import evolution

H = 0.6751            # Planck 2015
np.seterr(divide='ignore', invalid='ignore')   # lower bounds of 0 (consistent with no objects) give infinite log errors
plt.rcParams.update({'font.size': 9})
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
markers = ['o', 's', '^', 'v', 'D', 'p', 'h', '<', '>', '*', 'P', 'X', 'd', 'H', '8']
colors = plt.cm.tab20(np.linspace(0, 1, 20))


def number_density_panel(ax, feature, z, xlabel, ylabel, **kw):
    obs = number_density(feature=feature, z_target=z, h=H, quiet=True, **kw)
    n = 0
    for i, (name, row) in enumerate(obs.target_observation.iterrows()):
        d, ul = row['Data'], row['ULimit']
        if row['DataType'] != 'data':          # parametric fits and simulations are not drawn here
            continue
        c, m = colors[(2 * n) % 20 + (n // 10) % 2], markers[n % len(markers)]
        ax.errorbar(d[~ul, 0], np.log10(d[~ul, 1]), yerr=[np.log10(d[~ul, 1] / d[~ul, 3]), np.log10(d[~ul, 2] / d[~ul, 1])],
                    fmt=m, color=c, ms=3, mfc='white', elinewidth=0.6, capsize=0, lw=0)
        if ul.any():
            ax.errorbar(d[ul, 0], np.log10(d[ul, 1]), yerr=0.3, uplims=True, fmt=m, color=c, ms=3, mfc='white', elinewidth=0.6)
        n += 1
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('%s, z = %g: %d datasets' % (feature, z, n), fontsize=9)


number_density_panel(axes[0, 0], 'GSMF', 0.0, r'$\log_{10}(M_*/{\rm M}_\odot)$ (Chabrier IMF)',
                     r'$\log_{10}(\Phi/{\rm Mpc^{-3}\,dex^{-1}})$', IMF_out='Chabrier')
axes[0, 0].set_xlim(7.5, 12.3); axes[0, 0].set_ylim(-6, -0.5)
number_density_panel(axes[0, 1], 'GLF_UV', 7.0, r'$M_{\rm UV}$', r'$\log_{10}(\Phi/{\rm Mpc^{-3}\,mag^{-1}})$')
axes[0, 1].set_xlim(-24, -12); axes[0, 1].set_ylim(-8.5, 0.5)
number_density_panel(axes[1, 0], 'QLF_UV', 6.0, r'$M_{1450}$', r'$\log_{10}(\Phi/{\rm Mpc^{-3}\,mag^{-1}})$')
axes[1, 0].set_xlim(-30, -16); axes[1, 0].set_ylim(-12, -3)

ax = axes[1, 1]
obs = evolution('IGM_xHI', quiet=True, h=H)
n_obs = n_sim = 0
for name, row in obs.target_observation.iterrows():
    d = row['Data']
    z, zhi, zlo, y, yhi, ylo = d.T
    if row['DataType'] == 'data_sim':
        ax.plot(z, y, color='0.6', lw=0.6, zorder=1)
        n_sim += 1
        continue
    lim = row['ULimit'] | row['LLimit']
    c, m = colors[(2 * n_obs) % 20], markers[n_obs % len(markers)]
    ax.errorbar(z[~lim], y[~lim], xerr=[z[~lim] - zlo[~lim], zhi[~lim] - z[~lim]], yerr=[np.abs(y[~lim] - ylo[~lim]), np.abs(yhi[~lim] - y[~lim])],
                fmt=m, color=c, ms=3, mfc='white', elinewidth=0.6, lw=0, zorder=3)
    for flag, kw in ((row['ULimit'], dict(uplims=True)), (row['LLimit'], dict(lolims=True))):
        if flag.any():
            ax.errorbar(z[flag], y[flag], yerr=0.4 * y[flag], fmt=m, color=c, ms=3, mfc='white', elinewidth=0.6, zorder=3, **kw)
    n_obs += 1
ax.set_yscale('log')
ax.set_xlim(1.5, 14); ax.set_ylim(1e-6, 1.5)
ax.set_xlabel('redshift')
ax.set_ylabel(r'$\overline{x}_{\rm HI}$')
ax.set_title('IGM_xHI: %d observations, %d simulations (grey)' % (n_obs, n_sim), fontsize=9)

for ax, lab in zip(axes.flat, 'abcd'):
    ax.text(0.03, 0.95, '(%s)' % lab, transform=ax.transAxes, va='top')
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figure.png'), dpi=200)
print('written figure.png')
