"""Draw Figures 1-4 of the manuscript (PDF and PNG) into paper/figures/.

Usage: python scripts/make_figures.py
"""
import csv
import json
import math
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lrdwidth import stats as S
from lrdwidth.paths import DATA, PAPER, PHOTOSPHERE, RESULTS, WIDTHS

OUT = str(PAPER / "figures")
os.makedirs(OUT, exist_ok=True)

RED = '#B8323C'; BLUE = '#2F6FA7'; INK = '#1a1a1a'; MUTE = '#6b6b6b'; GOLD = '#B07A12'
mpl.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 9,
    'axes.labelsize': 9.5, 'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'legend.fontsize': 8.3, 'axes.linewidth': 0.7, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True, 'xtick.major.size': 3.5, 'ytick.major.size': 3.5,
    'xtick.minor.size': 2, 'ytick.minor.size': 2, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5, 'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'legend.frameon': False, 'axes.edgecolor': INK, 'text.color': INK, 'axes.labelcolor': INK,
    'xtick.color': INK, 'ytick.color': INK, 'pdf.fonttype': 42, 'savefig.dpi': 300, 'legend.handletextpad': 0.3,
})
W1 = 3.4; W2 = 7.1

def save(fig, name):
    fig.savefig(f'{OUT}/{name}.pdf')
    fig.savefig(f'{OUT}/{name}.png')
    plt.close(fig)

def kms_axis(ax, lo, hi, ticks=(1000, 1500, 2000, 3000, 4000, 6000), axis='y'):
    tk = [t for t in ticks if lo <= t <= hi]
    mn = [np.log10(t) for t in range(500, 9000, 100) if lo <= t <= hi]
    if axis == 'y':
        ax.set_ylim(np.log10(lo), np.log10(hi))
        ax.set_yticks([np.log10(t) for t in tk]); ax.set_yticklabels([str(t) for t in tk])
        ax.set_yticks(mn, minor=True)
    else:
        ax.set_xlim(np.log10(lo), np.log10(hi))
        ax.set_xticks([np.log10(t) for t in tk]); ax.set_xticklabels([str(t) for t in tk])
        ax.set_xticks(mn, minor=True)

def fit2(rows):
    x = np.array([r['x'] for r in rows]); y = np.array([r['y'] for r in rows]); C = np.array([r['cov'] for r in rows])
    b, s, p, ah = S.post2(x, y, C[:, 0, 0], C[:, 1, 1], C[:, 0, 1])
    return float((p*ah).sum()), S.quant(b, p.sum(1))[2], S.quant(s, p.sum(0))

# ------------------------------------------------------------------ Figure 1: the relation in both data sets
W1in = json.load(open(WIDTHS / 'regression_inputs.json'))
G = W1in['single_gaussian']; Q = json.load(open(DATA / 'bx26' / 'regression_inputs.json'))
qso = json.load(open(DATA / 'qso1' / 'qso1_fit.json'))
ov = {r['published_id'] for r in csv.DictReader(open(DATA / 'bx26' / 'overlap_with_this_work.tsv'), delimiter='\t')}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 3.2))
for cls, col, mk, z in (('non-LRD', BLUE, '^', 2), ('LRD', RED, 'o', 4)):
    r = [t for t in G if (t['cls'] == 'LRD') == (cls == 'LRD')]
    x = np.array([t['x'] for t in r])+43; y = np.array([t['y'] for t in r]); C = np.array([t['cov'] for t in r])
    a, b, sq = fit2(r); s = sq[2]; xx = np.linspace(x.min()-.06, x.max()+.06, 50); yy = a+b*(xx-43)
    a1.fill_between(xx, yy-s, yy+s, color=col, alpha=.14 if cls == 'LRD' else .09, lw=0, zorder=z-2)
    a1.plot(xx, yy, color=col, lw=1.3, ls='-' if cls == 'LRD' else (0, (4, 2)), zorder=z)
    a1.errorbar(x, y, xerr=np.sqrt(C[:, 0, 0]), yerr=np.sqrt(C[:, 1, 1]), fmt='none', ecolor=col, elinewidth=.6, alpha=.55, zorder=z)
    a1.scatter(x, y, s=24 if mk == 'o' else 28, marker=mk, facecolor=col, edgecolor='white', linewidth=.6, zorder=z+1)
xq = math.log10(qso['L_delensed_erg_s']); yq = math.log10(qso['FWHM_kms']); eq = qso['FWHM_err_conditional_kms']/qso['FWHM_kms']/math.log(10)
exq = math.hypot(math.sqrt(qso['log_cov_xy'][0][0]), qso['mu_error']/qso['mu']/math.log(10))
a1.errorbar(xq, yq, xerr=exq, yerr=eq, fmt='none', ecolor=INK, elinewidth=.7, zorder=6)
a1.scatter(xq, yq, s=34, marker='D', facecolor='white', edgecolor=INK, linewidth=.9, zorder=7)
a1.annotate('A2744-QSO1', (xq, yq), xytext=(6, 7), textcoords='offset points', fontsize=8, color=INK)
a1.legend(handles=[Line2D([], [], marker='o', color=RED, mec='white', ms=5.5, lw=1.3, label='LRDs (29)'),
                   Line2D([], [], marker='^', color=BLUE, mec='white', ms=6, lw=1.3, ls=(0, (4, 2)), label='other broad-line sources (15)')],
          loc='lower right', handlelength=2.6, borderaxespad=.6)
a1.set_xlim(41.45, 43.72); kms_axis(a1, 700, 4600)
a1.set_xlabel(r'$\log\,(L_{\rm H\alpha}\,/\,{\rm erg\,s^{-1}})$'); a1.set_ylabel(r'broad H$\alpha$ FWHM  (km s$^{-1}$)')
a1.text(.04, .95, '(a) this work', transform=a1.transAxes, va='top', fontsize=9.5)
def bxrow(t): return dict(id=t['id'], x=t['logLbol']-45, y=t['y'], cov=[[t['sx']**2, 0], [0, t['sy']**2]], z=t['z'], cls=t['cls'])
R = [bxrow(t) for t in Q]; XL = [r for r in R if r['cls'] == 'LRD']; lo = min(r['x'] for r in XL); hi = max(r['x'] for r in XL)
XC = [r for r in R if r['cls'] != 'LRD' and r['z'] > 3 and lo <= r['x'] <= hi]; assert len(XL) == 48 and len(XC) == 97
rng = np.random.default_rng(20261112)
for rows, col, mk, z, key in ((XC, BLUE, '^', 2, 'C'), (XL, RED, 'o', 4, 'LRD')):
    x = np.array([t['x'] for t in rows])+45; y = np.array([t['y'] for t in rows]); o = np.array([t['id'] in ov for t in rows])
    xj = x+rng.uniform(-.032, .032, len(x))
    a, b, sq = fit2(rows); s = sq[2]; xx = np.linspace(lo+45-.07, hi+45+.07, 50); yy = a+b*(xx-45)
    a2.fill_between(xx, yy-s, yy+s, color=col, alpha=.14 if key == 'LRD' else .08, lw=0, zorder=z-2)
    a2.plot(xx, yy, color=col, lw=1.3, ls='-' if key == 'LRD' else (0, (4, 2)), zorder=z)
    sz = 22 if mk == 'o' else 20
    a2.scatter(xj[~o], y[~o], s=sz, marker=mk, facecolor=col, edgecolor='white', linewidth=.5, zorder=z+1, alpha=1 if key == 'LRD' else .8)
    a2.scatter(xj[o], y[o], s=sz-4, marker=mk, facecolor='white', edgecolor=col, linewidth=.9, zorder=z+1)
a2.legend(handles=[Line2D([], [], marker='o', color=RED, mec='white', ms=5.5, lw=1.3, label='LRDs (48)'),
                   Line2D([], [], marker='^', color=BLUE, mec='white', ms=5.5, lw=1.3, ls=(0, (4, 2)), label='comparison sources (97)'),
                   Line2D([], [], marker='o', color='none', mec=MUTE, mfc='white', ms=5, label='also in our sample')],
          loc='lower right', handlelength=2.6, borderaxespad=.6)
a2.set_xlim(43.55, 45.35); kms_axis(a2, 540, 7600)
a2.set_xlabel(r'$\log\,(L_{\rm bol}\,/\,{\rm erg\,s^{-1}})$, from broad H$\alpha$'); a2.set_ylabel(r'broad H$\alpha$ FWHM  (km s$^{-1}$)')
a2.text(.04, .95, '(b) Baccus & Xu table', transform=a2.transAxes, va='top', fontsize=9.5)
fig.subplots_adjust(left=.075, right=.99, bottom=.135, top=.975, wspace=.22); save(fig, 'Figure1')

# ------------------------------------------------------------------ Figure 2: scatter under every width definition
r1_raw = json.load(open(RESULTS / 'width_luminosity.json'))
r1 = {('W1_' + k[6:] if k.startswith('width_') else k): v for k, v in r1_raw.items()}
r2 = json.load(open(RESULTS / 'robustness.json'))
rows = [('single-Gaussian FWHM', 'W1_single_gaussian'), (r'interquartile width $W_{50}$', 'W1_nonparametric_W50'),
        ('exponential FWHM', 'W1_exponential'), ('Lorentzian FWHM', 'W1_lorentzian'),
        ('nonparametric FWHM', 'W1_nonparametric_FWHM'), ('two-Gaussian FWHM', 'W1_double_gaussian')]
fig, ax = plt.subplots(figsize=(W1, 3.75)); y = 0; yt = []; yl = []
def head(t, yy): ax.text(-.03, yy, t, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=8.6, fontweight='bold')
def row(lab, qL, qO):
    global y
    for q, col, mk, dy in ((qL, RED, 'o', -.17), (qO, BLUE, '^', .17)):
        if q is None: continue
        ax.plot([q[0], q[4]], [y+dy]*2, color=col, lw=.9, solid_capstyle='butt')
        ax.plot([q[1], q[3]], [y+dy]*2, color=col, lw=2.6, solid_capstyle='butt')
        ax.plot(q[2], y+dy, marker=mk, color=col, mec='white', mew=.6, ms=5.2 if mk == 'o' else 5.6, zorder=5)
    yt.append(y); yl.append(lab); y += 1
head('This work', y); y += 1
for lab, k in rows: row(lab, r1[k]['LRD']['sigma'], r1[k]['other']['sigma'])
row(r'LSF width $\times0.5$', r2['frozen_LSFhalf']['LRD']['sigma'], r2['frozen_LSFhalf']['other']['sigma'])
y += .35; head('Baccus & Xu', y); y += 1
row('single-Gaussian FWHM', r2['bx_main']['LRD']['sigma'], r2['bx_main']['other']['sigma'])
sd = [0.112282, None, 0.123240, None, 0.136016]
ax.axvspan(sd[0], sd[4], color=INK, alpha=.16, lw=0, zorder=0); ax.axvline(sd[2], color=INK, lw=.8, zorder=1)
ax.text(sd[2]+.008, 8.35, 'SDSS quasars', va='center', ha='left', fontsize=8.3, color=INK)
ax.axvspan(.05, .08, color=RED, alpha=.07, lw=0, zorder=0)
ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_ylim(y-.35, -.7); ax.set_xlim(0, .47); ax.set_xlabel(r'intrinsic scatter $\sigma_{\rm int}$ (dex)')
ax.tick_params(axis='y', which='both', left=False, right=False)
ax.legend(handles=[Line2D([], [], marker='o', color=RED, mec='white', ms=5.2, lw=.9, label='LRDs'),
                   Line2D([], [], marker='^', color=BLUE, mec='white', ms=5.6, lw=.9, label='other broad-line sources')],
          loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=2, handlelength=1.8, borderaxespad=.2, columnspacing=1.0)
fig.subplots_adjust(left=.43, right=.965, top=.935, bottom=.115); save(fig, 'Figure2')

# ------------------------------------------------------------------ Figure 3: orbital scatter budget (unified 95% ceiling)
rep_cal = json.load(open(RESULTS / 'repeat_calibration.json'))
S0 = r1['W1_single_gaussian']['LRD']['sigma'][2]
S1 = rep_cal['posterior_with_floors']['formal_zero_floor']['LRD_p95_onesided']  # exact one-sided 95% upper limit: 0.07474 dex
ori = r2['orient']
ori_edge60 = ori['thin_60_90']
fig, a = plt.subplots(figsize=(W1, 3.1)); th = np.linspace(0, np.pi/2, 300)
def ell(s, sr=0.):
    e = math.sqrt(s*s-sr*sr/4); return 2*e*np.cos(th), e*np.sin(th)
x1, y1 = ell(S1); x0, y0 = ell(S0)
a.fill_between(x1[::-1], 0, y1[::-1], color=RED, alpha=.07, lw=0)
a.fill_between(x0[::-1], 0, y0[::-1], color=RED, alpha=.15, lw=0)
a.plot(x0, y0, color=RED, lw=1.5); a.plot(x1, y1, color=RED, lw=1.1, ls=(0, (4, 2)))
lab = [(r'thin disk, $i<45^\circ$', ori['thin_0_45'], 'bottom', '-'),
       (r'$H/R=0.5$', ori['h0.5_0_60'], 'bottom', (0, (5, 1.5, 1, 1.5))),
       (r'$H/R=0.7$', ori['h0.7_0_60'], 'bottom', (0, (5, 1.5, 1, 1.5))),
       (r'$H/R=1$', ori['h1_0_60'], 'bottom', (0, (5, 1.5, 1, 1.5))),
       (r'thin disk, $i>60^\circ$', ori_edge60, 'top', '-')]
for k, v, va, ls in lab:
    a.axhline(v, color=BLUE, lw=.8, ls=ls)
    a.text(.236, v+(.003 if va == 'bottom' else -.0025), k, color=BLUE, fontsize=8.3, ha='right', va=va)
a.text(.012, .026, 'allowed', color=RED, fontsize=9, va='center')
a.legend(handles=[Line2D([], [], color=RED, lw=1.5, label=r'$\sigma_{\rm int}=%.3f$ dex' % S0),
                  Line2D([], [], color=RED, lw=1.1, ls=(0, (4, 2)), label=r'$%.3f$ dex (95%% upper limit)' % S1)],
         loc='upper left', bbox_to_anchor=(0, .84), handlelength=2.4, borderaxespad=.7)
a.set_xlim(0, .24); a.set_ylim(0, .235)
a.set_xlabel(r'scatter of $\log(M/R)$ at fixed luminosity (dex)'); a.set_ylabel(r'orientation scatter of $\log$ FWHM (dex)')
fig.subplots_adjust(left=.155, right=.975, bottom=.14, top=.975); save(fig, 'Figure3')

# ------------------------------------------------------------------ Figure 4: 3-panel photospheric comparison & scaling tests
SXMAX = .3
_PHOT_MAP = {
    'A_Curtis_BX26_sources.tsv': PHOTOSPHERE / 'atmosphere_bx26_widths.tsv',
    'B_Curtis_frozen_sources.tsv': PHOTOSPHERE / 'atmosphere_our_widths.tsv',
    'C_deGraaff_unique_sources.tsv': PHOTOSPHERE / 'blackbody_sources.tsv',
}
def load(fn):
    rr = list(csv.DictReader(open(_PHOT_MAP[fn]), delimiter='\t')); f = lambda k: np.array([float(r[k]) for r in rr])
    return rr, f('x'), f('t'), f('y'), f('sx'), f('st'), f('sy')

fig, (pa, pb_ax, pc) = plt.subplots(1, 3, figsize=(W2, 2.82))

# Panel (a): FWHM vs v_K(R_ph)
rr, x, t, y, sx, st, sy = load('A_Curtis_BX26_sources.tsv'); k = sx < SXMAX
v = np.log10(S.VK0_KMS)+.25*x+t; ev = np.sqrt((.25*sx)**2+st**2)
pa.errorbar(v[k], y[k], xerr=ev[k], yerr=sy[k], fmt='none', ecolor=RED, elinewidth=.55, alpha=.6, zorder=3)
pa.scatter(v[k], y[k], s=20, marker='o', facecolor=RED, edgecolor='white', linewidth=.45, zorder=4, label='atmosphere fits (%d)' % k.sum())
rr, x, t, y, sx, st, sy = load('C_deGraaff_unique_sources.tsv'); own = np.array([r['width_source'] != 'BX26_dedup' for r in rr])
v = np.log10(S.VK0_KMS)+.25*x+t; ev = np.sqrt((.25*sx)**2+st**2)
pa.errorbar(v, y, xerr=ev, yerr=sy, fmt='none', ecolor=GOLD, elinewidth=.55, alpha=.6, zorder=3)
pa.scatter(v[own], y[own], s=22, marker='s', facecolor=GOLD, edgecolor='white', linewidth=.45, zorder=5, label='blackbody, our widths (%d)' % own.sum())
pa.scatter(v[~own], y[~own], s=18, marker='s', facecolor='white', edgecolor=GOLD, linewidth=.95, zorder=5, label='blackbody, B&X widths (%d)' % (~own).sum())
xx = np.array([2.55, 3.25])
for m, ls in ((1, (0, (1, 1.6))), (2, (0, (4, 2))), (3, '-')):
    pa.plot(xx, xx+math.log10(m), color=INK, lw=.85, ls=ls, zorder=1)
for m, xm in ((1, 3.12), (2, 3.12), (3, 3.01)):
    pa.text(xm, xm+math.log10(m)-.016, r'$%d\,v_{\rm K}$' % m if m > 1 else r'$v_{\rm K}$', fontsize=7.6, rotation=36, ha='center', va='top')
kms_axis(pa, 430, 1500, (500, 700, 1000, 1400), axis='x'); kms_axis(pa, 740, 5600, (1000, 1500, 2000, 3000, 5000))
pa.set_xlabel(r'$v_{\rm K}(R_{\rm ph})$ for $\lambda=1$  (km s$^{-1}$)'); pa.set_ylabel(r'broad H$\alpha$ FWHM  (km s$^{-1}$)')
pa.text(.04, .96, '(a)', transform=pa.transAxes, va='top', fontsize=9, fontweight='normal')
pa.legend(loc='upper left', bbox_to_anchor=(0.01, 0.90), handlelength=0.9, borderaxespad=.2, fontsize=6.8, handletextpad=.35, labelspacing=.26)

# Panel (b): joint posterior of (b, c) in FWHM ~ L_ph^b T_eff^c
LAB = dict(A=(.265, -.05), B=(.09, .95), C=(.44, -.12))
for tag, fn, col in (('A', 'A_Curtis_BX26_sources.tsv', RED), ('B', 'B_Curtis_frozen_sources.tsv', BLUE), ('C', 'C_deGraaff_unique_sources.tsv', GOLD)):
    rr, x, t, y, sx, st, sy = load(fn); k = sx < SXMAX
    g1, g2, s, p = S.post3(x[k], t[k], y[k], sx[k], st[k], sy[k]); pb = p.sum(2); l68, l95 = S.contour_levels(pb)
    pb_ax.contour(g1, g2, pb.T, levels=[l95], colors=[col], linewidths=.85, linestyles=[(0, (4, 2))])
    pb_ax.contour(g1, g2, pb.T, levels=[l68], colors=[col], linewidths=1.25)
    pb_ax.text(*LAB[tag], tag, color=col, fontsize=8.5, fontweight='bold', ha='center', va='center', zorder=7)
pb_ax.plot(.25, 1, marker='*', color=INK, ms=9.5, mec='white', mew=.55, zorder=6)
pb_ax.annotate('photosphere', (.25, 1), xytext=(.34, 1.62), fontsize=7.5, arrowprops=dict(arrowstyle='-', lw=.5, color=INK), ha='left')
FS = json.load(open(PHOTOSPHERE / 'simulation_results.json'))['results']
pts = []
for key in ('Curtis_A_own', 'Curtis_B_own'):
    h = FS[key]['scenarios']['measured']['methods']['flat']['Hline_bc']
    pts.append((h['b_distribution']['values'][2], h['c_distribution']['values'][2]))
    pb_ax.plot(*pts[-1], marker='P', color=INK, ms=6.0, mec='white', mew=.55, zorder=6)
    pb_ax.annotate(key[7], pts[-1], xytext=(5.0, 1), textcoords='offset points', fontsize=6.3, color=INK, va='center', zorder=7)
for pt in pts:
    pb_ax.annotate('', pt, xytext=(-.02, -1.25), arrowprops=dict(arrowstyle='-', lw=.45, color=INK))
pb_ax.text(-.04, -1.32, r'$L_{\rm H\alpha}$ only', fontsize=7.5, ha='left', va='top')
pb_ax.set_xlim(-.18, .72); pb_ax.set_ylim(-2.9, 2.35)
pb_ax.set_xlabel(r'exponent of $L_{\rm ph}$'); pb_ax.set_ylabel(r'exponent of $T_{\rm eff}$')
pb_ax.text(.04, .96, '(b)', transform=pb_ax.transAxes, va='top', fontsize=9)
pb_ax.legend(handles=[Line2D([], [], color=RED, lw=1.25, label='A: atm., B&X (18)'),
                      Line2D([], [], color=BLUE, lw=1.25, label='B: atm., ours (15)'),
                      Line2D([], [], color=GOLD, lw=1.25, label='C: blackbody (9)')],
             loc='lower right', handlelength=1.2, borderaxespad=.3, fontsize=6.9, labelspacing=.26)

# Panel (c): partial residuals at fixed dust-free L_Halpha (task 11a datasets.json)
F_ds = json.load(open(PHOTOSPHERE / 'samples.json'))
for ds_key, col, mk, lab_txt, z_ord in (('Curtis_A_own', RED, 'o', 'A: B&X widths (13)', 4),
                                        ('Curtis_B_own', BLUE, '^', 'B: our widths (15)', 5)):
    ds_rows = F_ds[ds_key]
    xl = np.array([float(r['line_x']) for r in ds_rows])
    xp = np.array([float(r['x']) for r in ds_rows])
    yy = np.array([float(r['y']) for r in ds_rows])
    sx_p = np.array([float(r['sx']) for r in ds_rows])
    sy_p = np.array([float(r['sy']) for r in ds_rows])
    rx = xp - np.polyval(np.polyfit(xl, xp, 1), xl)
    ry = yy - np.polyval(np.polyfit(xl, yy, 1), xl)
    pc.errorbar(rx, ry, xerr=sx_p, yerr=sy_p, fmt='none', ecolor=col, elinewidth=.55, alpha=.5, zorder=z_ord-1)
    pc.scatter(rx, ry, s=20 if mk == 'o' else 24, marker=mk, facecolor=col, edgecolor='white', linewidth=.45, zorder=z_ord, label=lab_txt)
xx_c = np.array([-.56, .34])
pc.plot(xx_c, .25*xx_c, color=INK, lw=.85, ls=(0, (4, 2)), zorder=1)
pc.axhline(0, color=INK, lw=.7, ls=(0, (1, 1.6)), zorder=1)
pc.text(-.50, -.132, 'slope 0.25', fontsize=7.5, rotation=30, ha='left', va='top', color=INK)
pc.text(-.53, -.008, r'$L_{\rm H\alpha}$ only', fontsize=7.5, ha='left', va='top', color=INK)
pc.set_xlim(-.56, .34); pc.set_ylim(-.21, .21)
pc.set_xlabel(r'$\Delta\log L_{\rm ph}$ at fixed $L_{\rm H\alpha}$ (dex)')
pc.set_ylabel(r'$\Delta\log$ FWHM at fixed $L_{\rm H\alpha}$ (dex)')
pc.text(.04, .96, '(c)', transform=pc.transAxes, va='top', fontsize=9)
pc.legend(loc='lower right', handlelength=1.0, borderaxespad=.3, fontsize=6.9, handletextpad=.35, labelspacing=.26)

fig.subplots_adjust(left=.075, right=.992, bottom=.165, top=.97, wspace=.35)
save(fig, 'Figure4')
print('figures written to', OUT)
