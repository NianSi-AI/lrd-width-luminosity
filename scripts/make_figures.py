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
# A2744-QSO1 with the narrow-line FWHM fixed at 150 km/s (narrow_degeneracy_refit.json)
q1 = [x_ for x_ in json.load(open(RESULTS / 'narrow_degeneracy_refit.json'))['critical_10_plus_QSO1'] if x_['source_id'] == 'A2744-QSO1'][0]['refit_W150_sensitivity']
xq2 = q1['log10_LHa_delensed_erg_s']; yq2 = math.log10(q1['FWHM_kms']); eq2 = q1['FWHM_err_cond_kms']/q1['FWHM_kms']/math.log(10)
a1.plot([xq, xq2], [yq, yq2], color=INK, lw=.6, ls=(0, (1.5, 1.5)), zorder=5)
a1.errorbar(xq2, yq2, xerr=exq, yerr=eq2, fmt='none', ecolor=INK, elinewidth=.7, zorder=6)
a1.scatter(xq2, yq2, s=30, marker='D', facecolor=INK, edgecolor='white', linewidth=.6, zorder=7)
a1.annotate('A2744-QSO1', (xq, yq), xytext=(6, 7), textcoords='offset points', fontsize=8, color=INK)
# label the three most luminous LRDs, which set the slope
_off = {'RUBIES-EGS-42046': (-3, 6, 'bottom'), 'RUBIES-EGS-55604': (-7, -2, 'top'), 'RUBIES-EGS-49140': (-6, 0, 'center')}
for t_ in sorted([t for t in G if t['cls'] == 'LRD'], key=lambda t: -t['x'])[:3]:
    dx_, dy_, va_ = _off.get(t_['source_id'], (-4, 5, 'bottom'))
    a1.annotate(t_['source_id'].replace('RUBIES-', ''), (t_['x']+43, t_['y']), xytext=(dx_, dy_), textcoords='offset points',
                fontsize=6.3, color=RED, ha='right', va=va_)
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
a2.set_xlim(43.55, 45.35); kms_axis(a2, 640, 5200, (1000, 1500, 2000, 3000, 4000, 5000))
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
rc = json.load(open(RESULTS / 'review_checks.json'))
RV = {v_['key']: v_ for v_ in rc['variants']}
S_CONS = rc['variant_extremes']['max_p95_onesided']
sd = json.load(open(DATA / 'sdss' / 'regression_results.json'))['regressions']['single_gaussian_QC']['sigma_dex']
fig, (ax, bx) = plt.subplots(1, 2, figsize=(W2, 3.9))
state = dict(y=0, yt=[], yl=[])
def head(a_, t, st):
    a_.text(-.03, st['y'], t, transform=a_.get_yaxis_transform(), ha='right', va='center', fontsize=8.2, fontweight='bold'); st['y'] += 1
def row(a_, st, lab, qL, qO=None, up=None):
    y_ = st['y']
    for q, col, mk, dy in ((qL, RED, 'o', -.17 if qO is not None else 0.), (qO, BLUE, '^', .17)):
        if q is None: continue
        a_.plot([q[0], q[4]], [y_+dy]*2, color=col, lw=.9, solid_capstyle='butt')
        a_.plot([q[1], q[3]], [y_+dy]*2, color=col, lw=2.6, solid_capstyle='butt')
        a_.plot(q[2], y_+dy, marker=mk, color=col, mec='white', mew=.6, ms=5.0 if mk == 'o' else 5.4, zorder=5)
    if up is not None:
        a_.plot(up, y_, marker='<', color=RED, mfc='white', mew=.8, ms=4.2, zorder=6)
    st['yt'].append(y_); st['yl'].append(lab); st['y'] += 1
def finish(a_, st, xmax, fs=8.0):
    a_.set_yticks(st['yt']); a_.set_yticklabels(st['yl'], fontsize=fs); a_.set_ylim(st['y']-.35, -.7); a_.set_xlim(0, xmax)
    a_.set_xlabel(r'intrinsic scatter $\sigma_{\rm int}$ (dex)'); a_.tick_params(axis='y', which='both', left=False, right=False)
    a_.axvspan(sd[0], sd[4], color=INK, alpha=.16, lw=0, zorder=0); a_.axvline(sd[2], color=INK, lw=.8, zorder=1)
# (a) width definitions
head(ax, 'This work', state)
for lab, k in rows: row(ax, state, lab, r1[k]['LRD']['sigma'], r1[k]['other']['sigma'])
state['y'] += .35; head(ax, 'Baccus & Xu', state)
row(ax, state, 'single-Gaussian FWHM', r2['bx_main']['LRD']['sigma'], r2['bx_main']['other']['sigma'])
finish(ax, state, .47, 8.3)
ax.text(sd[2]+.008, 0, 'SDSS quasars', va='center', ha='left', fontsize=8.0, color=INK)
ax.text(.97, .975, '(a)', transform=ax.transAxes, ha='right', va='top', fontsize=9)
ax.legend(handles=[Line2D([], [], marker='o', color=RED, mec='white', ms=5.0, lw=.9, label='LRDs'),
                   Line2D([], [], marker='^', color=BLUE, mec='white', ms=5.4, lw=.9, label='other broad-line sources')],
          loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=2, handlelength=1.8, borderaxespad=.2, columnspacing=1.0)
# (b) analysis variants of the single-Gaussian FWHM (LRDs only)
st2 = dict(y=0, yt=[], yl=[])
groups = (('Fit', ('baseline', 'lsf_half', 'win6800_noSII', 'win8000_maskSII', 'lowest_snr', 'swap42232', 'g395m_only')),
          ('Narrow line', ('reg_frozen_universal_W150_29LRD_15oth', 'reg_frozen_universal_W30_29LRD_15oth',
                           'reg_frozen_OIII_anchored_29LRD_15oth', 'oiii_qc', 'narrow_lt600')),
          ('Sample', ('no_top3', 'parent_sn3', 'any_label')),
          ('Repeat spectra', ('floor_same', 'floor_all', 're_add', 're_mult')))
for gi, (gname, keys) in enumerate(groups):
    if gi: st2['y'] += .35
    head(bx, gname, st2)
    for k in keys:
        v_ = RV[k]
        lbl = 'EGS-42232 on G395M' if k == 'swap42232' else v_['label']
        row(bx, st2, lbl + (f" ($N$={v_['N']})" if v_['N'] != 29 else ''), v_['sigma_q'], None, v_['sigma_p95_onesided'])
finish(bx, st2, .175, 7.4)
bx.axvline(S_CONS, color=RED, lw=.8, ls=(0, (1, 1.5)), zorder=1)
bx.text(S_CONS-.003, st2['y']-1.0, 'conservative ceiling', fontsize=7.0, color=RED, va='bottom', ha='right', rotation=90)
bx.text(.97, .975, '(b)', transform=bx.transAxes, ha='right', va='top', fontsize=9)
bx.legend(handles=[Line2D([], [], marker='<', color=RED, mfc='white', ls='none', ms=4.2, label='one-sided 95% limit')],
          loc='lower right', bbox_to_anchor=(1.0, 1.0), handlelength=1.0, borderaxespad=.2)
fig.subplots_adjust(left=.205, right=.99, top=.94, bottom=.11, wspace=.98); save(fig, 'Figure2')

# ------------------------------------------------------------------ Figure 3: orbital scatter budget (unified 95% ceiling)
rep_cal = json.load(open(RESULTS / 'repeat_calibration.json'))
S0 = r1['W1_single_gaussian']['LRD']['sigma'][2]
S1 = rep_cal['posterior_with_floors']['formal_zero_floor']['LRD_p95_onesided']  # exact one-sided 95% upper limit: 0.07474 dex
ori = r2['orient']
ori_edge60 = ori['thin_60_90']
fig, a = plt.subplots(figsize=(W1, 3.1)); th = np.linspace(0, np.pi/2, 300)
def ell(s, sr=0.):
    e = math.sqrt(s*s-sr*sr/4); return 2*e*np.cos(th), e*np.sin(th)
x1, y1 = ell(S1); x0, y0 = ell(S0); x2, y2 = ell(S_CONS)
a.fill_between(x1[::-1], 0, y1[::-1], color=RED, alpha=.07, lw=0)
a.fill_between(x0[::-1], 0, y0[::-1], color=RED, alpha=.15, lw=0)
a.plot(x0, y0, color=RED, lw=1.5); a.plot(x1, y1, color=RED, lw=1.1, ls=(0, (4, 2))); a.plot(x2, y2, color=RED, lw=1.0, ls=(0, (1, 1.5)))
lines_spec = [
    (r'thin disk, $i<45^\circ$', ori['thin_0_45'],  0.251, 'right', +0.0026, 'bottom', '-'),
    (r'thin disk, $i>45^\circ$', ori['thin_45_90'], 0.251, 'right', +0.0014, 'bottom', '-'),
    (r'thin disk, $i>60^\circ$', ori_edge60,        0.251, 'right', -0.0022, 'top',    '-'),
    (r'$H/R=0.5\ (i<60^\circ)$', ori['h0.5_0_60'],  0.132, 'left',  +0.0026, 'bottom', (0, (5, 1.5, 1, 1.5))),
    (r'$H/R=0.7$',               ori['h0.7_0_60'],  0.106, 'left',  +0.0022, 'bottom', (0, (5, 1.5, 1, 1.5))),
    (r'$H/R=1$',                 ori['h1_0_60'],    0.146, 'left',  -0.0022, 'top',    (0, (5, 1.5, 1, 1.5))),
]
for txt, v, xpos, ha_, dy_, va_, ls in lines_spec:
    a.axhline(v, color=BLUE, lw=.8, ls=ls, zorder=2)
    a.text(xpos, v + dy_, txt, color=BLUE, fontsize=7.9, ha=ha_, va=va_, zorder=4)
a.text(.012, .025, 'allowed', color=RED, fontsize=9, va='center')
a.legend(handles=[Line2D([], [], color=RED, lw=1.5, label=r'$\sigma_{\rm int}=%.3f$ dex' % S0),
                  Line2D([], [], color=RED, lw=1.1, ls=(0, (4, 2)), label=r'$%.3f$ dex (95%% upper limit)' % S1),
                  Line2D([], [], color=RED, lw=1.0, ls=(0, (1, 1.5)), label=r'$%.2f$ dex (conservative)' % S_CONS)],
         loc='upper left', bbox_to_anchor=(0, .84), handlelength=2.4, borderaxespad=.7)
a.set_xlim(0, .255); a.set_ylim(0, .235)
a.set_xlabel(r'scatter of $\log(M/R)$ at fixed luminosity (dex)'); a.set_ylabel(r'orientation scatter $\sigma_k$ of $\log$ FWHM (dex)')
fig.subplots_adjust(left=.16, right=.975, bottom=.14, top=.975); save(fig, 'Figure3')

# ------------------------------------------------------------------ Figure 4: 3-panel photospheric comparison & scaling tests
SXMAX = .3
MUTE4 = '#555555'
K_ISO = (3.0 / (8.0 * math.log(2.0))) ** (-0.5)  # 1.3596 (f_vir = 0.54, isotropic 3D Gaussian)
K_RING = 2.0                                     # 2.0000 (f_vir = 0.25, edge-on circular ring)
_PHOT_MAP = {
    'A_Curtis_BX26_sources.tsv': PHOTOSPHERE / 'atmosphere_bx26_widths.tsv',
    'B_Curtis_frozen_sources.tsv': PHOTOSPHERE / 'atmosphere_our_widths.tsv',
    'C_deGraaff_unique_sources.tsv': PHOTOSPHERE / 'blackbody_sources.tsv',
}
def load(fn):
    rr = list(csv.DictReader(open(_PHOT_MAP[fn]), delimiter='\t')); f = lambda k: np.array([float(r[k]) for r in rr])
    return rr, f('x'), f('t'), f('y'), f('sx'), f('st'), f('sy')

fig, (pa, pb_ax, pc) = plt.subplots(1, 3, figsize=(W2, 2.86))
fig.subplots_adjust(left=0.075, right=0.992, bottom=0.175, top=0.97, wspace=0.35)

# Panel (a): FWHM vs measured v_K(R_ph) at R = R_ph, lambda = 1
rr_a, x_a, t_a, y_a, sx_a, st_a, sy_a = load('A_Curtis_BX26_sources.tsv'); ka = sx_a < SXMAX
va = np.log10(S.VK0_KMS) + .25*x_a + t_a; eva = np.sqrt((.25*sx_a)**2 + st_a**2)
pa.errorbar(va[ka], y_a[ka], xerr=eva[ka], yerr=sy_a[ka], fmt='none', ecolor=RED, elinewidth=.5, alpha=.48, zorder=3)
pa.scatter(va[ka], y_a[ka], s=19, marker='o', facecolor=RED, edgecolor='white', linewidth=.45, zorder=4, label='A: atm., B&X (%d)' % ka.sum())

rr_b, x_b, t_b, y_b, sx_b, st_b, sy_b = load('B_Curtis_frozen_sources.tsv'); kb = sx_b < SXMAX
vb = np.log10(S.VK0_KMS) + .25*x_b + t_b; evb = np.sqrt((.25*sx_b)**2 + st_b**2)
pa.errorbar(vb[kb], y_b[kb], xerr=evb[kb], yerr=sy_b[kb], fmt='none', ecolor=BLUE, elinewidth=.5, alpha=.48, zorder=3)
pa.scatter(vb[kb], y_b[kb], s=22, marker='^', facecolor=BLUE, edgecolor='white', linewidth=.45, zorder=5, label='B: atm., ours (%d)' % kb.sum())

rr_c, x_c, t_c, y_c, sx_c, st_c, sy_c = load('C_deGraaff_unique_sources.tsv')
own = np.array([r['width_source'] != 'BX26_dedup' for r in rr_c])
vc = np.log10(S.VK0_KMS) + .25*x_c + t_c; evc = np.sqrt((.25*sx_c)**2 + st_c**2)
pa.errorbar(vc, y_c, xerr=evc, yerr=sy_c, fmt='none', ecolor=GOLD, elinewidth=.5, alpha=.52, zorder=3)
pa.scatter(vc[own], y_c[own], s=20, marker='s', facecolor=GOLD, edgecolor='white', linewidth=.45, zorder=5, label='C: BB, ours (%d)' % own.sum())
pa.scatter(vc[~own], y_c[~own], s=17, marker='s', facecolor='white', edgecolor=GOLD, linewidth=.9, zorder=5, label='C: BB, B&X (%d)' % (~own).sum())

kms_axis(pa, 430, 1500, (500, 700, 1000, 1400), axis='x')
kms_axis(pa, 740, 5800, (1000, 1500, 2000, 3000, 5000), axis='y')

xx = np.array([2.55, 3.25])
ref_lines = [
    (1.0,    (0, (1.1, 1.6)),           INK,   0.85),
    (K_ISO,  (0, (3.8, 1.5, 1.0, 1.5)), MUTE4, 0.80),
    (K_RING, (0, (4.2, 2.0)),           INK,   0.85),
    (3.0,    '-',                       INK,   0.85),
]
for mult, ls, col, lw in ref_lines:
    pa.plot(xx, xx + math.log10(mult), color=col, lw=lw, ls=ls, zorder=1)

p0 = pa.transData.transform((2.7, 2.7))
p1 = pa.transData.transform((3.0, 3.0))
rot_deg = math.degrees(math.atan2(p1[1] - p0[1], p1[0] - p0[0]))

line_labels = [
    (1.0,    3.085, r'$v_{\rm K}\ (f_{\rm vir}{\simeq}1)$',     INK,   -0.015, 'top'),
    (K_ISO,  2.955, r'$1.36\,v_{\rm K}\ (f_{\rm vir}{=}0.54)$', MUTE4, +0.013, 'bottom'),
    (K_RING, 3.070, r'$2\,v_{\rm K}\ (f_{\rm vir}{=}0.25)$',    INK,   -0.015, 'top'),
    (3.0,    2.965, r'$3\,v_{\rm K}$',                          INK,   -0.016, 'top'),
]
for mult, xm, txt, col, dy, va_ in line_labels:
    pa.text(xm, xm + math.log10(mult) + dy, txt, fontsize=6.55, color=col,
            rotation=rot_deg, rotation_mode='anchor', ha='center', va=va_, zorder=2)

pa.set_xlabel(r'$v_{\rm K}(R_{\rm ph})$ from measured $L_{\rm ph},\,T_{\rm eff}$ ($\mathrm{km\,s^{-1}}$)', fontsize=8.8)
pa.set_ylabel(r'broad H$\alpha$ FWHM  ($\mathrm{km\,s^{-1}}$)')
pa.text(0.04, 0.965, '(a)', transform=pa.transAxes, va='top', fontsize=9)
pa.text(0.965, 0.965, r'$R=R_{\rm ph},\ \lambda=1$', transform=pa.transAxes, va='top', ha='right', fontsize=7.2, color=MUTE4)
pa.legend(loc='upper left', bbox_to_anchor=(0.01, 0.905), handlelength=0.9, borderaxespad=0.18,
          fontsize=6.45, handletextpad=0.32, labelspacing=0.22)

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
pb_ax.set_xlabel(r'exponent of $L_{\rm ph}$', fontsize=9.3); pb_ax.set_ylabel(r'exponent of $T_{\rm eff}$', fontsize=9.3)
pb_ax.text(.04, .965, '(b)', transform=pb_ax.transAxes, va='top', fontsize=9)
pb_ax.legend(handles=[Line2D([], [], color=RED, lw=1.25, label='A: atm., B&X (18)'),
                      Line2D([], [], color=BLUE, lw=1.25, label='B: atm., ours (15)'),
                      Line2D([], [], color=GOLD, lw=1.25, label='C: blackbody (9)')],
             loc='lower right', handlelength=1.2, borderaxespad=.3, fontsize=6.8, labelspacing=.25)

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
    i29 = [r['temperature_id'] for r in ds_rows].index('RUBIES-UDS29813')
    if ds_key == 'Curtis_B_own':
        pc.annotate('UDS-29813', (rx[i29], ry[i29]), xytext=(4, -9), textcoords='offset points', fontsize=6.6, color=INK)
xx_c = np.array([-.56, .34])
pc.plot(xx_c, .25*xx_c, color=INK, lw=.85, ls=(0, (4, 2)), zorder=1)
PS = json.load(open(RESULTS / 'photosphere_simulations.json'))
for k_, col in (('Curtis_A_own', RED), ('Curtis_B_own', BLUE)):
    sm_ = PS[k_]['sim']['hierarchical_baseline']['Hph_p16_50_84'][1]
    pc.plot(xx_c, sm_*xx_c, color=col, lw=.75, ls=(0, (5, 1.5, 1, 1.5)), zorder=1, alpha=.85)
pc.axhline(0, color=INK, lw=.7, ls=(0, (1, 1.6)), zorder=1)
pc.text(-.50, -.132, 'slope 0.25', fontsize=7.5, rotation=30, ha='left', va='top', color=INK)
pc.text(-.53, -.008, r'$L_{\rm H\alpha}$ only', fontsize=7.5, ha='left', va='top', color=INK)
pc.set_xlim(-.56, .34); pc.set_ylim(-.21, .21)
pc.set_xlabel(r'$\Delta\log L_{\rm ph}$ at fixed $L_{\rm H\alpha}$ (dex)', fontsize=9.3)
pc.set_ylabel(r'$\Delta\log$ FWHM at fixed $L_{\rm H\alpha}$ (dex)', fontsize=9.3)
pc.text(.04, .965, '(c)', transform=pc.transAxes, va='top', fontsize=9)
pc.legend(loc='lower right', handlelength=1.0, borderaxespad=.3, fontsize=6.8, handletextpad=.35, labelspacing=.25)

save(fig, 'Figure4')
print('figures written to', OUT)
