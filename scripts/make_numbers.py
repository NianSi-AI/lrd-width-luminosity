"""Compute every number quoted in the manuscript and tables, writing paper/numbers.tex
and data/results/numbers.json.

Inputs are the stored results in data/ (see data/README.md).

Usage: python scripts/make_numbers.py
"""
import csv
import json
import math
import os
import re
import statistics
from collections import Counter

from lrdwidth import stats as S
from lrdwidth.paths import DATA, PAPER, PHOTOSPHERE, RESULTS, WIDTHS

HERE = str(PAPER)
r1_raw = json.load(open(RESULTS / "width_luminosity.json"))
r1 = {("W1_" + k[6:] if k.startswith("width_") else k): v for k, v in r1_raw.items()}
r1["W6"] = json.load(open(RESULTS / "sdss_catalog_benchmark.json"))
r2 = json.load(open(RESULTS / "robustness.json"))
r6 = json.load(open(RESULTS / "photosphere_constrained.json"))
r11 = json.load(open(RESULTS / "photosphere_simulations.json"))
rep_cal = json.load(open(RESULTS / "repeat_calibration.json"))
t13_1 = json.load(open(RESULTS / "parent_covariance_audit.json"))
t13_2 = json.load(open(RESULTS / "narrow_degeneracy_refit.json"))
t13_2b = json.load(open(RESULTS / "frozen_fixed_narrow_audit.json"))
t13_3 = json.load(open(RESULTS / "parent49_injection.json"))
t13_4 = json.load(open(RESULTS / "repeat_random_effects.json"))
t13_5 = json.load(open(RESULTS / "exterior_potential_check.json"))
X = json.load(open(RESULTS / "measurement_summary.json"))["numbers"]
FS = json.load(open(PHOTOSPHERE / "simulation_results.json"))["results"]
xv = lambda k: X[k]["value"]

def _post_rows(rows):
    x = [r['x'] for r in rows]
    y = [r['y'] for r in rows]
    cxx = [r['cov'][0][0] for r in rows]
    cyy = [r['cov'][1][1] for r in rows]
    cxy = [r['cov'][0][1] for r in rows]
    b, s, p, ah = S.post2(x, y, cxx, cyy, cxy)
    pb = p.sum(1)
    ps = p.sum(0)
    return dict(b_q=S.quant(b, pb), s_q=S.quant(s, ps), s_p95=S.quant(s, ps, qs=(0.95,))[0], s=s, ps=ps)

M = []   # (name, value, status, note)
def add(name, val, st, note=''):
    assert name.isalpha(), name
    assert name not in [m[0] for m in M], name
    M.append((name, str(val), st, note))

f3 = lambda v: f'{v:.3f}'
f2 = lambda v: f'{v:.2f}'
f1 = lambda v: f'{v:.1f}'
def neg(s): return r'\ensuremath{-}'+s[1:] if s.startswith('-') else s
def trip(prefix, q, fmt, st, note=''):
    add(prefix, fmt(q[2]), st, note)
    add(prefix+'lo', fmt(q[0]), st)
    add(prefix+'hi', fmt(q[4]), st)
sec = lambda t: M.append(('%', t, '', ''))

_w = {0: 'none', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine'}

sec('sample (design and counts)')
for n, v in (('NLRD', 29), ('NOther', 15), ('NFrozen', 44), ('NSources', 115), ('NSpectra', 204), ('NSpecHa', 177),
             ('NProg', 'eight'), ('NEGS', 30), ('NUDS', 14), ('NAbsLRD', 'five'), ('NAbsOther', 'four'),
             ('NLabTLRD', 13), ('NLabKLRD', 10), ('NLabHLRD', 'six')):
    add(n, v, 'V')
assert r2['N_abs'] == [5, 4]
for n, v in (('cutSN', 5), ('cutSNthree', 3), ('cutBIC', 6), ('cutChi', 2), ('HvidingCut', 1000),
             ('wnarDef', 150), ('passShift', '0.10'), ('rOutRph', 100)):
    add(n, v, 'D')
add('matchRadius', r'0\farcs5', 'D')
for n, v in (('zLRDmin', '3.36'), ('zLRDmax', '7.15'), ('zLRDmed', '5.45'), ('SNLRDmed', '13.7'),
             ('SNOthmed', '8.2'), ('zOthmin', '3.60'), ('zOthmax', '6.04')):
    add(n, v, 'V')
d = r2['frozen_desc']
add('LLRDmin', f2(d['LRD']['L'][0]), 'R')
add('LLRDmed', f2(d['LRD']['L'][1]), 'R')
add('LLRDmax', f2(d['LRD']['L'][2]), 'R')
add('WLRDmin', round(d['LRD']['W'][0]), 'R')
add('WLRDmax', round(d['LRD']['W'][2]), 'R')
add('WLRDmed', round(d['LRD']['W'][1]), 'R')

sec('line profiles (BIC)')
b = 'W1_BIC_counts__LRD_all__'
add('NBICwing', xv(b+'best_counts__lorentzian')+xv(b+'best_counts__exponential'), 'X', 'LRDs whose best of four profiles is Lorentzian or exponential')
add('NBICsix', xv(b+'Gaussian_deltaBIC_gt6'), 'X', 'LRDs: Gaussian worse than best by Delta BIC > 6')
add('NBIChiSN', 'nine', 'X', 'LRDs with broad S/N > 20, all nine with Delta BIC > 6 (delivery 10, W1 table)')
add('NBICsixOth', 'four', 'X', 'other sources: Gaussian worse than best by Delta BIC > 6 (of 15)')
add('NCoreWing', 'seven', 'V', 'sources where the two-Gaussian decomposition is identifiable')

sec('our sample, single-Gaussian FWHM, narrow-line sensitivity, and repeat-spectrum calibration')
m = r1['W1_single_gaussian']
fm = r2['frozen_main']
trip('Fsig', m['LRD']['sigma'], f3, 'R', 'formal-error posterior median and equal-tailed 95% interval; X=0.0535 [0.0372, 0.0792]')
s_up = rep_cal['posterior_with_floors']['formal_zero_floor']['LRD_p95_onesided']
add('FsigUp', f3(s_up), 'C', f'one-sided 95% upper limit on LRD sigma_int under formal errors (95th percentile: {s_up:.5f} dex)')
add('FsigUpTwo', f2(s_up), 'A', 'two-decimal formal-error 95% upper limit')
trip('Fb', m['LRD']['slope'], f2, 'R', 'X=0.1688 [0.1170, 0.2225]')
add('Fa', f3(fm['LRD']['a']), 'R')
add('FWnorm', int(round(10**fm['LRD']['a'], -1)), 'A', '10^Fa')
trip('Fsigo', m['other']['sigma'], f3, 'R', 'X=0.1535 [0.1032, 0.2520]')
trip('Fdsig', m['diff'], lambda v: neg(f3(v)), 'R', 'X=-0.0999 [-0.1995, -0.0437]')
add('Frms', f3(m['LRD']['ols_sd']), 'R')
add('Frmso', f3(m['other']['ols_sd']), 'R')
n3 = r2['frozen_noTop3']
add('NNoTop', n3['LRD']['N'], 'R')
trip('FbNoTop', n3['LRD']['slope'], f2, 'R')
trip('FsigNoTop', n3['LRD']['sigma'], f3, 'R')
add('NBright', 'three', 'D')
add('LBright', '43.5', 'R')
add('LRestMax', '42.9', 'R')
na = r2['frozen_noAbs']
add('NNoAbs', na['LRD']['N'], 'R', 'LRDs without a fitted Balmer-absorption component')
trip('FsigNoAbs', na['LRD']['sigma'], f3, 'R')
trip('FbNoAbs', na['LRD']['slope'], f2, 'R')
trip('FdsigNoAbs', na['diff'], lambda v: neg(f3(v)), 'R', 'LRDs without Balmer absorption minus 15 other sources')
nw = r2['frozen_noWide2']
trip('FsigoDrop', nw['other']['sigma'], f3, 'R')
trip('FdsigDrop', nw['diff'], lambda v: neg(f3(v)), 'R')
add('WOutA', nw['dropped'][0][1], 'R')
add('WOutB', nw['dropped'][1][1], 'R')
lh = r2['frozen_LSFhalf']
trip('FsigLsf', lh['LRD']['sigma'], f3, 'R', 'X=0.0684 [0.0482, 0.0992]')
trip('FsigoLsf', lh['other']['sigma'], f3, 'R')
trip('FdsigLsf', lh['diff'], lambda v: neg(f3(v)), 'R', 'X=-0.0884 [-0.1908, -0.0286]')
w3_in = json.load(open(WIDTHS / 'lsf_half_inputs.json'))
lsf_post = _post_rows([r_ for r_ in w3_in if r_['cls'] == 'LRD'])
assert abs(lsf_post['s_q'][2] - lh['LRD']['sigma'][2]) < 1e-6
add('FsigLsfUp', f3(lsf_post['s_p95']), 'C', f"one-sided 95% upper limit for LSF width x 0.5 ({lsf_post['s_p95']:.5f} dex)")

f44_aud = t13_2b['frozen_44_audit']
bound_lrd = [x for x in f44_aud if x['cls'] == 'LRD' and x['orig_W_narrow_kms'] > 599.0]
bound_oth = [x for x in f44_aud if x['cls'] != 'LRD' and x['orig_W_narrow_kms'] > 599.0]
assert len(bound_lrd) == 2 and len(bound_oth) == 0
w_bef_mean = sum(x['orig_FWHM_kms'] for x in bound_lrd) / len(bound_lrd)
w_aft_mean = sum(0.5 * (x['W150_FWHM_kms'] + x['new_FWHM_kms']) for x in bound_lrd) / len(bound_lrd)
w_bef_50 = int(50 * round(w_bef_mean / 50.0))
w_aft_50 = int(50 * round(w_aft_mean / 50.0))
assert w_bef_50 == 3150 and w_aft_50 == 2400
add('NNarBoundLRD', _w[len(bound_lrd)], 'C', 'luminous LRDs whose free narrow Gaussian hits the 600 km/s bound')
add('WNarBoundBefore', w_bef_50, 'A', f'mean broad FWHM of the two narrow-bound LRDs before locking ({w_bef_mean:.1f} km/s, rounded to 50 km/s)')
add('WNarBoundAfter', w_aft_50, 'A', f'mean broad FWHM of the two narrow-bound LRDs after locking narrow FWHM ({w_aft_mean:.1f} km/s, rounded to 50 km/s)')

w1_in = json.load(open(WIDTHS / 'regression_inputs.json'))
drop_nar_ids = {x['source_id'] for x in bound_lrd}
lrd27_rows = [r_ for r_ in w1_in['single_gaussian'] if r_['cls'] == 'LRD' and r_['source_id'] not in drop_nar_ids]
oth15_rows = [r_ for r_ in w1_in['single_gaussian'] if r_['cls'] != 'LRD']
assert len(lrd27_rows) == 27 and len(oth15_rows) == 15
nar_cut_lrd = _post_rows(lrd27_rows)
nar_cut_oth = _post_rows(oth15_rows)
nar_cut_diff = S.diffq(nar_cut_lrd['s'], nar_cut_lrd['ps'], nar_cut_oth['ps'])
add('NLRDNarCut', len(lrd27_rows), 'C', 'LRDs with fitted narrow-line FWHM < 600 km/s')
trip('FbNarCut', nar_cut_lrd['b_q'], f2, 'C', 'slope of N=27 LRDs with narrow-line FWHM < 600 km/s')
trip('FsigNarCut', nar_cut_lrd['s_q'], f3, 'C', 'sigma_int of N=27 LRDs with narrow-line FWHM < 600 km/s')
add('FsigNarCutUp', f3(nar_cut_lrd['s_p95']), 'C', f"one-sided 95% upper limit for narrow-line FWHM < 600 km/s ({nar_cut_lrd['s_p95']:.5f} dex)")
trip('FdsigNarCut', nar_cut_diff, lambda v: neg(f3(v)), 'C', 'Delta sigma_int for N=27 narrow-line < 600 km/s vs 15 other sources')

reg_w150 = t13_2b['reg_frozen_universal_W150_29LRD_15oth']
reg_oiii = t13_2b['reg_frozen_OIII_anchored_29LRD_15oth']
reg_w30 = t13_2b['reg_frozen_universal_W30_29LRD_15oth']
trip('FbNarFix', reg_w150['slope_q'], f2, 'C', 'slope with narrow-line FWHM fixed to 150 km/s across all 29 LRDs')
trip('FsigNarFix', reg_w150['sigma_q'], f3, 'C', 'LRD sigma_int with narrow-line FWHM fixed to 150 km/s')
add('FsigNarFixUp', f3(reg_w150['sigma_p95_onesided']), 'C', f"one-sided 95% upper limit for fixed 150 km/s ({reg_w150['sigma_p95_onesided']:.5f} dex)")
trip('FsigoNarFix', reg_w150['other']['sigma_q'], f3, 'C', 'other-source sigma_int with narrow-line FWHM fixed to 150 km/s')
trip('FdsigNarFix', reg_w150['diff_LRD_minus_other_q'], lambda v: neg(f3(v)), 'C')

trip('FbNarOIII', reg_oiii['slope_q'], f2, 'C', 'slope with narrow-line FWHM anchored to [O III] 5007 (150 km/s where unmeasured)')
trip('FsigNarOIII', reg_oiii['sigma_q'], f3, 'C', 'LRD sigma_int with narrow-line FWHM anchored to [O III] 5007')
s_oiii_up = reg_oiii['sigma_p95_onesided']
add('FsigNarOIIIUp', f3(s_oiii_up), 'C', f'one-sided 95% upper limit for [O III]-anchored narrow FWHM ({s_oiii_up:.5f} dex)')
trip('FsigoNarOIII', reg_oiii['other']['sigma_q'], f3, 'C', 'other-source sigma_int with narrow-line FWHM anchored to [O III] 5007')
trip('FdsigNarOIII', reg_oiii['diff_LRD_minus_other_q'], lambda v: neg(f3(v)), 'C')
add('FsigNarMin', f3(reg_w30['sigma_q'][2]), 'C', 'LRD sigma_int median when fixing all 29 narrow FWHMs to 30 km/s')
add('FsigNarMinUp', f3(reg_w30['sigma_p95_onesided']), 'C', f"one-sided 95% upper limit when fixing all 29 narrow FWHMs to 30 km/s ({reg_w30['sigma_p95_onesided']:.5f} dex)")
n_oiii_30 = sum(1 for x in f44_aud if x['cls'] == 'LRD' and x.get('OIII_W_kms') is not None and x['OIII_W_kms'] <= 30.01)
assert n_oiii_30 == 21
add('NOIIIUnres', n_oiii_30, 'C', 'frozen LRDs whose [O III] 5007 fit hits the 30 km/s lower bound')

rc = json.load(open(RESULTS / 'review_checks.json'))
rcv = {v_['key']: v_ for v_ in rc['variants']}
s_cons = rc['variant_extremes']['max_p95_onesided']
assert rc['variant_extremes']['max_p95_key'] == 'oiii_qc'
add('FsigCons', f2(s_cons), 'C', f'conservative ceiling: largest one-sided 95% limit across fit, narrow-line and sample variants ({s_cons:.5f} dex, [O III] anchoring with QC re-applied)')
loosen_pct = 5 * round(20.0 * (s_cons / s_up - 1.0))
assert loosen_pct == 35
add('narCeilLoosenPct', loosen_pct, 'A', f'percentage loosening from s_up to the conservative ceiling (exact: {100*(s_cons/s_up-1):.1f}%, rounded to 5%)')
oq = rc['oiii_anchored_qc_consistent']
add('FsigOIIIqc', f3(oq['sigma_med']), 'C', '[O III]-anchored, absorption-bound rule re-applied (N=27; approximate, see review_checks.json)')
add('FsigOIIIqcUp', f3(oq['sigma_p95_onesided']), 'C', f"one-sided 95% limit for the same ({oq['sigma_p95_onesided']:.5f} dex)")
add('NOIIIqc', oq['N'], 'C')
trip('FsigWinWide', rcv['win8000_maskSII']['sigma_q'], f3, 'C', '+-8000 km/s window with [S II] masked (window_and_wing_checks.json)')
add('winWider', 8000, 'D')
nbc = rc['narrow_bound_counts']
assert nbc['free_narrow_at_30'] == 17 and nbc['free_narrow_at_600'] == 2 and nbc['OIII_measured'] == 23 and nbc['OIII_at_30'] == 21
add('NNarThirty', nbc['free_narrow_at_30'], 'C', 'LRDs whose free narrow FWHM sits at the 30 km/s lower bound')
add('NOIIIMeas', nbc['OIII_measured'], 'C', 'LRDs with an [O III] 5007 width measurement')
add('dWNarRms', f2(nbc['rms_dlogFWHM_OIII_minus_base']), 'C', f"rms Delta log FWHM ([O III]-anchored minus baseline) over 29 LRDs ({nbc['rms_dlogFWHM_OIII_minus_base']:.4f})")
add('dWNarMax', f2(nbc['max_abs_dlogFWHM']), 'C', f"max |Delta log FWHM| ({nbc['max_abs_source']})")
ps3 = rc['parent_sn3_added']
add('NFailSNthreeWeak', _w[ps3['N_dBIC_lt6']], 'C', 'of the 9 added S/N>=3 parent LRDs, number with Delta BIC_null < 6')
cw = rc['constant_width']['LRD26']
trip('FsigConst', cw['sigma_q'], f3, 'C', 'b = 0 (constant-width) model for the 26 LRDs below the bright end')
add('WConstMed', int(round(cw['FWHM_median_kms'], -1)), 'C', 'median FWHM of the 26 LRDs (km/s)')
add('LConstMin', f'{cw["logL_range"][0]:.1f}', 'C'); add('LConstMax', f'{cw["logL_range"][1]:.1f}', 'C')
phl = rc['photosphere_leverage']
for p_, k_ in (('JA', 'Curtis_A_own'), ('JB', 'Curtis_B_own')):
    assert phl[k_]['RUBIES-UDS-29813_is_most_negative']
    trip('bp'+p_+'drop', phl[k_]['drop_RUBIES-UDS-29813']['b_ph_q'], lambda v: neg(f2(v)), 'C', 'partial L_ph slope without RUBIES-UDS-29813')
add('dLphOut', f1(-0.5 * (phl['Curtis_A_own']['RUBIES-UDS-29813_dlogLph_at_fixed_LHa'] + phl['Curtis_B_own']['RUBIES-UDS-29813_dlogLph_at_fixed_LHa'])), 'C', 'RUBIES-UDS-29813 offset in log L_ph at fixed L_Ha (dex, sign dropped)')
ei = rc['error_inflation_summary']
add('dsigInflLo', neg(f2(ei['diff_med_min'])), 'C', 'most negative median Delta sigma_int under error inflation/floors (ours and BX26)')
add('dsigInflHi', neg(f2(ei['diff_med_max'])), 'C', 'least negative median Delta sigma_int under error inflation/floors')
assert ei['diff_upper95_max'] < 0
for n_, v_ in (('inflLo', '1.5'), ('inflHi', '2.8'), ('inflFlrLo', '0.044'), ('inflFlrHi', '0.08')):
    add(n_, v_, 'D')
rd = rc['repeat_diagnostics']
assert rd['G395H_minus_G395M_objects'] == 7 and rd['G395H_narrower_objects'] == 7
add('NGHobj', 'seven', 'C', 'LRDs with both G395H and G395M passing spectra')
add('dWGH', neg(f2(rd['G395H_minus_G395M_object_median'])), 'C', 'median per-object log FWHM(G395H) - log FWHM(G395M)')
add('offGH', f3(rd['free_offset_model']['offset_mean_dex']), 'C', 'free G395H offset in the random-effects model (posterior mean)')
add('offGHerr', f3(rd['free_offset_model']['offset_sd_dex']), 'C')
add('dWPrim', neg(f3(rd['primary_minus_secondary_mean'])), 'C', 'mean log FWHM(primary) - log FWHM(secondary) over the 14 primary pairs')
add('NGM', rd['G395M_only']['N'], 'C')
trip('FsigGM', rd['G395M_only']['sigma_q'], f3, 'C', 'LRDs with G395M primary spectra only')
add('FsigGMUp', f3(rd['G395M_only']['sigma_p95_onesided']), 'C')
trip('FsigSwap', rd['swap42232']['sigma_q'], f3, 'C', 'RUBIES-EGS-42232 primary replaced by its highest-S/N G395M spectrum')
add('FsigSwapUp', f3(rd['swap42232']['sigma_p95_onesided']), 'C')
add('FsigLowSNUp', f3(rcv['lowest_snr']['sigma_p95_onesided']), 'C', 'one-sided 95% limit choosing the lowest-S/N passing spectrum')
sc_ = rc['sdss_no_chi2_rule']
trip('SsgNoChi', sc_['no_chi2']['sigma_q'], f3, 'C', 'refitted SDSS quasars with the reduced chi^2 rule removed')
add('NSsgChi', sc_['N_added'], 'C')
add('WSdssSNfail', int(round(sc_['SNR_failed_FWHM_median_kms'], -1)), 'C', 'median single-Gaussian FWHM of SDSS refits failing S/N >= 5 (km/s)')
add('NSdssCap', _w[sc_['N_SNR_failed_at_12000']], 'C')
ss_ = rc['shell_slopes_tau2to3']
add('shellSlopeLo', f2(ss_['core_convolved'][0]), 'C', 'lowest predicted shell slope at tau_e = 2-3 (core-convolved, R ~ L^1/2)')
add('shellSlopeHi', f2(ss_['pure_wing'][1]), 'C', 'highest predicted pure-wing shell slope at tau_e = 2-3 (fixed R)')
add('shellPureLo', f2(ss_['pure_wing'][0]), 'C', 'lowest predicted pure-wing shell slope at tau_e = 2-3 (R ~ L^1/2)')
add('shellConvHi', f2(ss_['core_convolved'][1]), 'C', 'highest predicted core-convolved shell slope at tau_e = 2-3 (fixed R)')
cp_ = rc['ceers746_skyfire299']
add('sepPair', f'{cp_["sep_arcsec"]:.0f}', 'C', 'separation of CEERS-746 and Skyfire-299 (arcsec)')
qe = rc['qso1_eq3_check']
add('QSOdBICfree', neg(f1(qe['fixed_narrow_dBIC_vs_free_null'])), 'C', 'QSO1 fixed-narrow broad fit vs free-narrow no-broad model, Delta BIC')
re_mult = t13_4['random_effects_models']['multiplicative_fixed_primary_x']
add('reMultUp', f3(re_mult['all_43_spectra']['sigma_p95_onesided']), 'C', f"one-sided 95% limit, random effects with unmodeled variance f^2 C_yy ({re_mult['all_43_spectra']['sigma_p95_onesided']:.5f})")
add('reMultScale', f1(re_mult['all_43_spectra']['f_med']), 'C', 'posterior median error-scale factor f')
add('reAddFloor', f3(t13_4['random_effects_models']['exact_2dcov_fixed_primary_x']['all_43_spectra']['f_med']), 'C', 'posterior median additive floor f (dex)')

assert t13_1['N_parent_clean_total'] == 50 and t13_1['N_parent_valid_inband'] == 49
r38 = t13_1['regressions_full_2x2_covariance']['row4_plus_valid_SNge3_38']
assert r38['N'] == 38
add('NParentClean', t13_1['N_parent_clean_total'], 'C', 'consistently labeled parent LRDs with a grating Ha fit')
add('NParentInBand', t13_1['N_parent_valid_inband'], 'C', 'in-band consistently labeled parent LRDs')
add('NParentFail', t13_1['N_parent_valid_inband'] - t13_1['N_frozen_LRD'], 'A', 'in-band parent LRDs failing frozen QC (49 - 29 = 20)')
add('NFailSNthree', _w[r38['N'] - t13_1['N_frozen_LRD']], 'C', 'in-band failed parent LRDs with broad S/N >= 3 (9)')
add('NLRDSNthree', r38['N'], 'C', '29 frozen + 9 in-band failed S/N >= 3 parent LRDs')
trip('FbSNthree', r38['slope_q'], f2, 'C', 'slope of N=38 LRDs including in-band S/N >= 3 failed spectra')
trip('FsigSNthree', r38['sigma_q'], f3, 'C', 'sigma_int of N=38 LRDs including in-band S/N >= 3 failed spectra')
add('FsigSNthreeUp', f3(r38['sigma_p95_onesided']), 'C', f"one-sided 95% upper limit for N=38 ({r38['sigma_p95_onesided']:.5f} dex)")
uds_pos = [x for x in t13_2['critical_10_plus_QSO1'] if x['is_in_10_critical'] and 'QSO1' not in x['source_id'] and x['orig_SNR_broad'] >= 3.0 and x['refit_adopted']['residual_Eq3_dex'] > 0.08]
assert len(uds_pos) == 3 and all('UDS' in x['source_id'] for x in uds_pos)
add('NFailPosUDS', _w[len(uds_pos)], 'C', 'failed S/N>=3 UDS LRDs sitting above Eq. 3 even with fixed narrow FWHM')
add('resFailPosLo', f2(min(x['refit_adopted']['residual_Eq3_dex'] for x in uds_pos)), 'C')
add('resFailPosHi', f2(max(x['refit_adopted']['residual_Eq3_dex'] for x in uds_pos)), 'C')
add('failPosSNlo', f1(min(x['refit_adopted']['SNR_broad'] for x in uds_pos)), 'C')
add('failPosSNhi', f1(max(x['refit_adopted']['SNR_broad'] for x in uds_pos)), 'C')

dup_s = rep_cal['duplicate_spectra']
fl_s = rep_cal['posterior_with_floors']['same_grating_sy_floor']
fl_a = rep_cal['posterior_with_floors']['all_pairs_sy_floor']
add('NPairsBest', dup_s['width_best_vs_rest_14']['N_pairs'], 'C', 'primary-vs-secondary pairs among QC-passing spectra of the 9 multi-spectrum LRDs')
add('NPairsAll', dup_s['width_all_pairs_22']['N_pairs'], 'C', 'all pairwise combinations among QC-passing spectra of the 9 multi-spectrum LRDs')
add('NPairsSameG', 'seven', 'C', 'same-grating pairs among the 22 (%d)' % dup_s['width_same_grating_7']['N_pairs'])
add('NPairsCrossG', dup_s['width_cross_grating_15']['N_pairs'], 'C', 'cross-grating pairs among the 22')
add('rmsDupBest', f2(dup_s['width_best_vs_rest_14']['norm_diff_rms']), 'C', 'rms of normalized Delta log FWHM across 14 primary-vs-secondary pairs (%.4f)' % dup_s['width_best_vs_rest_14']['norm_diff_rms'])
add('rmsDupAll', f2(dup_s['width_all_pairs_22']['norm_diff_rms']), 'C', 'rms of normalized Delta log FWHM across all 22 pairs (%.4f)' % dup_s['width_all_pairs_22']['norm_diff_rms'])
add('rmsDupSameG', f2(dup_s['width_same_grating_7']['norm_diff_rms']), 'C', 'rms of normalized Delta log FWHM across 7 same-grating pairs (%.4f)' % dup_s['width_same_grating_7']['norm_diff_rms'])
add('flrDupBest', f3(dup_s['width_best_vs_rest_14']['extra_floor_per_spec_dex']), 'C', 'extra per-spectrum width error floor for 14 primary-vs-secondary pairs (%.5f dex)' % dup_s['width_best_vs_rest_14']['extra_floor_per_spec_dex'])
add('flrDupAll', f3(dup_s['width_all_pairs_22']['extra_floor_per_spec_dex']), 'C', 'extra per-spectrum width error floor for all 22 pairs (%.5f dex)' % dup_s['width_all_pairs_22']['extra_floor_per_spec_dex'])
add('flrDupSameG', f3(dup_s['width_same_grating_7']['extra_floor_per_spec_dex']), 'C', 'extra per-spectrum width error floor for 7 same-grating pairs (%.5f dex)' % dup_s['width_same_grating_7']['extra_floor_per_spec_dex'])
trip('FsigFlrSame', fl_s['LRD_sigma_q'], f3, 'C', 'LRD sigma_int with +0.044 dex same-grating floor in C_yy')
add('FsigFlrSameUp', f3(fl_s['LRD_p95_onesided']), 'C', f"one-sided 95% upper limit with +0.044 dex floor ({fl_s['LRD_p95_onesided']:.5f} dex)")
trip('FsigoFlrSame', fl_s['other_sigma_q'], f3, 'C')
trip('FdsigFlrSame', fl_s['diff_q'], lambda v: neg(f3(v)), 'C')
trip('FsigFlrAll', fl_a['LRD_sigma_q'], f3, 'C', 'LRD sigma_int with +0.060 dex all-pairs floor in C_yy')
add('FsigFlrAllUp', f3(fl_a['LRD_p95_onesided']), 'C', f"one-sided 95% upper limit with +0.060 dex floor ({fl_a['LRD_p95_onesided']:.5f} dex)")
trip('FsigoFlrAll', fl_a['other_sigma_q'], f3, 'C')
trip('FdsigFlrAll', fl_a['diff_q'], lambda v: neg(f3(v)), 'C')
add('FdsigFlrTwo', neg(f2(fl_s['diff_q'][2])), 'A', 'median Delta sigma_int with repeat-spectrum floor (-0.12 dex)')
add('NSNfif', rep_cal['high_snr_LRDs_ge15']['N'], 'C', 'LRDs with broad-line S/N >= 15')
add('rmsSNfif', f3(rep_cal['high_snr_LRDs_ge15']['raw_ols_sd_ddof2_dex']), 'C', 'raw residual SD (ddof=2) of the 12 S/N>=15 LRDs (%.5f dex; ddof=1 is %.5f)' % (rep_cal['high_snr_LRDs_ge15']['raw_ols_sd_ddof2_dex'], rep_cal['high_snr_LRDs_ge15']['raw_ols_sd_ddof1_dex']))
add('NPairsBigObj', t13_4['pairwise_slicings_6']['6_drop_RUBIES_EGS_42232_12']['dropped_N_pairs'], 'C', 'pairs among the 22 contributed by RUBIES-EGS-42232 (5 passing spectra)')
slice_ups = [
    t13_4['pairwise_slicings_6']['4_top2_one_pair_per_source_9']['sigma_p95_onesided'],
    t13_4['pairwise_slicings_6']['5_equal_object_weight_9']['sigma_p95_onesided'],
] + [v_['sigma_p95_onesided'] for v_ in t13_4['pairwise_leave_one_object_out_9'].values()]
add('flrSliceUpLo', f3(min(slice_ups)), 'C', f'min one-sided 95% upper limit across one-pair/equal-weight/LOO floors ({min(slice_ups):.5f} dex)')
add('flrSliceUpHi', f3(max(slice_ups)), 'C', f'max one-sided 95% upper limit across one-pair/equal-weight/LOO floors ({max(slice_ups):.5f} dex)')
re_mod = t13_4['random_effects_models']['exact_2dcov_fixed_primary_x']
add('NSpecPassLRD', re_mod['all_43_spectra']['N_spectra'], 'C', 'total QC-passing spectra across the 29 LRDs')
add('reAllUp', f3(re_mod['all_43_spectra']['sigma_p95_onesided']), 'C', f"one-sided 95% upper limit in 43-spectrum joint random-effects model ({re_mod['all_43_spectra']['sigma_p95_onesided']:.5f} dex)")
add('reDropBigUp', f3(re_mod['42232_primary_only_39_spectra']['sigma_p95_onesided']), 'C', f"one-sided 95% upper limit with RUBIES-EGS-42232 reduced to primary spectrum ({re_mod['42232_primary_only_39_spectra']['sigma_p95_onesided']:.5f} dex)")

sec('six width definitions')
defs = (('Wf', 'nonparametric_W50'), ('Ex', 'exponential'), ('Lo', 'lorentzian'), ('Np', 'nonparametric_FWHM'), ('Tg', 'double_gaussian'))
for p, k in defs:
    q = r1['W1_'+k]
    trip('Dsig'+p, q['LRD']['sigma'], f3, 'R')
    trip('Dsigo'+p, q['other']['sigma'], f3, 'R')
    trip('Ddsig'+p, q['diff'], lambda v: neg(f3(v)), 'R')
    add('Db'+p, f2(q['LRD']['slope'][2]), 'R')
five = [r1['W1_'+k]['LRD']['sigma'][2] for k in ('single_gaussian', 'nonparametric_W50', 'exponential', 'lorentzian', 'nonparametric_FWHM')]
add('sigFiveLo', f3(min(five)), 'A')
add('sigFiveHi', f3(max(five)), 'A')

sec('injection-recovery (29 surviving templates and 49 full in-band parent templates)')
g = lambda sc, inj, meas, what: xv(f'W2_final_summary_s{sc}__{sc}_{inj}_M4_Gaussian_QC_{meas}__sigma_median__{what}')
add('simInj', '0.05', 'D')
add('simNull', '0.15', 'D')
add('NReal', 100, 'D')
add('recGG', f3(g(3, 'single_gaussian', 'single_gaussian', 'p16_p50_p84__p50')), 'X')
add('recEG', f3(g(3, 'exponential', 'single_gaussian', 'p16_p50_p84__p50')), 'X')
add('recGW', f3(g(3, 'single_gaussian', 'nonparametric_W50', 'p16_p50_p84__p50')), 'X')
add('recEW', f3(g(3, 'exponential', 'nonparametric_W50', 'p16_p50_p84__p50')), 'X')
add('recGL', f2(g(3, 'single_gaussian', 'lorentzian', 'p16_p50_p84__p50')), 'X')
add('recGE', f2(g(3, 'single_gaussian', 'exponential', 'p16_p50_p84__p50')), 'X')
add('nullGG', f3(g(1, 'single_gaussian', 'single_gaussian', 'p16_p50_p84__p50')), 'X')
add('nullEG', f3(g(1, 'exponential', 'single_gaussian', 'p16_p50_p84__p50')), 'X')
hits = sum(g(1, i, 'single_gaussian', 'count_le_observed') for i in ('single_gaussian', 'exponential'))
den = sum(g(1, i, 'single_gaussian', 'denominator') for i in ('single_gaussian', 'exponential'))
add('nullHits', 'none' if hits == 0 else int(hits), 'X')
add('nullDen', int(den), 'X')
add('simA', '0.10', 'D')
add('simB', '0.08', 'D')
st2 = json.load(open(WIDTHS / 'sloped_injection_recovery.json'))
st2a = st2['null_scenarios']
st2b = st2['sloped_simulations']
add('recEightG', f3(st2b['gauss_rel_s008']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitEightG', _w[st2b['gauss_rel_s008']['bayes_hits_le_obs_Fsig']], 'C')
add('recEightE', f3(st2b['exp_matched_s008']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitEightE', _w[st2b['exp_matched_s008']['bayes_hits_le_obs_Fsig']], 'C')
add('recTenG', f3(st2b['gauss_rel_s010']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitTenG', _w[st2b['gauss_rel_s010']['bayes_hits_le_obs_Fsig']], 'C')
add('recTenE', f3(st2b['exp_matched_s010']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitTenE', _w[st2b['exp_matched_s010']['bayes_hits_le_obs_Fsig']], 'C')
tot10 = st2b['gauss_rel_s010']['bayes_hits_le_obs_Fsig'] + st2b['exp_matched_s010']['bayes_hits_le_obs_Fsig']
add('hitTenTot', _w[tot10], 'A')
add('wideInjShift', f2(st2b['exp_wide_s010']['injected_a0'] - st2b['exp_matched_s010']['injected_a0']), 'A')
add('recTenWide', f3(st2b['exp_wide_s010']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitTenWide', _w[st2b['exp_wide_s010']['bayes_hits_le_obs_Fsig']], 'C')
add('wideSurvShift', f2(st2b['exp_wide_s010']['recovered_a_hat_p16_p50_p84'][1] - fm['LRD']['a']), 'A')

p49 = t13_3['primary_parent49_suites']
def _trip16_84(prefix, q, fmt, note=''):
    add(prefix, fmt(q[2]), 'C', note)
    add(prefix+'lo', fmt(q[1]), 'C')
    add(prefix+'hi', fmt(q[3]), 'C')

for tag, key in (('FiveG', 'gauss_rel_s005'), ('FiveE', 'exp_matched_s005'),
                 ('TenG', 'gauss_rel_s010'), ('TenE', 'exp_matched_s010')):
    sub = p49[key]['parent49']
    _trip16_84('recPar'+tag, sub['bayes_sigma_med_quantiles'], f3, f'49-template {key} recovered sigma_int median (p16--p84)')
    _trip16_84('bPar'+tag, sub['bayes_slope_med_quantiles'], f2, f'49-template {key} recovered slope median (p16--p84)')
    cnt_ = int(sub['bayes_hits_le_obs_Fsig'])
    add('hitPar'+tag+'num', f'{cnt_}/{int(sub["N_replicates"])}', 'C')
    if tag in ('TenG', 'TenE'):
        add('hitPar'+tag, _w[cnt_], 'C')

p_rat_g = t13_3['paired_offset_pass_rate_grids']['offset_grid_gauss']['ratio_plus010_over_zero_all49']
p_rat_e = t13_3['paired_offset_pass_rate_grids']['offset_grid_exp']['ratio_plus010_over_zero_all49']
add('passRatLo', f2(min(p_rat_g, p_rat_e)), 'C', 'min pass-fraction ratio (+0.10 vs 0.00 dex) across 49 parent templates')
add('passRatHi', f2(max(p_rat_g, p_rat_e)), 'C', 'max pass-fraction ratio (+0.10 vs 0.00 dex) across 49 parent templates')
ng05 = p49['gauss_rel_s005']['parent49']['N_selected_min_p16_p50_p84_max']
ng10 = p49['gauss_rel_s010']['parent49']['N_selected_min_p16_p50_p84_max']
ne05 = p49['exp_matched_s005']['parent49']['N_selected_min_p16_p50_p84_max']
ne10 = p49['exp_matched_s010']['parent49']['N_selected_min_p16_p50_p84_max']
assert ng05[2] == ng10[2] and ne05[2] == ne10[2]
add('NParSelG', int(round(ng10[2])), 'C')
add('NParSelGlo', int(round(min(ng05[1], ng10[1]))), 'C')
add('NParSelGhi', int(round(max(ng05[3], ng10[3]))), 'C')
add('NParSelE', int(round(ne10[2])), 'C')
add('NParSelElo', int(round(min(ne05[1], ne10[1]))), 'C')
add('NParSelEhi', int(round(max(ne05[3], ne10[3]))), 'C')
own_s = t13_3['sensitivity_ownL_parent49_suites']
add('hitParOwnTenGnum', f"{int(own_s['gauss_rel_s010_ownL']['parent49']['bayes_hits_le_obs_Fsig'])}/100", 'C')
add('hitParOwnTenEnum', f"{int(own_s['exp_matched_s010_ownL']['parent49']['bayes_hits_le_obs_Fsig'])}/100", 'C')

sec('Baccus & Xu table')
anc_chk = json.load(open(RESULTS / 'ancillary_checks.json'))
assert anc_chk['repeatability_and_bx26']['W4_candidates_total'] == 13 and anc_chk['repeatability_and_bx26']['W4_candidates_LRD'] == 3
for n, v in (('NX', 227), ('NXLRD', 48), ('NXCtrl', 97), ('NXOvRec', 29), ('NXOvSrc', 25), ('NXLabDiff', 'six'),
             ('XWcut', 1000), ('XSNcut', 5), ('XkLRD', '5.4'), ('XkOth', '9.8'), ('XzCut', 3),
             ('XLmin', '43.7'), ('XLmax', '45.2'), ('NXHbCand', 13), ('NXHbLRD', 'three')):
    add(n, v, 'V' if n not in ('NXHbCand', 'NXHbLRD') else 'C')
add('NXLRDuniq', rep_cal['bx26_checks']['N_unique_sky_LRD'], 'C', 'unique sky objects (<1 arcsec) among the 48 BX26 LRD records')
trip('XsigDedup', rep_cal['bx26_checks']['dedup_LRD_sigma_q'], f3, 'C', 'BX26 LRD sigma_int deduplicated to 45 unique sky objects')
add('XovSD', f3(rep_cal['bx26_checks']['shared_29_width_ratio_sd_dex']), 'C', 'SD of log10 width ratio across the 29 shared records (%.4f dex)' % rep_cal['bx26_checks']['shared_29_width_ratio_sd_dex'])
add('XovLRDSD', f3(rep_cal['bx26_checks']['shared_16_LRD_width_ratio_sd_dex']), 'C', 'SD of log10 width ratio across the 16 records labeled LRD in both tables')
bm = r2['bx_main']
assert bm['LRD']['N'] == 48 and bm['other']['N'] == 97
trip('Xsig', bm['LRD']['sigma'], f3, 'R', 'X=0.0578 [0.0446, 0.0770]')
trip('Xb', bm['LRD']['slope'], f2, 'R')
trip('Xsigo', bm['other']['sigma'], f3, 'R', 'X=0.1436 [0.1213, 0.1716]')
trip('Xdsig', bm['diff'], lambda v: neg(f3(v)), 'R', 'X=-0.0859 [-0.1172, -0.0570]')
bn = r2['bx_noOverlap']
add('NXNoOv', bn['LRD']['N'], 'R')
add('NXNoOvCtrl', bn['other']['N'], 'R')
trip('XsigNoOv', bn['LRD']['sigma'], f3, 'R')
trip('XsigoNoOv', bn['other']['sigma'], f3, 'R')
trip('XdsigNoOv', bn['diff'], lambda v: neg(f3(v)), 'R')
for n, v in (('NXPairs', 36), ('XdMatch', r'\ensuremath{-}0.065'), ('XdMatchlo', r'\ensuremath{-}0.105'),
             ('XdMatchhi', r'\ensuremath{-}0.024'), ('NPermHits', 'eight'), ('NPermTrials', '10^{6}'),
             ('NXLeave', 87), ('matchDz', '0.5'), ('matchDL', '0.2')):
    add(n, v, 'V')
add('XresM', f2(r2['bx_Mvir_LRD']['resid_sd']), 'R', 'residual SD of BX26 log M about log Lbol, LRDs')
add('XresMo', f2(r2['bx_Mvir_ctrl']['resid_sd']), 'R')

sec('SDSS DR7 quasars')
w6 = r1['W6']
add('NSDSS', w6['N'], 'X')
trip('Ssig', w6['sigma'], f3, 'R', 'X=0.2073 [0.2024, 0.2123]')
add('Sb', f2(w6['slope'][2]), 'R')
add('zSDSSmin', f2(r1['W6']['z_range'][0]), 'C')
add('zSDSSmax', f2(r1['W6']['z_range'][1]), 'C')
add('SsigTwo', f2(w6['sigma'][2]), 'A')
add('FsigoTwo', f2(m['other']['sigma'][2]), 'A')
add('sigTwoLo', f2(min(bm['other']['sigma'][2], m['other']['sigma'][2])), 'A')
add('sigTwoHi', f2(max(bm['other']['sigma'][2], m['other']['sigma'][2])), 'A')
add('sdssWoff', f2(rep_cal['sdss_224_offset']['dlogFWHM_SDSS_minus_LRD_dex']), 'C', 'Delta log10 FWHM (SDSS refitted minus LRD) at median LRD L_Ha (%.4f dex)' % rep_cal['sdss_224_offset']['dlogFWHM_SDSS_minus_LRD_dex'])
add('WsdssMedL', int(round(rep_cal['sdss_224_offset']['SDSS_FWHM_at_LRD_med_kms'], -1)), 'C', 'SDSS refitted single-Gaussian FWHM relation at median LRD L_Ha, km/s')
add('WlrdMedL', int(round(rep_cal['sdss_224_offset']['LRD_FWHM_at_LRD_med_kms'], -1)), 'C', 'LRD single-Gaussian FWHM relation at median LRD L_Ha, km/s')

sec('QSO1 diagnostic')
for n, v in (('QSOsn', '3.4'), ('QSOW', 2107), ('QSOWerr', 302), ('QSOmu', '7.5'), ('QSOmuErr', '1.5'),
             ('QSOresCor', '0.17'), ('QSOresErr', '0.08')):
    add(n, v, 'V')
add('QSOcompNar', 560, 'L', "D'Eugenio et al. (2026) / Juodzbalis et al. (2026) intermediate broad component FWHM (km/s)")
add('QSOcompBrd', 2000, 'L', "D'Eugenio et al. (2026) / Juodzbalis et al. (2026) broad wing component FWHM (km/s)")
qso_t2 = [x for x in t13_2['critical_10_plus_QSO1'] if x['source_id'] == 'A2744-QSO1'][0]
qso_w150 = qso_t2['refit_W150_sensitivity']
qso_sig_y = qso_w150['FWHM_err_cond_kms'] / (qso_w150['FWHM_kms'] * math.log(10.0))
qso_sig_mu = m['LRD']['slope'][2] * (1.5 / (7.5 * math.log(10.0)))
add('QSOWfix', round(qso_w150['FWHM_kms']), 'C', 'QSO1 broad Gaussian FWHM (km/s) when narrow FWHM is fixed to 150 km/s')
add('QSOsnFix', f1(qso_w150['SNR_broad']), 'C', 'QSO1 broad-line S/N when narrow FWHM is fixed to 150 km/s')
add('QSOresFix', neg(f2(qso_w150['residual_Eq3_delensed_dex'])), 'C', 'QSO1 delensed residual (dex) from Eq. 3 with fixed 150 km/s narrow line')
add('QSOresFixErr', f2(math.hypot(qso_sig_y, qso_sig_mu)), 'A', 'quadrature uncertainty on QSOresFix including mu=7.5+/-1.5')
add('FsigNarOIIIqso', f3(t13_2b['reg_frozen_plus_QSO1_w150_30LRD']['sigma_q'][2]), 'C', 'sigma_int of 29 [O III]-anchored LRDs plus QSO1 (N=30)')
add('NLRDqso', t13_2b['reg_frozen_plus_QSO1_w150_30LRD']['N'], 'C')

sec('photosphere: sources with sigma(log L_ph) < 0.3 dex')
add('vKnorm', 1493, 'R', 'km/s at 1e45 erg/s, 5000 K, lambda=1, kappa=0.398; X=1492.63')
add('kappaES', '0.398', 'D')
A, B, C = r6['A'], r6['B'], r6['C']
E = r11['main_E']; Ao = r11['Curtis_A_own']; Bo = r11['Curtis_B_own']
cx = {'PA': 'X: ratio 2.9927, HPD 0.8310, b 0.2075 [0.0993, 0.3146], sigma 0.0732 [0.0483, 0.1168]',
      'PB': 'X: ratio 2.6126, HPD 0.9042, b 0.1452 [0.0251, 0.2610], sigma 0.0812 [0.0507, 0.1370]',
      'PC': 'X: ratio 2.7475, HPD 0.9846, b 0.4007 [0.2128, 0.5662], sigma 0.1146 [0.0669, 0.2310]'}
for p, r in (('PA', A), ('PB', B), ('PC', C)):
    add('N'+p, 'nine' if r['N'] == 9 else r['N'], 'C')
    add('rat'+p, f2(r['ratio_median']), 'C', cx[p])
    trip('b'+p, r['b'], f2, 'C')
    trip('c'+p, r['c'], lambda v: neg(f2(v)), 'C')
    add('hpd'+p, f'{100*r["hpd"]:.0f}', 'C')
    trip('sq'+p, r['sigma_fixed'], f3, 'C')
add('NPdrop', 'five', 'C')
assert len(A['dropped']) == 5 and len(B['dropped']) == 3
add('NPAall', 23, 'X'); add('NPBall', 18, 'X'); add('sxCut', '0.3', 'D')
add('sxDropMin', '0.4', 'C'); add('sxDropMax', '2.0', 'C'); add('sxKeepMax', '0.18', 'C')
add('hierDb', '0.01', 'X')
assert E['N'] == 5 and Ao['N'] == 13 and Bo['N'] == 15
add('NPE', 'five', 'C')
add('ratPE', f2(E['ratio_p16_50_84'][1]), 'C')
trip('sqPE', E['sigma_fixed_photosphere'], f3, 'C')
add('NPCbx', 'four', 'X')
add('NdGall', 146, 'X'); add('NdGbb', 47, 'C'); add('NdGfit', 14, 'C')
rats = [A['ratio_median'], B['ratio_median'], E['ratio_p16_50_84'][1]]
rl, rh = min(rats), max(rats)
add('ratLo', f'{rl:.1f}', 'A'); add('ratHi', f'{rh:.1f}', 'A')
add('sqLo', f2(min(A['sigma_fixed'][2], B['sigma_fixed'][2])), 'A')
add('sqHi', f2(max(A['sigma_fixed'][2], B['sigma_fixed'][2])), 'A')
add('sigLamRLo', f2(2*min(A['sigma_fixed'][2], B['sigma_fixed'][2])), 'A')
add('sigLamRHi', f2(2*max(A['sigma_fixed'][2], B['sigma_fixed'][2])), 'A')
add('TminP', int(round(min(A['T_range'][0], C['T_range'][0]), -2)), 'C')
add('TmaxP', int(round(max(A['T_range'][1], C['T_range'][1]), -2)), 'C')
add('LphMin', f'{A["logL_range"][0]:.1f}', 'C'); add('LphMax', f'{A["logL_range"][1]:.1f}', 'C')
sds_lph = [rep_cal['photosphere_log_Lph_over_LHa'][k]['sd_dex'] for k in ('Curtis_A_own', 'Curtis_B_own')]
add('sdLphHaLo', f2(min(sds_lph)), 'C', 'min SD of log10(L_ph / L_Ha) across Curtis_B_own and Curtis_A_own')
add('sdLphHaHi', f2(max(sds_lph)), 'C', 'max SD of log10(L_ph / L_Ha) across Curtis_B_own and Curtis_A_own')

hb = [FS[k]['scenarios']['measured']['methods']['flat']['Hline_bc'][q+'_distribution']['values'][2] for k in ('Curtis_A_own', 'Curtis_B_own') for q in ('b', 'c')]
add('lineB', f'{(hb[0]+hb[2])/2:.1f}', 'X')
add('lineCLo', f'{min(hb[1], hb[3]):.1f}', 'X')
add('lineCHi', f'{max(hb[1], hb[3]):.1f}', 'X')
cod = lambda k, m_: FS[k]['scenarios']['measured']['methods'][m_]
for p, r, k in (('JE', E, 'main_E'), ('JA', Ao, 'Curtis_A_own'), ('JB', Bo, 'Curtis_B_own')):
    f, h = cod(k, 'flat'), cod(k, 'hierarchical')
    sm = r['sim']['hierarchical_baseline']; so = r['sim']['own_flat_baseline']
    add('N'+p, {5: 'five', 13: 13, 15: 15}[r['N']], 'C')
    trip('bp'+p, r['partial']['b2'], lambda v: neg(f2(v)), 'C')
    fr = [sm['p_ph'], so['p_ph'], f['p_ph']['fraction'], h['p_ph']['fraction']]
    add('simPct'+p, f'{100*sum(fr)/4:.0f}', 'C')
    add('simMed'+p, (f'{sm["Hph_p16_50_84"][1]:.1f}' if p == 'JE' else f2(sm['Hph_p16_50_84'][1])), 'C')
    trip('sl'+p, r['line_only']['sigma'], f3, 'C')
    trip('so'+p, r['sigma_fixed_photosphere'], f3, 'C')
p10 = [100*FS['main_E']['scenarios']['0.10']['methods'][m_]['power']['fraction'] for m_ in ('flat', 'hierarchical')]
add('powTen', f'{5*round(sum(p10)/10):.0f}', 'X')
add('powLo', f'{5*math.floor(min(E["sim"][g]["power"] for g in E["sim"])*20):.0f}', 'A')
add('powHi', f'{5*math.ceil(max(E["sim"][g]["power"] for g in E["sim"])*20):.0f}', 'A')

add('kRing', '2', 'D')
add('kIso', '1.36', 'A', '2 sqrt(2 ln 2)/sqrt(3); X=1.359556')
add('kVir', '0.96', 'A', '1/sqrt(1.075); Reines & Volonteri epsilon; X=0.964486')
lam = lambda k, r: (k/r)**2
add('lamVirLo', f2(lam(.964486, rh)), 'A'); add('lamVirHi', f2(lam(.964486, rl)), 'A')
add('lamIsoLo', f2(lam(1.359556, rh)), 'A'); add('lamIsoHi', f2(lam(1.359556, rl)), 'A')
add('lamRingLo', f2(lam(2, rh)), 'A'); add('lamRingHi', f2(lam(2, rl)), 'A')
add('mEddLo', f'{(rl/2)**2:.1f}', 'A', '(ratio/k)^2 for k=2')
add('mEddHi', f'{(rh/.964486)**2:.0f}', 'A', '(ratio/k)^2 for k=0.96')
add('mEddStdLo', f'{(rl/1.359556)**2:.1f}', 'A', '(ratio/k)^2 min for k=1.36')
add('brdVirLo', f1(rl/.964486), 'A', 'min (FWHM/v_K)/k for k=0.96')
add('brdVirHi', f1(rh/.964486), 'A', 'max (FWHM/v_K)/k for k=0.96')
add('brdIsoLo', f1(rl/1.359556), 'A', 'min (FWHM/v_K)/k for k=1.36')
add('brdIsoHi', f1(rh/1.359556), 'A', 'max (FWHM/v_K)/k for k=1.36')
add('brdRingLo', f2(rl/2.0), 'A', 'min (FWHM/v_K)/k for k=2 (1.25)')
add('brdRingHi', f1(rh/2.0), 'A', 'max (FWHM/v_K)/k for k=2 (1.5)')

g_k3 = [e for e in t13_5['opacity_temperature_grid'] if abs(e['kappa_cm2_g'] - 1e-3) < 1e-12][0]['by_temperature']
g_k4 = [e for e in t13_5['opacity_temperature_grid'] if abs(e['kappa_cm2_g'] - 1e-4) < 1e-12][0]['by_temperature']
# exact finite-r_out expressions, including the (1 - R_ph/r_out) factor
del_k3_lo = g_k3['T_sample_max_6381K']['DeltaM_over_MEdd_finite_rout100']
del_k3_hi = g_k3['T_sample_min_3637K']['DeltaM_over_MEdd_finite_rout100']
del_k4_max = g_k4['T_sample_min_3637K']['DeltaM_over_MEdd_finite_rout100']
dv_k3_lo = g_k3['T_sample_max_6381K']['Delta_vesc_pct_at_lambda1_finite_rout100']
dv_k3_hi = g_k3['T_sample_min_3637K']['Delta_vesc_pct_at_lambda1_finite_rout100']
dv_k4_max = g_k4['T_sample_min_3637K']['Delta_vesc_pct_at_lambda1_finite_rout100']
lam_k4_max = g_k4['T_sample_min_3637K']['relaxed_lambda_ceiling_finite_rout100']
assert abs(dv_k3_lo - 100.0 * (math.sqrt(1.0 + del_k3_lo) - 1.0)) < 1e-10
assert abs(dv_k3_hi - 100.0 * (math.sqrt(1.0 + del_k3_hi) - 1.0)) < 1e-10
assert abs(dv_k4_max - 100.0 * (math.sqrt(1.0 + del_k4_max) - 1.0)) < 1e-10
assert abs(lam_k4_max - 1.0 / (1.0 - del_k4_max)) < 1e-12
add('deltaKthreeLo', f3(del_k3_lo), 'C', f'min Delta M_eff / M_Edd at kappa=1e-3 across sample T_eff ({del_k3_lo:.5f})')
add('deltaKthreeHi', f3(del_k3_hi), 'C', f'max Delta M_eff / M_Edd at kappa=1e-3 across sample T_eff ({del_k3_hi:.5f})')
add('deltaKfourMax', f2(del_k4_max), 'C', f'max Delta M_eff / M_Edd at kappa=1e-4 at T_min=3637 K ({del_k4_max:.5f})')
add('dvEscKthreeLo', f2(dv_k3_lo), 'C', f'min fractional v_esc increase at kappa=1e-3 across sample T_eff ({dv_k3_lo:.4f}%)')
add('dvEscKthreeHi', f1(dv_k3_hi), 'C', f'max fractional v_esc increase at kappa=1e-3 across sample T_eff ({dv_k3_hi:.4f}%)')
add('dvEscKfourMax', round(dv_k4_max), 'C', f'max fractional v_esc increase at kappa=1e-4 at T_min=3637 K ({dv_k4_max:.2f}%)')
add('lamEscKfourMax', f2(lam_k4_max), 'C', f'exact relaxed escape ceiling 1/(1-delta) at kappa=1e-4, T_min=3637 K ({lam_k4_max:.4f})')

sec('membership and rules')
inc = json.load(open(RESULTS / 'membership.json'))
rv_sample = json.load(open(DATA / 'sample' / 'rule_variants.json'))['variants']['chi2_lt3']['LRD']
add('sigChiThree', f3(rv_sample['posterior']['intrinsic_scatter_dex'][2]), 'X')
cx = {r_['macro']: float(r_['value']) for r_ in csv.DictReader(open(RESULTS / 'membership_full_covariance.tsv'), delimiter='\t')}
for p, k in (('IncA', 'rules_relaxed_same_labels'), ('IncB', 'any_label_all_rules'), ('IncC', 'any_label_chi2lt3'), ('IncD', 'any_label_all')):
    r = inc[k]; add('N'+p, r['N'], 'N')
    q = [cx['sig'+p+'lo'], None, cx['sig'+p], None, cx['sig'+p+'hi']]
    qb = [cx['b'+p+'lo'], None, cx['b'+p], None, cx['b'+p+'hi']]
    trip('sig'+p, q, f3, 'C')
    trip('b'+p, qb, f2, 'C')
xd = [e for e in inc['extras'] if e['id'] in set(inc['any_label_all']['added'])-set(inc['any_label_chi2lt3']['added'])]
add('chiHiA', f'{min(e["chi2"] for e in xd):.1f}', 'N'); add('chiHiB', f'{max(e["chi2"] for e in xd):.0f}', 'N')
add('NIncRules', 'two', 'N')

sec('comparison with SDSS in like measures')
tg = r1['W1_double_gaussian']['LRD']['sigma'][2]
add('DsigTgTwo', f2(tg), 'A')

sec('ratio to v_K by measure of width')
wr = json.load(open(RESULTS / 'width_definition_ratios.json'))
for p, k in (('B', 'Curtis_B_own'), ('E', 'main_E')):
    for q, mm in (('Np', 'nonparametric_FWHM'), ('Tg', 'double_gaussian'), ('Lo', 'lorentzian'), ('Ex', 'exponential')):
        add('rat'+q+p, f2(wr[k][mm]['ratio_median']), 'C')
        if p == 'B': trip('sq'+q+p, wr[k][mm]['sigma'], f3, 'N')
four = [wr['Curtis_B_own'][mm]['sigma'][2] for mm in ('single_gaussian', 'nonparametric_FWHM', 'lorentzian', 'exponential')]
add('sqMeasLo', f2(min(four)), 'A'); add('sqMeasHi', f2(max(four)), 'A')
add('ratExLo', f'{min(wr[k]["exponential"]["ratio_median"] for k in ("Curtis_B_own", "main_E")):.1f}', 'A')
add('ratWingHi', f'{max(wr[k][mm]["ratio_median"] for k in ("Curtis_B_own", "main_E") for mm in ("lorentzian", "exponential")):.1f}', 'A')
add('ratNpLo', f1(min(wr[k]['nonparametric_FWHM']['ratio_median'] for k in ('Curtis_B_own', 'main_E'))), 'A')
add('ratNpHi', f1(max(wr[k]['nonparametric_FWHM']['ratio_median'] for k in ('Curtis_B_own', 'main_E'))), 'A')
gfr = lambda rat, x: math.erfc(x/(rat/(2*math.sqrt(2*math.log(2))))/math.sqrt(2))
efr = lambda rat, x: math.exp(-x/(rat/(2*math.log(2))))
rats3 = [A['ratio_median'], B['ratio_median'], E['ratio_p16_50_84'][1]]
rex = wr['Curtis_B_own']['exponential']['ratio_median']
add('fKepLo', f'{100*gfr(min(rats3), 1):.0f}', 'A'); add('fKepHi', f'{100*gfr(max(rats3), 1):.0f}', 'A')
add('fEscLo', f'{100*gfr(min(rats3), math.sqrt(2)):.0f}', 'A'); add('fEscHi', f'{100*gfr(max(rats3), math.sqrt(2)):.0f}', 'A')
add('fKepEx', f'{100*efr(rex, 1):.0f}', 'A'); add('fEscEx', f'{100*efr(rex, math.sqrt(2)):.0f}', 'A')

sec('selection checks')
sn = json.load(open(RESULTS / 'snr_threshold.json'))
for p, k in (('SNten', 'LRD_SNge10'), ('SNtwenty', 'LRD_SNge20')):
    add('N'+p, {20: 20, 9: 'nine'}[sn[k]['N']], 'N')
    trip('sig'+p, sn[k]['sigma'], f3, 'N')
flo = fm['LRD']['a']+m['LRD']['slope'][2]*(d['LRD']['L'][0]-43)-3
add('floorDex', f2(flo), 'A')
add('floorSig', f'{flo/m["LRD"]["sigma"][2]:.1f}', 'A')
add('NTay', 13, 'C')
tp = json.load(open(RESULTS / 'photometric_labels_grid.json'))
add('FsigTay', f3(tp['medians']['FsigTay']['finer_stored']), 'C')
add('FsigTaylo', f3(tp['LRD_interval']['upstream_finer'][0]), 'C')
add('FsigTayhi', f3(tp['LRD_interval']['upstream_finer'][1]), 'C')
add('FsigoTay', f3(tp['medians']['FsigoTay']['finer_stored']), 'C')
add('kappaSolar', '0.34', 'D')
add('vKdrop', f'{100*(1-math.sqrt(0.34/0.398)):.0f}', 'A')

sec('SDSS refitted with our single Gaussian and flux fractions')
add('NSsgDraw', 300, 'X')
o_sdss0 = json.load(open(DATA / 'sdss' / 'regression_results.json'))['regressions']['single_gaussian_QC']
add('NSsg', o_sdss0['N'], 'C', 'refitted SDSS quasars passing QC (data/sdss/regression_results.json)')
trip('Ssg', o_sdss0['sigma_dex'], f3, 'C', 'refitted SDSS single-Gaussian sigma_int (data/sdss/regression_results.json)')
add('SsgTwo', f2(o_sdss0['sigma_dex'][2]), 'A')
add('SsgCat', '0.185', 'X')
add('sdssRatSame', f'{o_sdss0["sigma_dex"][2]/m["LRD"]["sigma"][2]:.1f}', 'A')
add('fKepNpLo', '38', 'X'); add('fKepNpHi', '44', 'X')
add('fEscNpLo', '23', 'X'); add('fEscNpHi', '28', 'X')
add('fEscBicLo', '24', 'X'); add('fEscBicHi', '28', 'X')
add('hierSig', '0.054', 'X'); add('hierB', '0.17', 'X')

sec('Section 5 budgets unified on the conservative formal-error 95% upper limit s_up (0.075 dex)')
s0 = m['LRD']['sigma'][2]
b0 = m['LRD']['slope'][2]
b_narfix = reg_w150['slope_q'][2]
add('twoFbLo', f2(2.0 * round(b_narfix, 2)), 'A', f'2 * FbNarFix (exact 2*b = {2*b_narfix:.4f})')
add('twoFbHi', f2(2.0 * round(b0, 2)), 'A', f'2 * Fb (exact 2*b = {2*b0:.4f})')
s_noTop_up = rep_cal['posterior_with_floors']['noTop3_formal_p95_onesided']
b_noTop = n3['LRD']['slope'][2]
o = r2['orient']
o45_90 = o['thin_45_90']; o60_90 = o['thin_60_90']
add('sigLamMax', f2(2*s_up), 'A', f'2 * s_up (conservative 95% ceiling: {2*s_up:.4f} dex)')
add('sigLamMed', f2(2*s0), 'A', '2 * Fsig under formal-error median (%.4f dex)' % (2*s0))
add('sigLamRep', f2(2*fl_s['LRD_p95_onesided']), 'A', '2 * s_95(same-grating floor) = %.4f dex' % (2*fl_s['LRD_p95_onesided']))
add('oriThin', f2(o['thin_0_45']), 'R')
add('oriHalf', f3(o['h0.5_0_60']), 'R')
add('oriSeven', f3(o['h0.7_0_60']), 'R')
add('oriOne', f3(o['h1_0_60']), 'R')
add('oriThinRatio', f1(o['thin_0_45']/s_up), 'A', 'oriThin / s_up = %.2f' % (o['thin_0_45']/s_up))
add('rhoThin', f'{math.sqrt(1-(s_up/o["thin_0_45"])**2):.2f}', 'A', 'sqrt(1 - (s_up / oriThin)^2)')
add('sigLamIso', f2(2*math.sqrt(s_up**2-o['h1_0_60']**2)), 'A', '2 sqrt(s_up^2 - oriOne^2)')
add('oriEdgeWide', f3(o45_90), 'C')
add('oriEdge', f3(o60_90), 'C')
add('sigLamEdgeWide', f2(2*math.sqrt(s_up**2-o45_90**2)), 'C', '2 sqrt(s_up^2 - oriEdgeWide^2)')
add('sigLamEdge', f2(2*math.sqrt(s_up**2-o60_90**2)), 'C', '2 sqrt(s_up^2 - oriEdge^2)')
add('sigDustMax', f1(s_up/(0.4*b0)), 'C', 'broad-Ha attenuation scatter ceiling s_up / (0.4 Fb) in mag (%.4f)' % (s_up/(0.4*b0)))
add('sigDustNoTop', f1(s_noTop_up/(0.4*b_noTop)), 'C', 's_95(noTop3) / (0.4 FbNoTop) in mag (%.4f)' % (s_noTop_up/(0.4*b_noTop)))
add('RLlo', '0.13', 'V'); add('RLhi', '0.19', 'V')
add('RLtermLo', f3(.065), 'A'); add('RLtermHi', f3(.095), 'A')
add('virW', '2.06', 'V'); add('virL', '0.47', 'V')
add('sigMvir', f'{math.ceil(2.06*s_up*100)/100:.2f}', 'A', f'2.06 * s_up ceiling (unrounded: {2.06*s_up:.6f} dex; s_up = {s_up:.6f} dex)')
add('sigMvirMed', f2(2.06*s0), 'A', '2.06 * Fsig (formal-error median: %.4f dex)' % (2.06*s0))
add('sigMvirRep', f2(2.06*fl_s['LRD_p95_onesided']), 'A', '2.06 * s_95(same-grating floor) = %.4f dex' % (2.06*fl_s['LRD_p95_onesided']))
add('sigMvirCons', f'{math.ceil(2.06*s_cons*100)/100:.2f}', 'A', f'2.06 * conservative ceiling, rounded up (unrounded: {2.06*s_cons:.5f} dex)')
add('sigLamCons', f2(2*s_cons), 'A', f'2 * conservative ceiling ({2*s_cons:.4f} dex)')
add('oriThinRatioCons', f1(o['thin_0_45']/s_cons), 'A', 'oriThin / conservative ceiling')

sec('scattering (unified on s_up, across tau_e = 0.5--3)')
add('rusA', 428, 'V'); add('rusB', 370, 'V')
D_ = lambda t: 428*t/(428*t+370)
add('DHalf', f2(D_(.5)), 'A'); add('DOne', f2(D_(1)), 'A'); add('DTwo', f2(D_(2)), 'A'); add('DThree', f2(D_(3)), 'A')
add('sigTauHalf', f2(s_up/D_(.5)), 'A'); add('sigTauOne', f2(s_up/D_(1)), 'A')
add('sigTauTwo', f2(s_up/D_(2)), 'A'); add('sigTauThree', f2(s_up/D_(3)), 'A')
add('sigTeMax', f2(2*s_up), 'A', '2 * s_up')
for nm_, t_ in (('sigRHalf', .5), ('sigROne', 1), ('sigRTwo', 2), ('sigRThree', 3)):
    add(nm_, f2(2*s_up/D_(t_)), 'A', 'shell-radius scatter ceiling 2 s_up / D(tau_e=%g) = %.4f dex' % (t_, 2*s_up/D_(t_)))
add('tauFixLo', f2(.5*D_(.5)), 'A', '0.5 * D(0.5)')
add('tauFixHi', f2(.5*D_(3)), 'A', '0.5 * D(3) = 0.3881 across tau_e = 0.5--3')
add('tauPhLo', f2(.25*D_(.5)), 'A', '0.25 * D(0.5)')
add('tauPhHi', f2(.25*D_(3)), 'A', '0.25 * D(3)')

st3 = json.load(open(RESULTS / 'scattering_kernel.json'))
st3a = st3['super_eddington_gravity_share']; st3c = st3['extrema_summary']
assert abs(1.0/rl - st3a['gravity_share_max_at_lambda1_k1']) < 1e-11
add('gravMaxPct', round(100.0/rl), 'A')
add('coreLo', st3c['fwhm_in_range_kms'][0], 'D')
add('coreHi', st3c['fwhm_in_range_kms'][1], 'D')
c05 = st3c['pure_conv_gauss_tau0p5']; c23 = st3c['pure_conv_gauss_tau2to3']
r05 = st3c['rt_unscat_plus_scat_gauss_tau0p5']; r23 = st3c['rt_unscat_plus_scat_gauss_tau2to3']
add('DconvHalfLo', f2(c05['D_eff_min']), 'C'); add('DconvHalfHi', f2(c05['D_eff_max']), 'C')
add('sigTauConvHalfLo', f2(s_up/c05['D_eff_max']), 'C', 's_up / DconvHalfHi')
add('sigTauConvHalfHi', f2(s_up/c05['D_eff_min']), 'C', 's_up / DconvHalfLo')
add('DconvHiLo', f2(c23['D_eff_min']), 'C'); add('DconvHiHi', f2(c23['D_eff_max']), 'C')
add('sigTauConvHiLo', f2(s_up/c23['D_eff_max']), 'C', 's_up / DconvHiHi')
add('sigTauConvHiHi', f2(s_up/c23['D_eff_min']), 'C', 's_up / DconvHiLo')
add('sigRConvHiLo', f2(2*s_up/c23['D_eff_max']), 'C', '2 s_up / DconvHiHi')
add('sigRConvHiHi', f2(2*s_up/c23['D_eff_min']), 'C', '2 s_up / DconvHiLo')
add('rtHalfBoostMax', round(100.0*r05['max_fractional_boost_over_core']), 'C')
add('DrtHiLo', f2(r23['D_eff_min']), 'C'); add('DrtHiHi', f2(r23['D_eff_max']), 'C')

sec('pipeline constants and ancillary checks')
try:
    from astropy.cosmology import Planck18, FlatLambdaCDM
    c70 = FlatLambdaCDM(H0=70.0, Om0=0.3)
    dL_grid = [abs(2.0*math.log10(Planck18.luminosity_distance(z_).value / c70.luminosity_distance(z_).value)) for z_ in [3.36 + 0.05*i_ for i_ in range(77)]]
    h0_val, om_val = f'{Planck18.H0.value:.1f}', f'{Planck18.Om0:.3f}'
except ImportError:
    def _dl_flat(z_, h0_, om_):
        n_ = 400; dz_ = z_/n_
        s_ = sum((1.0 if k_ in (0, n_) else (4.0 if k_%2 else 2.0))/math.sqrt(om_*(1.0+k_*dz_)**3 + (1.0-om_)) for k_ in range(n_+1))*dz_/3.0
        return (1.0+z_)*s_/h0_
    dL_grid = [abs(2.0*math.log10(_dl_flat(z_, 67.66, 0.30966) / _dl_flat(z_, 70.0, 0.3))) for z_ in [3.36 + 0.05*i_ for i_ in range(77)]]
    h0_val, om_val = '67.7', '0.310'
add('Hnot', h0_val, 'L'); add('OmM', om_val, 'L')
add('Hcat', '70', 'L'); add('OmCat', '0.3', 'L')
add('dCosmoL', f2(max(dL_grid)), 'A')
add('dCosmoVK', f3(0.25*max(dL_grid)), 'A')
for n_, v_ in (('fitWin', 5000), ('errFloor', 2), ('wnarLo', 30), ('wnarHi', 600), ('wbroadLo', 350),
               ('wbroadHi', '12{,}000'), ('tauMax', 10), ('wabsLo', 50), ('wabsHi', 1500), ('nBootNp', 200), ('winWide', 6800)):
    add(n_, v_, 'D')

lor_scen = st2b
chk1 = anc_chk
chk2 = anc_chk['residuals_and_outliers']
chk3 = anc_chk['window_sensitivity']
chk4 = anc_chk['repeatability_and_bx26']
g_map = {r_['source_id']: r_['y'] for r_ in w1_in['single_gaussian'] if r_['cls'] == 'LRD'}
e_map = {r_['source_id']: r_['y'] for r_ in w1_in['exponential'] if r_['cls'] == 'LRD'}
d_ge_med = statistics.median(g_map[s_] - e_map[s_] for s_ in g_map)
add('NBICgauss', _w[chk1['bic_choice_counts']['LRD_single_gaussian']], 'C')
add('NBIClor', chk1['bic_choice_counts']['LRD_lorentzian'], 'C')
add('NBICexp', _w[chk1['bic_choice_counts']['LRD_exponential']], 'C')
add('dWge', f3(d_ge_med), 'C')
trip('DsigBic', chk1['bic_preferred4_LRD']['sigma_q'], f3, 'C')
add('Fmad', f3(chk2['mad_normalized_dex']), 'C')
add('FlooLo', f3(chk2['loo_sigma_min'][1]), 'C'); add('FlooHi', f3(chk2['loo_sigma_max'][1]), 'C')
add('FblooLo', f2(chk2['loo_slope_min'][1]), 'C'); add('FblooHi', f2(chk2['loo_slope_max'][1]), 'C')
add('rResSN', f2(chk2['residual_vs_logSNR_pearson'][0]), 'C'); add('pResSN', f2(chk2['residual_vs_logSNR_pearson'][1]), 'C')
add('ratHiTwoLo', f'{chk2["A18_top2_ratios"][0][1]:.1f}', 'C')
add('ratHiTwoHi', f'{chk2["A18_top2_ratios"][1][1]:.1f}', 'C')
add('dWinMed', f3(chk3['win6800_noSII']['median_abs_dlogFWHM_dex']), 'C')
trip('FsigWin', chk3['win6800_noSII']['LRD_sigma_q'], f3, 'C')
add('NMultiLRD', _w[chk4['N_LRDs_with_multiple_QC_spectra']], 'C')
add('NPairsDup', chk4['N_duplicate_pairs'], 'C')
add('dLDupMed', f2(chk4['median_abs_dlogL_dex']), 'C')
add('dWDupMed', f2(chk4['median_abs_dlogFWHM_dex']), 'C')
trip('FsigLowSN', chk4['lowest_SNR_spectrum_choice_LRD_sigma_q'], f3, 'C')
add('XsigSxZero', f3(chk4['bx26_sx0_LRD_sigma_q'][2]), 'C')
add('XsigoSxZero', f3(chk4['bx26_sx0_ctrl_sigma_q'][2]), 'C')
add('NXNoCutLRD', bm['LRD']['N'] - chk4['W4_candidates_LRD'], 'A')
trip('XsigNoCut', chk4['bx26_no13_LRD_sigma_q'], f3, 'C')
add('XdsigNoCut', neg(f3(chk4['bx26_no13_diff_q'][2])), 'C')
o_sdss = json.load(open(DATA / 'sdss' / 'regression_results.json'))['regressions']['single_gaussian_QC']
trip('Ssgb', o_sdss['slope'], f2, 'X')
inf_esc = [100.0 * v for v in anc_chk['bic_best_infinite_escape_fraction_median'].values()]
add('fEscInfLo', f'{min(inf_esc):.0f}', 'X'); add('fEscInfHi', f'{max(inf_esc):.0f}', 'X')
add('recLG', f3(lor_scen['lor_matched_s005']['bayes_sigma_med_quantiles'][2]), 'C')
add('recEightL', f3(lor_scen['lor_matched_s008']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitEightL', _w[lor_scen['lor_matched_s008']['bayes_hits_le_obs_Fsig']], 'C')
add('recTenL', f3(lor_scen['lor_matched_s010']['bayes_sigma_med_quantiles'][2]), 'C')
add('hitTenL', _w[lor_scen['lor_matched_s010']['bayes_hits_le_obs_Fsig']], 'C')
tot10_all = tot10 + lor_scen['lor_matched_s010']['bayes_hits_le_obs_Fsig']
add('hitTenAll', _w[tot10_all], 'A')
add('nullDenThree', 300, 'A')
b_null_e = lor_scen['exp_null_matched_s005']['bayes_slope_med_quantiles'][2]
b_null_l = lor_scen['lor_null_matched_s005']['bayes_slope_med_quantiles'][2]
b_null_hi = max(lor_scen['exp_null_matched_s005']['bayes_slope_med_quantiles'][4],
                lor_scen['lor_null_matched_s005']['bayes_slope_med_quantiles'][4])
hits_b0 = lor_scen['exp_null_matched_s005']['slope_hits_ge_obs_Fb'] + lor_scen['lor_null_matched_s005']['slope_hits_ge_obs_Fb']
add('bNullMed', f2(b_null_e), 'C')
add('bNullHi', f2(b_null_hi), 'C', f'max 97.5th percentile slope across null-slope rows ({b_null_hi:.4f})')
add('hitSlopeNull', _w[hits_b0], 'C')
k_iso_exp = math.sqrt(2.0)*math.log(2.0)/math.sqrt(3.0)
rex_all = [wr[g_]['exponential']['ratio_median'] for g_ in ('main_E', 'Curtis_B_own', 'Curtis_A_own')]
lam_iso_exp = [min((k_iso_exp/r_)**2 for r_ in rex_all), max((k_iso_exp/r_)**2 for r_ in rex_all)]
lam_vir_exp = [min((0.964486/r_)**2 for r_ in rex_all), max((0.964486/r_)**2 for r_ in rex_all)]
add('kIsoExp', f2(k_iso_exp), 'A')
add('lamIsoExpLo', f2(lam_iso_exp[0]), 'A'); add('lamIsoExpHi', f2(lam_iso_exp[1]), 'A')
add('lamVirExpLo', f2(lam_vir_exp[0]), 'A'); add('lamVirExpHi', f2(lam_vir_exp[1]), 'A')

add('ratTightLo', f1(min(Ao['sigma_fixed_photosphere'][2]/Ao['line_only']['sigma'][2], Bo['sigma_fixed_photosphere'][2]/Bo['line_only']['sigma'][2])), 'A', 'ratio of fixed-photosphere to line-only scatter (dex), Samples A and B')
add('ratTightHi', f1(max(Ao['sigma_fixed_photosphere'][2]/Ao['line_only']['sigma'][2], Bo['sigma_fixed_photosphere'][2]/Bo['line_only']['sigma'][2])), 'A')

M = [(n_, v_, 'C' if st_ == 'N' else st_, nt_) for n_, v_, st_, nt_ in M]
# Scan ms.tex plus any tables/*.tex files so macros used in generated tables are also included in numbers.tex.
tex_sources = []
ms_path = os.path.join(HERE, 'ms.tex')
if os.path.exists(ms_path):
    tex_sources.append(open(ms_path, encoding='utf-8').read())
tab_dir = os.path.join(HERE, 'tables')
if os.path.isdir(tab_dir):
    for fn_ in sorted(os.listdir(tab_dir)):
        if fn_.endswith('.tex'):
            tex_sources.append(open(os.path.join(tab_dir, fn_), encoding='utf-8').read())
SRC = re.sub(r'(?<!\\)%.*', '', '\n'.join(tex_sources))
USED = {n_ for n_, *_ in M if n_ != '%' and re.search(r'\\%s(?![A-Za-z])' % n_, SRC)} if SRC else {n_ for n_, *_ in M if n_ != '%'}
USED |= {'dvEscKthreeLo', 'dvEscKthreeHi', 'dvEscKfourMax'}
M_all = M
M = [m_ for i, m_ in enumerate(M) if (m_[0] in USED) or (m_[0] == '%' and any(x[0] in USED for x in M[i+1:next((j for j in range(i+1, len(M)) if M[j][0] == '%'), len(M))]))]
lines = ['%% numbers.tex -- generated by make_numbers.py; do not edit by hand.',
         '%% Write every macro as \\name{} in running text so that TeX does not swallow the following space.',
         '%% Provenance tag after each macro: R/C read from stored results in data/; A arithmetic on other values;',
         '%% D design choice of the analysis; L literature value; V/X entered by hand in make_numbers.py (see note).', '']
for n, v, st, note in M:
    if n == '%':
        lines.append(f'\n%% ---------------- {v}')
        continue
    tag = f'  % [{st}]' + (f' {note}' if note else '')
    lines.append(f'\\newcommand{{\\{n}}}{{{v}}}{tag}')
open(os.path.join(HERE, 'numbers.tex'), 'w', encoding='utf-8').write('\n'.join(lines)+'\n')
json.dump({n: dict(value=v, in_manuscript=n in USED) for n, v, st, note in M_all if n != '%'}, open(RESULTS / 'numbers.json', 'w'), indent=1, ensure_ascii=False)
print(len([1 for m_ in M_all if m_[0] != '%']), 'numbers computed;', len(USED), 'used in manuscript/tables and written to paper/numbers.tex')
