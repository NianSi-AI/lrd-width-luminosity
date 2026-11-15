"""Generate Tables 1, 2, 3, and Appendix Table B1 of the manuscript into paper/tables/.

Usage: python scripts/make_tables.py
"""
import csv
import json
import os

from lrdwidth.paths import DATA, PAPER, RESULTS, WIDTHS

T = str(PAPER / "tables")
os.makedirs(T, exist_ok=True)
r1_raw = json.load(open(RESULTS / "width_luminosity.json"))
r1 = {("W1_" + k[6:] if k.startswith("width_") else k): v for k, v in r1_raw.items()}
r1["W6"] = json.load(open(RESULTS / "sdss_catalog_benchmark.json"))
r2 = json.load(open(RESULTS / "robustness.json"))
rep_cal = json.load(open(RESULTS / "repeat_calibration.json"))
src = str(DATA / "sample" / "table1_sources.csv")

rows = list(csv.DictReader(open(src)))
assert len(rows) == 44
open(os.path.join(T, 'Table1_full.csv'), 'w').write(open(src).read())
key = ['RUBIES-EGS-42046', 'CEERS-746', 'RUBIES-UDS-119957', 'Skyfire-299', 'RUBIES-EGS-15825']
t1 = [
    r'\begin{deluxetable*}{llcccllcccc}',
    r'\tabletypesize{\footnotesize}',
    r'\tablewidth{0pt}',
    r'\tablecaption{The sample\label{tab:sample}}',
    r'\tablehead{\colhead{Source} & \colhead{Class} & \colhead{R.A.} & \colhead{Decl.} & \colhead{$z$} & \colhead{Program} & \colhead{Grating} & \colhead{Label} & \colhead{FWHM} & \colhead{$\log L_{\rm H\alpha}$} & \colhead{S/N} \\',
    r'\colhead{} & \colhead{} & \colhead{(deg)} & \colhead{(deg)} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{(km s$^{-1}$)} & \colhead{(erg s$^{-1}$)} & \colhead{}}',
    r'\startdata',
]
for k in key:
    r = [x for x in rows if x['name'] == k][0]
    t1.append(f"{r['name']} & {'LRD' if r['cls'] == 'LRD' else 'other'} & {float(r['ra']):.5f} & ${float(r['dec']):.5f}$ & {float(r['z']):.3f} & {r['program']} & {r['grating']} & {r['label_ref']} & ${r['fwhm']}\\pm{r['fwhm_err']}$ & ${float(r['logL']):.2f}\\pm{float(r['logL_err']):.2f}$ & {r['snr']} \\\\")
t1 += [
    r'\enddata',
    r'\tablecomments{FWHM is the intrinsic width of the single broad Gaussian and $L_{\rm H\alpha}$ the luminosity of that component before absorption, without dust correction. Uncertainties are $1\sigma$, conditional on the adopted model. Labels: T25 = \citet{Taylor2025}; K26 = \citet{Kocevski2026}; H25 = \citet{Hviding2025}. Five rows are shown. The full table of 44 sources, with the five other width measures and all covariances, is published in machine-readable form.}',
    r'\end{deluxetable*}',
]
open(os.path.join(T, 'Table1_stub.tex'), 'w').write('\n'.join(t1) + '\n')

# ------------------------------------------------------------------ Table 2: width-luminosity fits + variants
def P(q):
    return f'${q[2]:.3f}$ (${q[0]:.3f}$ to ${q[4]:.3f}$)' if q[0] < 0 else f'${q[2]:.3f}$ (${q[0]:.3f}$--${q[4]:.3f}$)'

ND = r'\nodata'

def line(lab, NL, bq, sL, NO, sO, dq):
    bs = (f'${bq[2]:.2f}$ (${bq[0]:.2f}$ to ${bq[4]:.2f}$)' if bq[0] < 0 else f'${bq[2]:.2f}$ (${bq[0]:.2f}$--${bq[4]:.2f}$)') if bq else ND
    return f"{lab} & {NL} & {bs} & {P(sL) if sL else ND} & {NO} & {P(sO) if sO else ND} & {P(dq) if dq else ND} \\\\"

W = lambda k: r1['W1_' + k]
lh = r2['frozen_LSFhalf']
n3 = r2['frozen_noTop3']
na = r2['frozen_noAbs']
nw = r2['frozen_noWide2']
bm = r2['bx_main']
bn = r2['bx_noOverlap']
w6 = r1['W6']
fl_same = rep_cal['posterior_with_floors']['same_grating_sy_floor']
fl_all = rep_cal['posterior_with_floors']['all_pairs_sy_floor']
o_sdss = json.load(open(DATA / 'sdss' / 'regression_results.json'))['regressions']['single_gaussian_QC']

t2 = [
    r'\centerwidetable',
    r'\begin{deluxetable*}{lclclcl}',
    r'\tabletypesize{\scriptsize}',
    r'\tabcolsep=1.7pt',
    r'\tablewidth{0pt}',
    r'\tablecaption{Width--luminosity fits\label{tab:fits}}',
    r'\tablehead{\colhead{} & \multicolumn{3}{c}{LRDs} & \multicolumn{2}{c}{Other broad-line sources} & \colhead{} \\ \cline{2-4} \cline{5-6}',
    r'\colhead{Width measure or variant} & \colhead{$N$} & \colhead{$b$} & \colhead{$\sigma_{\rm int}$ (dex)} & \colhead{$N$} & \colhead{$\sigma_{\rm int}$ (dex)} & \colhead{$\Delta\sigma_{\rm int}$ (dex)}}',
    r'\startdata',
    r'\cutinhead{This work: six measures of width}',
]
for lab, k in (
    ('Single-Gaussian FWHM\\tablenotemark{a}', 'single_gaussian'),
    (r'Interquartile width $W_{50}$', 'nonparametric_W50'),
    ('Exponential FWHM', 'exponential'),
    ('Lorentzian FWHM', 'lorentzian'),
    ('Nonparametric FWHM', 'nonparametric_FWHM'),
    ('Two-Gaussian FWHM', 'double_gaussian'),
):
    t2.append(line(lab, 29, W(k)['LRD']['slope'], W(k)['LRD']['sigma'], 15, W(k)['other']['sigma'], W(k)['diff']))

t2.append(r'\cutinhead{This work: single-Gaussian FWHM, error-floor and sample variants}')
t2.append(line(r'Repeat floor, same grating ($+0.044$ dex)\tablenotemark{a}', 29, fl_same['LRD_slope_q'], fl_same['LRD_sigma_q'], 15, fl_same['other_sigma_q'], fl_same['diff_q']))
t2.append(line(r'Repeat floor, all pairs ($+0.060$ dex)\tablenotemark{a}', 29, fl_all['LRD_slope_q'], fl_all['LRD_sigma_q'], 15, fl_all['other_sigma_q'], fl_all['diff_q']))
t2.append(line(r'LSF width $\times0.5$\tablenotemark{a}', 29, lh['LRD']['slope'], lh['LRD']['sigma'], 15, lh['other']['sigma'], lh['diff']))
t2.append(line('Three brightest LRDs removed', 26, n3['LRD']['slope'], n3['LRD']['sigma'], 15, n3['other']['sigma'], n3['diff']))
assert na['LRD']['N'] == 24 and na['other']['N'] == 15
t2.append(r'LRDs without Balmer absorption & \NNoAbs{} & $\FbNoAbs{}$ ($\FbNoAbslo{}$--$\FbNoAbshi{}$) & $\FsigNoAbs{}$ ($\FsigNoAbslo{}$--$\FsigNoAbshi{}$) & \NOther{} & $\Fsigo{}$ ($\Fsigolo{}$--$\Fsigohi{}$) & $\FdsigNoAbs{}$ ($\FdsigNoAbslo{}$ to $\FdsigNoAbshi{}$) \\')
t2.append(r'Narrow-line FWHM $<\wnarHi{}$~\kms{}\,\tablenotemark{a} & \NLRDNarCut{} & $\FbNarCut{}$ ($\FbNarCutlo{}$--$\FbNarCuthi{}$) & $\FsigNarCut{}$ ($\FsigNarCutlo{}$--$\FsigNarCuthi{}$) & \NOther{} & $\Fsigo{}$ ($\Fsigolo{}$--$\Fsigohi{}$) & $\FdsigNarCut{}$ ($\FdsigNarCutlo{}$ to $\FdsigNarCuthi{}$) \\')
t2.append(r'Narrow-line FWHM fixed to \wnarDef{}~\kms{}\,\tablenotemark{a} & \NLRD{} & $\FbNarFix{}$ ($\FbNarFixlo{}$--$\FbNarFixhi{}$) & $\FsigNarFix{}$ ($\FsigNarFixlo{}$--$\FsigNarFixhi{}$) & \NOther{} & $\FsigoNarFix{}$ ($\FsigoNarFixlo{}$--$\FsigoNarFixhi{}$) & $\FdsigNarFix{}$ ($\FdsigNarFixlo{}$ to $\FdsigNarFixhi{}$) \\')
t2.append(r'Narrow-line FWHM fixed to [\ion{O}{3}]\tablenotemark{a} & \NLRD{} & $\FbNarOIII{}$ ($\FbNarOIIIlo{}$--$\FbNarOIIIhi{}$) & $\FsigNarOIII{}$ ($\FsigNarOIIIlo{}$--$\FsigNarOIIIhi{}$) & \NOther{} & $\FsigoNarOIII{}$ ($\FsigoNarOIIIlo{}$--$\FsigoNarOIIIhi{}$) & $\FdsigNarOIII{}$ ($\FdsigNarOIIIlo{}$ to $\FdsigNarOIIIhi{}$) \\')
rcw = {v_['key']: v_ for v_ in json.load(open(RESULTS / 'review_checks.json'))['variants']}
ws_ = json.load(open(RESULTS / 'ancillary_checks.json'))['window_sensitivity']['win8000_maskSII']
t2.append(line(r'[\ion{O}{3}]-anchored, absorption rule re-applied\tablenotemark{a,e}', 27, None, rcw['oiii_qc']['sigma_q'], ND, None, None))
t2.append(line(r'Fit window $\pm\winWider{}$~\kms{}, [\ion{S}{2}] masked', 29, ws_['LRD_slope_q'], ws_['LRD_sigma_q'], 15, ws_['other_sigma_q'], ws_['diff_q']))
t2.append(line(r'Lowest-S/N passing spectrum\tablenotemark{a}', 29, None, rcw['lowest_snr']['sigma_q'], ND, None, None))
t2.append(r'In-band parent ${\rm S/N}\geq\cutSNthree{}$ included\tablenotemark{a} & \NLRDSNthree{} & $\FbSNthree{}$ ($\FbSNthreelo{}$--$\FbSNthreehi{}$) & $\FsigSNthree{}$ ($\FsigSNthreelo{}$--$\FsigSNthreehi{}$) & \nodata & \nodata & \nodata \\')
t2.append(line('Two widest other sources removed', 29, nw['LRD']['slope'], nw['LRD']['sigma'], 13, nw['other']['sigma'], nw['diff']))

t2.append(r'\cutinhead{Table of \citet{BaccusXu2026}}')
t2.append(line(r'All records\tablenotemark{b}', 48, bm['LRD']['slope'], bm['LRD']['sigma'], 97, bm['other']['sigma'], bm['diff']))
t2.append(line('Shared records removed', bn['LRD']['N'], bn['LRD']['slope'], bn['LRD']['sigma'], bn['other']['N'], bn['other']['sigma'], bn['diff']))
t2.append(r'Matched pairs\tablenotemark{c} & 36 & \nodata & \nodata & 36 & \nodata & $-0.065$ ($-0.105$ to $-0.024$) \\')

t2.append(r'\cutinhead{SDSS DR7 quasars \citep{Shen2011}}')
t2.append(f"Multi-Gaussian FWHM (catalog) & {ND} & {ND} & {ND} & {w6['N']} & {P(w6['sigma'])} & {ND} \\\\")
sdss_b = o_sdss['slope']
sdss_bs = f"${sdss_b[2]:.2f}$ (${sdss_b[0]:.2f}$--${sdss_b[4]:.2f}$)"
t2.append(f"Single-Gaussian FWHM, refitted\\tablenotemark{{d}} & {ND} & {sdss_bs} & {ND} & 224 & $0.123$ ($0.112$--$0.136$) & {ND} \\\\")

t2 += [
    r'\enddata',
    r'\tablecomments{Posterior medians, with equal-tailed 95\% intervals in parentheses. $b$ is the slope of $\log{\rm FWHM}$ on $\log L$ (shown in the LRDs column for the refitted SDSS row as well), and $\Delta\sigma_{\rm int}$ is the scatter of the LRDs minus that of the other sources. Each width measure is regressed on the luminosity and covariance of its own fit ($L_{\rm H\alpha}$ for our sample and SDSS; $L_{\rm bol}$ for \citealt{BaccusXu2026}).}',
    r'\tablenotetext{a}{One-sided 95\% upper limits (95th posterior percentiles) on LRD $\sigma_{\rm int}$ are \FsigUp{}~dex under formal errors, \FsigFlrSameUp{}~dex with the $+\flrDupSameG{}$~dex same-grating repeat-spectrum floor added in quadrature to $C_{yy}$, \FsigFlrAllUp{}~dex with the $+\flrDupAll{}$~dex all-pairs floor, \FsigLsfUp{}~dex for LSF width $\times0.5$, \FsigNarCutUp{}~dex for narrow-line FWHM $<\wnarHi{}$~\kms{}, \FsigNarFixUp{}~dex for narrow-line FWHM fixed to \wnarDef{}~\kms{}, \FsigNarOIIIUp{}~dex when anchored to [\ion{O}{3}]~$\lambda5007$ (\wnarDef{}~\kms{} where unmeasured; [\ion{O}{3}] hits the \wnarLo{}~\kms{} lower bound in \NOIIIUnres{} of the \NOIIIMeas{} LRDs where it is measured, and fixing all \NLRD{} to \wnarLo{}~\kms{} gives $\sigma_{\rm int}=\FsigNarMin{}$~dex, $<\FsigNarMinUp{}$~dex), \FsigOIIIqcUp{}~dex when the absorption-bound rule is re-applied to the [\ion{O}{3}]-anchored fits, \FsigLowSNUp{}~dex for the lowest-S/N passing spectra, and \FsigSNthreeUp{}~dex with in-band failed ${\rm S/N}\geq\cutSNthree{}$ parent LRD spectra included (Sections~\ref{sec:relation} and \ref{sec:tests}). Random-effects fits to all repeat spectra give \reAllUp{}~dex (additive floor) and \reMultUp{}~dex (multiplicative).}',
    r'\tablenotetext{e}{Two LRDs whose absorption optical depth reaches its bound in the [\ion{O}{3}]-anchored refit are removed. Refit covariances are not stored, so the change from the 29-object row is computed with the baseline covariances.}',
    r'\tablenotetext{b}{The \NXLRD{} LRD records correspond to \NXLRDuniq{} unique objects on the sky; deduplicating to one record per object gives $\sigma_{\rm int}=\XsigDedup{}$~dex (\XsigDeduplo{}--\XsigDeduphi{}).}',
    r'\tablenotetext{c}{Pairs matched in grating, redshift ($|\Delta z|<0.5$), and luminosity ($|\Delta\log L|<0.2$). The entry is the difference of the residual standard deviations, with its bootstrap interval.}',
    r'\tablenotetext{d}{A random 300 of the 3755 quasars refitted with the procedure of Section~\ref{sec:fit}; 224 pass the rules of Section~\ref{sec:sample}.}',
    r'\end{deluxetable*}',
]
open(os.path.join(T, 'Table2_fits.tex'), 'w').write('\n'.join(t2) + '\n')

# ------------------------------------------------------------------ Table 3: physical synthesis of the two readings
t3 = [
    r'\begin{deluxetable*}{p{0.25\textwidth}p{0.34\textwidth}p{0.34\textwidth}}',
    r'\tabletypesize{\scriptsize}',
    r'\tablewidth{0pt}',
    r'\tablecaption{Empirical constraints on the two readings of LRD broad \Ha{} widths\label{tab:readings}}',
    r'\tablehead{\colhead{Empirical property} & \colhead{Bound orbital motion (${\rm FWHM}=k\sqrt{GM/R}$)} & \colhead{Scattering envelope ($\lambda\geq1$, $R\geq\Rph$)}}',
    r'\startdata',
    r'Velocity scale \& wings: ${\rm FWHM}/\vK=\ratLo{}\text{--}\ratHi{}$; \fEscNpLo{}--\fEscNpHi{}\% of flux at $|v|>\sqrt{2}\,\vK$'
    r' & For a spherical photosphere with $R\geq\Rph$: $\lambda\leq\lamRingHi{}$ ($k=\kRing{}$, edge-on ring) and $\lambda\leq\lamVirLo{}\text{--}\lamIsoHi{}$ ($k=\kVir{}\text{--}\kIso{}$). Relaxes as $\Rph/R$ if $R<\Rph$, or drops out if no continuum photosphere exists.'
    r' & Escape speed $\leq\sqrt{2}\,\vK\,\lambda^{-1/2}$; core broadened $\brdVirLo{}\text{--}\brdVirHi{}\times$ ($k=\kVir{}$), $\brdIsoLo{}\text{--}\brdIsoHi{}\times$ ($k=\kIso{}$), or $\brdRingLo{}\text{--}\brdRingHi{}\times$ ($k=\kRing{}$), matched by $\tau_e\simeq2$--3 at $T_e\sim10^{4}$~K. \\[3pt]',
    r'Width scatter at fixed $\LHa$: $\sigint<\FsigUp{}$~dex (baseline; $<\FsigCons{}$~dex across analysis variants, $<\reAllUp{}\text{--}\reMultUp{}$~dex with repeat-spectrum variance)'
    r' & $\sigma_{M/R}\leq\sigLamMax{}$~dex; requires $H/R\simeq1$ ($\sigma_{M/R}\leq\sigLamIso{}$~dex) or inclinations near the disk plane ($i>45^\circ$, $\sigma_{M/R}\leq\sigLamEdgeWide{}\text{--}\sigLamEdge{}$~dex), and broad-\Ha{} dust scatter $\sigma_A\leq\sigDustMax{}$~mag at $b=\Fb{}$.'
    r' & $\sigma_{\log\tau_e}\leq\sigTauThree{}\text{--}\sigTauTwo{}$~dex across $\tau_e=2$--3 (\sigTauConvHiLo{}--\sigTauConvHiHi{}~dex with core convolution), or $\sigma_{\log T_e}\leq\sigTeMax{}$~dex; shell-radius budget $\sigma_{\log R}\lesssim\sigRThree{}\text{--}\sigRTwo{}$~dex at fixed geometry. \\[3pt]',
    r'Width--luminosity slope: $b=\Fb{}$ (\FbNarFix{}--\FbNarOIII{} with narrow-line bounds or \NBright{} brightest LRDs removed)'
    r' & Requires $M/R\propto\LHa^{2b}=\LHa^{\twoFbLo{}\text{--}\twoFbHi{}}$ across the population.'
    r' & Recombination shell ($\LHa\propto R\,\tau_e^{2}$) predicts ${\rm FWHM}\propto\LHa^{\shellSlopeLo{}\text{--}\shellSlopeHi{}}$ at $\tau_e\simeq2$--3 depending on $R(\LHa)$: matches $b=\Fb{}$, but the narrow-line-corrected $b$ lies at or below its lower edge. \\[3pt]',
    r'Null $\Lph^{1/4}\Teff$ scaling at fixed $\LHa$ despite $\sigma(\log\Lph/\LHa)=\sdLphHaLo{}\text{--}\sdLphHaHi{}$~dex'
    r' & Line-emitting radius $R$ must track ionizing/recombination luminosity (photoionized BLR) rather than the optical continuum photosphere $\Rph$.'
    r' & Recombination and Thomson scattering occur in the same ionized zone ($\LHa\propto R\,\tau_e^{2}$), coupling $\tau_e$ directly to $\LHa$ rather than $\Rph$. \\[3pt]',
    r'Single-epoch virial mass ($\Mvir\propto\LHa^{\virL{}}\,{\rm FWHM}^{\virW{}}$)'
    r' & Scatters intrinsically by $<\sigMvir{}$~dex (\sigMvirMed{}~dex at formal median, $<\sigMvirCons{}$~dex at the conservative ceiling) at fixed $\LHa$; enclosed mass is $(\mEddLo{}\text{--}\mEddHi{})\,(R/\Rph)\,M_{\rm Edd}(\Lph)$ if $R\geq\Rph$.'
    r' & Estimator still scatters by $<\sigMvir{}$~dex at fixed $\LHa$ (restating luminosity), while central mass is $\lambda^{-1}M_{\rm Edd}(\Lph)$. \\',
    r'\enddata',
    r'\tablecomments{All scatter ceilings in columns 2 and 3 are evaluated at the baseline one-sided 95\% upper limit $\sigint=\FsigUp{}$~dex; adopting the same-grating repeat-spectrum upper limit (\FsigFlrSameUp{}~dex) tightens each linear bound by a quarter, whereas the conservative ceiling (\FsigCons{}~dex) loosens them by \narCeilLoosenPct{}\%. $M_{\rm Edd}(\Lph)=\kappa_{\rm es}\Lph/(4\pi G c)=8\times10^{5}\,\Msun\,(\Lph/10^{44}\,\ergs)$.}',
    r'\end{deluxetable*}',
]
open(os.path.join(T, 'Table3_readings.tex'), 'w').write('\n'.join(t3) + '\n')

# ------------------------------------------------------------------ Table B1: injection-recovery summary (Appendix B)
X = json.load(open(RESULTS / 'measurement_summary.json'))['numbers']
st2 = json.load(open(WIDTHS / 'sloped_injection_recovery.json'))['sloped_simulations']
lor_scen = st2

def w2_entry(sc, inj, meas):
    base = f'W2_final_summary_s{sc}__{sc}_{inj}_M4_Gaussian_QC_{meas}__sigma_median__'
    p16 = X[base + 'p16_p50_p84__p16']['value']
    p50 = X[base + 'p16_p50_p84__p50']['value']
    p84 = X[base + 'p16_p50_p84__p84']['value']
    cnt = int(round(X[base + 'count_le_observed']['value']))
    den = int(round(X[base + 'denominator']['value']))
    return f'${p50:.3f}$ (${p16:.3f}$--${p84:.3f}$)', f'{cnt}/{den}'

def q_entry(d_obj):
    q = d_obj['bayes_sigma_med_quantiles']
    cnt = int(d_obj['bayes_hits_le_obs_Fsig'])
    return f'${q[2]:.3f}$ (${q[1]:.3f}$--${q[3]:.3f}$)', f'{cnt}/100'

def b_entry(d_obj):
    if 'bayes_slope_med_quantiles' not in d_obj:
        return ND
    qb = d_obj['bayes_slope_med_quantiles']
    return f'${qb[2]:.2f}$ (${qb[1]:.2f}$--${qb[3]:.2f}$)' if qb[1] >= 0 else f'${qb[2]:.2f}$ (${qb[1]:.2f}$ to ${qb[3]:.2f}$)'

tb1 = [
    r'\begin{deluxetable*}{llccccc}',
    r'\tabletypesize{\scriptsize}',
    r'\tabcolsep=2.2pt',
    r'\tablewidth{0pt}',
    r'\tablecaption{Injection--recovery tests on the LRD noise and LSF templates\label{tab:injections}}',
    r'\tablehead{\colhead{Injected profile \& relation} & \colhead{Width measure} & \colhead{$\sigma_{\rm parent}$} & \colhead{Injected $b$} & \colhead{Recovered $\sigma_{\rm int}$ (dex)} & \colhead{Trials $\leq\Fsig{}$~dex} & \colhead{Recovered $b$}}',
    r'\startdata',
    r'\cutinhead{Cross-profile matrix at fiducial scatter}',
]
for lab_i, inj in (('Gaussian, sloped', 'single_gaussian'), ('Exponential, sloped', 'exponential')):
    for lab_m, meas in (('Single-Gaussian', 'single_gaussian'), (r'Interquartile $W_{50}$', 'nonparametric_W50'),
                        ('Exponential', 'exponential'), ('Lorentzian', 'lorentzian')):
        s_str, h_str = w2_entry(3, inj, meas)
        tb1.append(f'{lab_i} & {lab_m} & $0.05$ & $0.17$ & {s_str} & {h_str} & {ND} \\\\')

s_lg, h_lg = q_entry(lor_scen['lor_matched_s005'])
tb1.append(f"Lorentzian, matched locus & Single-Gaussian & $0.05$ & $0.15$ & {s_lg} & {h_lg} & {b_entry(lor_scen['lor_matched_s005'])} \\\\")

tb1.append(r'\cutinhead{Parent-scatter and truncation tests}')
for lab_i, d_obj, sp, bi in (
    ('Gaussian, observed locus', st2['gauss_rel_s008'], '0.08', '0.17'),
    ('Exponential, matched locus', st2['exp_matched_s008'], '0.08', '0.18'),
    ('Lorentzian, matched locus', lor_scen['lor_matched_s008'], '0.08', '0.15'),
    ('Gaussian, observed locus', st2['gauss_rel_s010'], '0.10', '0.17'),
    ('Exponential, matched locus', st2['exp_matched_s010'], '0.10', '0.18'),
    ('Lorentzian, matched locus', lor_scen['lor_matched_s010'], '0.10', '0.15'),
    (r'Exponential, $+0.32$~dex wider', st2['exp_wide_s010'], '0.10', '0.18'),
):
    s_str, h_str = q_entry(d_obj)
    tb1.append(f'{lab_i} & Single-Gaussian & ${sp}$ & ${bi}$ & {s_str} & {h_str} & {b_entry(d_obj)} \\\\')

for lab_i, inj in (('Gaussian, lum.-independent', 'single_gaussian'), ('Exponential, lum.-independent', 'exponential')):
    s_str, h_str = w2_entry(1, inj, 'single_gaussian')
    tb1.append(f'{lab_i} & Single-Gaussian & $0.15$ & $0.00$ & {s_str} & {h_str} & {ND} \\\\')

tb1.append(r'\cutinhead{Null-slope tests ($b=0$)}')
for lab_i, d_obj in (('Exponential, matched width', lor_scen['exp_null_matched_s005']),
                     ('Lorentzian, matched width', lor_scen['lor_null_matched_s005'])):
    s_str, h_str = q_entry(d_obj)
    tb1.append(f'{lab_i} & Single-Gaussian & $0.05$ & $0.00$ & {s_str} & {h_str} & {b_entry(d_obj)} \\\\')

tb1.append(r'\cutinhead{Full in-band parent templates (\NParentInBand{} LRDs)\tablenotemark{a}}')
tb1.append(r'Gaussian, observed locus & Single-Gaussian & $0.05$ & $0.17$ & $\recParFiveG{}$ ($\recParFiveGlo{}$--$\recParFiveGhi{}$) & \hitParFiveGnum{} & $\bParFiveG{}$ ($\bParFiveGlo{}$--$\bParFiveGhi{}$) \\')
tb1.append(r'Exponential, matched locus & Single-Gaussian & $0.05$ & $0.18$ & $\recParFiveE{}$ ($\recParFiveElo{}$--$\recParFiveEhi{}$) & \hitParFiveEnum{} & $\bParFiveE{}$ ($\bParFiveElo{}$--$\bParFiveEhi{}$) \\')
tb1.append(r'Gaussian, observed locus & Single-Gaussian & $0.10$ & $0.17$ & $\recParTenG{}$ ($\recParTenGlo{}$--$\recParTenGhi{}$) & \hitParTenGnum{} & $\bParTenG{}$ ($\bParTenGlo{}$--$\bParTenGhi{}$) \\')
tb1.append(r'Exponential, matched locus & Single-Gaussian & $0.10$ & $0.18$ & $\recParTenE{}$ ($\recParTenElo{}$--$\recParTenEhi{}$) & \hitParTenEnum{} & $\bParTenE{}$ ($\bParTenElo{}$--$\bParTenEhi{}$) \\')

tb1 += [
    r'\enddata',
    r'\tablecomments{Each row summarizes \NReal{} simulated samples generated on the pixel grids, noise vectors, and grating LSFs of the observed LRD spectra (\NLRD{} surviving templates in the first three blocks; all \NParentInBand{} in-band parent templates in the bottom block), refitted with the pipeline of Section~\ref{sec:fit}, and filtered by the quality rules of Section~\ref{sec:sample}. Recovered $\sigma_{\rm int}$ and $b$ list the median (16th--84th percentiles) of the Bayesian posterior medians across the \NReal{} realizations. In the null-slope ($b=0$) rows, \hitSlopeNull{} of \nullDen{} realizations yield a posterior-median slope as steep as $b=\Fb{}$ (97.5th percentile $\leq\bNullHi{}$).}',
    r'\tablenotetext{a}{Across the \NParentInBand{}-template suites, a median of \NParSelG{} (\NParSelGlo{}--\NParSelGhi{}, Gaussian) and \NParSelE{} (\NParSelElo{}--\NParSelEhi{}, exponential) sources pass Section~\ref{sec:sample} per realization; replacing refitted/published $\LHa$ with pipeline-fitted $\LHa$ on all failed templates yields $\sigma_{\rm int}\leq\Fsig{}$~dex in \hitParOwnTenGnum{} and \hitParOwnTenEnum{} realizations at $\sigma_{\rm parent}=\simA{}$~dex.}',
    r'\end{deluxetable*}',
]
open(os.path.join(T, 'TableB1_injections.tex'), 'w').write('\n'.join(tb1) + '\n')

print('tables written:', sorted(os.listdir(T)))
