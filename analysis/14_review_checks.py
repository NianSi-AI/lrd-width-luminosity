"""Robustness checks added after internal review (round 2), all computed from stored per-source
measurements (no spectral refits):

  1. Analysis-variant summary of the LRD intrinsic scatter (median, 95% interval, one-sided 95%
     upper limit where the per-source inputs are stored), used in Section 4 and Figure 2;
  2. the [O III]-anchored narrow-line variant with the absorption-bound rule of Section 2.3
     re-applied (approximation: baseline covariances, see note in the output);
  3. a constant-width (b = 0) model for all 29 LRDs and for the 26 below the bright end;
  4. the photospheric partial slope with RUBIES-UDS-29813 removed (Samples A and B);
  5. LRD-minus-control scatter differences with all formal errors inflated or floored,
     in our sample and in the Baccus-Xu (BX26) table;
  6. repeat-spectrum diagnostics: G395H-minus-G395M widths, primary-spectrum selection offset,
     the G395M-only subsample, and RUBIES-EGS-42232 with its G235M primary replaced by G395M;
  7. SDSS refitted quasars with the reduced-chi^2 rule removed;
  8. the A2744-QSO1 fixed-narrow residual recomputed against Equation (3);
  9. electron-scattering shell slope predictions at tau_e = 2-3;
 10. the sky separation of CEERS-746 and Skyfire-299.

Writes data/results/review_checks.json.
"""

import csv
import json
import math
import os

for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np

from lrdwidth import stats as S
from lrdwidth.paths import DATA, PAPER, PHOTOSPHERE, RESULTS, WIDTHS
from lrdwidth.photosphere import arr, load_sample, partial

Q5 = (0.025, 0.16, 0.5, 0.84, 0.975)


def post(rows, cyy_scale=1.0, cxx_scale=1.0, floor=0.0, slope=None):
    """Grid posterior of Equation (1)-(2); `slope=0` fixes b = 0 (constant-width model)."""
    x = np.array([r["x"] for r in rows], float)
    y = np.array([r["y"] for r in rows], float)
    C = np.array([r["cov"] for r in rows], float)
    kw = {} if slope is None else dict(b_lo=slope - 1e-9, b_hi=slope + 1e-9, nb=1)
    b, s, p, _ = S.post2(x, y, C[:, 0, 0] * cxx_scale ** 2, C[:, 1, 1] * cyy_scale ** 2 + floor ** 2,
                         C[:, 0, 1] * cxx_scale * cyy_scale, **kw)
    ps = p.sum(0)
    out = dict(N=len(rows), sigma_q=[float(v) for v in S.quant(s, ps)],
               sigma_p95_onesided=float(S.quant(s, ps, qs=(0.95,))[0]))
    if slope is None:
        out["slope_q"] = [float(v) for v in S.quant(b, p.sum(1))]
    return out, s, ps


def summary(rows, **kw):
    return post(rows, **kw)[0]


def diff(rows_a, rows_b, **kw):
    ra, s, pa = post(rows_a, **kw)
    rb, _, pb = post(rows_b, **kw)
    return dict(LRD=ra, other=rb, diff_q=[float(v) for v in S.diffq(s, pa, pb)])


def main():
    out = {}
    w1 = json.loads((WIDTHS / "regression_inputs.json").read_text())["single_gaussian"]
    lrd = [r for r in w1 if r["cls"] == "LRD"]
    oth = [r for r in w1 if r["cls"] != "LRD"]
    assert len(lrd) == 29 and len(oth) == 15
    base = summary(lrd)
    assert abs(base["sigma_p95_onesided"] - 0.07472) < 2e-4

    # ---------------- 1-2. analysis variants ----------------
    rob = json.loads((RESULTS / "robustness.json").read_text())
    rep = json.loads((RESULTS / "repeat_calibration.json").read_text())
    nar = json.loads((RESULTS / "frozen_fixed_narrow_audit.json").read_text())
    anc = json.loads((RESULTS / "ancillary_checks.json").read_text())
    mem = json.loads((RESULTS / "membership.json").read_text())
    par = json.loads((RESULTS / "parent_covariance_audit.json").read_text())
    ree = json.loads((RESULTS / "repeat_random_effects.json").read_text())
    win = json.loads((WIDTHS / "window_and_wing_checks.json").read_text())
    lsf = json.loads((WIDTHS / "lsf_half_inputs.json").read_text())

    by_id = {r["source_id"]: r for r in lrd}
    f44 = {x["source_id"]: x for x in nar["frozen_44_audit"] if x["cls"] == "LRD"}
    assert len(f44) == 29

    def refit_rows(wkey, lkey, drop=()):
        # Variant widths/luminosities with the baseline covariance of each source (approximation).
        return [dict(source_id=s, x=f44[s][lkey] - 43.0, y=math.log10(f44[s][wkey]), cov=by_id[s]["cov"])
                for s in by_id if s not in drop]

    tau_bound = sorted(s for s, x in f44.items() if "tau_abs" in x["new_bound_parameters"])
    assert tau_bound == ["RUBIES-EGS-29489", "RUBIES-EGS-55604"]
    oiii_approx_all = summary(refit_rows("new_FWHM_kms", "new_log10_LHa"))
    oiii_approx_qc = summary(refit_rows("new_FWHM_kms", "new_log10_LHa", drop=tau_bound))
    stored_oiii = nar["reg_frozen_OIII_anchored_29LRD_15oth"]
    # Shift the stored (exact-covariance) result by the change found under the approximation.
    shift_med = oiii_approx_qc["sigma_q"][2] - oiii_approx_all["sigma_q"][2]
    shift_p95 = oiii_approx_qc["sigma_p95_onesided"] - oiii_approx_all["sigma_p95_onesided"]
    out["oiii_anchored_qc_consistent"] = dict(
        note="[O III]-anchored refit with the Section 2.3 absorption-bound rule re-applied; tau_abs reaches "
             "its bound in the two dropped LRDs. Variant covariances are not stored, so the effect of the drop "
             "is computed with baseline covariances and added to the stored exact-covariance result.",
        dropped=tau_bound, N=27,
        approx_all29=oiii_approx_all, approx_qc27=oiii_approx_qc,
        sigma_med=stored_oiii["sigma_q"][2] + shift_med,
        sigma_p95_onesided=stored_oiii["sigma_p95_onesided"] + shift_p95,
        slope_med=oiii_approx_qc["slope_q"][2],
    )

    low_rows = [dict(r) for r in lrd]
    for i, r in enumerate(low_rows):
        if r["source_id"] in win["lowest_snr_substitute_rows_9"]:
            low_rows[i] = win["lowest_snr_substitute_rows_9"][r["source_id"]]
    lowsn = summary(low_rows)
    assert abs(lowsn["sigma_q"][2] - anc["repeatability_and_bx26"]["lowest_SNR_spectrum_choice_LRD_sigma_q"][2]) < 2e-3

    bundle = json.loads((WIDTHS / "parent_and_repeat_covariances.json").read_text())
    groups = bundle["qc_spectra_by_group_29_43"]
    g42232 = [v for v in groups.values() if v[0]["source_id"] == "RUBIES-EGS-42232"][0]
    assert g42232[0]["grating"].startswith("G235M")
    m395 = max((s for s in g42232 if s["grating"].startswith("G395M")), key=lambda s: s["snr"])
    rows_42232 = [dict(r) if r["source_id"] != "RUBIES-EGS-42232" else dict(source_id=r["source_id"], x=m395["x"], y=m395["y"], cov=m395["cov"]) for r in lrd]
    swap42232 = summary(rows_42232)

    table1 = {r["name"]: r for r in csv.DictReader((PAPER / "tables" / "Table1_full.csv").open())}
    g395m_rows = [r for r in lrd if table1[r["source_id"]]["grating"] == "G395M"]
    g395m = summary(g395m_rows)

    V = []

    def add(key, label, group, q, p95=None, N=None, source=""):
        V.append(dict(key=key, label=label, group=group, N=N, sigma_q=[float(v) for v in q],
                      sigma_p95_onesided=None if p95 is None else float(p95), source=source))

    add("baseline", "Baseline (single Gaussian)", "baseline", base["sigma_q"], base["sigma_p95_onesided"], 29, "regression_inputs.json")
    lh = rob["frozen_LSFhalf"]["LRD"]
    add("lsf_half", r"LSF width $\times0.5$", "fit", lh["sigma"], summary([r for r in lsf if r["cls"] == "LRD"])["sigma_p95_onesided"], 29, "lsf_half_inputs.json")
    for k, lab in (("win6800_noSII", r"Window $\pm6800$ km s$^{-1}$"), ("win8000_maskSII", r"Window $\pm8000$ km s$^{-1}$")):
        add(k, lab, "fit", win["window_sensitivity"][k]["LRD_sigma_q"], None, 29, "window_and_wing_checks.json (quantiles only)")
    add("lowest_snr", "Lowest-S/N spectrum", "fit", lowsn["sigma_q"], lowsn["sigma_p95_onesided"], 29, "window_and_wing_checks.json rows")
    add("swap42232", "42232 on G395M", "fit", swap42232["sigma_q"], swap42232["sigma_p95_onesided"], 29, "parent_and_repeat_covariances.json")
    add("g395m_only", "G395M primaries only", "fit", g395m["sigma_q"], g395m["sigma_p95_onesided"], len(g395m_rows), "regression_inputs.json")
    for k, lab in (("reg_frozen_universal_W150_29LRD_15oth", r"Narrow fixed 150 km s$^{-1}$"),
                   ("reg_frozen_universal_W30_29LRD_15oth", r"Narrow fixed 30 km s$^{-1}$"),
                   ("reg_frozen_OIII_anchored_29LRD_15oth", "Narrow tied to [O III]")):
        add(k, lab, "narrow", nar[k]["sigma_q"], nar[k]["sigma_p95_onesided"], 29, "frozen_fixed_narrow_audit.json")
    oq = out["oiii_anchored_qc_consistent"]
    add("oiii_qc", "[O III], QC re-applied", "narrow",
        [v + shift_med for v in stored_oiii["sigma_q"]], oq["sigma_p95_onesided"], 27, "this script (approximate)")
    drop_nar = [s for s, x in f44.items() if x["orig_W_narrow_kms"] > 599.0]
    nar_cut = summary([r for r in lrd if r["source_id"] not in drop_nar])
    add("narrow_lt600", r"Free narrow $<600$ km s$^{-1}$", "narrow", nar_cut["sigma_q"], nar_cut["sigma_p95_onesided"], 27, "regression_inputs.json")
    add("no_top3", "Without 3 brightest", "sample", rob["frozen_noTop3"]["LRD"]["sigma"],
        rep["posterior_with_floors"]["noTop3_formal_p95_onesided"], 26, "robustness.json")
    r38 = par["regressions_full_2x2_covariance"]["row4_plus_valid_SNge3_38"]
    add("parent_sn3", r"+9 failed, S/N$\geq3$", "sample", r38["sigma_q"], r38["sigma_p95_onesided"], 38, "parent_covariance_audit.json")
    add("any_label", r"Any LRD label, no $\chi^2$/abs.", "sample", mem["any_label_all"]["sigma"], None, mem["any_label_all"]["N"], "membership.json (quantiles only)")
    fl = rep["posterior_with_floors"]
    add("floor_same", "Same-grating floor", "repeat", fl["same_grating_sy_floor"]["LRD_sigma_q"], fl["same_grating_sy_floor"]["LRD_p95_onesided"], 29, "repeat_calibration.json")
    add("floor_all", "All-pairs floor", "repeat", fl["all_pairs_sy_floor"]["LRD_sigma_q"], fl["all_pairs_sy_floor"]["LRD_p95_onesided"], 29, "repeat_calibration.json")
    rem = ree["random_effects_models"]
    add("re_add", "Random effects, additive", "repeat", rem["exact_2dcov_fixed_primary_x"]["all_43_spectra"]["sigma_q"],
        rem["exact_2dcov_fixed_primary_x"]["all_43_spectra"]["sigma_p95_onesided"], 29, "repeat_random_effects.json")
    add("re_mult", "Random effects, multiplicative", "repeat", rem["multiplicative_fixed_primary_x"]["all_43_spectra"]["sigma_q"],
        rem["multiplicative_fixed_primary_x"]["all_43_spectra"]["sigma_p95_onesided"], 29, "repeat_random_effects.json")

    analysis = [v for v in V if v["group"] in ("fit", "narrow", "sample")]
    p95s = [v["sigma_p95_onesided"] for v in analysis if v["sigma_p95_onesided"] is not None]
    out["variants"] = V
    out["narrow_bound_counts"] = dict(
        N_LRD=len(f44),
        free_narrow_at_30=sum(1 for x in f44.values() if x["orig_W_narrow_kms"] <= 30.01),
        free_narrow_at_600=sum(1 for x in f44.values() if x["orig_W_narrow_kms"] >= 599.0),
        OIII_measured=sum(1 for x in f44.values() if x.get("OIII_W_kms") is not None),
        OIII_at_30=sum(1 for x in f44.values() if x.get("OIII_W_kms") is not None and x["OIII_W_kms"] <= 30.01),
        rms_dlogFWHM_OIII_minus_base=float(np.sqrt(np.mean([x["dlogFWHM_new_minus_orig"] ** 2 for x in f44.values()]))),
        max_abs_dlogFWHM=max(((abs(x["dlogFWHM_new_minus_orig"]), s) for s, x in f44.items()))[0],
        max_abs_source=max(((abs(x["dlogFWHM_new_minus_orig"]), s) for s, x in f44.items()))[1],
    )
    fa = par["failed_21_audit"]
    add9 = [r for r in fa if r["in_passband"] and r["SNR_broad"] >= 3.0 and "W_broad1" not in r["bound_parameters"]]
    assert len(add9) == 9
    out["parent_sn3_added"] = dict(N=len(add9), N_dBIC_lt6=sum(1 for r in add9 if r["delta_BIC_null"] < 6.0),
                                   N_dBIC_negative=sum(1 for r in add9 if r["delta_BIC_null"] < 0.0))
    out["variant_extremes"] = dict(
        note="Over fitting, narrow-line, and sample variants (repeat-floor rows excluded).",
        max_p95_onesided=max(p95s),
        max_p95_key=[v["key"] for v in analysis if v["sigma_p95_onesided"] == max(p95s)][0],
        max_upper95=max(v["sigma_q"][4] for v in analysis),
        max_upper95_key=max(analysis, key=lambda v: v["sigma_q"][4])["key"],
        max_median=max(v["sigma_q"][2] for v in analysis),
        max_median_key=max(analysis, key=lambda v: v["sigma_q"][2])["key"],
    )

    # ---------------- 3. constant-width model ----------------
    top3 = sorted(lrd, key=lambda r: -r["x"])[:3]
    lrd26 = [r for r in lrd if r not in top3]
    out["constant_width"] = {
        "LRD29": dict(post(lrd, slope=0.0)[0], raw_sd_logFWHM=float(np.std([r["y"] for r in lrd], ddof=1))),
        "LRD26": dict(post(lrd26, slope=0.0)[0], raw_sd_logFWHM=float(np.std([r["y"] for r in lrd26], ddof=1)),
                      logL_range=[min(r["x"] for r in lrd26) + 43, max(r["x"] for r in lrd26) + 43],
                      FWHM_range_kms=[10 ** min(r["y"] for r in lrd26), 10 ** max(r["y"] for r in lrd26)],
                      FWHM_median_kms=float(10 ** np.median([r["y"] for r in lrd26]))),
        "LRD26_free_slope": summary(lrd26),
        "top3": [r["source_id"] for r in top3],
    }

    # ---------------- 4. photospheric partial slope without RUBIES-UDS-29813 ----------------
    ph = {}
    for key in ("Curtis_A_own", "Curtis_B_own"):
        rows = load_sample(PHOTOSPHERE / "samples.json", key)
        res = {}
        for tag, keep in (("all", rows), ("drop_RUBIES-UDS-29813", [r for r in rows if r["id"] != "RUBIES-UDS29813"])):
            xl, xp, y = arr(keep, "xl"), arr(keep, "xp"), arr(keep, "y")
            b1, b2, s, p = partial(xl, xp, y, arr(keep, "cll"), arr(keep, "sx"), arr(keep, "cyy"), arr(keep, "cly"))
            res[tag] = dict(N=len(keep), b_ph_q=[float(v) for v in S.quant(b2, p.sum((0, 2)))])
        assert res["drop_RUBIES-UDS-29813"]["N"] == res["all"]["N"] - 1
        xl, xp = arr(rows, "xl"), arr(rows, "xp")
        rx = xp - np.polyval(np.polyfit(xl, xp, 1), xl)
        i29813 = [r["id"] for r in rows].index("RUBIES-UDS29813")
        res["RUBIES-UDS-29813_dlogLph_at_fixed_LHa"] = float(rx[i29813])
        res["RUBIES-UDS-29813_is_most_negative"] = bool(int(np.argmin(rx)) == i29813)
        ph[key] = res
    out["photosphere_leverage"] = ph

    # ---------------- 5. error inflation and floors ----------------
    bx = json.loads((DATA / "bx26" / "regression_inputs.json").read_text())
    bl = [t for t in bx if t["cls"] == "LRD"]
    lo, hi = min(t["logLbol"] for t in bl), max(t["logLbol"] for t in bl)
    bc = [t for t in bx if t["cls"] != "LRD" and t["z"] > 3 and lo <= t["logLbol"] <= hi]
    assert len(bl) == 48 and len(bc) == 97
    bxrows = lambda R: [dict(x=t["logLbol"] - 45.0, y=t["y"], cov=[[t["sx"] ** 2, 0.0], [0.0, t["sy"] ** 2]]) for t in R]
    infl = {}
    for tag, kw in (("formal", {}), ("x1.5", dict(cyy_scale=1.5, cxx_scale=1.5)), ("x2.0", dict(cyy_scale=2.0, cxx_scale=2.0)),
                    ("x2.8", dict(cyy_scale=2.8, cxx_scale=2.8)), ("floor0.044", dict(floor=0.044)),
                    ("floor0.060", dict(floor=0.060)), ("floor0.080", dict(floor=0.080))):
        ours, bxd = diff(lrd, oth, **kw), diff(bxrows(bl), bxrows(bc), **kw)
        infl[tag] = dict(ours_LRD_sigma_med=ours["LRD"]["sigma_q"][2], ours_other_sigma_med=ours["other"]["sigma_q"][2],
                         ours_diff_q=ours["diff_q"], bx_LRD_sigma_med=bxd["LRD"]["sigma_q"][2],
                         bx_other_sigma_med=bxd["other"]["sigma_q"][2], bx_diff_q=bxd["diff_q"])
    out["error_inflation"] = infl
    out["error_inflation_summary"] = dict(
        diff_med_min=min(min(v["ours_diff_q"][2], v["bx_diff_q"][2]) for k, v in infl.items() if k != "formal"),
        diff_med_max=max(max(v["ours_diff_q"][2], v["bx_diff_q"][2]) for k, v in infl.items() if k != "formal"),
        diff_upper95_max=max(max(v["ours_diff_q"][4], v["bx_diff_q"][4]) for k, v in infl.items() if k != "formal"),
    )

    # ---------------- 6. repeat-spectrum diagnostics ----------------
    pairs = list(csv.DictReader((WIDTHS / "repeat_pairs_22.tsv").open(), delimiter="\t"))
    hm = []
    for r in pairs:
        ga, gb = r["grating_a"][:5], r["grating_b"][:5]
        if {ga, gb} == {"G395H", "G395M"}:
            d = float(r["dlogFWHM_dex"])
            hm.append((r["source_id"], d if ga == "G395H" else -d))
    per_src = {}
    for s_, d_ in hm:
        per_src.setdefault(s_, []).append(d_)
    per_src_med = {s_: float(np.median(v)) for s_, v in per_src.items()}
    prim = [float(r["dlogFWHM_dex"]) for r in pairs if r["is_primary_pair"] == "True"]
    lowminushigh = [(-float(r["dlogFWHM_dex"]) if float(r["snr_a"]) > float(r["snr_b"]) else float(r["dlogFWHM_dex"])) for r in pairs]
    a_base, b_base = rob["frozen_main"]["LRD"]["a"], base["slope_q"][2]
    resid = {r["source_id"]: r["y"] - (a_base + b_base * r["x"]) for r in lrd}
    gr = {}
    for s_, r in resid.items():
        gr.setdefault(table1[s_]["grating"], []).append(r)
    off = ree["grating_offset_model"]
    out["repeat_diagnostics"] = dict(
        G395H_minus_G395M_pairs=len(hm), G395H_minus_G395M_objects=len(per_src),
        G395H_narrower_objects=sum(1 for v in per_src_med.values() if v < 0),
        G395H_minus_G395M_object_median=float(np.median(list(per_src_med.values()))),
        G395H_minus_G395M_per_object=per_src_med,
        primary_minus_secondary_mean=float(np.mean(prim)), N_primary_pairs=len(prim),
        lower_snr_wider_fraction=float(np.mean(np.array(lowminushigh) > 0)),
        residual_mean_by_primary_grating={g: dict(N=len(v), mean=float(np.mean(v))) for g, v in gr.items()},
        free_offset_model=off,
        swap42232=dict(file_id=m395["file_id"], grating=m395["grating"], snr=m395["snr"], **swap42232),
        G395M_only=g395m,
    )

    # ---------------- 7. SDSS without the chi^2 rule ----------------
    sd_rows = list(csv.DictReader((DATA / "sdss" / "measurements.tsv").open(), delimiter="\t"))
    sd_in = json.loads((DATA / "sdss" / "regression_inputs.json").read_text())["single_gaussian_QC"]
    chi_only = [r["source_id"] for r in sd_rows if r["failures"] == "chi2_lt2"]
    add_rows = [json.loads((DATA / "sdss" / "fits" / f"{s}.json").read_text())["regression"] for s in chi_only]
    sd_base = summary(sd_in)
    sd_nochi = summary(sd_in + add_rows)
    assert abs(sd_base["sigma_q"][2] - 0.12324) < 2e-4
    sn_fail = [float(r["FWHM_single_kms"]) for r in sd_rows if "SNR_ge5" in r["failures"] and r["FWHM_single_kms"]]
    out["sdss_no_chi2_rule"] = dict(N_base=len(sd_in), N_added=len(add_rows), base=sd_base, no_chi2=sd_nochi,
                                    N_SNR_failed=len(sn_fail), SNR_failed_FWHM_median_kms=float(np.median(sn_fail)),
                                    N_SNR_failed_at_12000=sum(1 for w in sn_fail if w > 11990.0))

    # ---------------- 8. QSO1 residual against Equation (3) ----------------
    fm = rob["frozen_main"]["LRD"]
    a3, b3 = fm["a"], json.loads((RESULTS / "width_luminosity.json").read_text())["width_single_gaussian"]["LRD"]["slope"][2]
    q1 = [x for x in json.loads((RESULTS / "narrow_degeneracy_refit.json").read_text())["critical_10_plus_QSO1"] if x["source_id"] == "A2744-QSO1"][0]
    w150 = q1["refit_W150_sensitivity"]
    r_eq3 = math.log10(w150["FWHM_kms"]) - (a3 + b3 * (w150["log10_LHa_delensed_erg_s"] - 43.0))
    out["qso1_eq3_check"] = dict(Eq3_a=a3, Eq3_b=b3, recomputed_residual_dex=r_eq3,
                                 stored_residual_Eq3_delensed_dex=w150["residual_Eq3_delensed_dex"],
                                 stored_residual_vs_fixed_narrow_relation=nar["QSO1_residual_vs_fixed_narrow_29LRD_relation"]["QSO1_W150_resid_dex"],
                                 free_narrow_SNR=q1["orig_SNR_broad"], free_narrow_dBIC_null=q1["orig_delta_BIC_null"],
                                 fixed_narrow_dBIC_vs_free_null=w150["delta_BIC_vs_free_null"])
    assert abs(r_eq3 - w150["residual_Eq3_delensed_dex"]) < 5e-3

    # ---------------- 9. shell slope predictions at tau_e = 2-3 ----------------
    D = lambda t: 428.0 * t / (428.0 * t + 370.0)
    sk = json.loads((RESULTS / "scattering_kernel.json").read_text())["extrema_summary"]
    c23, r23 = sk["pure_conv_gauss_tau2to3"], sk["rt_unscat_plus_scat_gauss_tau2to3"]
    out["shell_slopes_tau2to3"] = dict(
        note="slope = D/4 for R ~ L^(1/2) (tau_e ~ L^(1/4)) up to D/2 at fixed R (tau_e ~ L^(1/2))",
        pure_wing=[0.25 * D(2.0), 0.5 * D(3.0)],
        core_convolved=[0.25 * c23["D_eff_min"], 0.5 * c23["D_eff_max"]],
        transmitted_core=[0.25 * r23["D_eff_min"], 0.5 * r23["D_eff_max"]],
        observed_narrow_variant_slopes=[nar["reg_frozen_universal_W150_29LRD_15oth"]["slope_q"], stored_oiii["slope_q"]],
    )

    # ---------------- 10. CEERS-746 and Skyfire-299 ----------------
    a_, b_ = table1["CEERS-746"], table1["Skyfire-299"]
    ra1, de1, ra2, de2 = (math.radians(float(v)) for v in (a_["ra"], a_["dec"], b_["ra"], b_["dec"]))
    sep = math.degrees(math.acos(min(1.0, math.sin(de1) * math.sin(de2) + math.cos(de1) * math.cos(de2) * math.cos(ra1 - ra2)))) * 3600.0
    out["ceers746_skyfire299"] = dict(sep_arcsec=sep, z=[float(a_["z"]), float(b_["z"])],
                                      fwhm=[float(a_["fwhm"]), float(b_["fwhm"])], logL=[float(a_["logL"]), float(b_["logL"])],
                                      files=[a_["file"], b_["file"]])

    (RESULTS / "review_checks.json").write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", RESULTS / "review_checks.json")


if __name__ == "__main__":
    main()
