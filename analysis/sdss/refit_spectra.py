"""Fit broad Halpha in SDSS Legacy spectra with the model used for the JWST sample (lrdwidth.linefit).

Only the spectrum loader and the wavelength-dependent line-spread function differ. Fluxes in cgs, velocities in
km/s, Planck18 distances. One JSON per source is written to data/sdss/fits/; existing files are never overwritten.

Usage: python analysis/sdss/refit_spectra.py [--selfcheck] [--first]
The spectra listed in data/sdss/sample_spectra.json must first be downloaded; data/sdss/download_successes.json
lists them (keys: source_id, plate, mjd, fiber, z, filename, parent_index).
"""

import concurrent.futures
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18
from astropy.io import fits
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.sparse import csr_matrix
from scipy.special import ndtr

from lrdwidth import linefit as c
from lrdwidth.paths import DATA

O = DATA / "sdss"


def native(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        return {k: native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [native(v) for v in x]
    return x


def save(p, x):
    with Path(p).open("x") as f:
        json.dump(native(x), f, ensure_ascii=False, indent=2, allow_nan=False)


def sdss_response(lo, hi, grid, loglam, wdisp, center_A, step):
    # SDSS resolution is RMS sigma in pixels of log10(lambda).  The extra lambda/
    # center converts a log-lambda sigma into the linear velocity coordinate.
    wave = 10**loglam
    sigpix = wdisp * np.log(10) * step * c.C * wave / center_A
    good = np.isfinite(sigpix) & (sigpix > 0)
    sigma = np.interp(center_A * (1 + grid / c.C), wave[good], sigpix[good])
    K = (
        (ndtr((hi[:, None] - grid) / sigma) - ndtr((lo[:, None] - grid) / sigma))
        * 10
        / (hi - lo)[:, None]
    )
    K[K < 1e-12] = 0
    sig0 = float(np.interp(center_A, wave[good], sigpix[good]))
    return (csr_matrix(K), sig0)


def load(row):
    path = O / "spectra" / row["filename"]
    with fits.open(path, memmap=False) as hd:
        a = hd["COADD"].data
        head = hd[0].header
        f = np.array(a["flux"], float)
        iv = np.array(a["ivar"], float)
        ll = np.array(a["loglam"], float)
        wd = np.array(a["wdisp"], float)
        assert int(head["PLATEID"]) == row["plate"] and int(head["MJD"]) == row["mjd"]
        fiber = int(head.get("FIBERID", hd["SPALL"].data["FIBERID"][0] if "SPALL" in hd else -1))
        assert fiber == row["fiber"]
        step = float(head.get("COEFF1", 0.0001))
        assert abs(step - 0.0001) < 1e-07 and abs(np.median(np.diff(ll)) - 0.0001) < 1e-06
        masks = {k: int(np.count_nonzero(a[k])) for k in ["and_mask", "or_mask"] if k in a.names}
        vers = {k: head.get(k) for k in ["RUN2D", "VERS2D", "VERSCOMB", "COEFF0", "COEFF1"]}
    wave = 10**ll
    z = row["z"]
    center = c.REST["Halpha"] * 10000.0 * (1 + z)
    v = c.C * (wave / center - 1)
    ve = c.bin_edges(v)
    e = np.full(len(iv), np.inf)
    ok = iv > 0
    e[ok] = 1 / np.sqrt(iv[ok])
    valid = np.isfinite(f + e + wd) & (e > 0) & (f != 0) & (wd > 0)
    if valid.sum() < 10:
        return {"status": "insufficient_valid_pixels", "Nvalid": int(valid.sum())}
    valid &= e > np.percentile(e[valid], 5) * 0.5
    ref = np.interp(wave, wave[valid], median_filter(e[valid], size=11))
    valid &= e > 0.2 * ref
    select = valid & (abs(v) <= 5000)
    ix = np.flatnonzero(select)
    center_valid = wave[0] < center < wave[-1] and valid[np.argmin(abs(wave - center))]
    if not (center_valid and len(ix) >= 25 and (v[ix] < -1000).any() and (v[ix] > 1000).any()):
        return {
            "status": "not_covered_or_insufficient_valid_pixels",
            "Nvalid": len(ix),
            "center_valid": bool(center_valid),
        }
    # SDSS native flux units are 1e-17 erg/s/cm2/Angstrom, while the JWST fits is
    # in 1e-18 erg/s/cm2 per 1000 km/s. d(lambda)/dv=center_A/c, no extra (1+z).
    factor = 1e-17 * center / c.C * 1000 / c.FUNIT
    y = f[ix] * factor
    err = np.sqrt(e[ix] ** 2 + (0.02 * f[ix]) ** 2) * factor
    grid = np.arange(-15000.0, 15000.01, 10.0)
    K, sig0 = sdss_response(ve[ix], ve[ix + 1], grid, ll, wd, center, step)
    return dict(
        status="ready",
        line="Halpha",
        z=z,
        center_um=center / 10000.0,
        v=v[ix],
        lo=ve[ix],
        hi=ve[ix + 1],
        y=y,
        err=err,
        quoted_err=e[ix] * factor,
        grid=grid,
        K=K,
        R_center=c.C / (c.SQ * sig0),
        lsf_scale=1.0,
        grating="SDSS_native_wdisp",
        window=5000.0,
        Nvalid=len(ix),
        valid_fraction=float(select.sum() / np.sum(abs(v) <= 5000)),
        wavelength_um=wave[ix] / 10000.0,
        LSF_sigma_center_kms=sig0,
        log10lambda_step=step,
        header_versions=vers,
        mask_nonzero_counts_full_spectrum=masks,
        pixel_selection="native ivar>0; finite f/err/wdisp; f!=0; wdisp>0; frozen low-error filters; no extra OR-mask veto",
        native_wdisp=wd[ix],
    )


def null_fit(d, absorption):
    names, xx, lo, hi = c.setup(d, "single_gaussian", absorption, 0)
    xx[5:8] = [0, 0, 1000]
    free = [i for i, n in enumerate(names) if n not in ["F_broad1", "v_broad1", "W_broad1"]]
    best = None
    for seed in range(3):
        _, start, _, _ = c.setup(d, "single_gaussian", absorption, seed)
        start = np.clip(start, lo + 1e-08, hi - 1e-08)

        def residual(p):
            pp = xx.copy()
            pp[free] = p
            return (c.components(pp, d, "single_gaussian", absorption)["total"] - d["y"]) / d["err"]

        f = least_squares(
            residual,
            start[free],
            bounds=(lo[free], hi[free]),
            x_scale="jac",
            max_nfev=1200,
            ftol=2e-07,
            xtol=2e-07,
            gtol=2e-07,
        )
        if best is None or np.sum(f.fun**2) < np.sum(best.fun**2):
            best = f
    pp = xx.copy()
    pp[free] = best.x
    chi = float(best.fun @ best.fun)
    return dict(
        absorption=absorption,
        chi2=chi,
        BIC=chi + len(free) * np.log(len(d["y"])),
        names=names,
        params=pp,
        success=bool(best.success),
    )


def one(row):
    t = time.monotonic()
    dest = O / "fits" / (row["source_id"] + ".json")
    if dest.exists():
        return {"source_id": row["source_id"], "status": "already_exists"}
    try:
        d = load(row)
        out = {"source_id": row["source_id"], "input": row, "status": d["status"]}
        if d["status"] != "ready":
            out["coverage"] = d
            save(dest, out)
            return out
        nulls = [null_fit(d, False), null_fit(d, True)]
        null = nulls[1] if nulls[0]["BIC"] - nulls[1]["BIC"] > 6 else nulls[0]
        candidates = [c.fit_model(d, "single_gaussian", a, seeds=3) for a in [False, True]]
        diff = candidates[0]["BIC"] - candidates[1]["BIC"]
        f = candidates[1] if diff > 6 else candidates[0]
        p = dict(zip(f["names"], f["params"]))
        w = p["W_broad1"]
        db = null["BIC"] - f["BIC"]
        lsf = c.C / d["R_center"]
        checks = {
            "SNR_ge5": (f["broad_SNR_conditional"] or 0) >= 5,
            "deltaBIC_ge6": db >= 6,
            "width_gt_LSF": w > lsf,
            "width_not_bound": "W_broad1" not in f["bound_parameters"],
            "chi2_lt2": f["reduced_chi2"] < 2,
            "optimizer": f["success"],
            "abs_not_bound": not (
                f["absorption"] and bool(set(f["bound_parameters"]) & {"tau_abs", "W_abs"})
            ),
        }
        ii = [f["names"].index(k) for k in ["F_broad1", "W_broad1"]]
        J = np.diag(1 / np.log(10) / np.array(f["params"])[ii])
        cov = J @ np.array(f["covariance_conditional"])[np.ix_(ii, ii)] @ J.T
        lum = (
            f["broad_flux_erg_s_cm2"]
            * 4
            * np.pi
            * Planck18.luminosity_distance(row["z"]).to_value("cm") ** 2
        )
        regression = {
            "source_id": row["source_id"],
            "x": float(np.log10(lum) - 43),
            "y": float(np.log10(w)),
            "cov": cov,
            "logLHa_ergs": float(np.log10(lum)),
            "FWHM_kms": w,
        }
        assert np.isfinite(cov).all() and np.linalg.eigvalsh(cov).min() > -1e-12
        out.update(
            status="fitted",
            main=f,
            candidates=candidates,
            null_candidates=nulls,
            absorption_deltaBIC=diff,
            deltaBIC_null=db,
            LSF_FWHM_kms=lsf,
            quality=checks,
            quality_pass=all(checks.values()),
            regression=regression,
            data={k: v for k, v in d.items() if k not in ["K", "grid"]},
        )
        if db >= 6 and (f["broad_SNR_conditional"] or 0) > 0:
            try:
                out["nonparametric"] = c.nonparametric(
                    d, f, nboot=0, seed=2026111301 + row["parent_index"]
                )
            except Exception as e:
                out["nonparametric_error"] = {"error": str(e), "traceback": traceback.format_exc()}
        out["elapsed_s"] = time.monotonic() - t
        save(dest, out)
        return dict(
            source_id=row["source_id"],
            status="fitted",
            quality_pass=all(checks.values()),
            chi2=f["reduced_chi2"],
            FWHM=w,
            seconds=round(time.monotonic() - t, 2),
        )
    except Exception as e:
        out = {
            "source_id": row["source_id"],
            "status": "fit_failed",
            "error": str(e),
            "traceback": traceback.format_exc(),
        }
        if not dest.exists():
            save(dest, out)
        return out


def checks():
    # Frozen physical/unit selfchecks plus an independent response fixture.
    chk = c.self_checks()
    ll = np.arange(np.log10(7000), np.log10(9000), 0.0001)
    w = np.ones(len(ll)) * 1.2
    center = 8000.0
    v = c.C * (10**ll / center - 1)
    mask = abs(v) < 5000
    ve = c.bin_edges(v)
    ix = np.flatnonzero(mask)
    g = np.arange(-15000, 15000.01, 10)
    K, sigma = sdss_response(ve[ix], ve[ix + 1], g, ll, w, center, 0.0001)
    out = K @ (1000 * c.gauss(g, 0, 2000))
    assert abs(np.sum(out * np.diff(ve)[ix]) / 1000 - 1) < 1e-06
    assert abs(sigma / (c.C * np.log(10) * 0.0001 * 1.2) - 1) < 1e-06
    # Known 1e-14 CGS integrated line maps to 1e4 internal flux units.
    assert abs(1e-14 / c.FUNIT - 10000.0) < 1e-09
    return {"linefit": chk, "sdss_flux_response": "passed", "wdisp_sigma_conversion": "passed"}


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        save(O / "selfcheck.json", checks())
        print("selfchecks passed")
        sys.exit()
    (O / "fits").mkdir(exist_ok=True)
    rows = json.loads((O / "download_successes.json").read_text())
    if "--first" in sys.argv:
        rows = rows[:1]
    pending = [r for r in rows if not (O / "fits" / (r["source_id"] + ".json")).exists()]
    print("pending", len(pending), flush=True)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=int(os.environ.get("LRD_SDSS_WORKERS", "4"))
    ) as pool:
        for i, res in enumerate(pool.map(one, pending), 1):
            print(i, json.dumps(native(res), ensure_ascii=False), flush=True)
