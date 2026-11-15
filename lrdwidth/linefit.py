"""Emission-line fits of NIRSpec grating spectra: fluxes in cgs, widths in km/s.

The line-spread function is Gaussian with the default R(lambda) of msaexp for each grating, integrated exactly over
the pixels. The intrinsic profile is forwarded through that response; widths are not corrected by subtracting the
instrumental FWHM in quadrature.
"""

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import median_filter
from scipy.optimize import least_squares, nnls
from scipy.sparse import csr_matrix
from scipy.special import ndtr, ndtri

from .paths import DATA

C = 299792.458
FUNIT = 1e-18
SQ = 2 * np.sqrt(2 * np.log(2))
CALIB = DATA / "calibration"
if not CALIB.exists():
    raise FileNotFoundError("NIRSpec dispersion tables missing: run scripts/fetch_calibration.sh")
# Published instrument tables and msaexp default degree=2 lambda/R interpolation.
RCOEFS = {}
for grating in ["G395M", "G395H", "G235M"]:
    tab = Table.read(CALIB / f"jwst_nirspec_{grating.lower()}_disp.fits")
    RCOEFS[grating] = np.polyfit(
        np.asarray(tab["WAVELENGTH"]), np.asarray(tab["WAVELENGTH"]) / np.asarray(tab["R"]), 2
    )
REST = {
    "Halpha": 0.6564697,
    "Hbeta": 0.4862738,
    "OI8446": 10000.0 / (88631.146 - 76794.978),
    "Pabeta": 1.28217,
    "Pagamma": 1.09412,
    "HeI10830": 1.08333,
}
# OI vacuum reference is computed from NIST upper/lower energies [cm^-1].
# Fine-structure components are unresolved here; the line centroid remains free.
MODELS = ["single_gaussian", "double_gaussian", "lorentzian", "exponential"]


def resolution(wave, grating="G395M"):
    return np.asarray(wave) / np.polyval(RCOEFS[grating], wave)


def gauss(v, center, width):
    s = width / SQ
    return np.exp(-0.5 * ((v - center) / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def profile(v, center, width, kind):
    if kind in ["single_gaussian", "double_gaussian"]:
        return gauss(v, center, width)
    if kind == "lorentzian":
        g = width / 2
        return g / np.pi / ((v - center) ** 2 + g * g)
    if kind == "exponential":
        ve = width / (2 * np.log(2))
        return np.exp(-abs(v - center) / ve) / (2 * ve)
    raise ValueError(kind)


def bin_edges(x):
    assert np.all(np.diff(x) > 0)
    return np.r_[x[0] - (x[1] - x[0]) / 2, (x[1:] + x[:-1]) / 2, x[-1] + (x[-1] - x[-2]) / 2]


def response(lo, hi, grid, center_um, lsf_scale=1.0, grating="G395M"):
    # For an intrinsic delta at v_j, integrate the Gaussian LSF into observed pixel i.
    sigma = C / resolution(center_um * (1 + grid / C), grating) / SQ * lsf_scale
    K = (
        (ndtr((hi[:, None] - grid) / sigma) - ndtr((lo[:, None] - grid) / sigma))
        * np.diff(grid)[0]
        / (hi - lo)[:, None]
    )
    K[K < 1e-12] = 0
    return csr_matrix(K)


def load_line(path, z, line, lsf_scale=1.0, window=None, mask_rest_um=(), mask_halfwidth_kms=0.0):
    with fits.open(path) as h:
        t = h["SPEC1D"]
        a = t.data
        wave = np.array(a["wave"], dtype=float)
        f = np.array(a["flux"], dtype=float)
        e = np.array(a["err"], dtype=float)
        meta = t.header
        assert (
            t.columns["wave"].unit == "um"
            and t.columns["flux"].unit == t.columns["err"].unit == "uJy"
        )
        grating = meta["GRATING"]
        assert grating in RCOEFS
        if "POLY0" in meta:
            coefficients = [meta[f"POLY{i}"] for i in range(10) if f"POLY{i}" in meta]
            corr = np.polyval(coefficients, np.log(wave * 10000.0))
            assert np.all(np.isfinite(corr)) and np.all(
                corr > 0
            ), "Invalid published polynomial calibration"
            f *= corr
            e *= corr
        escale = np.array(a["escale"], dtype=float) if "escale" in a.names else np.ones_like(e)
        assert np.all(np.isfinite(escale)) and np.all(
            escale > 0
        ), "Invalid published uncertainty scale"
        # read_spectrum uses published flux directly; path_corr is NOT multiplied a second time.
        valid = np.isfinite(f + e) & (e > 0) & (f != 0) & (np.asarray(a["npix"]) > 0)
        if valid.sum() < 10:
            return dict(
                status="not_covered_or_insufficient_valid_pixels",
                line=line,
                z=z,
                Nvalid=int(valid.sum()),
            )
        minerr = np.percentile(e[valid], 5) * 0.5
        valid &= e > minerr
        med = median_filter(e[valid], size=11)
        ref = np.interp(wave, wave[valid], med)
        valid &= e > 0.2 * ref
        err = np.sqrt((e * escale) ** 2 + (0.02 * f) ** 2)
    center = REST[line] * (1 + z)
    v = C * (wave / center - 1)
    ve = bin_edges(v)
    # Default half-widths are those of the paper; window/mask_rest_um only serve the fit-window sensitivity runs.
    window = (4500.0 if line == "Hbeta" else 5000.0) if window is None else float(window)
    masked = np.zeros(len(v), dtype=bool)
    for rest in mask_rest_um:  # contaminating lines, vacuum rest wavelengths in um
        masked |= abs(v - C * (rest / REST[line] - 1)) <= mask_halfwidth_kms
    sel = valid & (abs(v) <= window) & ~masked
    ix = np.flatnonzero(sel)
    center_valid = wave[0] < center < wave[-1] and valid[np.argmin(abs(wave - center))]
    if not (center_valid and len(ix) >= 25 and (v[ix] < -1000).any() and (v[ix] > 1000).any()):
        return dict(
            status="not_covered_or_insufficient_valid_pixels",
            line=line,
            z=z,
            center_um=center,
            Nvalid=int(sel.sum()),
            center_valid=bool(center_valid),
        )
    # f_nu -> f_lambda -> f_v. Y is flux density per1000km/s in units of1e-18 CGS.
    factor = (
        (np.ones(len(wave)) * u.uJy).to_value(
            u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(wave * u.um)
        )
        * center
        * 10000.0
        / C
        * 1000
        / FUNIT
    )
    y = f[ix] * factor[ix]
    er = err[ix] * factor[ix]
    grid = np.arange(-15000.0, 15000.01, 10.0)
    return dict(
        status="ready",
        line=line,
        z=z,
        center_um=center,
        v=v[ix],
        lo=ve[ix],
        hi=ve[ix + 1],
        y=y,
        err=er,
        quoted_err=e[ix] * factor[ix],
        grid=grid,
        K=response(ve[ix], ve[ix + 1], grid, center, lsf_scale, grating),
        lsf_scale=lsf_scale,
        R_center=float(resolution(center, grating)),
        grating=grating,
        valid_fraction=float(sel.sum() / np.sum((abs(v) <= window) & ~masked)),
        window=window,
        Nvalid=len(ix),
        wavelength_um=wave[ix],
    )


def layout(line, kind, absorption):
    names = [
        "continuum",
        "continuum_slope",
        "F_narrow",
        "v_narrow",
        "W_narrow",
        "F_broad1",
        "v_broad1",
        "W_broad1",
    ]
    if kind == "double_gaussian":
        names += ["F_broad2", "v_broad2", "W_broad2"]
    if line == "Halpha":
        names += ["F_NII6585"]
    if line == "OI8446":
        names += ["F_CaII8500", "F_CaII8544", "W_CaII"]
    if line in ["Pagamma", "HeI10830"]:
        names += ["F_neighbor", "v_neighbor", "W_neighbor"]
    if absorption:
        names += ["tau_abs", "v_abs", "W_abs"]
    return names


def components(p, d, kind, absorption):
    t = dict(zip(layout(d["line"], kind, absorption), p))
    v = d["grid"]
    K = d["K"]
    b = 1000 * t["F_broad1"] * profile(v, t["v_broad1"], t["W_broad1"], kind)
    if kind == "double_gaussian":
        b += 1000 * t["F_broad2"] * gauss(v, t["v_broad2"], t["W_broad2"])
    n = 1000 * t["F_narrow"] * gauss(v, t["v_narrow"], t["W_narrow"])
    others = np.zeros_like(v)
    if d["line"] == "Halpha":
        for wave, rat in [(6549.86, 1 / 3), (6585.27, 1.0)]:
            center = C * (wave / (REST["Halpha"] * 10000.0) - 1) + t["v_narrow"]
            others += 1000 * t["F_NII6585"] * rat * gauss(v, center, t["W_narrow"])
    if d["line"] == "OI8446":
        # Nearby CaII are nuisance lines: no CaII line-ratio constraint. They share a free width.
        for wave, key in [(8500.35, "F_CaII8500"), (8544.44, "F_CaII8544")]:
            center = C * (wave / (REST["OI8446"] * 10000.0) - 1) + t["v_narrow"]
            others += 1000 * t[key] * gauss(v, center, t["W_CaII"])
    if d["line"] in ["Pagamma", "HeI10830"]:
        neighbor = "HeI10830" if d["line"] == "Pagamma" else "Pagamma"
        center = C * (REST[neighbor] / REST[d["line"]] - 1) + t["v_neighbor"]
        others += 1000 * t["F_neighbor"] * gauss(v, center, t["W_neighbor"])
    # Multiplicative screen on target-line emission. It cannot create negative intrinsic line emission.
    if absorption:
        tau = t["tau_abs"] * np.exp(-0.5 * ((v - t["v_abs"]) / (t["W_abs"] / SQ)) ** 2)
        transmission = np.exp(-tau)
        a = (n + b) * (1 - transmission)
    else:
        a = np.zeros_like(v)
    assert np.all(n + b - a >= -1e-10)
    cont = t["continuum"] + t["continuum_slope"] * d["v"] / 5000
    return dict(
        total=cont + K @ (n + b + others - a),
        continuum=cont,
        narrow=K @ n,
        broad=K @ b,
        nuisance=K @ others,
        absorption=K @ a,
        intrinsic_broad=b,
    )


def setup(d, kind, absorption, seed):
    v, y, e = (d["v"], d["y"], d["err"])
    edge = abs(v) > 3500
    cont = np.median(y[edge]) if edge.sum() > 5 else np.median(y)
    area = max(float(np.sum(np.maximum(y - cont, 0) * (d["hi"] - d["lo"])) / 1000), 0.05)
    maxF = max(100.0, area * 20)
    yscale = max(np.max(abs(y)), 1.0)
    names = layout(d["line"], kind, absorption)
    pars = {}
    pars.update(
        continuum=(cont, -3 * yscale, 3 * yscale),
        continuum_slope=(0.0, -6 * yscale, 6 * yscale),
        F_narrow=(area * 0.5, 0, maxF),
        v_narrow=(float(v[abs(v) < 600][np.argmax(y[abs(v) < 600])]), -600, 600),
        W_narrow=(150.0, 30.0, 600.0),
        F_broad1=(area * 0.5, 0, maxF),
        v_broad1=(0.0, -1500.0, 1500.0),
        W_broad1=([700.0, 1800.0, 4000.0][seed % 3], 350.0, 12000.0),
    )
    if kind == "double_gaussian":
        pars.update(
            F_broad1=(area * 0.35, 0, maxF),
            F_broad2=(area * 0.15, 0, maxF),
            v_broad2=(200.0, -1500.0, 1500.0),
            W_broad2=([2500.0, 4500.0, 1000.0][seed % 3], 350.0, 12000.0),
        )
    if d["line"] == "Halpha":
        pars["F_NII6585"] = (area * 0.03, 0, maxF)
    if d["line"] == "OI8446":
        pars.update(
            F_CaII8500=(area * 0.1, 0, maxF),
            F_CaII8544=(area * 0.1, 0, maxF),
            W_CaII=(1000.0, 100.0, 8000.0),
        )
    if d["line"] in ["Pagamma", "HeI10830"]:
        pars.update(
            F_neighbor=(area * 0.1, 0, maxF),
            v_neighbor=(0.0, -1500.0, 1500.0),
            W_neighbor=(1000.0, 100.0, 8000.0),
        )
    if absorption:
        pars.update(
            tau_abs=(0.5, 0, 10.0),
            v_abs=([0.0, -300.0, 300.0][seed % 3], -1500.0, 1500.0),
            W_abs=(250.0, 50.0, 1500.0),
        )
    x = np.array([pars[k][0] for k in names])
    lo = np.array([pars[k][1] for k in names])
    hi = np.array([pars[k][2] for k in names])
    return (names, x, lo, hi)


def fit_model(d, kind, absorption=False, seeds=3, start=None, y_override=None):
    names, _, lo, hi = setup(d, kind, absorption, 0)
    y = d["y"] if y_override is None else y_override

    def fun(p):
        return (components(p, d, kind, absorption)["total"] - y) / d["err"]

    best = None
    attempts = []
    for seed in range(seeds):
        _, x, _, _ = setup(d, kind, absorption, seed)
        if start is not None:
            x = np.array(start).copy()
        x = np.clip(x, lo + 1e-08, hi - 1e-08)
        f = least_squares(
            fun,
            x,
            bounds=(lo, hi),
            x_scale="jac",
            max_nfev=1200,
            ftol=2e-07,
            xtol=2e-07,
            gtol=2e-07,
        )
        attempts.append(
            dict(
                seed=seed,
                chi2=float(np.sum(f.fun * f.fun)),
                success=bool(f.success),
                nfev=int(f.nfev),
            )
        )
        if best is None or np.sum(f.fun * f.fun) < np.sum(best.fun * best.fun):
            best = f
    p = best.x
    chi2 = float(np.sum(best.fun**2))
    k = len(p)
    N = len(y)
    cov = np.linalg.pinv(best.jac.T @ best.jac, rcond=1e-12)
    errors = np.sqrt(np.maximum(np.diag(cov), 0))
    broad_idx = [names.index(n) for n in names if n in ["F_broad1", "F_broad2"]]
    Fb = float(p[broad_idx].sum())
    Fb_err = float(np.sqrt(np.sum(cov[np.ix_(broad_idx, broad_idx)])))
    bound = [
        names[i]
        for i in range(k)
        if min((p[i] - lo[i]) / (hi[i] - lo[i]), (hi[i] - p[i]) / (hi[i] - lo[i])) < 0.0001
    ]
    return dict(
        kind=kind,
        absorption=absorption,
        names=names,
        params=p.tolist(),
        errors_conditional=errors.tolist(),
        covariance_conditional=cov.tolist(),
        chi2=chi2,
        N=N,
        k=k,
        BIC=chi2 + k * np.log(N),
        reduced_chi2=chi2 / (N - k),
        success=bool(best.success),
        jacobian_rank=int(np.linalg.matrix_rank(best.jac)),
        bound_parameters=bound,
        broad_flux_erg_s_cm2=Fb * FUNIT,
        broad_flux_err_conditional_erg_s_cm2=Fb_err * FUNIT,
        broad_SNR_conditional=Fb / Fb_err if Fb_err > 0 else None,
        total_emitted_target_flux_erg_s_cm2=(Fb + p[names.index("F_narrow")]) * FUNIT,
        absorbed_flux_erg_s_cm2=(
            float(
                np.trapezoid(components(p, d, kind, absorption)["absorption"], d["v"])
                / 1000
                * FUNIT
            )
            if absorption
            else 0.0
        ),
        optimizer_attempts=attempts,
    )


def widths(v, y):
    # All moments/quantiles refer to the explicitly supplied finite velocity window.
    y = np.asarray(y)
    v = np.asarray(v)
    if np.any(y < 0) or not np.isfinite(y).all() or np.trapezoid(y, v) <= 0:
        return None
    area = np.trapezoid(y, v)
    pdf = y / area
    cumulative = np.r_[0, np.cumsum((pdf[:-1] + pdf[1:]) / 2 * np.diff(v))]
    cumulative /= cumulative[-1]
    quant = np.interp([0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9], cumulative, v)
    mean = np.trapezoid(v * pdf, v)
    sig = np.sqrt(max(np.trapezoid((v - mean) ** 2 * pdf, v), 0))
    ip = int(np.argmax(y))
    out = {}
    for fraction, name in [(0.5, "FWHM_kms"), (0.2, "FW20peak_kms")]:
        above = np.where(y >= fraction * y[ip])[0]
        if above[0] == 0 or above[-1] == len(y) - 1:
            out[name] = None
            continue
        i, j = (above[0], above[-1])
        left = np.interp(fraction * y[ip], [y[i - 1], y[i]], [v[i - 1], v[i]])
        right = np.interp(fraction * y[ip], [y[j + 1], y[j]], [v[j + 1], v[j]])
        out[name] = float(right - left)
    out.update(
        W50_kms=float(quant[5] - quant[1]),
        W20_kms=float(quant[4] - quant[2]),
        sigma_kms=float(sig),
        v50_kms=float(quant[3]),
        vpeak_kms=float(v[ip]),
        flux_window_erg_s_cm2=float(area / 1000 * FUNIT),
        velocity_window_kms=[float(v[0]), float(v[-1])],
    )
    return out


def nonparametric(d, fit, nboot=80, seed=0):
    """Nonnegative piecewise-linear intrinsic profile, with 5fold CV smoothness.
    Shape does not assume Gaussian/Lorentz/exponential; narrow/continuum/absorption
    subtraction still depends on the selected fit. Bootstrap errors are conditional.
    """
    comp = components(fit["params"], d, fit["kind"], fit["absorption"])
    res = d["y"] - comp["continuum"] - comp["narrow"] - comp["nuisance"] + comp["absorption"]
    knots = np.arange(-5000.0, 5000.1, 125.0)
    B = np.maximum(1 - np.abs(d["grid"][:, None] - knots) / 125.0, 0)
    A = np.asarray(d["K"] @ B)
    D2 = np.diff(np.eye(len(knots)), n=2, axis=0)
    Aw = A / d["err"][:, None]
    b = res / d["err"]
    norm = np.sum(Aw**2) / np.sum(D2**2)
    alphas = np.array([0.001, 0.01, 0.1, 1, 10, 100.0, 1000.0, 10000.0]) * norm
    cv = []
    for alpha in alphas:
        losses = []
        for fold in range(5):
            train = np.arange(len(b)) % 5 != fold
            val = ~train
            mat = np.vstack([Aw[train], np.sqrt(alpha) * D2])
            rhs = np.r_[b[train], np.zeros(len(D2))]
            sol = nnls(mat, rhs, maxiter=2000)[0]
            losses.extend((Aw[val] @ sol - b[val]) ** 2)
        cv.append(float(np.mean(losses)))
    ai = int(np.argmin(cv))
    alpha = alphas[ai]
    mat = np.vstack([Aw, np.sqrt(alpha) * D2])
    rhs = np.r_[b, np.zeros(len(D2))]
    sol = nnls(mat, rhs, maxiter=3000)[0]
    fine = np.arange(-5000.0, 5000.1, 10.0)
    shape = np.interp(fine, knots, sol)
    m = widths(fine, shape)
    # Fixed-alpha, fixed subtraction bootstrap describes only conditional noise, not all systematics.
    rng = np.random.default_rng(seed)
    boot = []
    shapes = []
    pred = A @ sol
    for _ in range(nboot):
        yy = pred + rng.normal(0, d["err"])
        x = nnls(mat, np.r_[yy / d["err"], np.zeros(len(D2))], maxiter=3000)[0]
        s = np.interp(fine, knots, x)
        z = widths(fine, s)
        if z:
            boot.append(z)
            shapes.append(s)
    intervals = {}
    for key in [
        "FWHM_kms",
        "W50_kms",
        "W20_kms",
        "FW20peak_kms",
        "sigma_kms",
        "flux_window_erg_s_cm2",
    ]:
        vals = [b[key] for b in boot if b.get(key) is not None]
        intervals[key] = (
            dict(N=len(vals), p16_p50_p84=np.percentile(vals, [16, 50, 84]).tolist())
            if vals
            else None
        )
    # A descriptive exponential wing slope, not a model-free fundamental moment.
    stderr = (
        np.std(np.array(shapes), axis=0, ddof=1)
        if len(shapes) > 1
        else np.ones_like(shape) * np.inf
    )
    dv = fine - m["vpeak_kms"]
    inner = max(300.0, (m["FWHM_kms"] or 2000.0) / 2)
    tail = (abs(dv) >= inner) & (abs(dv) <= 4000) & (shape > 3 * stderr) & (shape > 0)
    ve = None
    if np.sum(tail & (dv < 0)) >= 4 and np.sum(tail & (dv > 0)) >= 4:
        xx = np.column_stack([np.ones(tail.sum()), (dv[tail] > 0).astype(float), -abs(dv[tail])])
        ww = np.maximum(shape[tail] / np.maximum(stderr[tail], 1e-12), 1)
        beta = np.linalg.lstsq(xx * ww[:, None], np.log(shape[tail]) * ww, rcond=None)[0]
        if beta[2] > 0:
            ve = float(1 / beta[2])
    m["ve_wing_regression_kms"] = ve
    observed_positive = np.maximum(res, 0)
    obs = widths(d["v"], observed_positive)
    # Alternative regularization values are retained; these are sensitivities, not CI endpoints.
    regularization = []
    for i in sorted(set([max(ai - 1, 0), ai, min(ai + 1, len(alphas) - 1)])):
        x = nnls(np.vstack([Aw, np.sqrt(alphas[i]) * D2]), rhs, maxiter=3000)[0]
        regularization.append(
            dict(alpha=float(alphas[i]), metrics=widths(fine, np.interp(fine, knots, x)))
        )
    return dict(
        metrics=m,
        conditional_noise_intervals=intervals,
        observed_positive_part_profile=obs,
        definition="Cumulative W50=v75-v25; W20=v60-v40; FW20peak is20%-peak full width. Intrinsic shape is nonnegative piecewise-linear, forward-convolved with assumed Gaussian LSF.",
        subtraction="Selected parametric continuum+narrow+nuisance removed; fitted absorption added back. Both assumptions propagate as systematics, not in conditional noise interval.",
        moment_window="-5000 to+5000 km/s; all moments and quantiles truncated to this window",
        cv_alpha_grid=alphas.tolist(),
        cv_loss=cv,
        selected_alpha_index=ai,
        regularization_sensitivity=regularization,
        Nbootstrap=nboot,
        ve_definition="Common exponential log-slope on both recovered wings, separate intercepts, |v-vpeak| from max(FWHM/2,300) to4000km/s; only bins>3conditional-noise. Null if insufficient bins or non-decreasing wings.",
        knots_kms=knots.tolist(),
        intrinsic_knots=sol.tolist(),
        v_fine=fine.tolist(),
        profile_fine=shape.tolist(),
        observed_residual=res.tolist(),
        nonparametric_model_on_pixels=(A @ sol).tolist(),
        noise_scope="Independent Gaussian pixel noise; fixed selected regularization, narrow/continuum and absorption. Does not include extraction covariance or their parameter errors.",
    )


def self_checks():
    v = np.arange(-10000.0, 10000.01, 1.0)
    sigma = 500.0
    y = 1000 * gauss(v, 0, SQ * sigma)
    m = widths(v, y)
    assert np.isclose(np.trapezoid(y, v) / 1000, 1.0, rtol=1e-08)
    assert abs(m["FWHM_kms"] / (SQ * sigma) - 1) < 1e-05
    assert abs(m["W50_kms"] / (2 * ndtri(0.75) * sigma) - 1) < 1e-05
    assert abs(m["W20_kms"] / (2 * ndtri(0.6) * sigma) - 1) < 1e-05
    assert abs(m["FW20peak_kms"] / (2 * np.sqrt(2 * np.log(5)) * sigma) - 1) < 1e-05
    assert abs(m["sigma_kms"] / sigma - 1) < 1e-06
    assert np.isclose((1 * u.uJy).to_value(u.erg / u.s / u.cm**2 / u.Hz), 1e-29)
    assert 0.8448 < REST["OI8446"] < 0.845
    grid = np.arange(-15000.0, 15000.01, 10.0)
    cent = np.arange(-5000.0, 5000.1, 80.0)
    edges = bin_edges(cent)
    K = response(edges[:-1], edges[1:], grid, 4.0)
    intrinsic = 1000 * gauss(grid, 0, 1000.0)
    out = K @ intrinsic
    assert abs(np.sum(out * np.diff(edges)) / 1000 - 1) < 1e-06
    # Variances add for a Gaussian line and Gaussian LSF; pixel averaging contributes dv^2/12.
    measured = widths(cent, out)
    expected = np.sqrt((1000 / SQ) ** 2 + (C / resolution(4.0) / SQ) ** 2 + 80**2 / 12)
    assert abs(measured["sigma_kms"] / expected - 1) < 0.002
    return dict(
        gaussian_metrics="passed",
        response_flux_and_variance="passed",
        unit_conversion="passed",
        OI_vacuum_um=REST["OI8446"],
    )
