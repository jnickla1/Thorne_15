"""Read-only copy of the 'hh' calibration in calibrate_historical_combinations.py.

The original writes Results2/comb_*_scales.npy, so it is not run here. This copy reads
the same Results2/historical_* inputs and writes only
data/supplementary/dashboard_combination_calibration.json:

- the 7-method PIVW and equal-weight scales, which must reproduce Results2/comb_real_scales.npy;
- an 8-method PIVW (the 7 plus CGWL10y_for_halfU), with its own constant scale;
- the 7-method CSCM stacking weights and sharpening factor (nsharp).

Since the 2026-10-19 fix, Results2/historical_newsamples.npy stores sampled methods as
c + dev*scale (matching the likelihoods and the dashboard's se*scale widths), so no
rescaling is applied here. nsharp uses the same bounded scalar search, and the infilling
uses the same seed (20261019 + 7), as calibrate_historical_combinations.py; the stacked
draw is also seeded here.

python source/calibrate_dashboard_combinations.py --source-root /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import norm, gaussian_kde

REAL7 = ['EBMKF_ta4', 'GAM_AR1', 'lowess1dt36wnc', 'Kal_flexLin', 'FaIR_comb_unB', 'GWI_tot_CGWL', 'CGWL10y_sfUKCP']
REAL8 = REAL7 + ['CGWL10y_for_halfU']


def kl_scale(standard, standard_se, centers, lhund=-100):
    def rescale_KL(scale):
        yearKL = np.log(scale/standard_se[lhund:-10]) + (standard_se[lhund:-10]**2 + (standard[lhund:-10]-centers[lhund:-10])**2)/2/scale**2 - 0.5
        return np.nansum(yearKL)
    return float(minimize(rescale_KL, x0=0.04, bounds=[(0.001, 1)]).x[0])


def sharpen_samples(blended_in, nsharp, seed):
    """Copied from calibrate_historical_combinations.py."""
    nyears, nsamples = blended_in.shape
    k = int(np.floor(nsharp)); frac = float(nsharp - k)
    draws = k + (1 if frac > 0 else 0)
    rng = np.random.default_rng(seed)
    out = np.empty_like(blended_in)
    for y in range(nyears):
        idx = rng.integers(0, nsamples, size=(nsamples, draws))
        picks = blended_in[y, idx]
        if frac > 0:
            num = picks[:, :k].sum(axis=1) + frac*picks[:, -1]; den = k + frac
        else:
            num = picks.sum(axis=1); den = k if k > 0 else 1.0
        out[y] = num/den
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path('/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats'))
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1]/'data/supplementary/dashboard_combination_calibration.json')
    args = parser.parse_args()
    root = args.source_root
    sys.path.insert(0, str(root))
    from custom_stack_utils import predict_non_nan, infill_log_likelihoods, fit_stacking_weights_logp

    df = pd.read_csv(root/'Results2/historical_names_var_scale.csv', index_col=0)
    all_methods = df['method_name']
    err_vars100 = df['err_var100'].to_numpy()
    central_est = np.load(root/'Results2/historical_central_est.npy')
    node_lls = np.load(root/'Results2/historical_hermguass_lls.npy')
    newsamples = np.load(root/'Results2/historical_newsamples.npy')
    standard = np.load(root/'Results2/historical_standard.npy')
    standard_se = np.load(root/'Results2/historical_standard_se.npy')
    stored = np.load(root/'Results2/comb_real_scales.npy')
    out = dict(source='Read-only copy of the hh branch of calibrate_historical_combinations.py, using Results2/historical_* (HadCRUT5.1 observations).',
               stored_comb_real_scales=dict(equal_weight=float(stored[0]), pivw=float(stored[1])))

    for label, chosen in [('real7', REAL7), ('real8', REAL8)]:
        # As in the original, methods keep the file order of historical_names_var_scale.csv.
        mask = all_methods.isin(chosen).to_numpy()
        names = all_methods[mask].tolist()
        iv = 1/err_vars100[mask]
        c = central_est[mask]
        ivcenters = np.nansum(c*iv[:, None], axis=0)/np.nansum(~np.isnan(c)*iv[:, None], axis=0)
        cwcenters = np.nansum(c, axis=0)/np.nansum(~np.isnan(c), axis=0)
        entry = dict(methods=names, pivw_weights=dict(zip(names, (iv/iv.sum()).tolist())),
                     pivw_scale=kl_scale(standard, standard_se, ivcenters),
                     equal_weight_scale=kl_scale(standard, standard_se, cwcenters))
        if label == 'real7':
            np.testing.assert_allclose([entry['equal_weight_scale'], entry['pivw_scale']], stored, rtol=1e-4,
                                       err_msg='7-method scales do not reproduce Results2/comb_real_scales.npy')
            x, w = hermgauss(100)
            z = np.sqrt(2.0)*x; alpha = w/np.sqrt(np.pi)
            np.random.seed(20261019 + len(names))  # as calibrate_historical_combinations.py
            infilled = infill_log_likelihoods(node_lls[mask][:, 10:-10, :], penalty=0)
            stack = fit_stacking_weights_logp(infilled, alpha)
            samples = newsamples[mask]

            def kl_total(sharp):
                total = 0.0
                for y in range(len(standard)-100, len(standard)-10):
                    if np.isnan(standard[y]):
                        continue
                    mu, sd = float(standard[y]), float(max(standard_se[y], 1e-12))
                    nodes = mu + sd*z
                    log_q = np.log(np.maximum(gaussian_kde(sharp[y]).evaluate(nodes), 1e-300))
                    total += float(np.sum(alpha*(norm.logpdf(nodes, loc=mu, scale=sd) - log_q)))
                return total
            entry['cscm_stack_weights'] = dict(zip(names, stack.tolist()))
            blended = predict_non_nan(stack, samples, seed=20261019)
            nsharp = float(minimize_scalar(lambda s: kl_total(sharpen_samples(blended, s, seed=42)), bounds=(1, 10), method='bounded', options={'xatol': 0.01}).x)
            sharp = sharpen_samples(blended, nsharp, seed=2)
            entry['cscm'] = dict(nsharp=nsharp, kl=kl_total(sharp),
                rmse_last100=float(np.sqrt(np.nanmean((sharp.mean(axis=1)[-100:-1]-standard[-100:-1])**2))),
                mean_sd_last100=float(np.nanmean(sharp.std(axis=1)[-100:-10])))
            entry['cscm_nsharp'] = nsharp
        entry['pivw_rmse_last100'] = float(np.sqrt(np.nanmean((ivcenters[-100:-1]-standard[-100:-1])**2)))
        out[label] = entry
        print(label, json.dumps({k: v for k, v in entry.items() if k != 'methods'}, indent=1))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
