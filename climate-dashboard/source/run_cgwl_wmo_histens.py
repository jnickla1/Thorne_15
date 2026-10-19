"""Run CGWL10y_for_halfU (CGWL-WMOLC-ADCP) on the 100 heads/tails observational members.

The method file is imported unchanged. It returns blanks for any experiment other
than 'historical', because the WMO LC-ADCP forecasts describe the real world. Each
heads/tails member is a version of that real-world record, so this wrapper calls it
with experiment_type='historical'. The forecasts are converted to 1850-1900 with the
method's fixed 0.88 °C offset. Like CGWL10y_sfUKCP on histens members, the member's
own 1850-1899 mean is removed before the call and added back afterwards, which keeps
the 1981-2010 (satcal) members on their own baseline.

Writes data/supplementary/CGWL10y_for_halfU_histens_{baseline}.npz containing the
current mean and SE arrays (members x years). The method has no retrospective output.

Run in cleanPy from the climate-dashboard directory:
python source/run_cgwl_wmo_histens.py --source-root /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
"""
import argparse
import importlib.util
import contextlib
import io
from pathlib import Path
import numpy as np
import pandas as pd

METHOD = 'Methods/44_EarthModel_CGWL/CGWL10y_for_halfU_method.py'


def load_method(root):
    spec = importlib.util.spec_from_file_location('cgwl_for_halfU', root/METHOD)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_method


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path('/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats'))
    parser.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parents[1]/'data/supplementary')
    args = parser.parse_args()
    root = args.source_root
    args.outdir.mkdir(parents=True, exist_ok=True)
    run_method = load_method(root)
    hadcrut = pd.read_csv(root/'Common_Data/HadCRUT.5.1.csv')
    obs, upper, lower = (hadcrut[c].to_numpy() for c in ['Anomaly', 'Upper', 'Lower'])
    for baseline, name in [('rpi', 'sst_pseudo.csv'), ('satcal', 'sst_pseudo_1981-2010.csv')]:
        a = np.loadtxt(root/'Hist_heads_tails'/name, delimiter=',')
        years = a[:, 0].astype(int)
        assert years[0] == 1850, 'the forecast lookup assumes the record starts in 1850'
        sims = a[:, 1:]
        means = np.full((sims.shape[1], len(years)), np.nan); ses = means.copy()
        for member in range(sims.shape[1]):
            sim = sims[:, member]
            # Same uncertainty construction as evaluate_heads_tails_ensemble.py.
            unc = (sim + lower[:len(years)] - obs[:len(years)], sim + upper[:len(years)] - obs[:len(years)])
            offset = np.mean(sim[:50])
            with contextlib.redirect_stdout(io.StringIO()):
                result = run_method(years, sim - offset, unc, member, 'historical')
            means[member] = result['mean'](years, 0) + offset
            ses[member] = result['se'](years, 0)
        np.savez_compressed(args.outdir/f'CGWL10y_for_halfU_histens_{baseline}.npz', years=years, current=means, current_se=ses)
        print(baseline, 'members', sims.shape[1], '2025 mean', np.nanmean(means[:, -1]).round(4),
              'first finite year', years[np.isfinite(means).any(axis=0)][0])


if __name__ == '__main__':
    main()
