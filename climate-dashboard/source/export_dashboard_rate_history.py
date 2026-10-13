"""Export optional rate histories from the existing validated dashboard JSON.

Run after export_dashboard_climate_metrics.py and before assemble_dashboard_site.py.
No method reruns or changes to headline estimates are made.
"""
import csv
import json
from pathlib import Path
import numpy as np


def fitted_rates(matrix, window, degree, endpoint_offset):
    """Require a complete window; return the ensemble mean endpoint derivative."""
    a = np.asarray(matrix, dtype=float)
    out = np.full(a.shape[-1], np.nan)
    t = np.arange(window, dtype=float) - endpoint_offset
    design = np.column_stack([t ** p for p in range(degree, -1, -1)])
    derivative = np.linalg.pinv(design)[-2] * 10
    for end in range(window - 1, a.shape[-1]):
        values = a[:, end-window+1:end+1]
        valid = np.isfinite(values).all(axis=1)
        if valid.any():
            out[end-int(window-1-endpoint_offset)] = np.mean(values[valid] @ derivative)
    return out


def clean(values):
    return [float(v) if np.isfinite(v) else None for v in values]


def export_rates(source, output):
    d = json.loads(source.read_text())
    years = d['years']
    benchmark = fitted_rates(d['centered10'], 11, 1, 5)
    result = dict(source_sha256=d['metadata']['source_sha256'], years=years,
                  benchmark=clean(benchmark), methods={},
                  definitions=dict(benchmark='OLS slope through 11 centered 10-year observational means, year-5 to year+5; realized observations also on anthropogenic views. Requires observations through year+10.',
                      retrospective='Endpoint derivative of a trailing 15-year quadratic fit to retrospective method estimates from the full observation record; GWI current outputs are the existing explicit retrospective proxy. These are revised estimates, not vintage forecasts.',
                      current='Endpoint derivative of a trailing 15-year quadratic fit to annual current method outputs. This alternative may differ from the retrospective headline rate.',
                      units='degrees Celsius per decade; ensemble mean; missing where a complete fitting window is unavailable.'))
    for method, md in d['methods'].items():
        ret = fitted_rates(md['recent'], 15, 2, 14)
        cur = fitted_rates(md['current'], 15, 2, 14)
        assert abs(ret[-1]-md['headline_summary']['rate']) < 5e-7, (source.name,method,ret[-1],md['headline_summary']['rate'])
        result['methods'][method] = dict(retrospective=clean(ret), current=clean(cur))
        path = output / f'{source.stem}__{method}__rate_history.csv'
        with path.open('w', newline='') as f:
            f.write('# Warming rates in degrees Celsius per decade; ensemble mean.\n')
            for key in ['benchmark','retrospective','current']:
                f.write('# '+key+': '+result['definitions'][key]+'\n')
            f.write('# Display starts in 1965; earlier rows retained. Pre-1965 anthropogenic smoothing is method-constrained; FaIR prescribes its early mean and GWI is unavailable before 1950.\n')
            writer = csv.writer(f)
            writer.writerow(['year','benchmark_rate_C_per_decade','method_retrospective_rate_C_per_decade','method_current_rate_C_per_decade'])
            for i, year in enumerate(years):
                writer.writerow([year, result['benchmark'][i], result['methods'][method]['retrospective'][i], result['methods'][method]['current'][i]])
    (output / f'{source.stem}_rate_history.json').write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
    return dict(dataset=source.stem, methods=len(result['methods']), benchmark_last_year=max(y for y,v in zip(years,result['benchmark']) if v is not None))


def main():
    directory = Path(__file__).resolve().parents[1] / 'data'
    for baseline in ['rpi','satcal']:
        for goal in ['real','anthro']:
            print(export_rates(directory / f'{baseline}_{goal}.json', directory))


if __name__ == '__main__':
    main()
