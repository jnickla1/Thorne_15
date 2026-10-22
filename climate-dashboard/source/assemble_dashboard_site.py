"""Copy validated figure/data exports into the portable website.

Run after export_dashboard_climate_metrics.py and render_dashboard_figures.py.
No source-project writes, deployment, or credentials are involved.
"""
import csv
import json
import shutil
from export_dashboard_rate_history import export_rates
from pathlib import Path

CSV_NOTES = """# Dashboard temperature series are displayed from 1965; earlier data rows are retained.
# CAUTION: pre-1965 smoothed anthropogenic estimates are not directly comparable to later estimates.
# FaIR applies a prescribed mean during 1930-1964 and contracts current samples; GWI is unavailable before 1950.
# Treat those early anthropogenic method outputs separately; the pinch does not demonstrate low physical uncertainty.
# This warning does not invalidate annual observational data or other realised methods.
# These metadata lines are comments: use pandas.read_csv(path, comment="#") or skip leading # lines.
"""

def copy_public_csv(source, destination):
    """Add download metadata without changing analysis CSVs or numeric rows."""
    destination.write_bytes(CSV_NOTES.encode() + source.read_bytes())

def main():
    base = Path(__file__).resolve().parents[1]
    for baseline in ['rpi', 'satcal']:
        for goal in ['real', 'anthro']:
            export_rates(base / 'data' / f'{baseline}_{goal}.json', base / 'data')
    assets = base / 'site/dist/assets'
    for name in ['data', 'figures']:
        (assets / name).mkdir(parents=True, exist_ok=True)
    for path in (base / 'data').iterdir():
        if path.suffix in {'.json', '.csv'} and 'legacy' not in path.name:
            copy_public_csv(path, assets / 'data' / path.name) if path.suffix == '.csv' else shutil.copy2(path, assets / 'data' / path.name)
    for path in (base / 'figures').iterdir():
        if path.suffix in {'.json', '.png', '.pdf', '.svg'}:
            shutil.copy2(path, assets / 'figures' / path.name)
    count = 0
    for baseline in ['rpi', 'satcal']:
        for goal in ['real', 'anthro']:
            prefix = f'{baseline}_{goal}'
            data = json.loads((base / 'data' / f'{prefix}.json').read_text())
            for method, values in data['methods'].items():
                path = base / 'data' / f'{prefix}__{method}__current.csv'
                with path.open('w', newline='') as output:
                    writer = csv.writer(output)
                    writer.writerow(['year'] + [f'member_{i}' for i in range(len(values['current']))])
                    for i, year in enumerate(data['years']):
                        writer.writerow([year] + [row[i] if row[i] is not None else '' for row in values['current']])
                copy_public_csv(path, assets / 'data' / path.name)
                count += 1
    manifest = json.loads((base / 'figures/manifest.json').read_text())
    # Three temperature-series options for every method in every dataset
    # (78 for the main options, 96 with the three supplementary realised options).
    expected = 3 * count
    if len(manifest) != expected or not all((assets / 'figures' / item['file']).is_file() for item in manifest):
        raise ValueError(f'The complete {expected}-figure library (three temperature-series options) is required.')
    print(f'Website ready: {expected} variants, default PNG/PDF/SVG, {count} method CSVs.')
    # Refresh the translation bundle (English page strings and data-derived method text).
    import sys
    from build_i18n import main as build_i18n
    argv, sys.argv = sys.argv, [sys.argv[0]]
    try:
        if build_i18n():
            raise ValueError('Translation check failed; see build_i18n.py output above.')
    finally:
        sys.argv = argv

if __name__ == '__main__':
    main()
