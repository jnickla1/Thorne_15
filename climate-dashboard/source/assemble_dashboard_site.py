"""Copy validated figure/data exports into the portable website.

Run after export_dashboard_climate_metrics.py and render_dashboard_figures.py.
No source-project writes, deployment, or credentials are involved.
"""
import csv
import json
import shutil
from pathlib import Path

def main():
    base = Path(__file__).resolve().parents[1]
    assets = base / 'site/dist/assets'
    for name in ['data', 'figures']:
        (assets / name).mkdir(parents=True, exist_ok=True)
    for path in (base / 'data').iterdir():
        if path.suffix in {'.json', '.csv'} and 'legacy' not in path.name:
            shutil.copy2(path, assets / 'data' / path.name)
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
                shutil.copy2(path, assets / 'data' / path.name)
                count += 1
    manifest = json.loads((base / 'figures/manifest.json').read_text())
    if len(manifest) != 52 or not all((assets / 'figures' / item['file']).is_file() for item in manifest):
        raise ValueError('The complete 52-figure library is required.')
    print(f'Website ready: 52 variants, default PNG/PDF/SVG, {count} method CSVs.')

if __name__ == '__main__':
    main()
