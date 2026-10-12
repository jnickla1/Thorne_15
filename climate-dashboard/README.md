# Climate state dashboard — 2025

The original mockup has been updated, with 52 additional figures and a portable interactive website. The source project's `export_heads_tails_rate_summary.py` now reports current warming levels with retrospective rate fits and console-only endpoint checks, following the 12 October 2026 reporting correction. `plot_heads_tails_paper_probabilities.py` remains unchanged; archived originals and dashboard working copies are retained in `source/`.

## Review the outputs

- `figures/updated_mockup_2025.png`: the main 16:9 figure, with matching PDF and SVG exports.
- `figures/manifest.json`: all 52 combinations of warming definition, calibration, method, and raw/smoothed background.
- `site/dist/`: the complete static website, including local chart software and data. It supports zooming, hover values, image switching, downloads, and method-reference hyperlinks.
- `data/`: complete JSON exports, current-method matrices, summary and member-level crossing/rate CSVs, and probability CSVs. Use `*_corrected_headlines.csv` and JSON `headline_summary` for the final headlines. The explicitly named legacy exports and `data/rates_fixed/` retain intermediate audit calculations, rather than final dashboard numbers.
- `review/`: regression results, correction diffs, original method backups, and quartic-fit review material.
- `CHANGELOG.md`: the scientific changes and whether they affect current outputs. It records validation with the transferred cut2021 trajectories and the subsequent display refinements.

The nine realized options are PIVW, equal weighting, CGWL, EBM–Kalman, FaIR, GAM, GWI, flexible Kalman, and LOWESS. Anthropogenic views include PIVW, equal weighting, FaIR, and GWI. Both definitions have both calibrations and both background choices. The main mockup is additional to these 52 variants. The Background lines menu switches only the grey/black curves between raw observations and the selected method's smoothed ensemble; all other calculations remain unchanged.

## Scientific conventions

Historical definitive crossing times use the realized observational centered 20-year mean. The centered 10-year mean supplies the rate, with an OLS fit at crossing−5 through crossing+5 inclusive. The symmetric running means average the two half-year-shifted windows, retaining half weight on each outer observation. The centered 20-year record ends in 2015 for observations through 2025.

For the seven preindustrial members without an observed 1.0°C crossing, the last 25 finite centered 10-year means (1996–2020) supply a quartic fit. The first upward root after the final centered 20-year year estimates the crossing. These seven are identified as `quartic_estimate` in the member CSVs. Their model uncertainty is excluded from the reported range. No 1.5°C definitive crossing is projected. Ensemble member numbers are zero-based, 0–99; Oscar runs for them are r1–r100.

Dark crossing segments show the definitive rate, centered on the definitive crossing. Red comparison segments extend 15 years before the selected method's predicted crossing, with the rate from a quadratic endpoint fit to its past 15 annual smoothed outputs. The red 2025 segment spans 2010–2025 and is anchored at the current-method level and fitted recent rate.

The 2025 headline level uses current outputs. The headline rate uses corrected retrospective outputs for 2011–2025. GWI total and anthropogenic use current outputs as an explicitly marked temporary retrospective proxy. PIVW rates include within-method fit variance, between-method dispersion, and between-member dispersion. Rate labels show ±1.96 SE. Level intervals come from a mixture of calibrated member distributions. Crossing ranges condition on members that cross.

The saved paper probabilities are used for 0.5°C in 1970–2005, 1.0°C in 1990–2025, and 1.5°C in 2005–2025. Outside those windows, and for 2.0°C, the display uses a Gaussian approximation from the calibrated current means and uncertainties. The bottom curves reproduce the paper's cumulative-maximum and one-year Gaussian-smoothed density proxy; they are not a joint first-passage distribution. Bars independently count member first crossings, divided by all 100 members. All three panels use the same 0–0.3 probability/year y-axis across every method, baseline, and warming definition. Bars that exceed 0.3 are clipped visibly and marked with `//` plus their full height; hover values retain that height. Their colors and opacity come directly from the paper plotting scripts. Smoothed density curves use red. An explicit magnify switch expands only the rightmost panel's y-axis to a finer labeled scale, preserving its probability mass. The top-left warming plot starts at 1960.

Satellite views use anomalies relative to 1981–2010, with thresholds shifted by 0.7°C according to the source pipeline. Anthropogenic views retain the realized observational definitive benchmark as a comparison, and identify it accordingly. LOWESS is labeled linear without asserting a Gaussian kernel that the local method does not use.

## FaIR correction status

All 101 transferred all-forcing cut2021 matrices are present under the user-selected `resliced_headstails25/retrospective_4/` folder. Their raw row 272 matches the existing combined array's 2021 row exactly, including the valid column prefix and NaN padding.

The corrected combined retrospective branch takes its lagged matrix from cut2021, applies the same component biases as the current branch, and uses normalized coefficients 0.6125, 0.25, and 0.1375. The current pathway is preserved. Across 800 FaIR method/member/calibration cases, complete current sample arrays, means, SEs, threshold probabilities, likelihoods, and seeded resamples match the original backups exactly. All selected current arrays and both current combinations are unchanged in the 400 result files. Current and retrospective combination weights normalize to one where data exist.

FaIR combined preindustrial 2025 mean is 1.367695715°C current versus 1.367686220°C retrospective, leaving approximately 0.000009495°C from the remaining distribution/noise definitions. The earlier approximately 0.0023°C lag discrepancy is resolved. CGWL retains its separate 10-observation/10-forecast current versus 11-observation/10-forecast retrospective window definition; FaIR alignment does not remove that distinction.

The dashboard reads 400 refreshed retrospective result copies in `data/corrected_results/`. The original source project's cached `Results` files remain byte-for-byte unchanged. No rerun of the current-output pipeline is required for these tested cases. Other scripts that read the original caches still require retrospective cache adoption or regeneration to see the corrected histories.

## Host the website

Unzip `climate-dashboard-website.zip`. Its `index.html`, `app.js`, `styles.css`, `vendor/`, and `assets/` are a complete static site. No build, Python server computation, external CDN, or account connector is required.

Recommended public host: **GitHub Pages**, in a dedicated dashboard repository containing only the website exports. The interactive charts and downloads already work without server-side Python. GitHub Pages supplies static hosting and HTTPS; the lab server is an option if future features need on-demand Python computation or institution-managed access. See [GitHub Pages documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages).

For GitHub Pages, put those files at the publishing root (or in the configured `docs/` directory) and enable Pages for that branch/directory. Relative asset URLs support repository subpaths. A `.nojekyll` file is included. For the lab server, copy the same files into its web-served directory. No hostname-specific code changes are needed.

For a local preview, serve the files rather than opening `index.html` directly:

```sh
cd '/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/climate-dashboard/site/dist'
python3 -m http.server 8787
```

Then open `http://localhost:8787`. The URL preserves the selected warming definition, calibration, method, background, view, and right-panel magnification. Python computes the scientific results locally; the browser renders those exports interactively.

## Reproduce the figures

Use the existing `cleanPy` environment. From this dashboard directory:

```sh
conda activate cleanPy
python3 source/audit_and_rebuild_retrospective_caches.py
python3 source/export_dashboard_climate_metrics.py all --results-dir data/corrected_results
python3 source/render_dashboard_figures.py
python3 source/assemble_dashboard_site.py
```

The default source root is `/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats`; both analysis scripts accept `--source-root` if the project moves. Rebuilding caches compares the edited methods against the retained original backups and stops if current outputs change. Rebuilding figures alone uses portable JSON and does not need the source project. Keep backups outside the source `Methods/` folder, where additional filenames could be discovered as new methods.

Browser QA covered the desktop and narrow layouts, dataset/method/background controls, interactive zoom, figure switching, persisted option URLs, reference hyperlinks, and downloadable-file availability. Detailed scientific verification is in `review/corrections_regression_report.json` and its companion CSVs.


The dashboard exporter accepts both the current-level rate CSVs and older retrospective-level CSVs. Its `saved_summary.level_source` records the matching convention; headline levels always use current outputs. Endpoint comparison diagnostics from the root rate script remain console-only.
