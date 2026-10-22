# Codex handoff: publish the 21 Oct 2026 dashboard refresh

## Status (already done, 21 Oct)

- `data/` was re-exported from the current `Results/` (commit 2d95aaec reruns) with a refitted `data/supplementary/dashboard_combination_calibration.json`.
- The figures were re-rendered (`render_dashboard_figures.py`) and the site assembled (`assemble_dashboard_site.py`): 96 variants and 32 method CSVs, and the i18n check passed.
- The cache tags were bumped to `?v=20261021-data`:
  - every tag in `site/dist/index.html`;
  - the `.json?v=` tags in `app.js`, `member-ribbon.js` and `rate-history.js`;
  - `rate-history.js?v=` in `app.js`.
- A preview at http://localhost:8788 passed: PIVW 1.36 °C, 0.294 °C/decade, 6.24 %; the member ribbon loads; the rate panel draws; Figure view loads; there are no console errors.

### Expected changes compared with the live site (commit 0472105), realised, 1850–1900

| Method | 2025 level (°C) | P(>1.5 °C) |
|---|---|---|
| PIVW | 1.366 → 1.365 | 6.42 → 6.24 % |
| FaIR_comb_unB | 1.368 → 1.365 | 6.79 → 6.42 % |
| GWI_tot_CGWL | 1.332 → 1.330 | 0.68 → 0.58 % |
| CSCM | 1.369 → 1.368 | 6.89 → 6.62 % |

Every other method moves by less than 0.003 °C and 0.3 percentage points, and none of the rates change.

## Remaining steps

Run these from `/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/climate-dashboard`.

1. **Re-render and assemble** only if anything in `data/` has changed since 21 Oct 11:39; otherwise skip this step. Use the cleanPy Python, `~/opt/anaconda3/envs/cleanPy/bin/python`:
   ```
   $PY source/render_dashboard_figures.py && $PY source/assemble_dashboard_site.py
   ```
   It must end with "Website ready: 96 variants" and print an i18n line with no missing keys. If you rerun it, bump the tags to a new value; never reuse `20261021-data`.
2. **Compare with the live site.** For each method in `site/dist/assets/data/{rpi,satcal}_{real,anthro}.json`, compare against `~/Downloads/climate-state-dashboard/assets/data/` on:
   - `headline_summary.level`;
   - `headline_summary.rate`;
   - `probabilities[2][-1]`.
   
   Publish automatically only if both of these hold:
   - no level or rate moves by more than 0.01 (°C, or °C/decade for rates);
   - no probability moves by more than 1.5 percentage points.
   
   Larger changes: stop and report to John.
3. **Preview.** Run `python3 -m http.server 8788 --directory site/dist`, then check:
   - the page loads with no console errors;
   - the metric cards show numbers, not "pending";
   - Smoothed ensemble shows the ribbon hint text and no "does not match" error;
   - "Explore rates over time" draws;
   - Figure view shows an image.
4. **Publish.**
   1. Back up first: `tar czf /tmp/csd_backup.tgz -C ~/Downloads climate-state-dashboard`.
   2. In `~/Downloads/climate-state-dashboard`, check that `git fetch && git status` is clean and up to date with origin/main.
   3. `cp -R site/dist/. ~/Downloads/climate-state-dashboard/`. The only files left out are `.git` and the repo's own `README.md`, which `site/dist` doesn't have.
   4. Make sure no `__capture.js` or test locale (`i18n/he.json`, `ar.json`, `zh-Hans.json`, `fr.json`) was copied in.
   5. `git diff --stat`: expect about 130 changed assets plus the tag edits in `index.html`, `app.js`, `member-ribbon.js` and `rate-history.js`, and no deletions.
   6. Commit as "Refresh dashboard data from 21 Oct reruns", ending the message with the line `Co-Authored-By: Codex`, then run `git push origin main`.

## Context

- The Thorne_15 repo's `climate-dashboard/` source changes are not committed there; John asked to leave that folder's git setup alone.
- Other Claude sessions may edit `site/dist`, so re-read files before editing them.
- The live site is English-only. `i18n/locales.js` must list only en-GB.
- The human-readable steps are in `UPDATE_WEBSITE.md`.
