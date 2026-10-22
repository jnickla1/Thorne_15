# Updating the climate state dashboard

The live site is GitHub Pages from `~/Downloads/climate-state-dashboard` (remote `git@github.com:jnickla1/climate-state-dashboard.git`, branch `main`). That repo only holds a copy of `climate-dashboard/site/dist/`. Every script below runs from `climate-dashboard/` with the cleanPy Python: `~/opt/anaconda3/envs/cleanPy/bin/python`.

## Steps

1. **Supplementary inputs** (needed only when the historical outputs in `Results2/` or the `Results/` pickles change):
   - `source/run_cgwl_wmo_histens.py`: runs CGWL-WMOLC-ADCP on the 100 heads/tails members.
   - `source/calibrate_dashboard_combinations.py`: CSCM and PIVW-8 weights. It is read-only and asserts that the 7-method scales match `Results2/comb_real_scales.npy`.
2. **Export the data**: `source/export_dashboard_climate_metrics.py all`. This reads `Results/` directly; don't pass `--results-dir`, because `data/corrected_results` is stale.
3. **Render the figures**: `source/render_dashboard_figures.py`, which draws 96 variants plus the default figure (about 5 minutes).
4. **Assemble the site**: `source/assemble_dashboard_site.py`. It copies data and figures into `site/dist/assets`, writes the rate histories, and runs `source/build_i18n.py`.
5. **Bump the cache tags** whenever data or scripts change. Each `?v=YYYYMMDD-tag` must change to a new value:
   - every tag in `site/dist/index.html`;
   - the `.json?v=` tags in `app.js`, `member-ribbon.js` and `rate-history.js`;
   - the `rate-history.js?v=` tag in `app.js`.
6. **Preview**: `python3 -m http.server 8788 --directory site/dist`, then open http://localhost:8788. Check:
   - the headline cards;
   - Smoothed ensemble (member ribbon);
   - "Explore rates over time";
   - Figure view;
   - that the browser console shows no errors.
7. **Publish**:
   1. Back up the publish repo.
   2. `cp -R site/dist/. ~/Downloads/climate-state-dashboard/`
   3. Review `git status` and `git diff --stat`, then commit and push.

## Text and translations

- English text: `site/dist/i18n/en-GB.json`, plus `data-i18n` keys in `index.html`. Rerun `source/build_i18n.py` after editing.
- Adding a language: see `site/dist/i18n/README.md`.
