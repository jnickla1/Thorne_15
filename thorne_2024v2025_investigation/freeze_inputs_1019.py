"""Write inputs_1019_manifest.csv: path, mtime, md5 of every 19-Oct-vintage input used by the *_1019 scripts."""
import hashlib, glob, os, time, pandas as pd
R = '../'
pats = ['Results/resultshistens_rpi_real*.pickle', 'Results/resultshistens_rpi_anthro*.pickle',
        'Results/resultshistens_satcal_real*.pickle', 'Results/resultshistens_satcal_anthro*.pickle',
        'Results/headstails_method_fineprobs_histens_*.npy', 'Results/headstails_statistics_histens_*.csv',
        'Results/resultsfut_ESM1-2-LR_SSP126_constVolc*.pickle', 'Results/resultsfut_ESM1-2-LR_SSP245_constVolc*.pickle',
        'Results/resultsfut_ESM1-2-LR_SSP370_constVolc*.pickle', 'Results/resultsfut_NorESM_RCP45_VolcConst*.pickle',
        'hist_results.pickle', '261018_HadCRUT5.1_var_adjusted/*.csv', '261018_HadCRUT5.1_var_adjusted/*.npy',
        'Results2/historical_names_var_scale.csv', 'Results2/comb_*_scales.npy',
        'rates_histens/rates_histens_histens_*.csv', 'climate-dashboard/data/*_corrected_headlines.csv']
rows = []
for p in pats:
    for f in sorted(glob.glob(R + p)):
        rows.append(dict(path=f[len(R):], mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f))),
                         bytes=os.path.getsize(f), md5=hashlib.md5(open(f, 'rb').read()).hexdigest()))
d = pd.DataFrame(rows); d.to_csv('inputs_1019_manifest.csv', index=False)
print(len(d), 'files; mtime range', d.mtime.min(), '->', d.mtime.max())
print(d.assign(day=d.mtime.str[:10]).groupby(d.path.str.extract(r'^([^0-9]*)')[0]).agg(n=('path', 'size'), first=('mtime', 'min'), last=('mtime', 'max')).to_string())
