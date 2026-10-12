"""Independent, read-only dashboard analysis of the saved heads/tails results.

Original scripts and source data are never modified. Run in cleanPy, for example:
python histens_rates_dashboard.py histens_rpi_real --source-root /path/to/Thorne_15_codefigurestats

Crossings use the first upward crossing of adjacent finite annual values, with
linear interpolation. Missing/unobserved crossings remain NaN, never a future
sentinel or the nearest endpoint. The observed centered 20-year mean provides
the definitive benchmark for both warming definitions. Anthropogenic panels
explicitly call this the realized observational benchmark, not attributed truth.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import pickle
from scipy.ndimage import gaussian_filter1d
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq

REAL = ['CGWL10y_sfUKCP', 'EBMKF_ta4', 'FaIR_comb_unB', 'GAM_AR1',
        'GWI_tot_CGWL', 'Kal_flexLin', 'lowess1dt36wnc']
ANTHRO = ['FaIR_nonat_unB', 'GWI_anthro_CGWL']
METHOD_LABELS = {
    'PIVW_comb': 'PIVW combination', 'equal_weight_comb': 'Equal weighting',
    'CGWL10y_sfUKCP': 'CGWL (UKCP)', 'EBMKF_ta4': 'Energy balance Kalman filter',
    'FaIR_comb_unB': 'FaIR realized', 'GAM_AR1': 'GAM with AR(1)',
    'GWI_tot_CGWL': 'Global Warming Index (total)', 'Kal_flexLin': 'Flexible Kalman filter',
    'lowess1dt36wnc': 'LOWESS (linear)', 'FaIR_nonat_unB': 'FaIR anthropogenic',
    'GWI_anthro_CGWL': 'Global Warming Index (anthropogenic)'}
PAPER = 'https://essd.copernicus.org/preprints/essd-2025-825/'
METHOD_LINKS = {
    'PIVW_comb': 'https://doi.org/10.2307/3001666',
    'equal_weight_comb': PAPER,
    'EBMKF_ta4': 'https://doi.org/10.1175/JCLI-D-23-0580.1',
    'FaIR_comb_unB': 'https://doi.org/10.5194/gmd-11-2273-2018',
    'FaIR_nonat_unB': 'https://doi.org/10.5194/gmd-11-2273-2018',
    'GWI_tot_CGWL': 'https://doi.org/10.5194/essd-16-2625-2024',
    'GWI_anthro_CGWL': 'https://doi.org/10.5194/essd-16-2625-2024',
    'CGWL10y_sfUKCP': PAPER, 'GAM_AR1': PAPER, 'Kal_flexLin': PAPER,
    'lowess1dt36wnc': PAPER}

def centered_mean(x, length):
    """Same symmetric convention as cent10y_method and cent20y_method.

    Mean of the two half-year-shifted L-year means: L+1 observations,
    half weight on the two endpoints. No extrapolation at either edge.
    """
    x = np.asarray(x, float)
    h = length // 2
    out = np.full_like(x, np.nan)
    for i in range(h, x.shape[-1]-h):
        out[..., i] = (x[..., i-h:i+h].mean(axis=-1) +
                       x[..., i-h+1:i+h+1].mean(axis=-1))/2
    return out

def first_crossing(years, y, threshold, earliest=1960):
    y = np.asarray(y)
    valid = (np.isfinite(y[:-1]) & np.isfinite(y[1:]) &
             (years[:-1] >= earliest) & (y[:-1] < threshold) & (y[1:] >= threshold))
    ii = np.flatnonzero(valid)
    if not len(ii):
        return np.nan
    i = ii[0]
    return float(years[i] + (threshold-y[i])/(y[i+1]-y[i]))

def linear_fit(t, y):
    mask = np.isfinite(t) & np.isfinite(y)
    t, y = np.asarray(t)[mask], np.asarray(y)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    tc = t-t.mean()
    slope = np.dot(tc, y-y.mean())/np.dot(tc, tc)
    resid = y - y.mean() - slope*tc
    se = np.sqrt(np.dot(resid, resid)/(len(t)-2)/np.dot(tc, tc))
    return float(slope*10), float(se*10)

def definitive_rate(years, mean10, crossing):
    if not np.isfinite(crossing):
        return np.nan, np.nan
    t = crossing + np.arange(-5, 6)
    # np.interp preserves NaNs at the finite-data boundary.
    y = np.interp(t, years, mean10, left=np.nan, right=np.nan)
    if not np.isfinite(y).all():
        return np.nan, np.nan
    return linear_fit(t, y)

def quartic_crossing(years, mean10, mean20, threshold):
    """User-requested estimate for a crossing not yet observed by centered20.

    Fit the last 25 finite centered10 annual values. Select the first real,
    upward root after the final centered20 year. Preserve the distinction
    between this estimate and an observed centered20 crossing.
    """
    ix=np.flatnonzero(np.isfinite(mean10))[-25:]
    last20=years[np.isfinite(mean20)][-1]
    if len(ix)!=25: return np.nan,np.nan,np.nan,None
    q=np.polynomial.Polynomial.fit(years[ix],mean10[ix],4)
    roots=(q-threshold).roots()
    roots=sorted(float(z.real) for z in roots if abs(z.imag)<1e-7 and
                 last20 < z.real <= years[-1]+25 and q.deriv()(z.real)>0)
    if not roots: return np.nan,np.nan,np.nan,None
    cross=roots[0]
    t=cross+np.arange(-5,6)
    y=np.interp(t,years,mean10,left=np.nan,right=np.nan)
    missing=~np.isfinite(y)
    y[missing]=q(t[missing])
    rate,se=linear_fit(t,y)
    meta=dict(fit_start=float(years[ix[0]]),fit_end=float(years[ix[-1]]),
              coefficients=q.convert().coef.tolist(),rate_extrapolated_points=int(missing.sum()),
              model_uncertainty_included=False)
    return cross,rate,se,meta

def quadratic_fit(t, y, endpoint):
    mask = np.isfinite(t) & np.isfinite(y)
    t, y = np.asarray(t)[mask]-endpoint, np.asarray(y)[mask]
    if len(t) < 4:
        return np.nan, np.nan, np.full(3, np.nan)
    X = np.column_stack([t*t, t, np.ones(len(t))])
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != 3:
        return np.nan, np.nan, np.full(3, np.nan)
    residual = y-X@beta
    cov = (residual@residual)/(len(t)-3)*np.linalg.inv(X.T@X)
    return float(beta[1]*10), float(np.sqrt(max(cov[1,1],0))*10), beta

def predicted_rate(years, current, crossing):
    """Endpoint derivative of a quadratic fit to 15 past annual outputs.

    The 15 samples end at the interpolated crossing; no later result enters.
    Displayed comparison segment spans crossing-15 through crossing.
    """
    if not np.isfinite(crossing):
        return np.nan, np.nan
    t = crossing + np.arange(-14, 1)
    y = np.interp(t, years, current, left=np.nan, right=np.nan)
    if not np.isfinite(y).all():
        return np.nan, np.nan
    return quadratic_fit(t, y, crossing)[:2]

def summarize(crossings, rates, ses):
    crossings, rates, ses = map(lambda x: np.asarray(x, float), (crossings, rates, ses))
    n = int(np.isfinite(crossings).sum())
    rmask = np.isfinite(rates) & np.isfinite(ses)
    cross = crossings[np.isfinite(crossings)]
    return dict(year=float(cross.mean()) if n else np.nan,
                year_sd=float(cross.std(ddof=1)) if n>1 else np.nan,
                year_lo=float(np.quantile(cross,.025)) if n else np.nan,
                year_hi=float(np.quantile(cross,.975)) if n else np.nan,
                n_crossed=n, n_members=len(crossings), n_rate=int(rmask.sum()),
                rate=float(rates[rmask].mean()) if rmask.any() else np.nan,
                rate_se=float(np.sqrt(np.mean(ses[rmask]**2) + np.var(rates[rmask],ddof=1))) if rmask.sum()>1 else np.nan)

def mixture_interval(centers, ses):
    mask = np.isfinite(centers) & np.isfinite(ses) & (ses>0)
    centers, ses = np.asarray(centers)[mask], np.asarray(ses)[mask]
    lo, hi = min(centers-10*ses), max(centers+10*ses)
    return [float(brentq(lambda x: norm.cdf(x,loc=centers,scale=ses).mean()-q,lo,hi)) for q in [.025,.975]]

def clean_json(x):
    if isinstance(x, dict): return {k: clean_json(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [clean_json(v) for v in x]
    if isinstance(x, np.ndarray): return clean_json(x.tolist())
    if isinstance(x, (float,np.floating)): return round(float(x),7) if np.isfinite(x) else None
    if isinstance(x, np.integer): return int(x)
    return x


def load_paper_probabilities(root, experiment, methods):
    """Reuse the paper's calibrated current probabilities, unchanged by the fixes."""
    sums = {m: [] for m in methods}
    starts = [1970, 1990, 2005]
    ends = [2006, 2026, 2026]
    for member in range(100):
        names = pd.read_csv(root/'Results'/f'headstails_statistics_{experiment}{member}.csv').method_name.tolist()
        a = np.load(root/'Results'/f'headstails_method_fineprobs_{experiment}{member}.npy')
        for m in methods:
            sums[m].append(a[names.index(m)])
    out = {}
    for m, samples in sums.items():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            cdf = np.nanmean(samples, axis=0)
        panels=[]
        for j in range(3):
            y=np.arange(starts[j],ends[j]); c=cdf[:len(y),j]
            mono=np.maximum.accumulate(np.where(np.isfinite(c),c,0))
            pdf=np.clip(np.diff(mono,prepend=mono[0]),0,None)
            nz=np.flatnonzero(pdf>0)
            if len(nz):pdf=pdf[:nz[-1]+1]
            fine=gaussian_filter1d(np.repeat(pdf/12,12),12,mode='nearest')*12
            t=starts[j]+np.arange(len(fine))/12+1/24
            panels.append(dict(years=y,probabilities=c,density_years=t[:-6],density=fine[:-6]))
        out[m]=panels
    return out

def analyze(source_root, experiment, outdir, results_dir=None):
    root, out = Path(source_root), Path(outdir)
    out.mkdir(parents=True,exist_ok=True)
    _, baseline, goal = experiment.split('_')
    methods = REAL if goal=='real' else ANTHRO
    rawfile = root/'Hist_heads_tails'/('sst_pseudo.csv' if baseline=='rpi' else 'sst_pseudo_1981-2010.csv')
    a = np.loadtxt(rawfile,delimiter=','); years=a[:,0]; raw=a[:,1:].T
    shift=0 if baseline=='rpi' else -.7
    thresholds=np.array([.5,1,1.5,2])+shift
    histvar=pd.read_csv(root/'Results2/historical_names_var_scale.csv').set_index('method_name')
    iv=1/histvar.loc[methods,'err_var100'].to_numpy()
    scales=np.load(root/'Results2'/f'comb_{goal}_scales.npy')
    mean20=centered_mean(raw,20); mean10=centered_mean(raw,10)
    curves={m: [] for m in methods}; recent={m: [] for m in methods}; stds={m: [] for m in methods}
    legacy_recent={m: [] for m in methods}
    proxy_methods=set()
    digest=hashlib.sha256()
    for member in range(raw.shape[0]):
        filename=f'results{experiment}{member}.pickle'
        path=(Path(results_dir) if results_dir else root/'Results')/filename
        b=path.read_bytes(); digest.update(b)
        data=pickle.loads(b)
        legacy=pickle.loads((root/'Results'/filename).read_bytes()) if results_dir else data
        for m in methods:
            trend=data[m]['LT_trend']
            if not isinstance(trend,tuple):
                raise TypeError(f'{m}: expected saved tuple, got {type(trend).__name__}')
            c,se,r,rse=map(np.asarray,trend)
            curves[m].append(c); stds[m].append(se*histvar.loc[m,'best_alter_scale'])
            recent[m].append(r if np.isfinite(r).any() else c)
            lc,_,lr,_=map(np.asarray,legacy[m]['LT_trend'])
            legacy_recent[m].append(lr if np.isfinite(lr).any() else lc)
            if data[m].get('retrospective_source')=='current_proxy_no_retrospective_in_reader' or not np.isfinite(r).any():
                proxy_methods.add(m)
        # Saved means must match the observational convention exactly.
        np.testing.assert_allclose(data['cent20y']['LT_trend'][2],mean20[member],equal_nan=True,atol=1e-10)
    curves={k:np.array(v) for k,v in curves.items()}; stds={k:np.array(v) for k,v in stds.items()}; recent={k:np.array(v) for k,v in recent.items()}
    legacy_recent={k:np.array(v) for k,v in legacy_recent.items()}
    for m,w,scale in [('equal_weight_comb',np.ones(len(methods)),scales[0]),('PIVW_comb',iv,scales[1])]:
        for dest in [curves,recent,legacy_recent]:
            stack=np.array([dest[k] for k in methods]); weights=w[:,None,None]
            denom=np.sum(np.isfinite(stack)*weights,axis=0)
            with np.errstate(invalid='ignore',divide='ignore'):
                dest[m]=np.nansum(stack*weights,axis=0)/denom
        stds[m]=np.full_like(curves[m],scale)
    paper_probs=load_paper_probabilities(root,experiment,list(curves))
    definition='realized' if goal=='real' else 'anthropogenic'
    result=dict(experiment=experiment,baseline=baseline,goal=goal,definition=definition,years=years,
                baseline_label='1850–1900' if baseline=='rpi' else '1981–2010',thresholds=thresholds,
                raw=raw,centered20=mean20,centered10=mean10,methods={},
                reference_last_year=int(years[np.isfinite(mean20).any(axis=0)][-1]),
                weights={m:float(w) for m,w in zip(methods,iv/iv.sum())},
                metadata=dict(source_sha256=digest.hexdigest(),raw_sha256=hashlib.sha256(rawfile.read_bytes()).hexdigest(),
                              probability_definition='Calibrated current probabilities saved by the paper pipeline in its three assessment windows. Outside those windows and for the fourth threshold, Gaussian approximations use the same calibrated means/SE. These are pointwise exceedances, not a first-passage CDF.',
                              crossing_definition='First upward crossing, linearly interpolated between adjacent finite annual values.',
                              benchmark='Realized observational centered 20-year mean, including on anthropogenic panels.',
                              comparison_fit='Endpoint derivative of a quadratic fitted to 15 past annual current-method outputs.',
                              recent_saved_fit='Saved rate CSV: 15 annual retrospective outputs when available, otherwise current outputs.',
                              headline_definition='Current-method warming level; rate from the corrected retrospective 2011–2025 histories, with current proxies where explicitly identified.',
                              retrospective_current_proxies=sorted(proxy_methods),
                              results_directory=str(Path(results_dir).resolve()) if results_dir else str(root/'Results'),
                              definitive_fit='OLS through centered 10-year means at crossing-5 through crossing+5 inclusive.',
                              uncertainty='95% level range: equal-member Gaussian mixture with calibrated method SE. Rates: 1.96 times sqrt(mean(within-fit SE squared)+sample variance of member rates).',
                              censoring='Comparison crossing-year ranges condition on members that cross. Definitive preindustrial 1.0 C includes the seven requested quartic estimates; their model uncertainty is not included.'))
    definitive=[]; rows=[]; memberrows=[]
    for th in thresholds[:3]:
        cr=np.array([first_crossing(years,y,th) for y in mean20])
        observed=cr.copy()
        rs=np.array([definitive_rate(years,y,c) for y,c in zip(mean10,cr)])
        estimated=np.zeros(len(cr),dtype=bool); quartics={}
        # Only the seven non-crossing preindustrial 1.0 C members receive
        # the specifically requested quartic treatment. No 1.5 C projection.
        if baseline=='rpi' and np.isclose(th,1.0):
            for i in np.flatnonzero(~np.isfinite(cr)):
                c,rate,se,meta=quartic_crossing(years,mean10[i],mean20[i],th)
                if np.isfinite(c):
                    cr[i]=c;rs[i]=[rate,se];estimated[i]=True;quartics[str(i)]=meta
        s=summarize(cr,rs[:,0],rs[:,1]); s.update(threshold=float(th),crossings=cr,
                    observed_crossings=observed,n_observed=int(np.isfinite(observed).sum()),
                    n_estimated=int(estimated.sum()),estimated_members=np.flatnonzero(estimated),quartic_fits=quartics)
        definitive.append(s); rows.append(dict(kind='definitive',method='cent10y_20y',**{k:v for k,v in s.items() if k not in ['crossings','observed_crossings','estimated_members','quartic_fits']}))
        for i,(c,(rate,se)) in enumerate(zip(cr,rs)):
            memberrows.append(dict(member=i,kind='definitive',method='cent10y_20y',threshold=th,crossing=c,rate=rate,rate_se=se,
                              crossing_status='quartic_estimate' if estimated[i] else 'observed' if np.isfinite(c) else 'pending'))
    result['definitive']=definitive
    saved=pd.read_csv(root/'rates_histens'/f'rates_histens_{experiment}.csv')
    recent_method_fits={m:np.array([quadratic_fit(years[-15:],y[-15:],years[-1])[:2] for y in recent[m]]) for m in methods}
    def recent_headline_rate(method):
        if method in methods:
            fits=recent_method_fits[method]
            return summarize(np.zeros(len(fits)),fits[:,0],fits[:,1])
        # Match histens_rates.py's within-method fit variance plus the
        # between-method dispersion, then add between-member dispersion.
        w=iv/iv.sum() if method=='PIVW_comb' else np.ones(len(methods))/len(methods)
        slopes=np.array([recent_method_fits[m][:,0] for m in methods])
        ses=np.array([recent_method_fits[m][:,1] for m in methods])
        rate=np.sum(w[:,None]*slopes,axis=0)
        variance=np.sum(w[:,None]*(ses**2+(slopes-rate)**2),axis=0)
        return summarize(np.zeros(len(rate)),rate,np.sqrt(variance))
    for m in ['PIVW_comb','equal_weight_comb']+methods:
        c=curves[m]; se=stds[m]; r=recent[m]
        comparisons=[]
        for th in thresholds[:3]:
            cr=np.array([first_crossing(years,y,th) for y in c]); rs=np.array([predicted_rate(years,y,x) for y,x in zip(c,cr)])
            s=summarize(cr,rs[:,0],rs[:,1]); s.update(threshold=float(th),crossings=cr)
            comparisons.append(s); rows.append(dict(kind='comparison',method=m,**{k:v for k,v in s.items() if k!='crossings'}))
            for i,(x,(rate,ser)) in enumerate(zip(cr,rs)):
                memberrows.append(dict(member=i,kind='comparison',method=m,threshold=th,crossing=x,rate=rate,rate_se=ser))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            probabilities=np.array([np.nanmean(norm.cdf(c,loc=th,scale=se),axis=0) for th in thresholds])
        for j,panel in enumerate(paper_probs[m]):
            idx=np.searchsorted(years,panel['years'])
            probabilities[j,idx]=panel['probabilities']
        # Preserve CSV numbers for the requested default. Explicit alternative
        # recomputes the recent rate solely from current outputs.
        sr=saved[(saved.method==m)&(saved.level.astype(str)==str(int(years[-1])))].iloc[0]
        currentfits=np.array([quadratic_fit(years[-15:],y[-15:],years[-1])[:2] for y in c])
        curr=summarize(c[:,-1],currentfits[:,0],currentfits[:,1])
        curr_summary=dict(level=float(np.mean(c[:,-1])),level_interval=mixture_interval(c[:,-1],se[:,-1]),
                          rate=curr['rate'],rate_se=curr['rate_se'])
        # Accept current-level rate CSVs and retained older retrospective-level CSVs.
        saved_level = float(sr.year_cross)
        if np.isclose(saved_level, curr_summary['level'], atol=1e-10, rtol=0):
            saved_centers = c[:, -1]
            saved_level_source = 'current'
        else:
            saved_centers = legacy_recent[m][:, -1]
            np.testing.assert_allclose(saved_centers.mean(), saved_level, atol=1e-10, rtol=0,
                err_msg=f'{m}: saved level matches neither current nor legacy retrospective outputs; refresh the rate CSV.')
            saved_level_source = 'legacy_retrospective'
        saved_summary=dict(level=saved_level,level_interval=mixture_interval(saved_centers,se[:,-1]),
                           rate=float(sr.rate),rate_se=float(sr.rate_uncert_1sd),
                           level_source=saved_level_source)
        corrected_rate=recent_headline_rate(m)
        headline_summary=dict(level=curr_summary['level'],level_interval=curr_summary['level_interval'],
                              rate=corrected_rate['rate'],rate_se=corrected_rate['rate_se'])
        def fit_curve(arr):
            # Direct fit of the ensemble mean recreates weighted endpoint rate
            # for combinations where method availability is constant.
            mean=np.nanmean(arr[:,-15:],axis=0)
            _,_,beta=quadratic_fit(years[-15:],mean,years[-1])
            x=np.linspace(years[-1]-15,years[-1],61); t=x-years[-1]
            return dict(years=x,values=beta[0]*t*t+beta[1]*t+beta[2])
        result['methods'][m]=dict(label=METHOD_LABELS[m],link=METHOD_LINKS[m],current=c,recent=r,
                                  probabilities=probabilities,comparisons=comparisons,paper_distributions=paper_probs[m],
                                  current_summary=curr_summary,saved_summary=saved_summary,
                                  headline_summary=headline_summary,
                                  retrospective_endpoint=float(r[:,-1].mean()),
                                  retrospective_is_current_proxy=m in proxy_methods,
                                  current_fit=fit_curve(c),saved_fit=fit_curve(legacy_recent[m]),recent_fit=fit_curve(r))
    prefix=experiment.removeprefix('histens_')
    refs_path=Path(__file__).resolve().parents[1]/'data/method_references.json'
    if refs_path.exists():
        refs=json.loads(refs_path.read_text())
        for m,entry in result['methods'].items():
            entry['link']=refs['methods'][m]['primary_url']
            entry['references']=refs['methods'][m]['references']
    colors_path=Path(__file__).resolve().parents[1]/'data/paper_histogram_colors.json'
    if colors_path.exists():
        paper_colors=json.loads(colors_path.read_text())['methods']
        for m,entry in result['methods'].items():
            entry['histogram_color']=paper_colors[m]['color']
            entry['histogram_alpha']=paper_colors[m]['alpha']
    pd.DataFrame(rows).to_csv(out/f'{prefix}_crossings_rates.csv',index=False)
    pd.DataFrame(memberrows).to_csv(out/f'{prefix}_member_crossings_rates.csv',index=False)
    saved.to_csv(out/f'{prefix}_legacy_saved_rates.csv',index=False)
    pd.DataFrame([dict(method=m,**v['headline_summary']) for m,v in result['methods'].items()]).to_csv(out/f'{prefix}_corrected_headlines.csv',index=False)
    for j in range(4):
        p=pd.DataFrame({m:v['probabilities'][j] for m,v in result['methods'].items()},index=years.astype(int))
        p.index.name='year'; p.to_csv(out/f'{prefix}_probability_{thresholds[j]:.1f}C.csv')
    pd.DataFrame(raw.T,index=years.astype(int)).rename_axis('year').to_csv(out/f'{prefix}_raw_ensemble.csv')
    (out/f'{prefix}.json').write_text(json.dumps(clean_json(result),separators=(',',':'),allow_nan=False))
    print(prefix, 'definitive:',[(x['year'],x['n_crossed'],x['rate']) for x in definitive[:2]],flush=True)
    print('PIVW saved/current',saved_summary if m=='PIVW_comb' else result['methods']['PIVW_comb']['saved_summary'],result['methods']['PIVW_comb']['current_summary'],flush=True)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('experiment',choices=['all','histens_rpi_real','histens_rpi_anthro','histens_satcal_real','histens_satcal_anthro'])
    parser.add_argument('--source-root',type=Path,default=Path('/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats'))
    parser.add_argument('--outdir',type=Path,default=Path(__file__).resolve().parents[1]/'data')
    parser.add_argument('--results-dir',type=Path,help='Use verified corrected result copies while retaining original CSVs for comparison.')
    args=parser.parse_args()
    exps=[f'histens_{b}_{g}' for b in ['rpi','satcal'] for g in ['real','anthro']] if args.experiment=='all' else [args.experiment]
    for exp in exps: analyze(args.source_root,exp,args.outdir,args.results_dir)

if __name__=='__main__': main()
