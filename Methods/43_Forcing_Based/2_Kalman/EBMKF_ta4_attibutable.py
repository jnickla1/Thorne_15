"""Figure 16 attribution experiment, using the finalized ta4 historical model.

The comparison ends in 2023 to match the supplied multi-method assessments.
The default components use native-grid reductions from CanESM5, IPSL-CM6A-LR
and MIROC6, with actual DAMIP/SSP245 extensions through 2025 where available.
No post-2014 temperature extrapolation is used in this mode. The old MPI route
remains explicitly available with --component-source legacy-mpi. These are
illustrative component filters, not a joint attribution posterior.
Run from any directory; --output-dir permits an isolated review run.
"""
from pathlib import Path
import argparse, importlib, json, os, sys, hashlib
import numpy as np
import pandas as pd
from scipy.stats import norm
from netCDF4 import Dataset
ROOT = Path(__file__).resolve().parents[3]
# A review copy can explicitly point at the research repository.
ROOT = Path(os.environ.get('THORNE_ROOT', str(ROOT)))
sys.path.insert(0, str(ROOT))
import config
wrapper = importlib.import_module('Methods.43_Forcing_Based.2_Kalman.EBMKF_ta4_method')
ekf = wrapper.ekf
from ebmkf_v11_support import nonnegative_sqrt, ECO2_ERF_SLOPE


def annual(values):
    values = np.asarray(values).squeeze()
    if values.ndim != 1 or len(values) % 12:
        raise ValueError('Expected complete monthly global-mean input')
    return values.reshape(-1, 12).mean(axis=1)


def extend(values, end_year):
    """Preserve the legacy component-extension recipe and seed, locally."""
    values = np.asarray(values)
    years = np.arange(1850, 1850 + len(values))
    fit = years >= 1980
    coeff = np.polyfit(years[fit], values[fit], 1)
    sigma = np.sqrt(np.mean((values[fit] - np.polyval(coeff, years[fit]))**2))
    future = np.arange(years[-1] + 1, end_year + 1)
    rng = np.random.RandomState(320)
    return np.r_[values, np.polyval(coeff, future) + rng.normal(0, sigma, len(future))]


def component_observations(experiment, n, preind_base):
    base = Path(config.CLIMATE_DATA_PATH) / 'attributable-warming_ESM1-2-LR'
    ts = base / f'ts/aave_ts_Amon_MPI-ESM1-2-LR_hist-{experiment}_r4i1p1f1_gn_185001-201412.nc'
    ohc = base / f'ohca_aave/ohca_thetaot_Emon_MPI-ESM1-2-LR_hist-{experiment}_r4i1p1f1_gn_185001-201412.nc'
    toa = base / f'rt/combined_global_rt/hist-{experiment}/global_rt_Amon_MPI-ESM1-2-LR_hist-{experiment}_r4i1p1f1_gn_185001-186912.nc-201412.nc'
    arrays = []
    for path, variable in [(ts, 'ts'), (ohc, '__xarray_dataarray_variable__'), (toa, 'rt')]:
        with Dataset(path) as ds:
            arrays.append(annual(ds.variables[variable][:]))
    temperature, ocean, flux = arrays
    temperature = temperature - np.mean(temperature[:51]) + preind_base + ekf.offset
    ocean = (ocean - ocean[0]) / ekf.zJ_from_W
    return np.stack([extend(x, 1849+n) for x in [temperature, ocean, flux]], axis=1)[:, :, None]


def solve(observations, anomaly_reference=None):
    mean, variance, _, _ = ekf.ekf_run(observations, len(observations), retPs=3, diagnostics=False)
    reference = ekf.Teq1850 if anomaly_reference is None else anomaly_reference
    return mean.copy() - reference, nonnegative_sqrt(variance).copy()


def selected_component_erf(years, columns, percentile):
    """Same Gaussian/independent-component approximation as historical eCO2.

    A subset percentile is evaluated for that subset; marginal subset percentiles
    are not additive quantiles of the joint full non-aerosol forcing sum.
    """
    frames=[]
    for tag in ['p05','best','p95']:
        frame=pd.read_csv(Path(ekf.__file__).with_name(f'ERF_{tag}_aggregates_1750-2024.csv'),index_col=0)
        frame.index=np.floor(frame.index.astype(float)).astype(int)
        frames.append(frame[columns])
    low,best,high=frames
    mean=best.sum(axis=1)
    sigma=pd.Series(np.sqrt(np.sum(((high-low)/(2*norm.ppf(.95))).to_numpy()**2,axis=1)),index=mean.index)
    for series in [mean,sigma]:series.loc[2025]=np.polyval(np.polyfit(np.arange(2015,2025),series.loc[2015:2024],1),2025)
    return (mean+norm.ppf(percentile/100)*sigma).loc[years].to_numpy()



def configure_component(label, n):
    """Historical tuned parameters; change only which forcing is active."""
    wrapper.restore_baselines(); ekf.n_iters = n
    # Model simulations are not CERES observations; retain pre-CERES variance.
    ekf.TOA_meas_artif_var[2001-1850:] = ekf._orig_TOA_meas_artif_var[1999-1850]
    ekf.Q *= .2; ekf.cM2 *= .6; ekf.dfrA_float_var *= .35/40
    nf = len(ekf.data)
    ekf.tsi = np.full(nf, ekf.sw_in)
    ekf.opt_depth = np.full(nf, ekf.raw_aod_from_inverse_denominator(ekf.involcavg))
    ekf.lCo2 = np.full(nf, ekf._orig_lCo2[0])
    ekf.anthro_clouds = np.full(nf, ekf._orig_anthro_clouds[0])
    if label == 'aer':
        ekf.anthro_clouds = ekf._orig_anthro_clouds.copy()
    elif label == 'ghg':
        erf = selected_component_erf(ekf.data[:,0].astype(int), ['CO2','CH4','N2O','halogen'], ekf.eCO2_FORCING_PERCENTILE)
        ekf.lCo2 = ekf._orig_lCo2[0] + (erf-erf[0])/ECO2_ERF_SLOPE
    elif label in ['nat','all']:
        ekf.tsi = ekf._orig_tsi.copy(); ekf.opt_depth = ekf.data[:,3]*.001
        if label == 'all':
            ekf.lCo2 = ekf._orig_lCo2.copy()
            ekf.anthro_clouds = ekf._orig_anthro_clouds.copy()
    else:
        raise ValueError(label)
    ekf.use_observational_forcing_percentiles = False
    ekf.precompute_coeffs(False)


def mixture_summary(pairs):
    """Equal model mixture moments, retaining between-model spread.

    Reported Gaussian percentiles moment-match this mixture; this is not the
    standard error of its mean or an exact joint attribution posterior.
    """
    means = np.array([p[0] for p in pairs]); ses = np.array([p[1] for p in pairs])
    valid = np.isfinite(means) & np.isfinite(ses)
    count = valid.sum(axis=0)
    center = np.nansum(np.where(valid,means,np.nan),axis=0)/np.where(count,count,np.nan)
    variance = np.nansum(np.where(valid,ses**2+(means-center)**2,np.nan),axis=0)/np.where(count,count,np.nan)
    return center, np.sqrt(variance), count


def processed_components(models, n, preind_base, out, toa_convention, expected_settings):
    """Independent component filters; same model/member fields within each run."""
    assert len(models)==len(set(models)), "Duplicate model names would double-weight a model"
    raw_results = {}
    provenance = []
    model_rows = []
    cases = {'ghg':'hist-GHG_ssp245-GHG', 'nat':'hist-nat_ssp245-nat',
             'all':'historical_ssp245', 'aer':'hist-aer_ssp245-aer'}
    for model in models:
        model_data = {}
        for label, prefix in cases.items():
            base = Path(config.CLIMATE_DATA_PATH)/('attributable-warming_'+model)/'processed'
            member = 'r4i1p1f1' if model == 'IPSL-CM6A-LR' else 'r1i1p1f1'
            files = sorted(base.glob(prefix+'_'+member+'_*_annual.csv'))
            if not files and model == 'IPSL-CM6A-LR' and label == 'aer':
                files = sorted(base.glob('hist-aer_r4i1p1f1_*_annual.csv'))
            if not files:
                raise FileNotFoundError(f'No processed {model} {label} inputs in {base}')
            if len(files)!=1: raise ValueError(f'Ambiguous processed inputs: {files}')
            path = files[0]
            source_metadata = json.loads(path.with_name(path.name.replace('_annual.csv','_provenance.json')).read_text())
            assert source_metadata['model']==model and source_metadata['member']==member
            frame = pd.read_csv(path)
            frame = frame.loc[frame.year <= 1849+n]
            assert np.array_equal(frame.year.to_numpy(), np.arange(1850,1850+len(frame)))
            assert np.isfinite(frame[['tas_K','net_toa_Wm2','ohc_ZJ']].to_numpy()).all()
            count = len(frame)
            configure_component(label,count)
            actual_settings = {key:float(getattr(ekf,key)) for key in expected_settings}
            assert actual_settings == expected_settings, (model,label,actual_settings,expected_settings)
            temp = frame.tas_K.to_numpy()
            temp = temp-temp[:51].mean()+preind_base+ekf.offset
            ocean = (frame.ohc_ZJ.to_numpy()-frame.ohc_ZJ.iloc[0])/ekf.zJ_from_W
            flux = frame.net_toa_Wm2.to_numpy().copy()
            toa_offset = float(flux[:50].mean()) if toa_convention == 'first50-anomaly' else 0.
            flux -= toa_offset
            observations = np.stack([temp,ocean,flux],axis=1)[:,:,None]
            mean,se = solve(observations)
            assert np.isfinite(mean).all() and np.isfinite(se).all()
            pair = (np.pad(mean,(0,n-count),constant_values=np.nan),
                    np.pad(se,(0,n-count),constant_values=np.nan))
            model_data[label] = pair
            for i,year in enumerate(frame.year):
                model_rows.append(dict(model=model,case=label,year=int(year),mean=float(mean[i]),se=float(se[i]),
                    input_temperature_anomaly=float(temp[i]-preind_base-ekf.offset),
                    input_ocean_ZJ=float(ocean[i]*ekf.zJ_from_W),input_TOA_Wm2=float(flux[i])))
            provenance.append(dict(model=model,case=label,annual_csv=str(path),annual_csv_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),end_year=int(frame.year.iloc[-1]),
                years=count,TOA_baseline_subtracted_Wm2=toa_offset,settings=actual_settings))
        # Derive broad OHF within each model before combining model spread.
        # Conditional filter covariances between experiments are unavailable.
        a,g,t = model_data['all'], model_data['ghg'], model_data['nat']
        model_data['model_residual'] = (a[0]-g[0]-t[0],np.sqrt(a[1]**2+g[1]**2+t[1]**2))
        model_data['model_anthropogenic'] = (a[0]-t[0],np.sqrt(a[1]**2+t[1]**2))
        for derived in ['model_residual','model_anthropogenic']:
            mean,se = model_data[derived]
            for i in range(n):
                model_rows.append(dict(model=model,case=derived,year=1850+i,mean=float(mean[i]),se=float(se[i])))
        raw_results[model] = model_data
    summaries = {}; coverage = pd.DataFrame({'Year':np.arange(1850,1850+n)})
    for label in list(cases)+['model_residual','model_anthropogenic']:
        mean,se,count = mixture_summary([raw_results[m][label] for m in models])
        summaries[label] = (mean,se)
        coverage[label+'_models'] = count
    pd.DataFrame(model_rows).to_csv(out/'Nicklas_per_model_components.csv',index=False)
    coverage.to_csv(out/'Nicklas_model_coverage.csv',index=False)
    return summaries, dict(models=models,inputs=provenance,TOA_convention=toa_convention,
        extrapolation='None. Native historical/scenario data only; missing model/cases remain missing.',
        mixture='Equal weights among available models; Gaussian moment-matched intervals include within-filter variance and between-model spread, not standard error of the model mean.',
        OHF_residual='Within each model: ALL minus NAT minus GHG. Residual includes non-GHG anthropogenic forcing responses and nonadditivity/internal variability; no cross-experiment covariance is available.',
        density='Same common 1026 kg/m3 and 3989.24495292815 J/kg/K as native reductions; CanESM diagnostic convention sensitivity remains documented in GHG_REDUCTION_REPORT.md.',
        scope='Review attribution experiment; no model-specific tuning, drift subtraction, or input extrapolation. Historical ANT and observed total are unchanged.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path, default=None, help='Optional isolated output folder; default updates the established assessment CSV')
    ap.add_argument('--end-year', type=int, default=2023, choices=[2023, 2024, 2025])
    ap.add_argument('--component-end-year', type=int, choices=[2014], default=None, help='Stop model-only components at available 2014 data; preserve ANT/total through end-year')
    ap.add_argument('--ohf', choices=['aerosol-proxy', 'ant-minus-ghg', 'model-residual'], default='aerosol-proxy')
    ap.add_argument('--component-source', choices=['processed', 'legacy-mpi'], default='processed')
    ap.add_argument('--models', nargs='+', default=['CanESM5','IPSL-CM6A-LR','MIROC6'], choices=['CanESM5','IPSL-CM6A-LR','MIROC6'])
    ap.add_argument('--model-toa', choices=['raw','first50-anomaly'], default='raw', help='Explicit TOA-offset sensitivity; raw retains legacy convention')
    args = ap.parse_args()
    os.chdir(ROOT)
    out = args.output_dir.resolve() if args.output_dir else ROOT/'Methods/43_Forcing_Based/3_Human_Induced/anthropogenic-warming-assessment/results'
    out.mkdir(parents=True, exist_ok=True)
    flat_path = out/'Nicklas_GMST_timeseries_flat.csv' if args.output_dir else Path(__file__).resolve().parent/'Nicklas_GMST_timeseries.csv'
    metadata_path = out/'attribution-metadata.json' if args.output_dir else Path(__file__).resolve().parent/'attribution-metadata.json'
    obs = pd.read_csv(ROOT/'Common_Data/HadCRUT.5.1.csv').query('Time <= @args.end_year')
    years = obs.Time.to_numpy(dtype=int); n = len(years)
    temperature = obs.Anomaly.to_numpy(); preind_base = float(np.mean(temperature[:51]))
    inputs = (obs.Lower.to_numpy(), obs.Upper.to_numpy())
    # This call is authoritative for ANT: it includes the finalized AOD/solar,
    # early-warming and delayed-cooling preprocessing, OHC/TOA corrections,
    # covariance tuning, forcing percentiles and current-output SE factor.
    ant_mean, ant_se, _, _ = wrapper.run_method(years, temperature, inputs, '', 'historical')
    data = {'ant': (ant_mean-preind_base, ant_se)}
    settings = {k:float(getattr(ekf,k)) for k in ['gad','fdbkA','AEROSOL_FORCING_PERCENTILE','eCO2_FORCING_PERCENTILE','volcanic_target_efficiency','heatCp']}
    # Total retains the original unfiltered-forcing interpretation.
    wrapper.restore_baselines(); ekf.n_iters = n
    raw = np.stack([temperature+ekf.offset, ekf.oc_meas[:n], ekf.TOA_meas_artif[:n]], axis=1)[:, :, None]
    data['tot'] = solve(raw, anomaly_reference=preind_base+ekf.offset)
    # This is an internal absolute-temperature gauge, not the reporting datum.
    # Keep equilibrium/radiative calibration fixed when rebasing model inputs.
    component_reference = float(ekf.Teq1850-ekf.offset)
    component_settings = {}
    component_n=min(n,args.component_end_year-1849) if args.component_end_year else n
    model_metadata = None
    if args.component_source == 'processed':
        if args.component_end_year:
            raise ValueError('--component-end-year is only for the legacy MPI route')
        combined,model_metadata = processed_components(args.models,n,component_reference,out,args.model_toa,settings)
        for label in ['aer','ghg','nat']:
            data[label] = combined[label]
        data['model_residual'] = combined['model_residual']
        data['model_anthropogenic'] = combined['model_anthropogenic']
        # No ozone-only experiment was downloaded for this three-model set.
        data['ozone'] = (np.full(n,np.nan),np.full(n,np.nan))
    else:
        for label, experiment in [('aer','aer'), ('ozone','totalO3'), ('ghg','GHG'), ('nat','nat')]:
            wrapper.restore_baselines(); ekf.n_iters = component_n
            # These are simulated inputs, not CERES observations. Retain the
            # pre-CERES model TOA variance after 2000; ANT/total are unchanged.
            ekf.TOA_meas_artif_var[2001-1850:] = ekf._orig_TOA_meas_artif_var[1999-1850]
            # Components use the finalized anthropogenic ta4 noise settings.
            ekf.Q *= .2; ekf.cM2 *= .6; ekf.dfrA_float_var *= .35/40
            nf = len(ekf.data)
            raw_aod = ekf.data[:,3]*.001
            ekf.tsi = np.full(nf, ekf.sw_in)
            ekf.opt_depth = np.full(nf, ekf.raw_aod_from_inverse_denominator(ekf.involcavg))
            ekf.lCo2 = np.full(nf, ekf._orig_lCo2[0])
            ekf.anthro_clouds = np.full(nf, ekf._orig_anthro_clouds[0])
            if label == 'aer':ekf.anthro_clouds = ekf._orig_anthro_clouds.copy()
            elif label in ['ozone','ghg']:
                columns = ['O3'] if label=='ozone' else ['CO2','CH4','N2O','halogen']
                erf = selected_component_erf(ekf.data[:,0].astype(int),columns,ekf.eCO2_FORCING_PERCENTILE)
                # Apply the forcing anomaly to the same preindustrial carbon level.
                ekf.lCo2 = ekf._orig_lCo2[0] + (erf-erf[0])/ECO2_ERF_SLOPE
            elif label=='nat':
                ekf.tsi = ekf._orig_tsi.copy(); ekf.opt_depth = raw_aod
            # Recompute beta/radiative reference for the selected forcing, without
            # reapplying observational percentile helpers over this experiment.
            ekf.use_observational_forcing_percentiles = False
            ekf.precompute_coeffs(False)
            component_settings[label] = {key:float(getattr(ekf,key)) for key in settings}
            for key in settings:assert component_settings[label][key]==settings[key],(label,key)
            mean,se=solve(component_observations(experiment,component_n,component_reference))
            data[label]=(np.pad(mean,(0,n-component_n),constant_values=np.nan),np.pad(se,(0,n-component_n),constant_values=np.nan))
    # Keep the legacy, unused aerosol/ozone average as OHFa for schema parity.
    data['ohfa'] = ((data['aer'][0]+data['ozone'][0])/2, np.sqrt((data['aer'][1]**2+data['ozone'][1]**2)/2))
    if args.ohf=='ant-minus-ghg':
        # Covariance is unavailable: the quadrature width is explicitly approximate.
        data['ohf'] = (data['ant'][0]-data['ghg'][0],np.sqrt(data['ant'][1]**2+data['ghg'][1]**2))
    elif args.ohf=='model-residual':
        if args.component_source != 'processed': raise ValueError('model-residual requires processed inputs')
        data['ohf'] = data['model_residual']
    else:data['ohf'] = data['aer']
    flat = pd.DataFrame({'year':years}); assessment = pd.DataFrame(index=pd.Index(years,name='Year'))
    for key,label in [('ohf','OHF'),('ghg','GHG'),('ohfa','OHFa'),('ant','Ant'),('nat','Nat'),('tot','Tot')]:
        mean,se = data[key]
        assert mean.shape==(n,) and se.shape==(n,)
        valid_n=0 if key=='ohfa' and args.component_source=='processed' else n if key in ['ant','tot'] else component_n
        assert np.isfinite(mean[:valid_n]).all() and np.isfinite(se[:valid_n]).all() and (se[:valid_n]>=0).all()
        assert np.isnan(mean[valid_n:]).all() and np.isnan(se[valid_n:]).all()
        for p in [5,17,83,95,50]:
            values=mean+norm.ppf(p/100)*se
            flat[f'{key}_p{p}']=values; assessment[(label,str(p))]=values
    assessment.columns=pd.MultiIndex.from_tuples(assessment.columns,names=['variable','percentile'])
    flat.to_csv(flat_path,index=False)
    assessment.to_csv(out/'Nicklas_GMST_timeseries.csv')
    # Explicit alternatives do not match the plotter's canonical filename.
    variants = {
        'aerosol_proxy': data['aer'],
        'ant_minus_GHG': (data['ant'][0]-data['ghg'][0], np.sqrt(data['ant'][1]**2+data['ghg'][1]**2)),
        'aerosols_plus_ozone': (data['aer'][0]+data['ozone'][0], np.sqrt(data['aer'][1]**2+data['ozone'][1]**2)),
    }
    if args.component_source=='processed':
        variants.pop('aerosols_plus_ozone')
        variants['model_residual'] = data['model_residual']
    for name,(mean,se) in variants.items():
        alternative=assessment.copy()
        for p in [5,17,83,95,50]:alternative[('OHF',str(p))]=mean+norm.ppf(p/100)*se
        alternative.to_csv(out/f'Nicklas_GMST_{name}.csv')
    components=pd.DataFrame({'Year':years})
    for name,(mean,se) in data.items():components[name+'_mean']=mean;components[name+'_se']=se
    component_path=out/'Nicklas_components.csv'
    components.to_csv(component_path,index=False)

    metadata={'end_year':args.end_year,'component_end_year':1849+component_n,'model_TOA_uncertainty':'Pre-CERES 1999 variance held fixed from 2001 onward in model-only components; historical ANT/total unchanged','years':n,'historical_settings':settings,'component_settings':component_settings,'component_percentile_convention':'Aerosol uses historical aggregate percentile 85; GHG/ozone use the historical Gaussian-subset convention at percentile 17 with the same updated ERF release. Marginal subset percentiles are not additive joint quantiles.','OHF':args.ohf,'component_extension':'1980-2014 linear fit plus seeded residual noise through selected end year; retained from original script','ozone_correction':'o3 forcing replaces erroneous solar-column read','component_SE':'Native filter SE; historical ANT includes production 1.2 current-output multiplier','uncertainty':'Gaussian analytic marginal percentiles; not heads/tails ensemble uncertainty','OHF_alternatives':'aerosol_proxy; ant_minus_GHG; aerosols_plus_ozone (partial non-GHG forcing). Combined widths assume independent components; covariance is unavailable.', 'legacy_OHFa':'Unused legacy mean of aerosol/ozone runs; not a sum of forcing contributions'}
    metadata['temperature_anomaly_reference_years'] = [1850,1900]
    metadata['temperature_anomaly_reference_count'] = 51
    metadata['physical_equilibrium_reference_years'] = list(ekf.equilibrium_reference_years)
    metadata['temperature_reference_note'] = ('1850–1900 inclusive temperature anomalies; internal equilibrium and radiative calibration retained. Model anomalies are mapped to the unchanged equilibrium gauge. The separate first50-anomaly option is a TOA-flux sensitivity, not a temperature datum.')
    metadata['component_source'] = args.component_source
    metadata['processed_model_inputs'] = model_metadata
    if args.component_source=='processed':
        metadata.update(component_extension='None: real downloaded historical/scenario outputs.',
            component_end_year=args.end_year,
            legacy_OHFa='Unavailable: no ozone-only run in new model set; unused compatibility columns contain NaN.',
            OHF_alternatives='aerosol_proxy; ant_minus_GHG; model_residual (ALL-NAT-GHG). Quadrature conditional errors assume independence; across-model spread is retained.')
    metadata_path.write_text(json.dumps(metadata,indent=2))
    wrapper.restore_baselines()
    print('Wrote',out/'Nicklas_GMST_timeseries.csv')

if __name__=='__main__':main()
