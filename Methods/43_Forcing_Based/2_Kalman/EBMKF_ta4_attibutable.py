"""Figure 16 attribution experiment, using the finalized ta4 historical model.

The comparison ends in 2023 to match the supplied multi-method assessments.
Single-forcing MPI experiments end in 2014; 2015--2023 retain the original
1980--2014 linear-trend plus seeded-residual-noise extension. These component
experiments are illustrative attribution runs, not the heads/tails ensemble.
Run from any directory; --output-dir permits an isolated review run.
"""
from pathlib import Path
import argparse, importlib, json, os, sys
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
    temperature = temperature - np.mean(temperature[:50]) + preind_base + ekf.offset
    ocean = (ocean - ocean[0]) / ekf.zJ_from_W
    return np.stack([extend(x, 1849+n) for x in [temperature, ocean, flux]], axis=1)[:, :, None]


def solve(observations):
    mean, variance, _, _ = ekf.ekf_run(observations, len(observations), retPs=3, diagnostics=False)
    return mean.copy() - ekf.Teq1850, nonnegative_sqrt(variance).copy()


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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path, default=None, help='Optional isolated output folder; default updates the established assessment CSV')
    ap.add_argument('--end-year', type=int, default=2023, choices=[2023, 2024, 2025])
    ap.add_argument('--component-end-year', type=int, choices=[2014], default=None, help='Stop model-only components at available 2014 data; preserve ANT/total through end-year')
    ap.add_argument('--ohf', choices=['aerosol-proxy', 'ant-minus-ghg'], default='aerosol-proxy')
    args = ap.parse_args()
    os.chdir(ROOT)
    out = args.output_dir.resolve() if args.output_dir else ROOT/'Methods/43_Forcing_Based/3_Human_Induced/anthropogenic-warming-assessment/results'
    out.mkdir(parents=True, exist_ok=True)
    flat_path = out/'Nicklas_GMST_timeseries_flat.csv' if args.output_dir else Path(__file__).resolve().parent/'Nicklas_GMST_timeseries.csv'
    metadata_path = out/'attribution-metadata.json' if args.output_dir else Path(__file__).resolve().parent/'attribution-metadata.json'
    obs = pd.read_csv(ROOT/'Common_Data/HadCRUT.5.1.csv').query('Time <= @args.end_year')
    years = obs.Time.to_numpy(dtype=int); n = len(years)
    temperature = obs.Anomaly.to_numpy(); preind_base = float(np.mean(temperature[:50]))
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
    data['tot'] = solve(raw)
    component_settings = {}
    component_n=min(n,args.component_end_year-1849) if args.component_end_year else n
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
        mean,se=solve(component_observations(experiment,component_n,preind_base))
        data[label]=(np.pad(mean,(0,n-component_n),constant_values=np.nan),np.pad(se,(0,n-component_n),constant_values=np.nan))
    # Keep the legacy, unused aerosol/ozone average as OHFa for schema parity.
    data['ohfa'] = ((data['aer'][0]+data['ozone'][0])/2, np.sqrt((data['aer'][1]**2+data['ozone'][1]**2)/2))
    if args.ohf=='ant-minus-ghg':
        # Covariance is unavailable: the quadrature width is explicitly approximate.
        data['ohf'] = (data['ant'][0]-data['ghg'][0],np.sqrt(data['ant'][1]**2+data['ghg'][1]**2))
    else:data['ohf'] = data['aer']
    flat = pd.DataFrame({'year':years}); assessment = pd.DataFrame(index=pd.Index(years,name='Year'))
    for key,label in [('ohf','OHF'),('ghg','GHG'),('ohfa','OHFa'),('ant','Ant'),('nat','Nat'),('tot','Tot')]:
        mean,se = data[key]
        assert mean.shape==(n,) and se.shape==(n,)
        valid_n=n if key in ['ant','tot'] else component_n
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
    for name,(mean,se) in variants.items():
        alternative=assessment.copy()
        for p in [5,17,83,95,50]:alternative[('OHF',str(p))]=mean+norm.ppf(p/100)*se
        alternative.to_csv(out/f'Nicklas_GMST_{name}.csv')
    components=pd.DataFrame({'Year':years})
    for name,(mean,se) in data.items():components[name+'_mean']=mean;components[name+'_se']=se
    component_path=out/'Nicklas_components.csv'
    components.to_csv(component_path,index=False)

    metadata={'end_year':args.end_year,'component_end_year':1849+component_n,'model_TOA_uncertainty':'Pre-CERES 1999 variance held fixed from 2001 onward in model-only components; historical ANT/total unchanged','years':n,'historical_settings':settings,'component_settings':component_settings,'component_percentile_convention':'Aerosol uses historical aggregate percentile 85; GHG/ozone use the historical Gaussian-subset convention at percentile 17 with the same updated ERF release. Marginal subset percentiles are not additive joint quantiles.','OHF':args.ohf,'component_extension':'1980-2014 linear fit plus seeded residual noise through selected end year; retained from original script','ozone_correction':'o3 forcing replaces erroneous solar-column read','component_SE':'Native filter SE; historical ANT includes production 1.2 current-output multiplier','uncertainty':'Gaussian analytic marginal percentiles; not heads/tails ensemble uncertainty','OHF_alternatives':'aerosol_proxy; ant_minus_GHG; aerosols_plus_ozone (partial non-GHG forcing). Combined widths assume independent components; covariance is unavailable.', 'legacy_OHFa':'Unused legacy mean of aerosol/ozone runs; not a sum of forcing contributions'}
    metadata_path.write_text(json.dumps(metadata,indent=2))
    wrapper.restore_baselines()
    print('Wrote',out/'Nicklas_GMST_timeseries.csv')

if __name__=='__main__':main()
