"""Shared v11 repairs; variant state equations and noise choices remain local."""
from pathlib import Path
import numpy as np

GLOBAL_FLUX_TO_ZJ = 5.1006447 * 3.154

def load_inputs(module_path):
    root = Path(module_path).resolve().parents[3]
    data = np.genfromtxt(root/'Common_Data/toyKFmodelData9c.csv', delimiter=',')
    if data.shape[1] != 9 or not np.isfinite(data).all():
        raise ValueError('Invalid frozen Kalman input')
    if not np.array_equal(data[:,0],np.arange(1850,1850+len(data))):
        raise ValueError('Frozen Kalman input must be annual from 1850')
    return data

def precompute(namespace, print_them=False):
    g=namespace
    temps=g['temps'];data=g['data'];years=np.arange(1850,1850+len(temps))
    # The model's absolute temperature series sets these references.
    # The separate HadCRUT5.0.2 capacity file below NEVER sets Teq1850/T02.
    g['equilibrium_reference_years']=g.get('equilibrium_reference_years',(1850,1899))
    g['radiation_reference_years']=g.get('radiation_reference_years',(1999,2005))
    eq_start,eq_end=g['equilibrium_reference_years']
    rad_start,rad_end=g['radiation_reference_years']
    eq_ref=(years>=eq_start)&(years<=eq_end)
    ref=(years>=rad_start)&(years<=rad_end)
    if eq_ref.sum()!=eq_end-eq_start+1 or ref.sum()!=rad_end-rad_start+1:
        raise ValueError('Temperature reference period is incomplete')
    g['T02']=float(np.mean(temps[ref]))
    g['Teq1850']=float(np.mean(temps[eq_ref]))
    g['pindavg']=g['Teq1850']
    g['aerosol_reference']=float(np.mean(data[(data[:,0]>=rad_start)&(data[:,0]<=rad_end),7]))
    g['volcanic_aod_reference']=float(np.mean(data[(data[:,0]>=rad_start)&(data[:,0]<=rad_end),3]*.001))
    g['inbndf']=g['a_refl']*g['g_refl']*9.068
    swref=g['sw_in']*g['inbndf']/(9.7279+g['volcanic_aod_reference'])
    g['reference_absorbed_sw']=swref
    if 'volcanic_target_efficiency' in g:
        # v10 reference-preserving mapping. Do not scale raw AOD directly:
        # its 1999–2005 reference transmission and the SW calibration stay fixed.
        target=float(g['volcanic_target_efficiency'])
        if not np.isfinite(target) or target<=0:
            raise ValueError('Volcanic target efficiency must be positive')
        data_ref=(data[:,0]>=rad_start)&(data[:,0]<=rad_end)
        aerosol_ref_anomaly=float(np.mean(g['anthro_clouds'][:len(data)][data_ref]))
        g['volcanic_reference_absorbed_sw']=swref*(1+aerosol_ref_anomaly/swref)
        g['volcanic_denominator_offset']=9.7279
        g['volcanic_aod_scale']=(target*(9.7279+g['volcanic_aod_reference'])
                                 /g['volcanic_reference_absorbed_sw'])
        # Refresh the historical reciprocal mean if the target is edited live.
        raw_aod=data[:,3]*.001
        g['involcavg']=float(np.mean(1/(9.7279+g['volcanic_aod_reference']
                            +g['volcanic_aod_scale']*(raw_aod-g['volcanic_aod_reference']))))

    # Manuscript beta2 = dfaA (cloud), beta3 = dfaS (surface).
    # Both divide by the COMPLETE reference absorbed SW, including cloud,
    # surface and volcanic transmission, rather than incident solar alone.
    g['dfaA'] = g['fdbkA'] / g['reference_absorbed_sw']
    g['dfaS'] = g['fdbkS'] / g['reference_absorbed_sw']
    g['powp1']=g['fdbkW']*4/3.22
    g['B1B0']=12.74/g['sig']/g['T02']**(4-g['powp1'])
    g['rad1850']=swref*(1+g['dfaA']*(g['Teq1850']-g['T02'])+g['anthro_clouds'][0]/swref)*(1+g['dfaS']*(g['Teq1850']-g['T02']))
    carbon_reference=g['lCo2'][0] if g.get('v11_carbon_reference_actual',False) else 2.444
    g['B1']=g['rad1850']/5.670e-8/g['Teq1850']**(4-g['powp1'])+g['B1B0']*carbon_reference
    g['B0']=g['B1B0']/g['B1'];g['outbndf']=g['sig']*g['B1']
    g['capacity_per_m_ocean_depth']=.71*1030*4180/3.154e7
    g['Cs']=g.get('upper_ocean_capacity_reference',8.2);g['Cd']=g['capacity_per_m_ocean_depth']*g['deepdepth']
    g['shaldepth']=g['Cs']/g['capacity_per_m_ocean_depth']
    root=Path(g['__file__']).resolve().parents[3]
    refdata=np.genfromtxt(root/'Common_Data/HadCRUT5_capacity_reference_5.0.2.0.csv',delimiter=',',skip_header=1)
    warming=np.mean(refdata[(refdata[:,0]>=2015)&(refdata[:,0]<=2024),1])-np.mean(refdata[(refdata[:,0]>=1961)&(refdata[:,0]<=1980),1])
    if not np.isfinite(warming) or warming<=0:raise ValueError('Invalid heat-capacity reference warming')
    g['non_ocean_heat_capacity']=42/(GLOBAL_FLUX_TO_ZJ*warming)
    # Default total storage capacity; preserve the explicitly fitted MPI variant.
    g['heatCp'] = g.get('v11_fitted_heat_capacity',
                         g['Cs'] + g['non_ocean_heat_capacity'])
    if print_them:print('v11 configuration:',{k:g[k] for k in ['Cs','Cd','heatCp','fdbkS','fdbkA','reference_absorbed_sw','aerosol_reference']})

def complete_hidden_jacobian(update, x, year, *forcing):
    """Complex-step derivative of the variant's complete next-state map.

    This differentiates the actual augmented feedback/flux/exchange equations,
    including baseline recalibration in the five-state variants, without
    subtractive cancellation. Independent central differences test this map.
    """
    state=np.asarray(x);flat=state.reshape(-1)
    result=np.empty((len(flat),len(flat)))
    for j in range(len(flat)):
        perturbed=flat.astype(complex);perturbed[j]+=1e-20j
        result[:,j]=np.imag(update(perturbed.reshape(state.shape),year,*forcing)).reshape(-1)/1e-20
    return result

def nonnegative_sqrt(values):
    values=np.asarray(values)
    if not np.isfinite(values).all() or np.any(values < -1e-9):
        raise ValueError('Invalid variance; absolute values must not hide covariance errors')
    return np.sqrt(np.maximum(values,0))

def check_covariances(g,n):
    for key in ['P','Pminus','Phat']:
        matrices=np.asarray(g[key][:n])
        if not np.isfinite(matrices).all():raise ValueError(f'{key}: nonfinite covariance')
        symmetric=(matrices+matrices.swapaxes(-1,-2))/2
        eig=np.linalg.eigvalsh(symmetric)
        if eig.min() < -1e-8:
            k=int(np.unravel_index(np.argmin(eig),eig.shape)[0])
            raise ValueError(f'{key}: covariance not PSD at {1850+k}; minimum eigenvalue {eig.min()}')

def prepare_observed_two_state(model, years, temperature, uncertainty):
    years=np.asarray(years);n=len(years)
    if not np.array_equal(years,np.arange(1850,1850+n)) or n>len(model.data):
        raise ValueError('Kalman input years exceed the frozen forcing data; rebuild toyKFmodelData9c.csv first')
    model.n_iters=n
    reference=float(np.mean(model.data[:50,1]))
    model.temps=np.asarray(temperature)-np.mean(temperature[:50])+reference+model.offset
    model.R_tvar=((np.asarray(uncertainty[1])-np.asarray(uncertainty[0]))/4)**2
    model.oc_meas=model.data[:n,5]/model.zJ_from_W
    model.precompute_coeffs(False)
    model.observ=np.column_stack([model.temps,model.oc_meas])
    for name in ['xhat','P','F','xhatminus','Pminus','K','xhathat','Phat','Khat','Shat','S','xblind','lml','y','qqy']:
        old=getattr(model,name);setattr(model,name,np.zeros((n,)+old.shape[1:]))
    model.qqykerall.clear()

def scenario_aerosols(model, years, experiment):
    """Scenario ARI+ACI at the core's aerosol percentile, centred consistently.

    Match the historical piecewise-linear p05 / best-at-50 / p95 rule.
    SSP files provide component curves, so sum ARI and ACI in each file before
    interpolation. This is an approximate combined-percentile trajectory,
    not a reconstructed joint aerosol probability distribution.
    Cores without a percentile control retain the previous central estimate.
    """
    import pandas as pd
    import config
    percentile=float(getattr(model, 'AEROSOL_FORCING_PERCENTILE', 50.0))
    if not np.isfinite(percentile) or not 5<=percentile<=95:
        raise ValueError('AEROSOL_FORCING_PERCENTILE must be between 5 and 95')
    requested=np.asarray(years)
    if not np.isfinite(requested).all() or not np.array_equal(requested,requested.astype(int)):
        raise ValueError('Scenario aerosol forcing requires integer calendar years')
    attrs=experiment.split('_')
    scenario=attrs[2].lower() if attrs[1]=='ESM1-2-LR' else 'ssp245'
    root=Path(config.CLIMATE_DATA_PATH)/'SSP_inputdata/ERFs-Smith-ar6'
    def read_curve(suffix):
        source=root/f'ERF_{scenario}_1750-2500{suffix}.csv'
        frame=pd.read_csv(source).set_index('year')
        if frame.index.has_duplicates:
            raise ValueError(f'Duplicate years in {source.name}')
        selected=frame.loc[requested.astype(int)]
        values=(selected['aerosol-radiation_interactions'].to_numpy()
                +selected['aerosol-cloud_interactions'].to_numpy())
        if not np.isfinite(values).all():
            raise ValueError(f'Nonfinite aerosol forcing in {source.name}')
        return values
    central=read_curve('')
    if percentile==50:
        chosen=central
    elif percentile<50:
        low=read_curve('_pc05')
        if np.any(low>central+1e-10):
            raise ValueError('Unordered SSP p05/best aerosol curves')
        chosen=low+(central-low)*((percentile-5)/45)
    else:
        high=read_curve('_pc95')
        if np.any(central>high+1e-10):
            raise ValueError('Unordered SSP best/p95 aerosol curves')
        chosen=central+(high-central)*((percentile-50)/45)
    return chosen-model.aerosol_reference

def volcanic_erf(model, raw_aod, solar):
    """Reference-centred absorbed-SW anomaly, using the v11 AOD mapping."""
    denominator = getattr(model, 'volcanic_denominator', lambda raw: 9.7279+np.asarray(raw))
    return np.asarray(solar)*model.inbndf*(
        1/denominator(raw_aod)-1/denominator(model.volcanic_aod_reference))

def covariance_inverse(covariance):
    """Symmetric eigensolver inverse for resolved covariance modes.

    Preserve valid semidefinite/null modes without borrowing another time's
    gain. Genuine indefinite covariances remain errors; no eigenvalue clipping.
    NumPy's default pseudoinverse cutoff is 1e-15 times the largest mode.
    """
    covariance = np.asarray(covariance)
    if not np.isfinite(covariance).all():
        raise ValueError('Nonfinite covariance in gain calculation')
    symmetric = (covariance + covariance.T) / 2
    if np.linalg.eigvalsh(symmetric).min() < -1e-8:
        raise ValueError('Indefinite covariance in gain calculation')
    return np.linalg.pinv(symmetric, hermitian=True)

def aerosol_forcing_at_percentile(years, percentile):
    """Piecewise-linear interpolation of the supplied total-aerosol curves.

    The aggregate 'aerosol' column already represents ARI+ACI, including the
    source's combined uncertainty. Do not sum marginal component quantiles.
    Annual labels such as 2024.5 denote calendar year 2024, not a shifted year.
    Continue each 2015–2024 linear trend by one year for 2025, as for frozen 9c.
    """
    import pandas as pd
    percentile=float(percentile)
    if not np.isfinite(percentile) or not 5<=percentile<=95:
        raise ValueError('AEROSOL_FORCING_PERCENTILE must be between 5 and 95')
    requested=np.asarray(years)
    if not np.isfinite(requested).all() or not np.array_equal(requested,requested.astype(int)):
        raise ValueError('Aerosol forcing requires integer calendar years')
    curves=[]
    for tag in ['p05','best','p95']:
        path=Path(__file__).with_name(f'ERF_{tag}_aggregates_1750-2024.csv')
        frame=pd.read_csv(path,index_col=0)
        frame.index=np.floor(frame.index.to_numpy(dtype=float)).astype(int)
        if frame.index.has_duplicates or not np.array_equal(frame.index,np.arange(1750,2025)):
            raise ValueError(f'Unexpected annual forcing coverage in {path.name}')
        values=frame['aerosol'].astype(float).copy()
        if not np.isfinite(values).all():raise ValueError(f'Nonfinite forcing in {path.name}')
        fit_years=np.arange(2015,2025)
        slope,intercept=np.polyfit(fit_years,values.loc[fit_years].to_numpy(),1)
        values.loc[2025]=slope*2025+intercept
        if np.any(requested<1750) or np.any(requested>2025):
            raise ValueError('Observed aerosol-percentile curves cover 1750–2025 only')
        curves.append(values.loc[requested.astype(int)].to_numpy())
    low,central,high=curves
    if np.any(low>central+1e-10) or np.any(central>high+1e-10):
        raise ValueError('Unordered aggregate aerosol uncertainty curves')
    if percentile<=50:
        return low+(central-low)*((percentile-5)/45)
    return central+(high-central)*((percentile-50)/45)

def _apply_historical_aerosol_percentile(model):
    """Apply a core's own setting, accepting a module or its globals dict.

    Called before coefficient calibration so standalone runs and wrappers
    use the same selected aerosol history and reference-centred coefficients.
    """
    g=model if isinstance(model,dict) else vars(model)
    percentile=g['AEROSOL_FORCING_PERCENTILE']
    if '_ta2_frozen_aerosol' not in g:
        g['_ta2_frozen_aerosol']=g['data'][:,7].copy()
    chosen=aerosol_forcing_at_percentile(g['data'][:,0],percentile)
    if percentile==50:
        if not np.allclose(chosen,g['_ta2_frozen_aerosol'],atol=1e-8,rtol=0):
            raise ValueError('Best ERF CSV and frozen 9c aerosol series disagree; check source releases')
        chosen=g['_ta2_frozen_aerosol'].copy()
    g['data']=g['data'].copy()
    g['data'][:,7]=chosen
    start,end=g.get('radiation_reference_years',(1999,2005))
    reference=(g['data'][:,0]>=start)&(g['data'][:,0]<=end)
    g['aerosol_reference']=float(chosen[reference].mean())
    g['anthro_clouds']=chosen-g['aerosol_reference']


# Workbook ERF_best_aggregates_1750-2024.xlsx: Z sums these columns, AA is
# 10**((Z+V100)/V99). No aerosol, solar, volcanic or duplicate minor subtotal.
ECO2_ERF_COMPONENTS = ('CO2','CH4','N2O','O3','contrails','land_use',
                       'BC_on_snow','H2O_stratospheric','halogen')
ECO2_ERF_SLOPE = 12.743331707472763  # workbook V99
ECO2_ERF_INTERCEPT = 31.14922366816877  # workbook V100

def non_aerosol_forcing_gaussian(years):
    """Approximate independent-component Gaussian sum, in W/m2.

    Sigma_i=(p95_i-p05_i)/(2*Phi^-1(.95)); sigma_sum is quadrature.
    This is not a source-provided percentile of the joint non-aerosol sum.
    """
    import pandas as pd
    from scipy.stats import norm
    years=np.asarray(years)
    if not np.isfinite(years).all() or not np.array_equal(years,years.astype(int)) or np.any(years<1750) or np.any(years>2025):
        raise ValueError('Non-aerosol forcing years must be integer 1750–2025')
    frames=[]
    for tag in ['p05','best','p95']:
        frame=pd.read_csv(Path(__file__).with_name(f'ERF_{tag}_aggregates_1750-2024.csv'),index_col=0)
        frame.index=np.floor(frame.index.to_numpy(dtype=float)).astype(int)
        if not np.array_equal(frame.index,np.arange(1750,2025)):
            raise ValueError('Unexpected forcing-file years')
        frames.append(frame.loc[:,list(ECO2_ERF_COMPONENTS)])
    low,best,high=frames
    mean=best.sum(axis=1)
    widths=high.to_numpy()-low.to_numpy()
    if not np.isfinite(widths).all() or np.any(widths<0):raise ValueError('Invalid component forcing ranges')
    sigma=pd.Series(np.sqrt(np.sum((widths/(2*norm.ppf(.95)))**2,axis=1)),index=best.index)
    fit_years=np.arange(2015,2025)
    for values in [mean,sigma]:
        values.loc[2025]=np.polyval(np.polyfit(fit_years,values.loc[fit_years],1),2025)
    if np.any(sigma<0):raise ValueError('Negative extrapolated forcing sigma')
    return mean.loc[years.astype(int)].to_numpy(),sigma.loc[years.astype(int)].to_numpy()

def _apply_historical_eco2_percentile(model):
    """Keep frozen best eCO2 exactly; perturb log10(eCO2) by deltaERF/V99."""
    from scipy.stats import norm
    g=model if isinstance(model,dict) else vars(model)
    p=float(g['eCO2_FORCING_PERCENTILE'])
    if not np.isfinite(p) or not 5<=p<=95:raise ValueError('eCO2_FORCING_PERCENTILE must be between 5 and 95')
    if '_ta2_frozen_eco2' not in g:g['_ta2_frozen_eco2']=g['data'][:,2].copy()
    years=g['data'][:,0];mean,sigma=non_aerosol_forcing_gaussian(years)
    available=years<=2024
    expected=10**((mean[available]+ECO2_ERF_INTERCEPT)/ECO2_ERF_SLOPE)
    if not np.allclose(expected,g['_ta2_frozen_eco2'][available],atol=1e-6,rtol=0):
        raise ValueError('Frozen eCO2 and workbook forcing formula/source disagree')
    # The workbook extends 2025 via the eCO2/observed-CO2 trend, not total ERF.
    # Preserve that central continuation; extend the uncertainty width only.
    chosen=g['_ta2_frozen_eco2']*10**(norm.ppf(p/100)*sigma/ECO2_ERF_SLOPE)
    g['data']=g['data'].copy();g['data'][:,2]=chosen;g['lCo2']=np.log10(chosen)

def parameter_choice_density(model, reference_year=2024, print_report=True):
    """Local density of chosen forcing amplitudes/cloud feedback, not probability.

    Approximate independent Gaussians; aerosol sigma from the aggregate range,
    non-aerosol sigma from independent component widths. Actual selected aerosol
    amplitude is evaluated, including the earlier piecewise-percentile mapping.
    Units are defined in forcing coordinates, NOT in ppm or percentile ranks.
    """
    import json
    from scipy.stats import norm
    g=model if isinstance(model,dict) else vars(model)
    year=np.array([reference_year])
    aerosol_mean=float(aerosol_forcing_at_percentile(year,50)[0])
    aerosol_sigma=float((aerosol_forcing_at_percentile(year,95)[0]-aerosol_forcing_at_percentile(year,5)[0])/(2*norm.ppf(.95)))
    aerosol_value=float(aerosol_forcing_at_percentile(year,g['AEROSOL_FORCING_PERCENTILE'])[0])
    mu,sigma=non_aerosol_forcing_gaussian(year)
    co2_value=float(mu[0]+norm.ppf(g['eCO2_FORCING_PERCENTILE']/100)*sigma[0])
    # AR6 Table 7.10: cloud feedback 0.42, very likely (5–95%) [-0.10,0.94].
    cloud_sigma=(.94-(-.10))/(2*norm.ppf(.95))
    entries={}
    for name,value,mean,sd in [('aerosol_ERF',aerosol_value,aerosol_mean,aerosol_sigma),
                             ('non_aerosol_ERF',co2_value,float(mu[0]),float(sigma[0])),
                             ('cloud_feedback',float(g['fdbkA']),.42,cloud_sigma)]:
        if sd<=0:raise ValueError('Density requires positive uncertainty widths')
        entries[name]={'chosen':value,'mean':mean,'sigma':sd,'z':(value-mean)/sd,
                       'local_pdf':float(norm.pdf(value,loc=mean,scale=sd))}
    logpdf=sum(float(norm.logpdf(e['chosen'],loc=e['mean'],scale=e['sigma'])) for e in entries.values())
    z2=sum(e['z']**2 for e in entries.values())
    report={'reference_year':reference_year,'assumption':'independent approximate Gaussian forcing amplitudes and cloud feedback',
            'AEROSOL_FORCING_PERCENTILE':float(g['AEROSOL_FORCING_PERCENTILE']),
            'eCO2_FORCING_PERCENTILE':float(g['eCO2_FORCING_PERCENTILE']),
            'marginals':entries,'joint_local_density':float(np.exp(logpdf)),
            'joint_log_density':logpdf,'density_units':'(W/m2)^-2 * (W/m2/K)^-1',
            'density_relative_to_joint_peak':float(np.exp(-.5*z2)),
            'standardized_joint_density':float(np.exp(-.5*z2)/(2*np.pi)**1.5),
            'note':'Density, not probability mass or trajectory likelihood. Independence is an assumption, not an assessment of correlations.',
            'cloud_source':'IPCC AR6 WGI Table 7.10: https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7/'}
    report['relative_parameter_likelihood']=report['density_relative_to_joint_peak']
    # Match the density-to-peak ratio exp(-z^2/2) in one dimension.
    equivalent_z=float(np.sqrt(z2))
    report['gaussian_1d_equivalent_abs_z']=equivalent_z
    report['gaussian_1d_equivalent_lower_percentile']=float(100*norm.cdf(-equivalent_z))
    report['gaussian_1d_equivalent_upper_percentile']=float(100*norm.cdf(equivalent_z))
    if print_report:
        ratio=report['relative_parameter_likelihood']
        print(f"Relative parameter likelihood (peak = 1): {ratio:.6f} ({100*ratio:.2f}% of peak)")
        print(f"  1D Gaussian equivalent (same density/peak): |z|={equivalent_z:.3f}; lower-tail percentile={report['gaussian_1d_equivalent_lower_percentile']:.2f}% (upper-tail counterpart={report['gaussian_1d_equivalent_upper_percentile']:.2f}%)")
        print(f"  Aerosol percentile={g['AEROSOL_FORCING_PERCENTILE']:g}; eCO2 ERF percentile={g['eCO2_FORCING_PERCENTILE']:g}; fdbkA={g['fdbkA']:g}")
        print(f"  Gaussian z scores at {reference_year}: aerosol={entries['aerosol_ERF']['z']:.3f}, non-aerosol ERF={entries['non_aerosol_ERF']['z']:.3f}, cloud feedback={entries['cloud_feedback']['z']:.3f}")
        print("  Independent Gaussian approximation; ratio to the best-ERF / AR6-cloud peak, not a probability of model correctness.")
    return report

def reference_preserving_volcanic_denominator(model, raw_aod):
    """v10 attenuation in raw-AOD coordinates; same reference transmission."""
    g=model if isinstance(model,dict) else vars(model)
    return (g['volcanic_denominator_offset']+g['volcanic_aod_reference']
            +g['volcanic_aod_scale']*(np.asarray(raw_aod)-g['volcanic_aod_reference']))

def inverse_volcanic_denominator(model, inverse_denominator):
    """Return RAW AOD after averaging the mapped reciprocal denominator."""
    g=model if isinstance(model,dict) else vars(model)
    return (g['volcanic_aod_reference']
            +(1/np.asarray(inverse_denominator)-g['volcanic_denominator_offset']
              -g['volcanic_aod_reference'])/g['volcanic_aod_scale'])

def plot_two_state_components(model, years, observations, title):
    """Diagnostic views of the actual two-state run, including future inputs.

    Uses the active forcing/variance arrays rather than the 176-year frozen
    input table, so future runs retain their full time axis and correct drivers.
    """
    import matplotlib.pyplot as plt
    years=np.asarray(years); n=len(years)
    obs=np.asarray(observations)
    fig,axes=plt.subplots(4,2,figsize=(13,13),sharex=True,layout='constrained')
    ax=axes.flat
    temp_offset=model.offset
    ax[0].plot(years,obs[:,0]-temp_offset,color='0.7',label='Input temperature')
    ax[0].plot(years,model.xhat[:n,0]-temp_offset,label='Filtered')
    ax[0].plot(years,model.xhathat[:n,0]-temp_offset,label='Retrospective')
    ax[0].plot(years,model.xblind[:n,0]-temp_offset,label='Blind EBM')
    benchmark=np.full(n,np.nan)
    for i in range(10,n-10):
        benchmark[i]=.5*(obs[i-10:i+10,0].mean()+obs[i-9:i+11,0].mean())-temp_offset
    ax[0].plot(years,benchmark,'k',lw=1.5,label='Centered 20-year benchmark')
    ax[0].set_ylabel('Temperature anomaly (°C)');ax[0].legend(fontsize=8)
    ax[1].plot(years,obs[:,1]*model.zJ_from_W,color='0.7',label='Input OHC')
    ax[1].plot(years,model.xhat[:n,1]*model.zJ_from_W,label='Filtered')
    ax[1].plot(years,model.xhathat[:n,1]*model.zJ_from_W,label='Retrospective')
    ax[1].plot(years,model.xblind[:n,1]*model.zJ_from_W,label='Blind EBM')
    ax[1].set_ylabel('OHC (ZJ)');ax[1].legend(fontsize=8)
    ax[2].plot(years,model.opt_depth[:n]);ax[2].set_ylabel('Active raw AOD (after smoothing)')
    ax[3].plot(years,model.tsi[:n]);ax[3].set_ylabel('Incident solar (W m$^{-2}$)')
    ax[4].plot(years,10**model.lCo2[:n]);ax[4].set_ylabel('Effective CO₂ (ppm)')
    ax[5].plot(years,model.anthro_clouds[:n]);ax[5].set_ylabel('Centred ARI + ACI (W m$^{-2}$)')
    state=model.xhat[:n];t=state[:,0]
    sw=(model.tsi[:n]*model.inbndf/model.volcanic_denominator(model.opt_depth[:n])
        *(1+model.dfaA*(t-model.T02)+model.anthro_clouds[:n]/model.reference_absorbed_sw)
        *(1+model.dfaS*(t-model.T02)))
    lw=model.outbndf*t**(4-model.powp1)*(1-model.B0*model.lCo2[:n])
    ax[6].plot(years,sw,label='Absorbed SW');ax[6].plot(years,lw,label='Outgoing LW')
    ax[6].set_ylabel('Radiative flux (W m$^{-2}$)');ax[6].legend(fontsize=8)
    ax[7].plot(years,sw-lw,label='Net radiation')
    # Use compute_update to avoid inventing a different ocean-exchange convention.
    updates=np.array([model.compute_update(state[i],i) for i in range(n)])
    ax[7].plot(years,updates[:,1],label='OHC storage rate')
    ax[7].plot(years,updates[:,0]*(model.heatCp-model.Cs),label='Non-ocean storage rate')
    ax[7].set_ylabel('Energy balance (W m$^{-2}$)');ax[7].legend(fontsize=8)
    for item in axes.flat:item.grid(alpha=.2)
    for item in axes[-1]:item.set_xlabel('Year')
    fig.suptitle(title+'\nAerosol p'+str(model.AEROSOL_FORCING_PERCENTILE)
                 +', observational eCO₂ p'+str(model.eCO2_FORCING_PERCENTILE)+', fdbkA='+str(model.fdbkA)
                 +'; aerosol percentile applies to historical and SSP inputs')
    fig2,aa=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    from scipy.stats import norm
    for j,label in enumerate(['Temperature','OHC']):
        aa[j].hist(np.asarray(model.qqy[1:n,j]).ravel(),bins=24,density=True,color='0.7',label='Innovations')
        aa[j].plot(model.xnorm,model.qqyker[j],color='royalblue',label='Observation mixture')
        aa[j].plot(model.xnorm,norm.pdf(model.xnorm),'k--',label='Standard normal')
        aa[j].set(title=label,xlabel='Standardised innovation',ylabel='Density');aa[j].legend(fontsize=8)
    fig2.suptitle(title+' — innovation diagnostics')
    return [fig,fig2]


def four_state_diagnostic_series(model, years, observations):
    """Reconstruct ta4's physical fluxes without changing active state/inputs.

    The prediction for year k uses posterior state k-1 and forcing k. The
    correction is a same-state counterfactual: (SW multiplier - 1)*base SW,
    not a difference between two independently evolved temperature paths.
    """
    years=np.asarray(years);n=len(years)
    obs=np.asarray(observations).reshape(n,3)
    filtered=np.asarray(model.xhat[:n]).reshape(n,4).copy()
    retro=np.asarray(model.xhathat[:n]).reshape(n,4).copy()
    def fluxes(states,start):
        t=states[:,0];sl=slice(start,start+len(states))
        sw0=(model.tsi[sl]*model.inbndf/model.volcanic_denominator(model.opt_depth[sl])
             *(1+model.dfaA*(t-model.T02)+model.anthro_clouds[sl]/model.reference_absorbed_sw)
             *(1+model.dfaS*(t-model.T02)))
        lw=model.outbndf*t**(4-model.powp1)*(1-model.B0*model.lCo2[sl])
        return sw0,lw
    result={'year':years.copy(),'observations':obs.copy(),'filtered':filtered,'retrospective':retro}
    # Match the actual forward observation covariance, not prediction covariance S.
    # The filter uses 9*TOA_var plus N*cM2[2,2] as assigned TOA observation noise.
    result['TOA_nominal_sigma']=np.sqrt(9*np.asarray(model.TOA_meas_artif_var[:n]))
    result['TOA_effective_sigma']=np.sqrt(9*np.asarray(model.TOA_meas_artif_var[:n])+model.N*model.cM2[2,2])
    result['TOA_interval_definition']='Approximate Gaussian 95% assigned observation interval: input +/- 1.96*sqrt(9*TOA_var + N*cM2[2,2]); excludes forecast-state uncertainty.'

    for label,states in [('filtered',filtered),('retrospective',retro)]:
        sw0,lw=fluxes(states[:-1],1)
        correction=(states[:-1,3]-1)*sw0
        sw=sw0+correction
        for key,value in [('base_SW',sw0),('SW_correction',correction),('absorbed_SW',sw),('outgoing_LW',lw),('net_flux',sw-lw)]:
            result[label+'_'+key]=np.r_[np.nan,value]
    # Exact no-correction replay with the active post-wrapper parameter values.
    blind=np.zeros((n,4));blind[0]=[model.Teq1850,model.oc_meas[0],0,1]
    for i in range(1,n):blind[i]=model.compute_update_hidden(blind[i-1,:,None],i)[:,0]
    result['blind']=blind
    sw0,lw=fluxes(filtered[:-1],1)
    states=filtered[:-1]
    ocean_temp=(states[:,1]-model.Cs*(states[:,0]-model.Teq1850))/model.Cd+model.oc1850
    exchange=model.gad*(states[:,0]-model.Teq1850-ocean_temp+model.oc1850)
    net=result['filtered_net_flux'][1:]
    dt=(net-exchange)/model.heatCp
    result['ocean_storage_rate']=np.r_[np.nan,model.Cs*dt+exchange]
    result['non_ocean_storage_rate']=np.r_[np.nan,(model.heatCp-model.Cs)*dt]
    np.testing.assert_allclose(net,result['ocean_storage_rate'][1:]+result['non_ocean_storage_rate'][1:],rtol=0,atol=1e-12)
    np.testing.assert_allclose(net,np.asarray(model.xhatminus[1:n,2,0]),rtol=0,atol=1e-10)
    result['correction_definition']='(previous inferred SW multiplier - 1) * unadjusted absorbed SW at the same previous state and current-year forcing'
    return result


def plot_four_state_correction(series, title):
    """Separate inferred extra SW in flux units and multiplier departure."""
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(12,7),sharex=True,layout='constrained')
    years=series['year']
    axes[0].plot(years,series['filtered_SW_correction'],label='Applied prediction correction (filtered)')
    axes[0].plot(years,series['retrospective_SW_correction'],label='Retrospective reconstruction')
    axes[0].set_ylabel('Additional absorbed SW (W m$^{-2}$)')
    axes[0].set_title('Prediction at year k uses the inferred state from year k−1; positive adds heating')
    axes[1].plot(years,series['filtered'][:,3]-1,label='Filtered multiplier − 1')
    axes[1].plot(years,series['retrospective'][:,3]-1,label='Retrospective multiplier − 1')
    axes[1].set(ylabel='Multiplier departure (dimensionless)',xlabel='Year (end-of-year inferred state)')
    for ax in axes:ax.axhline(0,color='0.4',ls=':');ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.suptitle(title+' — inferred ta4 correction only')
    return fig


def plot_four_state_components(model, years, observations, title, reference_temperature=None):
    """Active ta4 wrapper inputs, posterior states, fluxes and correction."""
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    d=four_state_diagnostic_series(model,years,observations)
    years=d['year'];n=len(years);obs=d['observations']
    fig,axes=plt.subplots(4,2,figsize=(13,13),sharex=True,layout='constrained');a=axes.ravel()
    for j,scale,label in [(0,1,'Temperature anomaly (°C)'),(1,model.zJ_from_W,'OHC (ZJ)')]:
        offset=model.offset if j==0 else 0
        a[j].plot(years,obs[:,j]*scale-offset,color='0.7',label='Adjusted input')
        for key,caption in [('filtered','Filtered'),('retrospective','Retrospective'),('blind','Blind replay, multiplier = 1')]:
            a[j].plot(years,d[key][:,j]*scale-offset,label=caption)
        if j==0:
            reference=(obs[:,0]-offset if reference_temperature is None
                       else np.asarray(reference_temperature))
            if reference.shape!=(n,):raise ValueError('Reference temperature has wrong shape')
            d['reference_temperature_anomaly']=reference.copy()
            if reference_temperature is not None:
                a[j].plot(years,reference,color='0.85',lw=.8,label='Original temperature input')
            benchmark=np.full(n,np.nan)
            for i in range(10,n-10):benchmark[i]=.5*(reference[i-10:i+10].mean()+reference[i-9:i+11].mean())
            a[j].plot(years,benchmark,'k',lw=1.5,label='Centered 20-year original-input benchmark')
        a[j].set_ylabel(label);a[j].legend(fontsize=7)
    for j,value,label in [(2,model.opt_depth[:n],'Active AOD (after wrapper adjustment)'),
                          (3,model.tsi[:n],'Incident solar (W m$^{-2}$)'),
                          (4,10**model.lCo2[:n],'Effective CO₂ (ppm)'),
                          (5,model.anthro_clouds[:n],'Centred ARI + ACI (W m$^{-2}$)')]:
        a[j].plot(years,value);a[j].set_ylabel(label)
    a[6].plot(years,d['filtered_absorbed_SW'],label='Absorbed SW, with inferred correction')
    a[6].plot(years,d['filtered_base_SW'],'--',label='Same state, multiplier fixed to 1')
    a[6].plot(years,d['filtered_outgoing_LW'],label='Outgoing LW');a[6].set_ylabel('Prediction flux (W m$^{-2}$)');a[6].legend(fontsize=7)
    for key,label in [('filtered_net_flux','Net radiation'),('ocean_storage_rate','OHC storage rate'),('non_ocean_storage_rate','Non-ocean storage rate')]:a[7].plot(years,d[key],label=label)
    toa_sigma=d['TOA_effective_sigma'];toa_nominal=d['TOA_nominal_sigma']
    a[7].fill_between(years,obs[:,2]-1.96*toa_sigma,obs[:,2]+1.96*toa_sigma,
                      color='0.7',alpha=.22,label='TOA input: effective 95% interval',zorder=0)
    a[7].fill_between(years,obs[:,2]-1.96*toa_nominal,obs[:,2]+1.96*toa_nominal,
                      color='0.5',alpha=.18,label='TOA input: nominal 95% interval',zorder=0)
    a[7].plot(years,obs[:,2],color='maroon',lw=1,label='Adjusted TOA input')
    a[7].set_ylabel('Energy balance / TOA (W m$^{-2}$)');a[7].legend(fontsize=6.5)
    for ax in a:ax.grid(alpha=.2)
    for ax in axes[-1]:ax.set_xlabel('Year')
    fig.suptitle(title+'\nAerosol p'+str(model.AEROSOL_FORCING_PERCENTILE)+', observational eCO₂ p'+str(model.eCO2_FORCING_PERCENTILE)
                 +', fdbkA='+str(model.fdbkA)+', volcanic target='+str(model.volcanic_target_efficiency)+' W/m²/AOD')
    fig2,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for j,label in enumerate(['Temperature','OHC','TOA']):
        axes[j].hist(model.qqy[1:n,j],bins=24,density=True,color='0.7',label='Innovations')
        if j<2:mixture=model.qqyker[j]
        else:
            width=np.sqrt(9*np.asarray(model.TOA_meas_artif_var[1:n]))/np.sqrt(model.S[1:n,j,j])
            mixture=norm.pdf(model.xnorm[None,:],loc=model.qqy[1:n,j,None],scale=width[:,None]).mean(axis=0)
        axes[j].plot(model.xnorm,mixture,color='royalblue',label='Observation mixture')
        axes[j].plot(model.xnorm,norm.pdf(model.xnorm),'k--',label='Standard normal')
        axes[j].set(title=label,xlabel='Standardised innovation',ylabel='Density');axes[j].legend(fontsize=7)
    fig2.suptitle(title+' — innovation diagnostics')
    return [fig,fig2,plot_four_state_correction(d,title)],d
