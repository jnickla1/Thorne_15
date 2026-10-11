"""Rebuild only corrected retrospective caches, guarding every current series.

Reads original inputs/results. Writes new caches and regression tables outside
Methods and Results. Before adoption, --stage-root can validate staged edits.
"""
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import pickle
import sys
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore',category=RuntimeWarning)
REAL=['CGWL10y_sfUKCP','EBMKF_ta4','FaIR_comb_unB','GAM_AR1','GWI_tot_CGWL','Kal_flexLin','lowess1dt36wnc']
ANTHRO=['FaIR_nonat_unB','GWI_anthro_CGWL']
FAIR_DIR=Path('Methods/43_Forcing_Based/1_ERF_FaIR')
GWI_DIR=Path('Methods/43_Forcing_Based/3_Human_Induced')
FAIR=['FaIR_comb_unB','FaIR_nonat_unB','FaIR_all_unB','FaIR_anthro_unB']
GWI=['GWI_tot_CGWL','GWI_anthro_CGWL']

def module_from_source(path,source_file):
    rel=path.parts[path.parts.index('Methods'):]
    name='.'.join(rel[:-1]+(Path(rel[-1]).stem,))
    ns={'__file__':str(path),'__name__':name,'__package__':name.rsplit('.',1)[0]}
    exec(compile(source_file.read_text(),str(path),'exec'),ns)
    return ns['run_method']

def closure_array(result,key):
    fun=result['mean']
    return dict(zip(fun.__code__.co_freevars,[x.cell_contents for x in fun.__closure__]))[key]

def combine(series,weights):
    series=np.asarray(series);w=np.asarray(weights)[:,None]
    valid=np.isfinite(series)
    denominator=np.sum(valid*w,axis=0)
    with np.errstate(invalid='ignore',divide='ignore'):
        normw=np.where(valid,w,0)/denominator
        center=np.nansum(series*normw,axis=0)
    center[denominator==0]=np.nan
    sums=np.nansum(normw,axis=0)
    np.testing.assert_allclose(sums[denominator>0],1,atol=2e-15)
    return center,sums

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,default=Path('/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats'))
    p.add_argument('--backup-root',type=Path,default=Path(__file__).resolve().parents[1]/'review/backups')
    p.add_argument('--stage-root',type=Path)
    p.add_argument('--outdir',type=Path,default=Path(__file__).resolve().parents[1]/'data/corrected_results')
    p.add_argument('--review-dir',type=Path,default=Path(__file__).resolve().parents[1]/'review')
    args=p.parse_args();root=args.source_root.resolve();out=args.outdir;review=args.review_dir
    out.mkdir(parents=True,exist_ok=True);review.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(root))
    old={};new={}
    for directory,names in [(FAIR_DIR,FAIR),(GWI_DIR,GWI)]:
        for name in names:
            rel=directory/(name+'_method.py')
            backup=args.backup_root/rel
            if not backup.exists(): backup=Path(str(backup)+'.original')
            old[name]=module_from_source(root/rel,backup)
            new[name]=module_from_source(root/rel,(args.stage_root or root)/rel)
    had=pd.read_csv(root/'Common_Data/HadCRUT.5.1.csv')
    var=pd.read_csv(root/'Results2/historical_names_var_scale.csv').set_index('method_name')
    changes=[];weights_rows=[];fallback_rows=[];source_hashes={}
    for baseline in ['rpi','satcal']:
        rawfile='sst_pseudo.csv' if baseline=='rpi' else 'sst_pseudo_1981-2010.csv'
        raw=np.loadtxt(root/'Hist_heads_tails'/rawfile,delimiter=',');years=raw[:,0].astype(int)
        for member in range(100):
            t=raw[:,member+1]
            uncert=(t+had.Lower.to_numpy()-had.Anomaly.to_numpy(),t+had.Upper.to_numpy()-had.Anomaly.to_numpy())
            corrected={}
            for method in FAIR:
                exp=f'histens_{baseline}_real'
                a=old[method](years,t,uncert,member,exp);b=new[method](years,t,uncert,member,exp)
                # Complete samples include NaN padding and added noise; equality
                # guarantees the unchanged current distribution, not just means.
                np.testing.assert_array_equal(closure_array(a,'samp_cur'),closure_array(b,'samp_cur'))
                oldmeans=[a['mean'](years,k) for k in [0,1]]
                newmeans=[b['mean'](years,k) for k in [0,1]]
                oldses=[a['se'](years,k) for k in [0,1]]
                newses=[b['se'](years,k) for k in [0,1]]
                np.testing.assert_array_equal(oldmeans[0],newmeans[0]);np.testing.assert_array_equal(oldses[0],newses[0])
                corrected[method]=(newmeans[0],newses[0],newmeans[1],newses[1])
                for threshold in [0.5,1.,1.5]:
                    points=np.full(len(years),threshold if baseline=='rpi' else threshold-.7)
                    for two_sided in [True,False]:
                        np.testing.assert_array_equal(a['pvalue'](years,points,0,two_sided=two_sided),
                                                      b['pvalue'](years,points,0,two_sided=two_sided))
                # Current likelihood and resampling consume the same sample
                # arrays and unchanged function bodies; representative checks.
                yi=np.array([2025]);pt=np.array([1.5 if baseline=='rpi' else .8])
                np.testing.assert_array_equal(a['log_likelihood'](yi,pt,0),b['log_likelihood'](yi,pt,0))
                np.random.seed(7301);asamp=a['resample'](yi,12,0)
                np.random.seed(7301);bsamp=b['resample'](yi,12,0)
                np.testing.assert_array_equal(asamp,bsamp)
                changes.append(dict(baseline=baseline,member=member,method=method,
                    current_2025_before=oldmeans[0][-1],current_2025_after=newmeans[0][-1],
                    retrospective_2025_before=oldmeans[1][-1],retrospective_2025_after=newmeans[1][-1],
                    current_mean_max_change=float(np.nanmax(np.abs(newmeans[0]-oldmeans[0]))),
                    current_se_max_change=float(np.nanmax(np.abs(newses[0]-oldses[0]))),
                    retrospective_rate_before=float(np.polyfit(np.arange(15),oldmeans[1][-15:],2)[0]*28*10+np.polyfit(np.arange(15),oldmeans[1][-15:],2)[1]*10),
                    retrospective_rate_after=float(np.polyfit(np.arange(15),newmeans[1][-15:],2)[0]*28*10+np.polyfit(np.arange(15),newmeans[1][-15:],2)[1]*10)))
            for goal,selected in [('real',REAL),('anthro',ANTHRO)]:
                exp=f'histens_{baseline}_{goal}'
                source=root/'Results'/f'results{exp}{member}.pickle'
                blob=source.read_bytes();source_hashes[str(source)]=hashlib.sha256(blob).hexdigest()
                original=pickle.loads(blob);data=pickle.loads(blob)
                for method in FAIR:
                    if method in data:
                        np.testing.assert_array_equal(original[method]['LT_trend'][0],corrected[method][0])
                        np.testing.assert_array_equal(original[method]['LT_trend'][1],corrected[method][1])
                        data[method]['LT_trend']=corrected[method]
                        data[method]['retrospective_source']='corrected_FaIR_recipe_and_boundary'
                for method in selected:
                    c,se,r,rse=data[method]['LT_trend']
                    if not np.isfinite(r).any():
                        data[method]['LT_trend']=(c,se,c.copy(),se.copy())
                        data[method]['retrospective_source']='current_proxy_no_retrospective_in_reader'
                        fallback_rows.append(dict(experiment=exp,member=member,method=method,current_used_as_retrospective=True))
                for method in original:
                    before=original[method]['LT_trend'];after=data[method]['LT_trend']
                    np.testing.assert_array_equal(before[0],after[0]);np.testing.assert_array_equal(before[1],after[1])
                # Neither current method curves nor combination curves change.
                iv=1/var.loc[selected,'err_var100'].to_numpy()
                for label,w in [('PIVW',iv),('equal',np.ones(len(selected)))]:
                    before,_=combine([original[m]['LT_trend'][0] for m in selected],w)
                    after,current_sums=combine([data[m]['LT_trend'][0] for m in selected],w)
                    _,retro_sums=combine([data[m]['LT_trend'][2] for m in selected],w)
                    np.testing.assert_array_equal(before,after)
                    for year,cs,rs in zip(years,current_sums,retro_sums):
                        weights_rows.append(dict(experiment=exp,member=member,combination=label,year=year,
                                                 current_weight_sum=cs,retrospective_weight_sum=rs))
                with (out/source.name).open('wb') as fp:pickle.dump(data,fp,protocol=pickle.HIGHEST_PROTOCOL)
        # Directly check both GWI APIs, including their temporary proxy branches.
        t=raw[:,1];uncert=(t+had.Lower.to_numpy()-had.Anomaly.to_numpy(),t+had.Upper.to_numpy()-had.Anomaly.to_numpy())
        for method in GWI:
            a=old[method](years,t,uncert,0,f'histens_{baseline}_real')
            b=new[method](years,t,uncert,0,f'histens_{baseline}_real')
            for key in ['mean','se']:
                np.testing.assert_array_equal(a[key](years,0),b[key](years,0))
                np.testing.assert_array_equal(b[key](years,0),b[key](years,1))
            points=np.full(len(years),1.5 if baseline=='rpi' else .8)
            for key in ['pvalue','log_likelihood']:
                np.testing.assert_array_equal(a[key](years,points,0),b[key](years,points,0))
                np.testing.assert_array_equal(b[key](years,points,0),b[key](years,points,1))
            yi=np.array([2025])
            np.random.seed(7301);asamp=a['resample'](yi,12,0)
            np.random.seed(7301);bsamp=b['resample'](yi,12,0)
            np.random.seed(7301);rsamp=b['resample'](yi,12,1)
            np.testing.assert_array_equal(asamp,bsamp);np.testing.assert_array_equal(bsamp,rsamp)
        print('Passed current-output and fallback checks:',baseline,flush=True)
    # Verify original cached results remained byte-for-byte unchanged.
    for name,digest in source_hashes.items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
    frame=pd.DataFrame(changes);frame.to_csv(review/'corrections_current_output_regression.csv',index=False)
    summary=frame.drop(columns='member').groupby(['baseline','method']).mean(numeric_only=True).reset_index()
    summary.to_csv(review/'corrections_endpoint_summary.csv',index=False)
    pd.DataFrame(fallback_rows).to_csv(review/'retrospective_fallback_members.csv',index=False)
    pd.DataFrame(weights_rows).to_csv(review/'combination_weight_sums.csv',index=False)
    report=dict(current_output_changed=False,max_current_mean_change=float(frame.current_mean_max_change.max()),
                max_current_se_change=float(frame.current_se_max_change.max()),fair_member_calibration_method_cases=len(frame),
                corrected_result_files=len(source_hashes),source_results_unchanged=True,
                fallbacks=sorted(set(r['method'] for r in fallback_rows)),normalized_weights=True,
                corrected_results_dir=str(out.resolve()),original_result_sha256=source_hashes)
    (review/'corrections_regression_report.json').write_text(json.dumps(report,indent=2))
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':main()
