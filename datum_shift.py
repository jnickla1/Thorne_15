"""1850–1900 datum correction for methods whose output is anchored to a fixed reference (2026-10-21).

The preindustrial datum moved from 1850–1899 to 1850–1900 (inclusive), which lowers every input temperature series
by c = mean(1850–1900) - mean(1850–1899). Methods that follow their input are unaffected, but methods whose output is
tied to a fixed reference (hand-set bias constants: FaIR, GWI, eROF, KCC_human) would gain a bias of +c. Each of those
methods calls `correction` when run_method receives its temperature series and passes its result through `shift`,
which lowers the whole output distribution by c (equivalent to adding c to the method's bias constant).

c values:
  historical / heads-and-tails (HadCRUT5.1):            0.00254 °C
  MPI-ESM1-2-LR futures (mean over the 50 members):     0.00113 °C
  NorESM futures (aligned to HadCRUT5 observations):    0.00244 °C
"""
import numpy as np

HISTORICAL = 0.00254
MPI_FUTURE = 0.00113
NORESM_FUTURE = 0.00244


def correction(experiment_type, splits_futures=True):
    """Datum correction (°C) for this experiment; futures use the model values only if the method branches for futures."""
    parts = experiment_type.split('_')
    if splits_futures and parts[0] == 'fut':
        return MPI_FUTURE if parts[1] == 'ESM1-2-LR' else NORESM_FUTURE
    return HISTORICAL


def shift(result, dc):
    """Lower a method's output by dc: central estimates for tuple results; mean/pvalue/log_likelihood/resample for dicts."""
    if dc == 0:
        return result
    if isinstance(result, dict):
        out = dict(result)
        f = dict(result)
        if 'mean' in f:
            out['mean'] = lambda year_idx, k, _f=f['mean']: _f(year_idx, k) - dc
        if 'pvalue' in f:
            out['pvalue'] = lambda year_idx, point, k, *a, _f=f['pvalue'], **kw: _f(year_idx, np.asarray(point) + dc, k, *a, **kw)
        if 'log_likelihood' in f:
            out['log_likelihood'] = lambda year_idx, point, k, *a, _f=f['log_likelihood'], **kw: _f(year_idx, np.asarray(point) + dc, k, *a, **kw)
        if 'resample' in f:
            out['resample'] = lambda year_idx, n, k, *a, _f=f['resample'], **kw: _f(year_idx, n, k, *a, **kw) - dc
        return out
    if isinstance(result, (tuple, list)) and len(result) == 4:
        cur, cur_se, ret, ret_se = result
        return (np.asarray(cur) - dc, cur_se, np.asarray(ret) - dc, ret_se)
    return result
