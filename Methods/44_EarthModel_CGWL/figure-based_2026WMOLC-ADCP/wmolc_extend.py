"""
Extend GlobalT_WMOLC-ADCPforecast_1991-2020.csv with pseudo-ensembles digitized
from the WMO LC-ADCP 'decadal_exchange_timeseries_GlobalT_yrs1-10_<init>.png' figures.

What the figures show (verified against the CSV rows for start years 2023 & 2024):
  per-centre panel: square marker = centre ensemble mean, whisker caps = member min / max
  (RMS error ~0.4 px ~ 0.007 C).  The marker is a fixed-size symbol, not an IQR box.
  top 'Mean' panel marker ~ mean over centres WITHOUT an asterisk.

Method for each centre in a new figure:
  take that centre's members from the most recent CSV row (its 'shape template'),
  apply a monotone piecewise-linear map  t_min->lo, knot->m*, t_max->hi,
  with m* solved so the mapped members have exactly the plotted mean.
  -> exact plotted (min, mean, max), same member count and within-centre shape as last year.
  If the upper cap is clipped by the axis, hi is set from the template's asymmetry.
"""
import numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage, optimize, stats
import os
UP =os.getcwd()
INCLUDE_IITM = True   # False -> drop the clipped IITM-MoES* panel from 2026-2035
CSV = f'{UP}/../GlobalT_WMOLC-ADCPforecast_1991-2020.csv'
IMG = lambda y: f'{UP}/decadal_exchange_timeseries_GlobalT_yrs1-10_{y}.png'
# panel titles in reading order (after the top 'Mean' panel); init year -> forecast start = init+1
NAMES = {'2022': ['BCCR','BSC','CCCMA','CMCC','DWD/MPI','GFDL','LabESM','MIROC','MOHC'],
         '2023': ['BCCR','BSC','CCCMA','CMCC','DWD/MPI','GFDL','IPSL*','LabESM','MIROC','MOHC','NCAR*'],
         '2024': ['BSC','CCCMA','CMCC','DWD/MPI','GFDL','IPSL*','LabESM','MIROC','MOHC/KMA'],
         '2025': ['BCCR','BSC','CCCMA','CMCC','DWD/MPI','GFDL','IITM-MoES*','MIROC','MOHC/KMA']}
TEMPLATE_NAME = {'MOHC/KMA': 'MOHC'}          # which past centre's members to use as a template
# CSV column blocks (identified by DP segmentation against the digitized figures)
BLOCKS = {2023: [('CCCMA',0,10),('CMCC',10,30),('DWD/MPI',30,46),('MOHC',46,56),('BSC',56,66),
                 ('LabESM',66,75),('GFDL',75,85),('MIROC',85,105),('BCCR',105,125)],
          2024: [('CCCMA',0,10),('CMCC',10,30),('DWD/MPI',30,46),('MOHC',46,56),('BSC',56,66),
                 ('LabESM',66,75),('GFDL',75,85),('IPSL*',85,95),('MIROC',95,106),('NCAR*',106,126),
                 ('BCCR',126,149)]}   # 75:149 split is ambiguous at pixel precision (similar centres)

# ---------------------------------------------------------------- digitizing
def _panels(im):
    blk = im.sum(-1) < 60
    grey = np.all(np.abs(im - 128) < 6, -1)
    sp = []
    for r in range(im.shape[0]):
        xs = np.where(blk[r])[0]
        if len(xs) < 100: continue
        for ru in np.split(xs, np.where(np.diff(xs) > 1)[0] + 1):
            if len(ru) > 100: sp.append((r, ru[0], ru[-1]))
    sp.sort(key=lambda t: (t[1] // 40, t[0])); used = set(); out = []
    for i, (r, a, b) in enumerate(sp):
        if i in used: continue
        for j in range(i + 1, len(sp)):
            r2, a2, b2 = sp[j]
            if j not in used and abs(a2 - a) <= 6 and r2 > r:
                used |= {i, j}
                z = [k for k in range(r + 1, r2) if grey[k, a + 5:b - 5].sum() > 30][0]
                out.append(dict(top=r, bot=r2, zero=z, x1=max(b, b2), x0=min(a, a2))); break
    return sorted(out, key=lambda q: (q['top'], q['x0']))

def _measure(im, p):
    x0 = p['x1'] - 25; off = p['top'] - 15
    sub = im[off:p['bot'] + 15, x0:p['x1'] - 1]
    r, g, b = sub[..., 0], sub[..., 1], sub[..., 2]
    fill = ((r < 10) & (g < 10) & (b > 240)) | ((r > 240) & (g < 10) & (b < 10))
    lab, n = ndimage.label(fill); objs = ndimage.find_objects(lab)
    k = max((o[1].start, i) for i, o in enumerate(objs) if (o[1].stop - o[1].start) >= 4)[1]
    fill = lab == k + 1
    cols = np.where(fill.any(0))[0]; xc = (cols[0] + cols[-1]) // 2
    nonwhite = sub.min(-1) < 245; c = nonwhite[:, xc]; fr = np.where(fill[:, xc])[0]
    t = fr[0]
    while t > 0 and c[t - 1]: t -= 1
    bb = fr[-1]
    while c[bb + 1]: bb += 1
    d = (765 - sub[:, xc].sum(-1)).astype(float)
    pk = lambda lo, hi: (np.arange(lo, hi + 1) * d[lo:hi + 1]).sum() / d[lo:hi + 1].sum()
    s = p['zero'] - p['top']; val = lambda row: (p['zero'] - (row + off)) / s
    hi = val(pk(t, t + 2)); lo = val(pk(bb - 2, bb)); mk = val((fr[0] + fr[-1]) / 2)
    return dict(mk=mk, lo=lo, hi=hi, px=1 / s, clip_hi=hi > 1 - 1.5 / s)

def digitize(init):
    im = np.array(Image.open(IMG(init)).convert('RGB')).astype(int)
    P = _panels(im); names = ['MEAN'] + NAMES[init]
    assert len(P) == len(names), (init, len(P))
    return {n: _measure(im, p) for n, p in zip(names, P)}

# ---------------------------------------------------------------- pseudo-ensemble
def _resample(t, n):
    if len(t) == n: return np.sort(t)
    q = (np.arange(1, n + 1) - 0.5) / n
    return np.quantile(t, q)

def pseudo_members(template, n, lo, mk, hi, clip_hi=False):
    t = _resample(np.asarray(template, float), n)
    tmin, tmax, tmean = t.min(), t.max(), t.mean()
    if clip_hi:   # upper cap hidden by axis: borrow the template's upper/lower asymmetry
        hi = max(hi, mk + (mk - lo) * (tmax - tmean) / (tmean - tmin))
    def mapped(m):
        return np.where(t <= tmean, lo + (t - tmin) / (tmean - tmin) * (m - lo),
                        m + (t - tmean) / (tmax - tmean) * (hi - m))
    f = lambda m: mapped(m).mean() - mk
    if f(lo + 1e-6) * f(hi - 1e-6) < 0:
        return mapped(optimize.brentq(f, lo + 1e-6, hi - 1e-6))
    # fallback (plotted mean too skewed for the template): keep min/max members,
    # bend the interior with u**gamma until the mean matches
    u = (t - tmin) / (tmax - tmin)
    g = lambda lg: (lo + (hi - lo) * u ** np.exp(lg)).mean() - mk
    lg = optimize.brentq(g, -6, 6)
    return lo + (hi - lo) * u ** np.exp(lg)

def csv_blocks(sy, df):
    x = df[df.start_year == sy].iloc[0, 2:].values.astype(float); x = x[~np.isnan(x)]
    return {c: x[a:b] for c, a, b in BLOCKS[sy]}

def build(init, templates, counts=None, skip=(), override=None):
    dig = digitize(init); out = {}; info = []
    for c, o in (override or {}).items(): dig[c] = {**dig[c], **o}
    for c, v in dig.items():
        if c == 'MEAN' or c in skip: continue
        tn = TEMPLATE_NAME.get(c, c)
        tmpl = templates.get(tn)
        if tmpl is None:
            n = (counts or {}).get(c, 10)
            tmpl = stats.norm.ppf((np.arange(1, n + 1) - 0.375) / (n + 0.25))  # Blom scores
        n = (counts or {}).get(c, len(tmpl))
        out[c] = pseudo_members(tmpl, n, v['lo'], v['mk'], v['hi'], v['clip_hi'])
        info.append((c, n, v['mk'], v['lo'], v['hi'], v['clip_hi']))
    return out, dig, info

if __name__ == '__main__':
    df = pd.read_csv(CSV, na_values='NA')
    t23, t24 = csv_blocks(2023, df), csv_blocks(2024, df)

    # ---------- validation: rebuild the known 2024-2033 row from its figure + 2023 templates
    truth = np.concatenate(list(t24.values()))
    counts = {c: len(v) for c, v in t24.items()}          # oracle counts
    rec, dig, _ = build('2023', t23, counts)
    recon = np.concatenate(list(rec.values()))
    counts_prev = {c: len(v) for c, v in t23.items()}      # realistic: last year's counts, 10 for new
    rec2, _, _ = build('2023', t23, counts_prev)
    recon2 = np.concatenate(list(rec2.values()))
    qs = [5, 17, 50, 83, 95]
    print('VALIDATION  2024-2033  (n, mean, sd, q5,q17,q50,q83,q95, W1 vs truth)')
    for lab, z in [('truth', truth), ('recon, true counts', recon), ('recon, prev-yr counts', recon2)]:
        print(f'  {lab:22s} n={len(z):3d} mean={z.mean():.3f} sd={z.std(ddof=1):.3f} q=',
              np.round(np.percentile(z, qs), 3), f'W1={stats.wasserstein_distance(z, truth):.4f}')

    # ---------- new rows: templates = latest row where the centre appears
    templ = {**t23, **t24}
    newrows = {}
    for init, sy in [('2024', 2025), ('2025', 2026)]:
        ov = {}
        if init == '2025' and INCLUDE_IITM:
            # IITM panel is clipped at the axis top: only its lower cap (~0.97) is visible.
            # Top 'Mean' panel = equal-weight mean of plotted centres (checked on 2023/2024 figs),
            # so invert it for IITM's mean; assume symmetric spread and n=10.
            d0 = digitize(init)
            others = [v['mk'] for c, v in d0.items() if c not in ('MEAN', 'IITM-MoES*')]
            mk = d0['MEAN']['mk'] * (len(others) + 1) - sum(others)
            lo = d0['IITM-MoES*']['lo']
            ov = {'IITM-MoES*': dict(mk=mk, lo=lo, hi=2 * mk - lo, clip_hi=False)}
        mem, dig, info = build(init, templ, counts={'IITM-MoES*': 10}, override=ov,
                               skip=() if INCLUDE_IITM else ('IITM-MoES*',))
        z = np.concatenate(list(mem.values()))
        newrows[sy] = z
        noast = np.concatenate([v for c, v in mem.items() if '*' not in c])
        eqw = np.mean([v.mean() for c, v in mem.items() if '*' not in c])
        print(f'\n{sy}-{sy+9}: n={len(z)} mean={z.mean():.3f} sd={z.std(ddof=1):.3f} q=',
              np.round(np.percentile(z, qs), 3))
        allm = np.mean([v.mean() for v in mem.values()])
        print(f'  check vs top-panel marker {dig["MEAN"]["mk"]:.3f}: member-wt no-asterisk {noast.mean():.3f}, '
              f'equal-wt no-asterisk {eqw:.3f}, equal-wt all plotted {allm:.3f}')
        for c, n, mk, lo, hi, cl in info:
            print(f'   {c:11s} n={n:2d} mean={mk:.3f} min={lo:.3f} max={hi:.3f}{"  (max clipped -> extrapolated %.3f)" % mem[c].max() if cl else ""}')

    # ---------- write extended CSV
    ncol = max(df.shape[1] - 2, max(len(v) for v in newrows.values()))
    rows = []
    for sy, z in newrows.items():
        rows.append([sy, sy + 9] + list(np.round(z, 7)) + [np.nan] * (ncol - len(z)))
    header = open(CSV).readline().rstrip('\r\n').split(',')
    header += [''] * (ncol + 2 - len(header))
    ext = pd.concat([df.reindex(columns=df.columns), pd.DataFrame(rows, columns=df.columns)], ignore_index=True)
    ext.columns = header
    ext.to_csv('/mnt/user-data/outputs/GlobalT_WMOLC-ADCPforecast_1991-2020_ext2026.csv',
               index=False, na_rep='NA', float_format='%.8g', lineterminator='\r\n')

    # ---------- figure
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for z, lab in [(truth, 'CSV truth'), (recon, 'reconstructed (true counts)'), (recon2, 'reconstructed (prev-yr counts)')]:
        ax[0].plot(np.sort(z), np.linspace(0, 1, len(z)), label=lab)
    ax[0].set(title='Validation: 2024-2033 pooled members', xlabel='anomaly wrt 1991-2020 (°C)', ylabel='ECDF'); ax[0].legend(fontsize=8)
    for c, v in t24.items():
        r = rec[c]; ax[1].plot([v.mean()], [r.mean()], 'o'); ax[1].plot([v.min(), v.max()], [r.min(), r.max()], 'k+', ms=4)
        ax[1].annotate(c, (v.mean(), r.mean()), fontsize=7)
    ax[1].plot([0.1, 1.05], [0.1, 1.05], 'k:', lw=0.8)
    ax[1].set(title='Per-centre: mean (o), min/max (+)', xlabel='CSV', ylabel='digitized/reconstructed')
    fig.tight_layout(); fig.savefig('/mnt/user-data/outputs/wmolc_validation.png', dpi=130)
