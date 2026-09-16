#!/usr/bin/env python3
# P2-02 Stage B1: cis-MVMR — ITPKA vs RTF1 blood eQTL (eQTLGen) as joint exposures,
# outcomes = FinnGen R12 dermatitis/eczema + atopic dermatitis.
# Instruments: union of cis SNPs with p<5e-8 for either gene, LD-clumped r2<0.1 (1000G EUR).
# Method: MVMR-IVW = precision-weighted least squares of beta_GWAS on [beta_ITPKA, beta_RTF1].
import pandas as pd, numpy as np
from scipy import stats

DIR = '/workspace/coloc_panel'

# ---- load eQTLGen (hg19) -> universe (hg38) by rsid ----
univ = pd.read_csv(f'{DIR}/snp_universe.csv')
univ['key'] = univ.chromosome.astype(str) + ':' + univ.position.astype(str)

def load_eqtl(gene):
    df = pd.read_csv(f'{DIR}/eqtlgen_{gene}.csv').sort_values('p').drop_duplicates('rsid')
    df = df.merge(univ[['rsid','key','position']].rename(columns={'position':'pos38'}), on='rsid', how='inner')
    return df[['rsid','key','pos38','ea','nea','eaf','beta','se','p','n']].rename(
        columns={'beta':f'beta_{gene}','se':f'se_{gene}','p':f'p_{gene}','eaf':f'eaf_{gene}'})

ei = load_eqtl('ITPKA'); er = load_eqtl('RTF1')
m = ei.merge(er, on='rsid', suffixes=('_i','_r'))
# verify same effect allele coding across the two eQTLGen files
same = (m.ea_i == m.ea_r) & (m.nea_i == m.nea_r)
print(f'eQTLGen allele coding identical for {same.sum()}/{len(m)} shared SNPs; dropping mismatches')
m = m[same].copy()
m['key'] = m.key_i

# ---- load FinnGen GWAS, harmonize to eQTL ea ----
def load_gwas(ep):
    gw = pd.read_csv(f'{DIR}/gwas_finngen_{ep}.tsv', sep='\t', header=None,
        names=['chromosome','position','ref','alt','rsids','nearest_genes','pval','mlogp',
               'beta','sebeta','af_alt','af_alt_cases','af_alt_controls'])
    gw['key'] = gw.chromosome.astype(str) + ':' + gw.position.astype(str)
    return gw

def harmonize(mm, gw):
    g = mm.merge(gw[['key','ref','alt','beta','sebeta','pval']], on='key')
    exact = (g.ea_i == g.alt) & (g.nea_i == g.ref)
    swap  = (g.ea_i == g.ref) & (g.nea_i == g.alt)
    g = g[exact | swap].copy()
    sw = (g.ea_i == g.ref) & (g.nea_i == g.alt)  # recompute on filtered frame
    g.loc[sw, 'beta'] = -g.loc[sw, 'beta']
    return g

# ---- LD panel for clumping ----
ldsnps = pd.read_csv(f'{DIR}/ld_eur_snps.csv')
ld = np.load(f'{DIR}/ld_eur.npy')
ldidx = {r: i for i, r in enumerate(ldsnps.rsid)}

def clump(df, r2thr=0.1):
    df = df.sort_values('pmin')
    keep, blocked = [], np.zeros(len(df), bool)
    arr = df.reset_index(drop=True)
    for j in range(len(arr)):
        if blocked[j]:
            continue
        keep.append(j)
        rj = arr.rsid[j]
        if rj not in ldidx:
            continue
        for k2 in range(j+1, len(arr)):
            if blocked[k2]:
                continue
            rk = arr.rsid[k2]
            if rk in ldidx and abs(ld[ldidx[rj], ldidx[rk]])**2 > r2thr:
                blocked[k2] = True
    return arr.iloc[keep]

# ---- MVMR-IVW ----
def mvmr_ivw(d):
    X = d[['beta_ITPKA','beta_RTF1']].values
    y = d['beta'].values
    w = 1.0 / d.sebeta.values**2
    WX = X * np.sqrt(w)[:, None]; Wy = y * np.sqrt(w)
    beta, res, rank, sv = np.linalg.lstsq(WX, Wy, rcond=None)
    XtWX = WX.T @ WX
    cov_fe = np.linalg.inv(XtWX)
    resid = Wy - WX @ beta
    q = float(resid @ resid); dfree = len(d) - 2
    scale = max(q / dfree, 1.0)  # random-effects scale (>=1)
    cov_re = cov_fe * scale
    out = {}
    for i, name in enumerate(['ITPKA','RTF1']):
        se_fe = np.sqrt(cov_fe[i,i]); se_re = np.sqrt(cov_re[i,i])
        out[name] = dict(b=beta[i], se_fe=se_fe, p_fe=2*stats.norm.sf(abs(beta[i]/se_fe)),
                         se_re=se_re, p_re=2*stats.norm.sf(abs(beta[i]/se_re)))
    out['Q_resid'] = q; out['df'] = dfree; out['Q_p'] = 1-stats.chi2.cdf(q, dfree)
    return out

def ivw_uni(d, bcol, secoll_se_g='sebeta'):
    x = d[bcol].values; y = d['beta'].values; w = 1.0/d[secoll_se_g].values**2
    b = (w*x*y).sum() / (w*x*x).sum()
    se = np.sqrt(1.0/(w*x*x).sum())
    return b, se, 2*stats.norm.sf(abs(b/se))

def mvmr_gls_ld(d):
    """MVMR-IVW with LD-aware covariance: Sigma_jk = se_j*se_k*LD_jk (1000G EUR)."""
    d = d[d.rsid.isin(ldidx)].copy()
    X = d[['beta_ITPKA','beta_RTF1']].values
    y = d['beta'].values
    ii = [ldidx[r] for r in d.rsid]
    R = ld[np.ix_(ii, ii)]
    S = np.outer(d.sebeta.values, d.sebeta.values) * R
    # regularize slightly for numerical stability
    S += np.eye(len(d)) * 1e-10
    Si = np.linalg.inv(S)
    XtSiX = X.T @ Si @ X
    beta = np.linalg.solve(XtSiX, X.T @ Si @ y)
    cov = np.linalg.inv(XtSiX)
    resid = y - X @ beta
    q = float(resid @ Si @ resid); dfree = len(d) - 2
    scale = max(q/dfree, 1.0)
    cov = cov * scale
    out = {}
    for i, name in enumerate(['ITPKA','RTF1']):
        se = np.sqrt(cov[i,i])
        out[name] = dict(b=beta[i], se=se, p=2*stats.norm.sf(abs(beta[i]/se)))
    out['n'] = len(d); out['Q_p'] = 1-stats.chi2.cdf(q, dfree)
    return out

results = []
for ep, label in [('L12_DERMATITISECZEMA','Eczema'), ('L12_ATOPIC','AD')]:
    g = harmonize(m, load_gwas(ep))
    g['pmin'] = g[['p_ITPKA','p_RTF1']].min(axis=1)
    inst = g[(g.p_ITPKA < 5e-8) | (g.p_RTF1 < 5e-8)].copy()
    n_pre = len(inst)
    inst = clump(inst, 0.1)
    # instrument strength diagnostics
    F_i = (inst.beta_ITPKA/inst.se_ITPKA)**2
    F_r = (inst.beta_RTF1/inst.se_RTF1)**2
    # conditional Q (Sanderson-Windmeijer style): residual SS of precision-weighted
    # regression of one exposure's associations on the other's
    def condQ(bk, se_k, bo, se_o):
        wk = 1.0/se_k**2
        Xo = bo.values; Xk = bk.values
        W = np.diag(wk.values)
        XtWX = Xo @ W @ Xo
        beta_reg = (Xo @ W @ Xk) / XtWX
        resid = Xk - beta_reg*Xo
        return float(resid @ W @ resid)
    L = len(inst)
    Qc_i = condQ(inst.beta_ITPKA, inst.se_ITPKA, inst.beta_RTF1, inst.se_RTF1)
    Qc_r = condQ(inst.beta_RTF1, inst.se_RTF1, inst.beta_ITPKA, inst.se_ITPKA)
    res = mvmr_ivw(inst)
    bi, si, pi = ivw_uni(inst, 'beta_ITPKA')
    br, sr, pr = ivw_uni(inst, 'beta_RTF1')
    # single-instrument Wald ratios at the two leads
    def wald(rsid):
        row = inst[inst.rsid == rsid]
        if len(row) == 0:
            row = g[g.rsid == rsid]
        if len(row) == 0: return None
        row = row.iloc[0]
        out = {}
        for gg in ['ITPKA','RTF1']:
            b = row['beta']/row[f'beta_{gg}']; se = abs(row['sebeta']/row[f'beta_{gg}'])
            out[gg] = (b, se, 2*stats.norm.sf(abs(b/se)))
        return out
    w116 = wald('rs11635906'); w1942 = wald('rs1942')
    results.append(dict(endpoint=label, n_inst_pre=n_pre, n_inst=L,
        cor_beta=np.corrcoef(inst.beta_ITPKA, inst.beta_RTF1)[0,1],
        meanF_ITPKA=F_i.mean(), meanF_RTF1=F_r.mean(),
        condF_ITPKA=Qc_i/(L-1), condF_RTF1=Qc_r/(L-1),
        mvmr_b_ITPKA=res['ITPKA']['b'], mvmr_se_ITPKA=res['ITPKA']['se_re'], mvmr_p_ITPKA=res['ITPKA']['p_re'],
        mvmr_b_RTF1=res['RTF1']['b'], mvmr_se_RTF1=res['RTF1']['se_re'], mvmr_p_RTF1=res['RTF1']['p_re'],
        Q_resid=res['Q_resid'], Q_df=res['df'], Q_p=res['Q_p'],
        ivw_b_ITPKA=bi, ivw_p_ITPKA=pi, ivw_b_RTF1=br, ivw_p_RTF1=pr,
        wald116_ITPKA=w116['ITPKA'][0] if w116 else np.nan, wald116_p=w116['ITPKA'][2] if w116 else np.nan,
        wald1942_RTF1=w1942['RTF1'][0] if w1942 else np.nan, wald1942_p=w1942['RTF1'][2] if w1942 else np.nan))
    print(f"\n=== {label}: {n_pre} union instruments -> {L} after clump (r2<0.1) ===")
    print(f"cor(beta_ITPKA, beta_RTF1) across instruments = {results[-1]['cor_beta']:.3f}")
    print(f"mean F: ITPKA {F_i.mean():.1f}, RTF1 {F_r.mean():.1f} | cond F: ITPKA {Qc_i/(L-1):.1f}, RTF1 {Qc_r/(L-1):.1f}")
    print(f"MVMR-IVW: ITPKA b={res['ITPKA']['b']:+.4f} (SE {res['ITPKA']['se_re']:.4f}, p={res['ITPKA']['p_re']:.2e}) | "
          f"RTF1 b={res['RTF1']['b']:+.4f} (SE {res['RTF1']['se_re']:.4f}, p={res['RTF1']['p_re']:.2e})")
    print(f"Univariable IVW: ITPKA b={bi:+.4f} p={pi:.2e} | RTF1 b={br:+.4f} p={pr:.2e}")
    print(f"Residual Q={res['Q_resid']:.1f} (df {res['df']}, p={res['Q_p']:.3f})")
    if w116: print(f"Wald rs11635906 via ITPKA: b={w116['ITPKA'][0]:+.4f} p={w116['ITPKA'][2]:.2e}")
    if w1942: print(f"Wald rs1942 via RTF1: b={w1942['RTF1'][0]:+.4f} p={w1942['RTF1'][2]:.2e}")
    # sensitivity: looser clump r2<0.2 (WLS) and LD-aware GLS on r2<0.5 set
    inst2 = clump(g[(g.p_ITPKA < 5e-8) | (g.p_RTF1 < 5e-8)].copy(), 0.2)
    r2 = mvmr_ivw(inst2)
    print(f"[sens r2<0.2, L={len(inst2)}] MVMR: ITPKA b={r2['ITPKA']['b']:+.4f} p={r2['ITPKA']['p_re']:.2e} | "
          f"RTF1 b={r2['RTF1']['b']:+.4f} p={r2['RTF1']['p_re']:.2e}")
    inst5 = clump(g[(g.p_ITPKA < 5e-8) | (g.p_RTF1 < 5e-8)].copy(), 0.5)
    rg = mvmr_gls_ld(inst5)
    print(f"[sens LD-GLS, L={rg['n']}] MVMR: ITPKA b={rg['ITPKA']['b']:+.4f} p={rg['ITPKA']['p']:.2e} | "
          f"RTF1 b={rg['RTF1']['b']:+.4f} p={rg['RTF1']['p']:.2e}")
    results[-1].update(dict(
        sens02_b_ITPKA=r2['ITPKA']['b'], sens02_p_ITPKA=r2['ITPKA']['p_re'],
        sens02_b_RTF1=r2['RTF1']['b'], sens02_p_RTF1=r2['RTF1']['p_re'], sens02_L=len(inst2),
        gls_b_ITPKA=rg['ITPKA']['b'], gls_p_ITPKA=rg['ITPKA']['p'],
        gls_b_RTF1=rg['RTF1']['b'], gls_p_RTF1=rg['RTF1']['p'], gls_L=rg['n']))

pd.DataFrame(results).to_csv(f'{DIR}/cis_mvmr_results.csv', index=False)
print('\nsaved cis_mvmr_results.csv')
