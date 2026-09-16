#!/usr/bin/env python3
"""Step 2: reconstruct the 498-pair SASP x skin cis-MR scan (Supplementary Table S3).
Design faithful to manuscript Methods 2.2: lead cis-eQTL variant per gene (eQTLGen,
p<5e-8, clumped), single-instrument Wald ratio vs FinnGen R12 endpoints, first-order SE
(z_MR = z_GWAS). Harmonization by hg38 position + alleles (palindromes resolved by EAF,
ambiguous dropped), mirroring TwoSampleMR harmonise action=2."""
import csv, subprocess, math, os
from concurrent.futures import ThreadPoolExecutor
from pyliftover import LiftOver
from scipy import stats

BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release'
ENDPOINTS = ['L12_DERMATITISECZEMA', 'L12_ATOPIC', 'L12_ACTINKERA',
             'L12_ALOPECAREATA', 'L12_CHRONICULCEROFSKIN', 'L12_DECUBITANSULCERANDPRESSURE']
EP_LABEL = {'L12_DERMATITISECZEMA': 'dermatitis_eczema', 'L12_ATOPIC': 'atopic_dermatitis',
            'L12_ACTINKERA': 'actinic_keratosis', 'L12_ALOPECAREATA': 'alopecia_areata',
            'L12_CHRONICULCEROFSKIN': 'chronic_ulcer', 'L12_DECUBITANSULCERANDPRESSURE': 'decubitus_ulcer'}
OUT = '/workspace/scan498'

# ---- load lead variants, liftOver ----
lo = LiftOver(f'{OUT}/hg19ToHg38.over.chain.gz')
leads = []
with open(f'{OUT}/lead_variants_cis_hg19.csv') as f:
    for r in csv.DictReader(f):
        if r['status'] != 'ok':
            continue
        c = lo.convert_coordinate('chr' + r['chr_hg19'], int(r['pos_hg19']))
        r['chr_hg38'] = c[0][0].replace('chr', '') if c else ''
        r['pos_hg38'] = c[0][1] if c else ''
        leads.append(r)
print(f'{len(leads)} lead variants lifted')

# ---- tabix fetch ----
def fetch_one(args):
    ep, chrom, pos = args
    url = f'{BASE}/finngen_R12_{ep}.gz'
    try:
        out = subprocess.run(['tabix', url, f'{chrom}:{pos-500}-{pos+500}'],
                             capture_output=True, text=True, timeout=120).stdout
    except Exception:
        return (ep, chrom, pos, None)
    rows = []
    for ln in out.splitlines():
        if ln.startswith('#') or not ln.strip():
            continue
        f = ln.split('\t')
        if len(f) >= 13 and f[1] == str(pos):
            rows.append(f)
    return (ep, chrom, pos, rows)

tasks = [(ep, r['chr_hg38'], int(r['pos_hg38'])) for ep in ENDPOINTS for r in leads if r['pos_hg38'] != '']
print(f'{len(tasks)} tabix queries')
results = {}
with ThreadPoolExecutor(max_workers=8) as ex:
    for i, (ep, chrom, pos, rows) in enumerate(ex.map(fetch_one, tasks)):
        results[(ep, chrom, pos)] = rows
        if (i + 1) % 100 == 0:
            print(f'  {i+1}/{len(tasks)}')

# ---- harmonize + Wald ----
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def palindromic(a, b):
    return COMP.get(a) == b

def harmonize(ea, nea, eaf, ref, alt, beta, af_alt):
    """return (b_gwas_for_ea, af_for_ea, status)"""
    ea, nea, ref, alt = ea.upper(), nea.upper(), ref.upper(), alt.upper()
    if len(ea) > 1 or len(nea) > 1 or len(ref) > 1 or len(alt) > 1:
        return None, None, 'indel_dropped'
    pal = palindromic(ea, nea)
    if ea == alt and nea == ref:
        if pal and eaf != '' and af_alt != '':
            e, a = float(eaf), float(af_alt)
            if 0.42 < e < 0.58 and 0.42 < a < 0.58:
                return None, None, 'palindrome_ambiguous'
            if abs(e - a) > abs(e - (1 - a)):  # freq says flipped
                return -beta, 1 - a, 'ok_flipped_pal'
        return beta, af_alt, 'ok'
    if ea == ref and nea == alt:
        if pal and eaf != '' and af_alt != '':
            e, a = float(eaf), float(af_alt)
            if 0.42 < e < 0.58 and 0.42 < a < 0.58:
                return None, None, 'palindrome_ambiguous'
            if abs(e - (1 - a)) > abs(e - a):
                return beta, a, 'ok_pal'
        return -beta, (1 - af_alt) if af_alt != '' else '', 'ok_flipped'
    return None, None, 'allele_mismatch'

pairs = []
for r in leads:
    b_e, se_e = float(r['beta_eqtl']), float(r['se_eqtl'])
    for ep in ENDPOINTS:
        key = (ep, r['chr_hg38'], int(r['pos_hg38']))
        rows = results.get(key)
        rec = {'pair_id': f"sasp_{r['gene']}_x_{EP_LABEL[ep]}", 'gene': r['gene'],
               'exposure_id': r['exposure_id'], 'outcome_id': ep, 'outcome': EP_LABEL[ep],
               'lead_rsid': r['rsid'], 'chr_hg38': r['chr_hg38'], 'pos_hg38': r['pos_hg38'],
               'ea': r['ea'], 'nea': r['nea'], 'eaf_eqtl': r['eaf'],
               'beta_eqtl': b_e, 'se_eqtl': se_e, 'p_eqtl': float(r['p_eqtl']),
               'beta_gwas': '', 'se_gwas': '', 'p_gwas': '', 'af_outcome': '',
               'b_mr': '', 'se_mr_fo': '', 'p_mr': '', 'or_mr': '', 'or_lo': '', 'or_hi': '',
               'se_mr_delta': '', 'p_mr_delta': '', 'status': ''}
        if not rows:
            rec['status'] = 'variant_not_in_outcome'
        else:
            # prefer allele-matched row at multiallelic sites
            chosen = None
            for f in rows:
                ref, alt = f[2].upper(), f[3].upper()
                if {r['ea'], r['nea']} == {ref, alt}:
                    chosen = f; break
            if chosen is None and len(rows) == 1:
                chosen = rows[0]
            if chosen is None:
                rec['status'] = 'allele_mismatch_multiallelic'
            else:
                ref, alt = chosen[2].upper(), chosen[3].upper()
                p_g, b_g, se_g = float(chosen[6]), float(chosen[8]), float(chosen[9])
                af_alt = chosen[11]
                b_h, af_ea, st = harmonize(r['ea'], r['nea'], r['eaf'], ref, alt, b_g, af_alt)
                if b_h is None:
                    rec['status'] = st
                else:
                    b_mr = b_h / b_e
                    se_fo = se_g / abs(b_e)
                    z = b_mr / se_fo
                    p_mr = 2 * stats.norm.sf(abs(z))
                    se_d = math.sqrt(se_g**2 / b_e**2 + b_h**2 * se_e**2 / b_e**4)
                    p_d = 2 * stats.norm.sf(abs(b_mr / se_d))
                    rec.update({'beta_gwas': b_h, 'se_gwas': se_g, 'p_gwas': p_g,
                                'af_outcome': af_ea, 'b_mr': b_mr, 'se_mr_fo': se_fo,
                                'p_mr': p_mr, 'or_mr': math.exp(b_mr),
                                'or_lo': math.exp(b_mr - 1.959964 * se_fo),
                                'or_hi': math.exp(b_mr + 1.959964 * se_fo),
                                'se_mr_delta': se_d, 'p_mr_delta': p_d,
                                'status': 'ok' if st.startswith('ok') else st})
        pairs.append(rec)

# genes without instrument
n_noinst = 0
with open(f'{OUT}/lead_variants_cis_hg19.csv') as f:
    for r in csv.DictReader(f):
        if r['status'] != 'ok':
            n_noinst += 1
            for ep in ENDPOINTS:
                pairs.append({'pair_id': f"sasp_{r['gene']}_x_{EP_LABEL[ep]}", 'gene': r['gene'],
                              'exposure_id': r['exposure_id'], 'outcome_id': ep,
                              'outcome': EP_LABEL[ep], 'lead_rsid': '', 'chr_hg38': '', 'pos_hg38': '',
                              'ea': '', 'nea': '', 'eaf_eqtl': '', 'beta_eqtl': '', 'se_eqtl': '',
                              'p_eqtl': '', 'beta_gwas': '', 'se_gwas': '', 'p_gwas': '',
                              'af_outcome': '', 'b_mr': '', 'se_mr_fo': '', 'p_mr': '', 'or_mr': '',
                              'or_lo': '', 'or_hi': '', 'se_mr_delta': '', 'p_mr_delta': '',
                              'status': 'no_instrument_p5e-8'})

# ---- multiple testing ----
ok = [p for p in pairs if p['status'].startswith('ok')]
ps = sorted(float(p['p_mr']) for p in ok)
m = len(ps)
# BH q-values
import numpy as np
order = np.argsort([float(p['p_mr']) for p in ok])
ranks = np.empty(m, int); ranks[order] = np.arange(1, m + 1)
qvals = np.minimum.accumulate(
    (np.array([float(p['p_mr']) for p in ok])[order][::-1] * m / np.arange(m, 0, -1)))[::-1]
q = np.empty(m); q[order] = qvals
for p_, q_ in zip(ok, q):
    p_['q_mr_BH'] = float(q_)
    p_['pass_bonferroni'] = float(p_['p_mr']) < 0.05 / 498
for p_ in pairs:
    if 'q_mr_BH' not in p_:
        p_['q_mr_BH'] = ''; p_['pass_bonferroni'] = ''

bonf_hits = [(p_['pair_id'], p_['p_mr']) for p_ in ok if p_['pass_bonferroni']]
bh_hits = [(p_['pair_id'], p_['q_mr_BH']) for p_ in ok if p_['q_mr_BH'] != '' and p_['q_mr_BH'] < 0.05]
print(f'\ntested pairs (harmonized): {m}/498; no-instrument genes: {n_noinst}')
print(f'Bonferroni hits ({len(bonf_hits)}):')
for h in sorted(bonf_hits, key=lambda x: x[1]): print('  ', h[0], f'{h[1]:.2e}')
print(f'BH-FDR<0.05 hits ({len(bh_hits)}):')
for h in sorted(bh_hits, key=lambda x: x[1])[:15]: print('  ', h[0], f'{h[1]:.2e}')

flds = list(pairs[0].keys()) + ['q_mr_BH', 'pass_bonferroni']
with open(f'{OUT}/scan498_results_cis.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=flds)
    w.writeheader()
    for p_ in sorted(pairs, key=lambda x: (x['gene'], x['outcome'])):
        w.writerow(p_)
print(f'wrote {OUT}/scan498_results_cis.csv ({len(pairs)} rows)')
