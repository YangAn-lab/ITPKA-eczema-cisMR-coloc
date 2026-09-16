import pypdf
r = pypdf.PdfReader('/workspace/P2-02_manuscript_review_v1.pdf')
targets = {'Figure 1. Study design': None, 'Table 1. Mendelian': None, 'Figure 2. Main MR': None,
           'Table 3. Locus-wide': None, 'Figure S2. Signal structure': None, 'Figure 4. Tissue': None,
           'Figure 5. Mechanistic': None, 'Supplementary Table S1': None, 'Supplementary Table S2': None,
           'References': None}
for i, pg in enumerate(r.pages, 1):
    t = pg.extract_text() or ''
    for k in targets:
        if targets[k] is None and k in t:
            targets[k] = i
for k, v in targets.items():
    print(f'{k:35s} -> page {v}')