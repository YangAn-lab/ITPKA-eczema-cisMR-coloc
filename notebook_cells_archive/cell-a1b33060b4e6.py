for i, pg in enumerate(r.pages, 1):
    t = pg.extract_text() or ''
    hits = []
    if 'dilution' in t.lower() and i > 21: hits.append('S1-table?')
    if 'trait' in t.lower() and ('category' in t.lower() or 'opengwas' in t.lower()) and i > 21: hits.append('S2-table?')
    if hits: print(i, hits)
# also show last page number content head
print('last page head:', (r.pages[-1].extract_text() or '')[:120])