"""Verify all manuscript references against Crossref.

Part 1: the 44 existing references (DOI lookup + metadata match).
Part 2: candidate new references (bibliographic search -> resolve DOI -> verify).
Writes a verification table to /workspace/ref_verification.csv
"""
import json
import time
import urllib.request
import urllib.parse

HDR = {"User-Agent": "ref-check/1.0 (mailto:research@example.org)"}


def crossref_doi(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers=HDR)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)["message"]
    except Exception as e:
        return {"_error": str(e)}


def crossref_query(title, first_author=None, rows=3):
    q = urllib.parse.urlencode({"query.bibliographic": title, "rows": rows})
    req = urllib.request.Request("https://api.crossref.org/works?" + q, headers=HDR)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)["message"]["items"]
    except Exception as e:
        return [{"_error": str(e)}]


# ---- Part 1: existing 44 references: (n, doi, expected_journal, expected_year)
refs = [
    (1, "10.1038/s41467-025-58310-7", "Nat Commun", 2025),
    (2, "10.1038/s41467-022-32552-1", "Nat Commun", 2022),
    (3, "10.1038/s41588-021-00913-z", "Nat Genet", 2021),
    (4, "10.1038/s41586-022-05473-8", "Nature", 2023),
    (5, "10.1016/j.cels.2015.12.004", "Cell Syst", 2015),
    (6, "10.1074/jbc.M109.047050", "J Biol Chem", 2010),
    (7, "10.1016/j.cellsig.2011.11.010", "Cell Signal", 2012),
    (8, "10.1016/j.jbior.2015.09.004", "Adv Biol Regul", 2016),
    (9, "10.1007/s00018-009-0238-5", "Cell Mol Life Sci", 2010),
    (10, "10.1101/2020.08.10.244293", "bioRxiv", 2020),
    (11, "10.1126/science.aaz1776", "Science", 2020),
    (12, "10.1038/s41588-021-00924-w", "Nat Genet", 2021),
    (13, "10.1016/j.cell.2016.10.026", "Cell", 2016),
    (14, "10.1016/j.cell.2022.08.004", "Cell", 2022),
    (15, "10.1371/journal.pgen.1004383", "PLoS Genet", 2014),
    (16, "10.1371/journal.pgen.1008720", "PLoS Genet", 2020),
    (17, "10.1111/rssb.12388", "J R Stat Soc Series B", 2020),
    (18, "10.1016/j.jid.2018.12.018", "J Invest Dermatol", 2019),
    (19, "10.1016/j.jaci.2020.06.012", "J Allergy Clin Immunol", 2021),
    (20, "10.1126/science.1260419", "Science", 2015),
    (21, "10.1126/sciadv.abh2169", "Sci Adv", 2021),
    (22, "10.7554/eLife.34408", "eLife", 2018),
    (23, "10.1093/hmg/ddu328", "Hum Mol Genet", 2014),
    (24, "10.1186/s12979-026-00580-w", "Immun Ageing", 2026),
    (25, "10.1016/j.jdermsci.2024.04.002", "J Dermatol Sci", 2024),
    (26, "10.1038/ng.3424", "Nat Genet", 2015),
    (27, "10.1093/nar/gkae1142", "Nucleic Acids Res", 2025),
    (28, "10.1016/j.jaci.2020.01.042", "J Allergy Clin Immunol", 2020),
    (29, "10.1038/jid.2012.456", "J Invest Dermatol", 2013),
    (30, "10.1016/j.jaci.2015.01.020", "J Allergy Clin Immunol", 2015),
    (31, "10.1016/j.jdermsci.2025.09.006", "J Dermatol Sci", 2025),
    (32, "10.1111/all.13223", "Allergy", 2017),
    (33, "10.3390/ijms27052371", "Int J Mol Sci", 2026),
    (34, "10.1038/s41590-021-00927-z", "Nat Immunol", 2021),
    (35, "10.3389/fcell.2022.835675", "Front Cell Dev Biol", 2022),
    (36, "10.2147/jaa.s579390", "J Asthma Allergy", 2026),
    (37, "10.3390/ijms22115729", "Int J Mol Sci", 2021),
    (38, "10.1111/acel.70527", "Aging Cell", 2026),
    (39, "10.1073/pnas.0306907101", "Proc Natl Acad Sci USA", 2004),
    (40, "10.1038/ni980", "Nat Immunol", 2003),
    (41, "10.1371/journal.pone.0131071", "PLoS One", 2015),
    (42, "10.1016/j.jbior.2012.08.001", "Adv Biol Regul", 2013),
    (43, "10.1038/ng.2213", "Nat Genet", 2012),
    (44, "10.1007/s12016-024-09004-3", "Clin Rev Allergy Immunol", 2024),
]

rows = []
for n, doi, ej, ey in refs:
    m = crossref_doi(doi)
    if "_error" in m:
        rows.append({"n": n, "doi": doi, "status": "DOI_FAIL", "detail": m["_error"]})
    else:
        title = (m.get("title") or [""])[0]
        journal = (m.get("container-title") or [""])[0]
        year = None
        for k in ("published-print", "published-online", "issued"):
            if k in m and m[k].get("date-parts"):
                year = m[k]["date-parts"][0][0]
                break
        vol = m.get("volume", "")
        page = m.get("page", m.get("article-number", ""))
        auth = m.get("author", [])
        first = (auth[0].get("family", "") if auth else "")
        ok = (ey is None or year == ey or (ey and year and abs(year - ey) <= 1))
        rows.append({"n": n, "doi": doi, "status": "OK" if ok else "YEAR_MISMATCH",
                     "first_author": first, "year": year, "journal": journal,
                     "volume": vol, "page": page, "title": title[:110]})
    time.sleep(0.3)

# ---- Part 2: candidate new references (search + verify)
candidates = [
    ("epidemiology", "Atopic dermatitis Weidinger Novak Lancet 2016"),
    ("twin_heritability_1", "Importance of genetic factors in the etiology of atopic dermatitis: a twin study Thomsen"),
    ("twin_heritability_2", "Twin Studies of Atopic Dermatitis Interpretations and Applications in the Filaggrin Era Elmose Thomsen"),
    ("NFAT", "NFAT proteins: key regulators of T-cell development and function Macian"),
    ("GWAS_Catalog", "The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource Sollis"),
    ("weak_instruments", "Avoiding bias from weak instruments in Mendelian randomization studies Burgess Thompson"),
    ("heterogeneity_I2", "Quantifying heterogeneity in a meta-analysis Higgins Thompson"),
]

cand_rows = []
for tag, q in candidates:
    items = crossref_query(q)
    for it in items[:2]:
        if "_error" in it:
            cand_rows.append({"tag": tag, "status": "QUERY_FAIL", "detail": it["_error"]})
            continue
        title = (it.get("title") or [""])[0]
        journal = (it.get("container-title") or [""])[0]
        year = None
        for k in ("published-print", "published-online", "issued"):
            if k in it and it[k].get("date-parts"):
                year = it[k]["date-parts"][0][0]
                break
        auth = it.get("author", [])
        first = (auth[0].get("family", "") if auth else "")
        cand_rows.append({"tag": tag, "status": "FOUND", "doi": it.get("DOI", ""),
                          "first_author": first, "year": year, "journal": journal,
                          "volume": it.get("volume", ""), "page": it.get("page", ""),
                          "title": title[:120]})
    time.sleep(0.3)

import csv
with open("/workspace/ref_verification.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["n", "doi", "status", "first_author", "year",
                                      "journal", "volume", "page", "title", "detail"],
                       extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)

with open("/workspace/ref_candidates.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tag", "status", "doi", "first_author", "year",
                                      "journal", "volume", "page", "title", "detail"],
                       extrasaction="ignore")
    w.writeheader()
    for r in cand_rows:
        w.writerow(r)

print("=== EXISTING 44 ===")
for r in rows:
    print(r.get("n"), r.get("status"), "|", r.get("first_author", ""), r.get("year", ""),
          "|", r.get("journal", ""), r.get("volume", ""), r.get("page", ""),
          "|", (r.get("title") or r.get("detail", ""))[:80])
print()
print("=== CANDIDATES ===")
for r in cand_rows:
    print(r.get("tag"), "|", r.get("status"), "|", r.get("doi", ""), "|",
          r.get("first_author", ""), r.get("year", ""), "|", r.get("journal", ""),
          r.get("volume", ""), r.get("page", ""), "|", (r.get("title") or "")[:80])
