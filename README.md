# Code repository: Blood ITPKA and the risk of dermatitis and eczema (P2-02)

Analysis code for the manuscript:

> An Y, Li Y, Zhong X, Tang M. *Genetically proxied up-regulation of blood ITPKA and the risk of dermatitis and eczema: a cis-Mendelian randomization and colocalization study with locus-wide multi-gene and protein-level discrimination.* (2026, under review)

All analyses use **publicly available summary-level data**; no individual-level data are required. Software: Python 3.11 (pandas, numpy, scipy, matplotlib 3.10.9, pysam, pyliftover) and R 4.4.3 (coloc 5.2.3, susieR, data.table).

## Repository layout

| Directory | Contents | Main outputs |
|---|---|---|
| `01_scan_and_triage/` | cis-MR scan of 498 SASP-gene × skin-disease pairs (eQTLGen blood × FinnGen R12), Bonferroni/FDR triage, triage colocalization | Suppl. Table S3 |
| `02_main_coloc_panel/` | Main ITPKA colocalization (coloc.abf / coloc.susie, EUR + FIN LD panels), twelve-gene panel colocalization across eQTL resources and GWAS endpoints, window/prior sensitivity, cis-MVMR (LD-aware GLS + clumped IVW) | Tables 2, 3, 5; Suppl. Tables S4, S6, S13 |
| `03_pqtl_and_signal_attribution/` | Plasma pQTL colocalization (ITPKA, TYRO3; INTERVAL / Icelandic AGES / HELIC-MANOLIS), cystatin C signal attribution (region exclusion, conditioning), EAGLE full-summary-statistic colocalization, FDR-hit colocalization (FGF2/CSF2RB/EDN1), FIN-LD sensitivity, IL6ST/ANKRD55 positive-control discrimination, ±1 Mb panel extension | Suppl. Tables S9–S13, S17 |
| `04_phewas_immune_census/` | OpenGWAS phenome-wide scan processing, trait-level colocalization (9 non-atopic/blood-cell traits), purified immune-cell eQTL survey (33 datasets), CZ CELLxGENE Census ITPKA/ITPKB expression queries | Suppl. Tables S2, S7, S8, S16 |
| `05_figures/` | Final figure scripts (matplotlib; Figures 1, 2, 3, 4, 5, S2) | Figures (PNG 300 dpi + SVG) |
| `06_tables_manuscript/` | Supplementary-table builder, manuscript docx/review-PDF builders, reference renumbering/verification | Suppl. Tables S1–S17; manuscript builds |
| `notebook_cells_archive/` | Chronological dump of all executed analysis notebook cells (including SMR-HEIDI runs via the SMR CLI, single-cell queries, and one-off verification probes) | — |

## Data sources (all public)

- **eQTL**: eQTLGen phase I cis-eQTL (via OpenGWAS, dataset `eqtl-a-*`); GTEx v8 (GTEx Portal) and GTEx v10 / BLUEPRINT / DICE / OneK1K / GENCORD / CEDAR / Kasela / Bossini-Castillo via the eQTL Catalogue (https://www.ebi.ac.uk/eqtl/)
- **GWAS**: FinnGen R12 (L12_DERMATITISECZEMA, L12_ATOPIC, and extended skin endpoints; https://www.finngen.fi/); Oliva et al. 2025 AD meta-analysis (GWAS Catalog GCST90503108–90503111); EAGLE Eczema Consortium (public release)
- **pQTL**: INTERVAL (Sun 2018), Icelandic AGES (Gudjonsson 2022), HELIC-MANOLIS (Gilly 2020), via GWAS Catalog
- **LD**: 1000 Genomes Project high-coverage hg38 (EUR n=525; FIN n=99)
- **PheWAS**: OpenGWAS API
- **Single-cell**: CZ CELLxGENE Census (2025-01-30); He et al. GSE147424
- **Regulatory annotation**: Open Targets Platform (VEP, L2G), ENCODE-rE2G, Javierre 2016 PCHi-C (OSF peak matrix)
- **Expression localization**: Human Protein Atlas; GTEx median TPM

## Reproducibility notes

- Scripts read local copies of the summary statistics listed above; file paths at the top of each script indicate the expected inputs.
- The tri-allelic index variant rs11635906 (A/C/G) requires allele-aware harmonization; see `02_main_coloc_panel/run_coloc_eqtlgen_fixed.R` for the harmonization logic used throughout.
- All random-number-free analyses (deterministic); coloc.susie used default settings with 1000 Genomes LD.

## License

MIT (code). The manuscript text and figures are © the authors.
