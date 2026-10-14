# Broad but not specific: a negative-control-calibrated phenome-wide atlas of the triglyceride–glucose index and other insulin resistance surrogates in UK Biobank and NHANES

Analysis code for the study above.

The study compares 12 insulin resistance surrogate indices in 374,814 UK Biobank participants without diabetes:

- TyG, TyG-BMI, TyG-WC, TyG-WHtR, TyG-ABSI
- TG/HDL-C, AIP, METS-IR, CHG, CTI
- eGDR (sign-inverted), with and without its hypertension term

It also analyses their six components (ln triglycerides, ln glucose, HDL-C, BMI, waist, HbA1c), ApoB and a fasting-standardised TyG. Each is related to 540 incident diseases (3-character ICD-10 codes from first-occurrence fields and the cancer registry) and to all-cause and cardiovascular death.

Main steps:

1. **Phenome-wide Cox models.** Model 1 adjusts for demographics and lifestyle and is stratified by assessment centre; model 2 adds BMI. In total 17,886 models are fitted.
2. **Empirical calibration.** Systematic error N(μ, τ²) is estimated from 60 pre-specified negative-control outcomes (and from a strict 40-outcome subset), giving calibrated *P* values and FDR. A calibrated FDR below 0.05 defines a specific association.
3. **Further analyses:**
   - component decomposition
   - incremental C-index over demographic and clinical models
   - spline dose–response curves and subgroups
   - abdominal and cardiac MRI
   - TyG–ApoB discordance
4. **External anchoring in NHANES 2005–2018.** Construct validity against HOMA-IR, and linked mortality with and without HOMA-IR adjustment.
5. **Mendelian randomisation.** Triglycerides (GLGC 2013), BMI and WHR adjusted for BMI (GIANT 2015) against 10 FinnGen R12 endpoints; univariable IVW, MR-Egger and weighted median, and multivariable MR.

## Data

This research was conducted using the UK Biobank Resource under application number 1099827. No individual-level data, results or manuscript files are included in this repository. UK Biobank data are available to approved researchers through https://www.ukbiobank.ac.uk.

Public resources used:

| Resource | Use | Source |
|---|---|---|
| NHANES 2005–2018 (DEMO, GLU, INS, TRIGLY, HDL, TCHOL, GHB, BMX, BPX, BPQ, DIQ, CRP/HSCRP) | HOMA-IR anchoring | https://wwwn.cdc.gov/nchs/nhanes/ |
| NHANES public-use linked mortality files (2019) | Mortality follow-up | https://www.cdc.gov/nchs/data-linkage/mortality-public.htm |
| UK Biobank coding 19 (ICD-10) | Cancer-registry code titles | https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19 |
| GLGC 2013 triglycerides (`jointGwasMc_TG.txt.gz`) | MR exposure | http://csg.sph.umich.edu/willer/public/lipids2013/ |
| GIANT 2015 BMI, European (`SNP_gwas_mc_merge_nogc.tbl.uniq.gz`) and WHR adjusted for BMI (`GIANT_2015_WHRadjBMI_COMBINED_EUR.txt`) | MR exposures | https://portals.broadinstitute.org/collaboration/giant/ |
| FinnGen R12 summary statistics and tabix indexes | MR outcomes (I9_CHD, I9_HEARTFAIL, I9_AF, I9_STR, T2D, M13_GOUT, K11_CHOLELITH, NAFLD, N14_CHRONKIDNEYDIS, K11_ACUTPANC) | https://storage.googleapis.com/finngen-public-data-r12/summary_stats/ |

## Setup

Paths are written as placeholders. Before running, replace them as follows:

| Placeholder | Folder |
|---|---|
| `/path/to/ukb_data` | The UK Biobank extract |
| `/path/to/project` | This project (with `cache/`, `results/`, `figures/`, `external/`) |
| `/path/to/upstream` | Intermediate files listed below |

Put the NHANES `.xpt` files and mortality `.dat` files in `external/nhanes/`, coding 19 in `external/coding19.tsv`, the GWAS files in `external/gwas/glgc2013_tg/`, `external/gwas/giant2015_bmi/` and `external/gwas/giant2015_whradjbmi/`, and the FinnGen `.tbi` indexes in `external/finngen/`. `19_mr_finngen_fetch.py` reads only the BGZF blocks around each instrument with HTTP range requests, so the full FinnGen files do not have to be downloaded. Run the scripts from `scripts/` in numerical order. Run `14_nhanes.py` and `18`–`20` with `python3 -P`.

The UK Biobank extract is expected as comma-separated tables keyed by `eid`, with columns named by field (`p<field>_i<instance>_a<array>`), as exported from the Research Analysis Platform:

| File | Contents |
|---|---|
| `ukb_assessment_center.csv` | Baseline and imaging-visit assessment fields (height, diastolic blood pressure, abdominal and cardiac MRI) |
| `ukb_health_outcomes.csv` | First-occurrence dates, cancer registry (fields 40005, 40006) and death registry |
| `ukb_data_dictionary.csv` | Showcase data dictionary |

Intermediate files expected under `/path/to/upstream` (instance 0 unless stated):

| File | Contents |
|---|---|
| `ddw_ppdm/cache/cohort_full.pkl` | Baseline participants with baseline date, sex, birth date, BMI, smoking, Townsend index and ethnicity |
| `cd_screening/cache/bs_lipgluc.pkl` | Fields 30630, 30640, 30690, 30740, 30750, 30760, 30780, 30870 |
| `diet_barrier/cache/biochem.pkl` | Fields 30620, 30650, 30700, 30710, 30720, 30750 |
| `gut_hormone/cache/ac_cov.pkl` | Fields 48, 49, 54, 55, 74, 93, 1558, 2443, 2976, 3166, 4080, 21000, 22040 |
| `idea2_feasibility/cache/extra_cov.pkl` | Fields 74, 845, 2724, 2814, 6153, 6177, 30700, 30710, 30720, 30900, 30901, 30903 |

## Scripts

| Script | Step |
|---|---|
| `01`–`02` | Extract first-occurrence dates, height and blood pressure; base cohort |
| `03_phenome_outcomes.py` | All first-occurrence and cancer-registry outcomes as days from baseline; outcome metadata (≥200 incident events, chapters A–N) |
| `04_analysis_table.py` | Indices, components, ApoB, fasting-standardised TyG, covariates |
| `05_participant_flow.py` | Exclusion counts |
| `06_phewas_cox.R` | Phenome-wide Cox models (R `survival`, parallel) |
| `07`–`08` | Negative-control calibration (60 and 40 negative controls) |
| `09_component_decomposition.R` | Mutually adjusted components for every outcome |
| `10_incremental_prediction.R` | Harrell's C increments over models A and B |
| `11_dose_response_subgroups.R` | Natural-spline dose–response; subgroups by sex, age, ethnicity, BMI and fasting |
| `12_mri.py` | Indices and imaging phenotypes |
| `13_tyg_apob.R` | TyG–ApoB discordance |
| `14_nhanes.py` | NHANES HOMA-IR anchoring and mortality |
| `15_tables.py` | Tables 1 and 2, follow-up and event counts |
| `16`–`17` | Main and supplementary figures |
| `18_mr_instruments.py` | Genome-wide significant SNPs, GRCh38 positions from Ensembl REST, distance clumping (1 Mb), MHC excluded |
| `19_mr_finngen_fetch.py` | FinnGen R12 lookups by HTTP range reads |
| `20_mr_analysis.py` | Harmonisation, IVW, MR-Egger, weighted median, multivariable MR, conditional F |
| `21_mr_figure.py` | MR figure |

## Requirements

- Python 3 with pandas, numpy, scipy, statsmodels, lifelines, matplotlib and Pillow
- R with data.table, survival, splines and parallel

## License

MIT
