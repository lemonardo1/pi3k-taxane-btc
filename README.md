# Derived data and analysis code

**PI3K pathway alterations and attenuated taxane benefit in advanced biliary tract cancer**

Kim D, Chon HJ, Kim CG, Choi HJ. Submitted to *npj Precision Oncology* (2026).

This record contains every derived dataset and all analysis code required to reproduce the
results reported in the manuscript. It does **not** contain individual-level clinical or
sequencing data from the Yonsei index cohort (n = 287) or the CHA validation cohort (n = 167):
the consent and institutional approvals under which those cohorts were accrued
(Yonsei 4-2023-0485; CHA 2021-01-010, 2021-03-045, 2022-11-002) do not permit open
redistribution. De-identified data supporting the findings are available from the corresponding
author on reasonable request, subject to an institutional data-sharing agreement and approval by
the respective institutional review boards.

The ten public cohorts were obtained from cBioPortal (https://www.cbioportal.org/datasets);
mutation calls were retrieved through the cBioPortal REST API. `data/harmonized_public_btc.csv`
is the harmonised patient-level table produced from them.

## Contents

| Path | Description |
|---|---|
| `data/harmonized_public_btc.csv` | Harmonised patient-level table for the ten public biliary tract cancer cohorts: treatment annotation, PI3K pathway status, survival. |
| `data/pi3k_dms_scan.csv.gz` | Complete in-silico saturation mutagenesis background: 283,442 ESM-2 log-likelihood-ratio scores across the fourteen pathway genes. |
| `data/gene_level_weights.csv` | Empirical-Bayes shrunk gene-level weights (k = 4) used to build the continuous PI3K score. |
| `data/variant_function_scores.csv` | Variant-level mode-aware functional scores for the observed PI3K pathway variants, with CIViC status. |
| `data/evo2_variant_scores.csv` | DNA-level Evo 2 scores for the 56 observed missense variants (8,193-bp window). |
| `data/pathway_alteration_events.csv` | Per-event pathway alteration calls underlying the public-cohort prevalence and prognostic analyses. |
| `structure/pi3k_p110a_p85_complex.pdb` | Boltz-2 predicted p110alpha-p85 niSH2 complex coordinates. |
| `structure/pi3k_interface_metrics.csv` | Per-residue interface geometry computed on the predicted complex. |
| `structure/pik3ca_structure_variants.csv` | Per-variant distance to the p85 regulatory interface for the PIK3CA variants analysed. |
| `supplementary_tables/TableS1_public_cohort_inventory.csv` | Supplementary Table S1 source data: public cohort inventory. |
| `supplementary_tables/TableS2_pathway_and_gene_screens.csv` | Supplementary Table S2 source data. |
| `supplementary_tables/TableS3_pik3ca_structural_metrics.csv` | Supplementary Table S3 source data. |
| `supplementary_tables/TableS7_baseline_by_treatment_arm.csv` | Supplementary Table S7 source data. |
| `supplementary_tables/permutation_pathway_screen.csv` | Permutation pathway screen (10,000 replicates), full output. |
| `supplementary_tables/permutation_gene_screen.csv` | Permutation gene-level screen, full output. |
| `supplementary_tables/sensitivity_analyses.csv` | Sensitivity and robustness analyses (Supplementary Table S6 source data). |
| `supplementary_tables/trial_design_simulation.csv` | Trial-design event-requirement simulation output (Figure 4b). |
| `code/reproduce_analysis.py` | End-to-end reproduction script for the reported estimates. |
| `code/build_manuscript.py` | Figure and PDF build script. |
| `code/build_si.py` | Figure and PDF build script. |
| `code/build_cover.py` | Figure and PDF build script. |
| `code/build_figure1.py` | Figure and PDF build script. |
| `code/npjpdf.py` | Figure and PDF build script. |
| `code/environment.yml` | Conda environment specification with version pins. |

`MANIFEST.csv` lists the size and SHA-256 checksum of every file.

## Reproducing the reported estimates

    conda env create -f code/environment.yml
    conda activate pi3k-taxane-btc
    python code/reproduce_analysis.py

The public-cohort analyses (prevalence, prognostic pooling, variant scoring, structural and
DNA-level axes) run end to end from the files in this record. The index- and validation-cohort
interaction models require the restricted clinical tables described above.

## Licence

Data files: CC BY 4.0. Code: MIT.

## Citation

Please cite the manuscript above and this record by its Zenodo DOI.
