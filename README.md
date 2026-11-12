# ZF-Temporal-RNAseq

Gene-level RNA-seq data, analysis scripts and processed results for:

> Berkhout JH, Caiment F, Piersma AH, Legler J, Heusinkveld HJ. *Time-resolved transcriptomic signatures of zebrafish neurulation.* Manuscript under review.
<!-- TODO (Job): update title, journal and DOI on acceptance. -->

The study profiles whole zebrafish embryos at five stages of neurulation, in controls and after exposure to eight chemicals. It provides a stage-resolved reference for neural tube closure (NTC) and maps how each exposure shifts NTC-related pathway activity over time.

**Raw sequencing reads (FASTQ)** are deposited in BioStudies/ArrayExpress under accession **`E-MTAB-XXXXX`** <!-- TODO (Job): insert accession -->. This repository holds everything downstream of read quantification.

---

## Study design

| | |
|---|---|
| Organism | *Danio rerio*, wild-type AB strain |
| Sample | Pool of 30 whole embryos per sample |
| Exposure | Waterborne, in HEPES-buffered E3, from ≤ 3 hpf until sampling |
| Stages | **T1** bud (~10 hpf) · **T2** 5–8 somites · **T3** 11–14 somites · **T4** 18–21 somites · **T5** 24 hpf |
| Treatments | 8 chemicals × 3 concentrations × 5 stages, 2 biological replicates (independent clutches) |
| Controls | DMSO (0.01% v/v) and E3 vehicle controls, 4 replicates per vehicle per stage |
| Libraries | 280 (240 treated + 40 controls) |
| Sequencing | Utrecht Sequencing Facility (USEQ): TruSeq Stranded mRNA (poly(A)), Illumina NovaSeq X, 2 × 150 bp |
| Processing | fastp (v0.19.8) trimming and QC, STAR alignment and RSEM (v1.3.1) quantification against GRCz11 (Ensembl release 114) |

`control_pair_id` in the metadata links each treated sample to the vehicle control collected from the same clutch on the same day.

### Chemicals and dose codes

**S1 is the highest concentration and S3 the lowest.**

| Code | Chemical | Vehicle | S3 (low, µM) | S2 (mid, µM) | S1 (high, µM) |
|---|---|---|---:|---:|---:|
| `2-4D` | 2,4-Dichlorophenoxyacetic acid | E3 | 10 | 30 | 100 |
| `CAFF` | Caffeine | E3 | 100 | 300 | 1000 |
| `CYCL` | Cyclopamine | DMSO | 1 | 5 | 10 |
| `GLYP` | Glyphosate | E3 | 20 | 60 | 200 |
| `IMID` | Imidacloprid | E3 | 150 | 500 | 1500 |
| `RACI` | All-trans retinoic acid (ATRA) | DMSO | 0.003 | 0.01 | 0.03 |
| `TEBU` | Tebuconazole | DMSO | 6 | 20 | 60 |
| `THIR` | Thiram | DMSO | 0.01 | 0.03 | 0.1 |

### Sample IDs

`<tube>-<date>-<stage>-<chemical>-<dose>`, for example `112-0402-T1-CYCL-S1`. Controls end in `CTRL-DM` (DMSO) or `E3-CTRL` (E3).

The date code is MMDD for February–April and DDMM for May (`0705`, `1305`, `1505`, `2005`). Use the `experiment_date` column in the metadata rather than parsing it from the ID.

---

## Repository layout

```
ZF_RNA_Raw_Count.txt          gene x sample count matrix, all 280 samples
norm_data_controls_ALL.txt    size-factor-normalised counts, 40 control samples
metadata/                     sample metadata
DEGs_NTC/                     DESeq2 results, one file per chemical-dose-stage contrast
Scripts/                      R analysis scripts
Scripts/figures/              scripts that draw the manuscript figures
Output/Control/               DMSO vs E3 vehicle comparison (DE + GO:BP)
Output/PCA/                   control PCA plots (earlier versions)
Output/GSEA_NTC/              NTC gene-set collection and preranked GSEA results
Output/GSVA_NTC/              ssGSEA scores and limma results
Output/Figures/               manuscript figures, drawn by Scripts/figures/
```

### Data files

- **`ZF_RNA_Raw_Count.txt`**: 32,520 genes × 280 samples, tab-separated. Rows are Ensembl gene IDs (`ENSDARG…`). Values are RSEM expected counts, which are non-integer; scripts round them before DESeq2.
- **`norm_data_controls_ALL.txt`**: the 40 control samples, as counts divided by per-sample size factors. Input for `Pipeline_compare_controls.R`.
- **`metadata/20250520_ZF_RNAseq_metadata_JB.csv`**: one row per tube (315 rows). The 280 rows with `will_sequence = yes` match the columns of the count matrix. Key columns are `sample_id`, `compound_code`, `dose_code`, `dose_uM`, `timepoint_code`, `dev_stage_name`, `vehicle_type`, `replicate`, `experiment_date` and `control_pair_id`.
- **`metadata/norm_data_controls_{ALL,DMSO,E3}_metadata.tsv`**: metadata for the 40 control samples, all together and split by vehicle.
- **`DEGs_NTC/<CHEM>_<dose>_<stage>_vs_<vehicle>_C_<stage>_FDR_0.01_all_genes.txt`**: DESeq2 results for one treated group (2 samples) against the 4 controls of the same vehicle at the same stage. Columns are `ENSEMBL`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and `padj`. Each file lists every gene with an adjusted p-value (15,053–17,842 per contrast); `FDR_0.01` in the name is a label, not a filter. The Wald `stat` column is the ranking metric for GSEA.

There are 119 contrasts rather than 120. The library for sample `219-1305-T2-2-4D-S2` failed (almost no reads assigned), so the 2,4-D mid-dose (S2) contrast at T2 has no DESeq2 or GSEA output.

---

## Scripts

R scripts use paths relative to the repository root. Open `ZF-Temporal-RNAseq.Rproj` in RStudio, or `setwd()` to the repository root, before running them. The Python scripts find the repository root themselves.

| Step | Script | What it does | Writes to | Used in manuscript |
|---|---|---|---|---|
| 1 | `Scripts/Pipeline_compare_controls.R` | DMSO vs E3 vehicle effect per stage (DESeq2), GO:BP over-representation of vehicle-effect genes, control VST/PCA | `Output/Control/` | Vehicle-effect comparison |
| 2 | `Scripts/figures/Create_control_PCA_pastel.R` | PCA of the 40 control samples, coloured by stage (self-contained) | `Output/Figures/PCA_controls_pastel_notitle.pdf` | Control PCA figure |
| 3 | *DESeq2 contrast script (to be added)* <!-- TODO (Job): add the script that wrote DEGs_NTC/ --> | Each chemical–dose–stage group vs its vehicle control | `DEGs_NTC/` | Input for GSEA |
| 4 | `Scripts/GSEA_NTD.R` | Builds the NTC-focused GO:BP collection; preranked GSEA (clusterProfiler) on every `DEGs_NTC/` file | `Output/GSEA_NTC/` | GSEA results and heatmaps |
| 5 | `Scripts/Focused_GSVA_NTC_Chemical_perturbation.R` | ssGSEA (GSVA) scores per group and stage, including DMSO and E3 controls; limma tests of treated vs control scores | `Output/GSVA_NTC/` | GSVA results |
| 6 | `Scripts/figures/make_gsea_heatmaps.py` | GSEA heatmap: top 25 terms, plus a version with all FDR < 0.05 terms | `Output/Figures/` | Main-text GSEA heatmap (`THESIS_A5_top25_heatmap`) |
| 7 | `Scripts/figures/make_gsea_heatmaps_themed.py` | Same heatmaps, rows grouped by biological theme | `Output/Figures/` | Supplementary GSEA heatmap (`ARTICLE_heatmap_themed`) |
| 8 | `Scripts/figures/make_gsva_control_heatmaps.py` | Heatmaps of control ssGSEA scores across T1–T5 (DMSO and E3) | `Output/Figures/` | Supplementary GSVA heatmap |

Other scripts:

- **`Scripts/Create_control_PCAs.R`** and **`Scripts/Create_control_PCAs_poster.R`** draw earlier PCA versions (`Output/PCA/`). Run step 1 first, because they use the `vsd` object it creates.
- **`Scripts/Focused_GSVA_NTC_Controls.R`** is an earlier ssGSEA script that was not used for any figure. It groups samples by vehicle without restricting to controls, and its outputs are not part of this repository; the control GSVA results come from step 5. <!-- TODO (Job): confirm, or remove this script -->

### Run order

```bash
# from the repository root
Rscript Scripts/Pipeline_compare_controls.R
Rscript Scripts/figures/Create_control_PCA_pastel.R
# step 3: DESeq2 contrasts -> DEGs_NTC/ (script to be added)
Rscript Scripts/GSEA_NTD.R
Rscript Scripts/Focused_GSVA_NTC_Chemical_perturbation.R
python3 Scripts/figures/make_gsea_heatmaps.py
python3 Scripts/figures/make_gsea_heatmaps_themed.py
python3 Scripts/figures/make_gsva_control_heatmaps.py
```

Every output is already in `Output/`, so you can start at any step. GSEA uses `set.seed(42)`.

---

## Software

Analyses were run in **R 4.5.1** with DESeq2 1.48.1, clusterProfiler 4.16.0, enrichplot 1.28.4, GSVA 2.2.0, limma 3.64.3, org.Dr.eg.db 3.21.0 and GO.db 3.21.0.

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "DOSE", "GSVA",
                       "limma", "org.Dr.eg.db", "GO.db", "AnnotationDbi"))
install.packages(c("ggplot2", "ggrepel", "dplyr", "tidyr", "stringr", "readr", "tibble",
                   "purrr", "matrixStats", "pheatmap", "patchwork", "ggraph", "igraph",
                   "ggridges"))
```

The figure scripts need **Python 3** with `pandas`, `numpy`, `matplotlib` and `scipy`:

```bash
pip install pandas numpy matplotlib scipy
```

---

## Notes

- One metadata row has `compound_code` `CyCL` instead of `CYCL` (sample `73-0227-T1-CYCL-S1`). The GSVA script upper-cases compound codes, so this sample is analysed as cyclopamine.
- The metadata file contains `vehicle_type` twice; the first column holds the values.
- The metadata columns for RNA integrity, library preparation, lane and FASTQ file names are not filled in. <!-- TODO (Job): fill these from the USEQ QC report, or drop them -->

## License

Code is released under the MIT License (see `LICENSE`). <!-- TODO (Job): decide on a data license, e.g. CC BY 4.0 for the count matrix, metadata and results. -->

## Citation

If you use these data or scripts, please cite the paper above. <!-- TODO (Job): add the Zenodo DOI of the release that matches the published paper. -->

## Acknowledgements

Sequencing was performed by the Utrecht Sequencing Facility (USEQ). This work was carried out within the ONTOX project, funded by the European Union's Horizon 2020 Research and Innovation programme under grant agreement No 963845.

## Contact

Questions about the data or code: open an issue on this repository, or contact the corresponding author listed in the paper.
