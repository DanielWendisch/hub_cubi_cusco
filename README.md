# hub_cubi_cusco

Single-cell RNA-seq analysis of multiplexed organoid datasets (`hub_01`, `hub_02`, `hub_10`, `hub_11`): quality control, cell-line and doublet assignment, integration, and bulk deconvolution case studies. Written in R (Seurat, Quarto) with a Snakemake entry point.

![workflow](workflow.PNG)

> **Status: research work in progress.** This is an analysis log rather than a packaged tool. Raw data and rendered reports are not included, so the code will not run end-to-end without the original data.

## Layout

| Path | Contents |
|---|---|
| `scripts/` | Main pipeline: dataset setup (`QC_dataset_setup_hub_*.R`, `1_QC_dataset_setup.R`), QC dashboard (`QC_dashboard.qmd`), per-dataset processing and master notebooks (`hub_01_master.qmd`, `hub_10_master.qmd`), export to cellxgene / AnnData |
| `case_studies/` | Exploratory notebooks: BayesPrism deconvolution, distance measures, organoid differentiation assessment, CellBender output checks |
| `side_scripts/` | Small one-off analyses (e.g. IL1RAP) and AnnData conversion |
| `integration_by_cell_type.qmd`, `read_in_03_04_05_05.qmd`, `Cibersortx prep.qmd` | Integration and deconvolution preparation |
| `Snakefile`, `envs/` | Snakemake entry point (experimental) and conda environment for the QC dashboard |
| `renv.lock` | R package versions |
| `docs/metadata_notes.md` | Notes on how the cell metadata columns are derived (cellranger / CellBender / vireo) |

## Expected directory layout

Scripts use paths relative to the project root, with raw data in a sibling folder two levels up:

```
<parent>/
├── raw_data/cubi_cusco/        # cellranger multi + CellBender outputs (not included)
└── <projects>/hub_cubi_cusco/  # this repository
```

Generated files go to `output/` and `intermediate_data/`, which are git-ignored.

## Setup

```r
renv::restore()   # R packages from renv.lock
```

Quarto is needed to render the `.qmd` notebooks.

## Known issues

- The `Snakefile` is unfinished and does not run as is.
- Some notebooks were written interactively on Windows and may need path adjustments.

## License

See [LICENSE](LICENSE).
