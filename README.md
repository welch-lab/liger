<img src="https://github.com/mvfki/liger/raw/newObj/inst/extdata/logo.png" width="120">

[![R](https://github.com/mvfki/liger/actions/workflows/r.yml/badge.svg?branch=newObj)](https://github.com/mvfki/liger/actions/workflows/r.yml)[![codecov](https://codecov.io/gh/mvfki/liger/branch/newObj/graph/badge.svg?token=77TTU4GY8A)](https://codecov.io/gh/mvfki/liger)

# LIGER (Linked Inference of Genomic Experimental Relationships)

>Now we have a comprehensive documentation site for the latest development version of [rliger 2.0](https://mvfki.github.io/liger/index.html)!

LIGER (installed as `rliger` ) is a package for integrating and analyzing multiple single-cell datasets, developed by the Macosko lab and maintained/extended by the Welch lab. It relies on integrative non-negative matrix factorization to identify shared and dataset-specific factors.

Check out our [Cell paper](https://doi.org/10.1016/j.cell.2019.05.006) for a more complete description of the methods and analyses. To access data used in our SN and BNST analyses, visit our study "SCP466" on the
[Single Cell Portal](https://singlecell.broadinstitute.org/single_cell). 

LIGER can be used to compare and contrast experimental datasets in a variety of contexts, for instance:

* Across experimental batches
* Across individuals
* Across sex
* Across tissues
* Across species (e.g., mouse and human)
* Across modalities (e.g., scRNAseq and spatial transcriptomics data, scMethylation, or scATAC-seq)

Once multiple datasets are integrated, the package provides functionality for further data exploration,
analysis, and visualization. Users can:

* Identify clusters
* Find significant shared (and dataset-specific) gene markers
* Compare clusters with previously identified cell types
* Visualize clusters and gene expression using t-SNE and UMAP

We have also designed LIGER to interface with existing single-cell analysis packages, including
[Seurat](https://satijalab.org/seurat/).

## Feedback

Consider filling out our [feedback form](https://forms.gle/bhvp3K6tiHwf976r8) to help us improve the functionality and accessibility of LIGER.

## Usage

For usage examples and guided walkthroughs, check the `vignettes` directory of the repo.

* [Iterative Single-Cell Multi-Omic Integration Using Online iNMF](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/online_iNMF_tutorial.html)
* [Integrating unshared features with UINMF](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html)
* [scATAC and scRNA Integration using unshared features (UINMF)](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/SNAREseq_walkthrough.html)
* [Cross-species Analysis with UINMF](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/cross_species_vig.html)
* [Performing Parameter Selection](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Parameter_selection.html)
* [Integrating spatial transcriptomic and transcriptomic datasets using UINMF](https://www.dropbox.com/s/ho62az1kshj59d2/STARmap_dropviz_vig.html?dl=0) (Click to Download)
* [Integrating Multiple Single-Cell RNA-seq Datasets](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html)

* [Jointly Defining Cell Types from scRNA-seq and scATAC-seq](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)

* [Jointly Defining Cell Types from Single-Cell RNA-seq and DNA Methylation](https://welch-lab.github.io/liger/rna-methylation.html)

* [Running Liger directly on Seurat objects using Seurat wrappers](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html)


## Sample Datasets

The `rliger` package provides a small sample dataset for basic demos of the functions. You can find it in folder `liger/tests/testdata/small_pbmc_data.RDS`.

We also provide a set of scRNA-seq and scATAC-seq datasets for real-world style demos. These datasets are as follows:

* scRNA data from control and interferon-stimulated PBMCs. Raw data provided by [Kang, et.al., 2017](https://www.nature.com/articles/nbt.4042); The datasets were downsampled by applying the sample function without replacement yield 3000 cells for each matrix. You can download downsampled data from [here](https://www.dropbox.com/sh/u94ib3dkf9pb6nd/AABemvnxDgKDGRs8Ek5QGlXWa?dl=0):
     + `PBMC_control.RDS`;
     + `PBMC_interferon-stimulated.RDS`.

* scRNA data composed of two datasets of interneurons and oligodendrocytes from the mouse frontal cortex, two distinct cell types that should not align if integrated. Provided by [Saunders, A. et.al., 2018](https://doi.org/10.1016/j.cell.2018.07.028); you can access the pre-processed data from [here](https://www.dropbox.com/sh/u94ib3dkf9pb6nd/AABemvnxDgKDGRs8Ek5QGlXWa?dl=0):
     + `interneurons_and_oligo.RDS`;

* scATAC and scRNA data provided by [GreenleafLab](https://github.com/GreenleafLab/MPAL-Single-Cell-2019); you can access the pre-processed data from [here](https://www.dropbox.com/sh/5e9cy4qabs89kc8/AADy3jDxHx94j6A57t7A42u3a?dl=0):
     + `GSM4138872_scRNA_BMMC_D1T1.RDS`;
     + `GSM4138873_scRNA_BMMC_D1T2.RDS`;
     + `GSM4138888_scATAC_BMMC_D5T1_peak_counts.RDS`;
     + `GSM4138888_scATAC_BMMC_D5T1.RDS`.

<!-- * scATAC and scRNA data provided by 10X Genomics, access the pre-processed data from [here](https://umich.app.box.com/s/5hkhpou4inrulo570yhc38ay8g3ont5b). The data sources are: -->
* scATAC and scRNA data provided by 10X Genomics. The data sources are:
    + `pbmc.atac.expression.mat.RDS`: raw data can be accessed [here](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k), created by Cell Ranger ATAC 1.1.0;
    + `pbmc.rna.expression.mat.RDS`: raw data can be accessed [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3), created by Cell Ranger 3.0.0.
