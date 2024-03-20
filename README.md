<img src="https://github.com/welch-lab/liger/raw/newObj/inst/extdata/logo.png" width="120" style="display: inline;">

<a href="https://github.com/welch-lab/liger/actions/workflows/r.yml"><img src="https://github.com/welch-lab/liger/actions/workflows/r.yml/badge.svg?branch=newObj" alt="R" style="display: inline;"></a>
<a href="https://app.codecov.io/gh/mvfki/liger"><img src="https://codecov.io/gh/mvfki/liger/graph/badge.svg?token=77TTU4GY8" alt="codecov" style="display: inline;"></a>
<a href="https://cran.r-project.org/package=rliger"><img src="https://cranlogs.r-pkg.org/badges/rliger" alt="cran" style="display: inline;"></a>

# LIGER (Linked Inference of Genomic Experimental Relationships)

>**Now we have a comprehensive documentation site for the latest version of [rliger (2.0)](https://welch-lab.github.io/liger/index.html)!**

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

If you have any questions, comments, or suggestions, you are welcomed to [open an Issue](https://github.com/welch-lab/liger/issues)!

## Usage

For usage examples and guided walkthroughs, check the `vignettes` directory of the repo.

* [Integrating Multiple Single-Cell RNA-seq Datasets](https://welch-lab.github.io/liger/articles/Integrating_multi_scRNA_data.html)
* [Jointly Defining Cell Types from scRNA-seq and scATAC-seq](https://welch-lab.github.io/liger/articles/Integrating_scRNA_and_scATAC_data.html)
* [Iterative Single-Cell Multi-Omic Integration Using Online iNMF](https://welch-lab.github.io/liger/articles/online_iNMF_tutorial.html)
* [Integrating unshared features with UINMF](https://welch-lab.github.io/liger/articles/UINMF_vignette.html)
* [Integrating spatial transcriptomic and transcriptomic datasets using UINMF](https://welch-lab.github.io/liger/articles/STARmap_dropviz_vig.html)
* [scATAC and scRNA Integration using unshared features (UINMF)](https://welch-lab.github.io/liger/articles/SNAREseq_walkthrough.html)
* [Cross-species Analysis with UINMF](https://welch-lab.github.io/liger/articles/cross_species_vig.html)
* [Jointly Defining Cell Types from Single-Cell RNA-seq and DNA Methylation](https://welch-lab.github.io/liger/articles/rna_methylation.html)

Meanwhile, since version 2.0.0, LIGER is massively updated for usability and interoperability with other packages. Below are links to the introduction of new features.

* [Introduction to new liger object and other related classes](https://welch-lab.github.io/liger/articles/liger_object.html)
* [Running Liger directly on Seurat objects](https://welch-lab.github.io/liger/articles/liger_with_seurat.html)

If you need to refer to the tutorials for the old version of rliger, please check the [GitHub archive v1.0.1](https://github.com/welch-lab/liger/tree/v1.0.1/vignettes), download the desired rendered HTML files and open them in your browser.

## Sample Datasets

The `rliger` package provides different types of small toy dataset for basic demos of the functions. After attaching the package in an R session, you can load them with:

```R
data("pbmc")
data("pbmcPlot")
data("bmmc")
```

We also provide a set of datasets for real-world style demos, including scRNAseq, scATACseq, spatial transcriptomics and DNA methylation data.
They are described in detail in the articles that make use of them. Please check them out from the links above.
