<img src="https://github.com/welch-lab/liger/raw/newObj/inst/extdata/logo.png" width="120" style="display: inline;">

<a href="https://github.com/welch-lab/liger/actions/workflows/r.yml"><img src="https://github.com/welch-lab/liger/actions/workflows/r.yml/badge.svg?branch=master" alt="R" style="display: inline;"></a>
<a href="https://app.codecov.io/gh/welch-lab/liger"><img src="https://codecov.io/gh/welch-lab/liger/graph/badge.svg?token=chxwVaVsGp" alt="codecov" style="display: inline;"></a>
<a href="https://cran.r-project.org/package=rliger"><img src="https://cranlogs.r-pkg.org/badges/rliger" alt="cran" style="display: inline;"></a>

# LIGER (Linked Inference of Genomic Experimental Relationships)

<div style="background: #dddddd;">

>**NEWS** Oct., 2024
>
>- Checkout new cell factor alignment method (function [`centroidAlign()`](https://welch-lab.github.io/liger/reference/centroidAlign.html)), which aligns cell factor loading by moving soft clustering centroids. Its overall performance, in terms of batch effect removal and especially biological information conservation, out performs many public well-known methods. [**See benchmarking article here**](https://welch-lab.github.io/liger/articles/benchmark.html).
>- Checkout Consensus iNMF method (function [`runCINMF()`](https://welch-lab.github.io/liger/reference/runCINMF.html)), which runs regular iNMF multiple times with different random initialization and summarizes a consensus result with better confidence.
>- Please visit [*rliger* website](https://welch-lab.github.io/liger/index.html) for comprehensive documentation and [revised tutorial](https://welch-lab.github.io/liger/articles/Integrating_multi_scRNA_data.html) that walks through scRNAseq integration and analysis in detail
>- More [changelogs](https://welch-lab.github.io/liger/news/index.html)

</div>

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

## Citation

If you use LIGER in your research please cite our paper correspondingly:

* Generally the *Cell* paper should be cited:

>Joshua D. Welch and et al., Single-Cell Multi-omic Integration Compares and Contrasts Features of Brain Cell Identity, Cell, VOLUME 177, ISSUE 7, P1873-1887.E17 (2019), [https://doi.org/10.1016/j.cell.2019.05.006](https://doi.org/10.1016/j.cell.2019.05.006)

* For the *rliger* package:

>Liu, J., Gao, C., Sodicoff, J. et al. Jointly defining cell types from multiple single-cell datasets using LIGER. Nat Protoc 15, 3632–3662 (2020), [https://doi.org/10.1038/s41596-020-0391-8](https://doi.org/10.1038/s41596-020-0391-8)

* For online iNMF integration method:

>Gao, C., Liu, J., Kriebel, A.R. et al. Iterative single-cell multi-omic integration using online learning. Nat Biotechnol 39, 1000–1007 (2021), [https://doi.org/10.1038/s41587-021-00867-x](https://doi.org/10.1038/s41587-021-00867-x)

* For UINMF integration method:

>Kriebel, A.R., Welch, J.D. UINMF performs mosaic integration of single-cell multi-omic datasets using nonnegative matrix factorization. Nat Commun 13, 780 (2022), [https://doi.org/10.1038/s41467-022-28431-4](https://doi.org/10.1038/s41467-022-28431-4)

## Feedback

If you have any questions, comments, or suggestions, you are welcomed to [open an Issue](https://github.com/welch-lab/liger/issues)!

## Usage

For usage examples and guided walkthroughs of specific use cases, please check our articles below:

* [Integrating Multiple Single-Cell RNA-seq Datasets](https://welch-lab.github.io/liger/articles/Integrating_multi_scRNA_data.html)
* [Jointly Defining Cell Types from scRNA-seq and scATAC-seq](https://welch-lab.github.io/liger/articles/Integrating_scRNA_and_scATAC_data.html)
* [Iterative Single-Cell Multi-Omic Integration Using Online iNMF](https://welch-lab.github.io/liger/articles/online_iNMF_tutorial.html)
* [Integrating datasets using unshared features with UINMF](https://welch-lab.github.io/liger/articles/UINMF_vignette.html)
* [Integrating spatial transcriptomic and transcriptomic datasets using UINMF](https://welch-lab.github.io/liger/articles/STARmap_dropviz_vig.html)
* [Integrating scATAC and scRNA using unshared features (UINMF)](https://welch-lab.github.io/liger/articles/SNAREseq_walkthrough.html)
* [Cross-species Analysis with UINMF](https://welch-lab.github.io/liger/articles/cross_species_vig.html)
* [Jointly Defining Cell Types from Single-Cell RNA-seq and DNA Methylation](https://welch-lab.github.io/liger/articles/rna_methylation.html)

Meanwhile, since version 2.0.0, LIGER is massively updated for usability and interoperability with other packages. Below are links to the introduction of new features.

* [Introduction to new liger object and other related classes](https://welch-lab.github.io/liger/articles/liger_object.html)
* [Running LIGER directly on Seurat objects](https://welch-lab.github.io/liger/articles/liger_with_seurat.html)

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
