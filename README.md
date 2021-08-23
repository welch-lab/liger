<img src="docs/img/liger_cropped.png" width="150">

[![Build Status](https://api.travis-ci.com/welch-lab/liger.svg?branch=master)](https://api.travis-ci.com/welch-lab/liger.svg?branch=master)

# LIGER (Linked Inference of Genomic Experimental Relationships)

LIGER (installed as `rliger` ) is a package for integrating and analyzing multiple single-cell datasets, developed by the Macosko lab and maintained/extended by the Welch lab. It relies on integrative non-negative matrix factorization to identify shared and dataset-specific factors. 

Check out our [Cell paper](https://www.cell.com/cell/fulltext/S0092-8674%2819%2930504-5) for a more complete description of the methods and analyses. To access data used in our SN and BNST analyses, visit our [study](https://portals.broadinstitute.org/single_cell/study/SCP466) on the
Single Cell Portal. 

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

* NEW: [Iterative Single-Cell Multi-Omic Integration Using Online iNMF](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/online_iNMF_tutorial.html)

* NEW: [Integrating unshared features with UINMF](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html)
* NEW: [scATAC and scRNA Integration using unshared features (UINMF)](http://htmlpreview.github.io/?https://www.dropbox.com/s/ujsypyfpebvdoye/STARmap_dropviz_vig.html?dl=0)
* NEW: [Cross-species Analysis with UINMF](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/cross_species_vig.html)
* NEW: [Performing Parameter Selection](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Parameter_selection.html) 

* [Integrating Multiple Single-Cell RNA-seq Datasets](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html)

* [Jointly Defining Cell Types from scRNA-seq and scATAC-seq](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)

* [Jointly Defining Cell Types from Single-Cell RNA-seq and DNA Methylation](https://welch-lab.github.io/liger/rna-methylation.html)

* [Running Liger directly on Seurat objects using Seurat wrappers](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html)

## System Requirements

### Hardware requirements
The `rliger` package requires only a standard computer with enough RAM to support the in-memory operations. For minimal performance, please make sure that the computer has at least about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

* RAM: 16+ GB
* CPU: 4+ cores, 2.3 GHz/core

### Software requirements

The package development version is tested on *Linux* operating systems and *Mac OSX*.

* Linux: CentOS 7, Manjaro 5.3.18
* Mac OSX: Mojave (10.14.1), Catalina (10.15.2)

The `rliger` package should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the `rliger` package, users should have R version 3.4.0 or higher, and several packages set up from CRAN and other repositories. The user can check the dependencies in `DESCRIPTION`.

## Installation

LIGER is written in R and is also available on the Comprehensive R Archive Network (CRAN). Note that the package name is `rliger` to avoid a naming conflict with an unrelated package. To install the version on CRAN, follow these instructions:

1. Install [R](https://www.r-project.org/)  (>= 3.4)
2. Install [Rstudio](https://rstudio.com/products/rstudio/download/) (recommended)
3. Type the following R command:
```
install.packages('rliger')
```
To install the latest development version directly from GitHub, type the following commands instead of step 3:
```
install.packages('devtools')
library(devtools)
install_github('welch-lab/liger')
```
Note that the GitHub version requires installing from source, which may involve additional installation steps on MacOS (see below).

### Additional Steps for Installing LIGER from Source (recommended before step 3)
Installation from CRAN is easy because pre-compiled binaries are available for Windows and MacOS. However, a few additional steps are required to install from source on MacOS/Windows (e.g. Install RcppArmadillo).
(MacOS) Installing RcppArmadillo on R>=3.4 requires Clang >= 4 and gfortran-6.1. For newer versions of R (R>=3.5), it's recommended to follow the instructions in this [post](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/). Follow the instructions below if you have R version 3.4.0-3.4.4.

1. Install gfortran as suggested [here](https://gcc.gnu.org/wiki/GFortranBinaries)
2. Download clang4 from this [page](https://mac.R-project.org/libs/clang-4.0.0-darwin15.6-Release.tar.gz)
3. Uncompress the resulting zip file and type into Terminal (`sudo` if needed): 
```
mv /path/to/clang4/ /usr/local/ 
```
4. Create `.R/Makevars` file containing following:
```
# The following statements are required to use the clang4 binary
CC=/usr/local/clang4/bin/clang
CXX=/usr/local/clang4/bin/clang++
CXX11=/usr/local/clang4/bin/clang++
CXX14=/usr/local/clang4/bin/clang++
CXX17=/usr/local/clang4/bin/clang++
CXX1X=/usr/local/clang4/bin/clang++
LDFLAGS=-L/usr/local/clang4/lib
```
For example, use the following Terminal commands:
```
cd ~
mkdir .R
cd .R 
nano Makevars
``` 
Paste in the required text above and save with `Ctrl-X`.

### Additional Installation Steps for Online Learning using LIGER
The HDF5 library is required for implementing online learning in LIGER on data files in HDF5 format. It can be installed via one of the following commands:

| System                                    | Command
|:------------------------------------------|:---------------------------------|
|**OS X (using Homebrew or Conda)**                  | `brew install hdf5` or `conda install -c anaconda hdf5`
|**Debian-based systems (including Ubuntu)**| `sudo apt-get install libhdf5-dev` 
|**Systems supporting yum and RPMs**        | `sudo yum install hdf5-devel`

For Windows, the latest HDF5 1.12.0 is available at https://www.hdfgroup.org/downloads/hdf5/.

### Detailed Instructions for FIt-SNE Installation
Note that the runUMAP function (which calls the `uwot` package) also scales to large datasets and does not require additional installation steps.
However, using FIt-SNE is recommended for computational efficiency if you want to perform t-SNE on very large datasets.
Installing and compiling the necessary software requires the use of git, FIt-SNE, and FFTW. For a 
basic overview of installation, visit this [page](https://github.com/KlugerLab/FIt-SNE).

Basic installation for most Unix machines can be achieved with the following commands after downloading 
the latest version of FFTW from [here](http://www.fftw.org/). In the fftw directory, run:
```
./configure
make
make install
```
(Additional [instructions](http://www.fftw.org/fftw3_doc/Installation-and-Customization.html) if 
necessary).
Then in desired directory:
```
git clone https://github.com/KlugerLab/FIt-SNE.git
cd FIt-SNE
g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm
pwd
```
Use the output of `pwd` as the `fitsne.path` parameter in runTSNE. 

Note that the above instructions require root access. To install into a specified folder (such as your home directory) on a server, use the `--prefix` option:
```
./configure --prefix=<install_dir>
make
make install
git clone https://github.com/KlugerLab/FIt-SNE.git
cd FIt-SNE
g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -I<install_dir>/include/ -L<install_dir>/lib/ -o bin/fast_tsne -pthread -lfftw3 -lm
pwd
```

### Install Time and Expected Run Time

The installation process of `rliger` should take less than 30 minutes.

The expected run time is 1 - 4 hours depending on dataset size and downstream analysis of the user’s choice.

## Sample Datasets
The `rliger` package provides a small sample dataset for basic demos of the functions. You can find it in folder `liger/tests/testdata/small_pbmc_data.RDS`.

We also provide a set of scRNA-seq and scATAC-seq datasets for real-world style demos. These datasets are as follows:

* scATAC and scRNA data provided by 10X Genomics, access the pre-processed data from [here](https://umich.app.box.com/s/5hkhpou4inrulo570yhc38ay8g3ont5b). The data sources are:
    + `pbmc.atac.expression.mat.RDS`: raw data can be accessed [here](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k), created by Cell Ranger ATAC 1.1.0;
    + `pbmc.rna.expression.mat.RDS`: raw data can be accessed [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3), created by Cell Ranger 3.0.0.
  
* scATAC and scRNA data provided by [GreenleafLab](https://github.com/GreenleafLab/MPAL-Single-Cell-2019); you can access the pre-processed data from [here](https://umich.box.com/s/wip2nzpktn6fdnlpc83o1u7anjn4ue2c):
     + `GSM4138872_scRNA_BMMC_D1T1.RDS`;
     + `GSM4138873_scRNA_BMMC_D1T2.RDS`;
     + `GSM4138888_scATAC_BMMC_D5T1_peak_counts.RDS`;
     + `GSM4138888_scATAC_BMMC_D5T1.RDS`.

* scRNA data composed of two datasets of interneurons and oligodendrocytes from the mouse frontal cortex, two distinct cell types that should not align if integrated. Provided by [Saunders, A. et.al., 2018](https://doi.org/10.1016/j.cell.2018.07.028); you can access the pre-processed data from [here](https://umich.box.com/s/n1xfpu9hplrknu6to6u9kvbktsqcql8t):
     + `interneurons_and_oligo.RDS`;
     
* scRNA data from control and interferon-stimulated PBMCs. Raw data provided by [Kang, et.al., 2017](https://www.nature.com/articles/nbt.4042); The datasets were downsampled by applying the sample function without replacement yield 3000 cells for each matrix. You can download downsampled data from [here](https://umich.box.com/s/n1xfpu9hplrknu6to6u9kvbktsqcql8t):
     + `PBMC_control.RDS`;
     + `PBMC_interferon-stimulated.RDS`.

Corresponding tutorials can be found in section **Usage** above.

## License
This project is covered under the **GNU General Public License 3.0**.
