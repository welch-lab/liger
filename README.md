<img src="docs/img/liger_cropped.png" width="150">

[![Build Status](https://travis-ci.org/MacoskoLab/liger.svg?branch=master)](https://travis-ci.org/MacoskoLab/liger.svg?branch=master)

# LIGER (Linked Inference of Genomic Experimental Relationships)

LIGER (`liger`) is a package for integrating and analyzing multiple single-cell datasets, developed and maintained by the Macosko lab. It relies on integrative non-negative matrix factorization to identify shared and dataset-specific factors. 

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

## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Basic commands tutorial](http://htmlpreview.github.io/?https://github.com/MacoskoLab/liger/blob/master/vignettes/liger-vignette.html)

* [Integrating Multiple Single-Cell RNA-seq Datasets](https://macoskolab.github.io/liger/walkthrough_pbmc.html)

* [Jointly Defining Cell Types from scRNA-seq and snATAC-seq](https://macoskolab.github.io/liger/walkthrough_rna_atac.html)

* [Jointly Defining Cell Types from Single-Cell RNA-seq and DNA Methylation](https://macoskolab.github.io/liger/rna-methylation.html)

* [Running Liger directly on Seurat objects using Seurat wrappers](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html)

## System Requirements

### Hardware requirements
The `liger` package requires only a standard computer with enough RAM to support the in-memory operations. For minimal performance, please make sure that the computer has at least about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

* RAM: 16+ GB
* CPU: 4+ cores, 2.3 GHz/core

### Software requirements

The package development version is tested on *Linux* operating systems and *Mac OSX*.

* Linux: CentOS 7, Manjaro 5.3.18
* Mac OSX: Mojave (10.14.1), Catalina (10.15.2)

The `liger` package should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the `liger` package, users should have R version 3.4.0 or higher, and several packages set up from CRAN and other repositories. The user can check the dependencies in `DESCRIPTION`.

## Installation

`liger` is written in R and has a few other system requirements (Java) and recommended packages (umap in Python). To install the most recent development version, follow these instructions:

1. Install [R](https://www.r-project.org/)  (>= 3.4)
2. Install [Rstudio](https://www.rstudio.com/products/rstudio/download/) (recommended)
3. Make sure you have Java installed in your machine. Check by typing `java -version` into Terminal or Command Prompt. 
4. Use the following R commands.
```
install.packages('devtools')
library(devtools)
install_github('MacoskoLab/liger')
```

### Additional Installation Steps for MacOS (recommended before step 4)
Installing RcppArmadillo on R>=3.4 requires Clang >= 4 and gfortran-6.1. Follow the instructions below if you have R version 3.4.0-3.4.4. These instructions (using clang4) may also be sufficient for R>=3.5 but for newer versions of R, it's recommended to follow the instructions in this [post](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/). 

1. Install gfortran as suggested [here](https://gcc.gnu.org/wiki/GFortranBinaries)
2. Download clang4 from this [page](http://r.research.att.com/libs/clang-4.0.0-darwin15.6-Release.tar.gz)
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

### Running `liger` with Docker
If installing natively is difficult, you can run `liger` through our Docker image (available
publically), which also comes with Rstudio and Seurat (v2) installed.

1. Install [Docker](https://docs.docker.com/install/). 
2. Run the following in terminal:
```
docker run -d -p 8787:8787 docker.io/vkozareva/sc-liger:latest
```
3. Type [http://localhost:8787](http://localhost:8787) in any browser and enter "rstudio" as the 
username and password when prompted. `liger` and all of its dependencies are already installed in 
this environment.

If you wish to access local files in this container (mounting to `/data`) modify the command as follows:
```
docker run -d -v /path/to/local/directory:/data -p 8787:8787 docker.io/vkozareva/sc-liger:latest
```
Note that you will have to stop the container if you wish to allocate port `8787` to another application
later on. Further Docker documentation can be found [here](https://docs.docker.com/get-started/).

### Detailed Instructions for FIt-SNE Installation for use in runTSNE (recommended for large datasets)
Using FIt-SNE is recommended for computational efficiency when using runTSNE on very large datasets.
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

The installation process of `liger` should take less than 30 minutes.

The expected run time is 1 - 4 hours depending on dataset size and downstream analysis of the user’s choice.

## Sample Datasets
The `liger` package provides a small simulated dataset for basic demos of the functions, you can find it in folder `liger/tests/testdata/small_pbmc_data.RDS`.

We also provide a set of scRNA-seq and scATAC-seq datasets for real-world style demos. These datasets are as follows:

* scATAC and scRNA data provided by 10X Genomics, access the pre-processed data from [here](https://umich.app.box.com/s/5hkhpou4inrulo570yhc38ay8g3ont5b). The data sources are:
    + `pbmc.atac.expression.mat.RDS`: raw data can be accessed [here](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k), created by Cell Ranger ATAC 1.1.0;
    + `pbmc.rna.expression.mat.RDS`: raw data can be accessed [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3), created by Cell Ranger 3.0.0.
  
* scATAC and scRNA data provided by [GreenleafLab](https://github.com/GreenleafLab/MPAL-Single-Cell-2019); you can access the pre-processed data from [here](https://umich.box.com/s/wip2nzpktn6fdnlpc83o1u7anjn4ue2c):
     + `GSM4138872_scRNA_BMMC_D1T1.RDS`;
     + `GSM4138873_scRNA_BMMC_D1T2.RDS`;
     + `GSM4138888_atac_gene_counts_D5T1.RDS`;
     + `GSM4138888_atac_promoter_counts_D5T1.RDS`.

* scRNA data composed of two datasets of interneurons and oligodendrocytes from the mouse frontal cortex, two distinct cell types that should not align if integrated. Provided by [Saunders, A. et.al., 2018](https://doi.org/10.1016/j.cell.2018.07.028); you can access the pre-processed data from [here](https://umich.box.com/s/n1xfpu9hplrknu6to6u9kvbktsqcql8t):
     + `interneurons_and_oligo.RDS`;
     
* scRNA data from control and interferon-stimulated PBMCs. Raw data provided by [Kang, et.al., 2017](https://www.nature.com/articles/nbt.4042); The datasets were downsampled by applying the sample function without replacement yield 3000 cells for each matrix. You can download downsampled data from [here](https://umich.box.com/s/n1xfpu9hplrknu6to6u9kvbktsqcql8t):
     + `PBMC_control.RDS`;
     + `PBMC_interferon-stimulated.RDS`.

Corresponding tutorials can be found in section **Usage** above.

## License
This project is covered under the **GNU General Public License 3.0**.
