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


