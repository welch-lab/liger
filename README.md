# Analogizer

Analogizer is a package for integrating and analyzing multiple single-cell datasets, developed and maintained by the Macosko lab. It relies on integrative non-negative matrix factorization to identify shared and dataset-specific factors. 

## Installation

Analogizer is written in R and has a few other system requirements (Java) and recommended packages (umap in Python). To install the most recent development version, follow these instructions:

1. Install [R](https://www.r-project.org/)  (>= 3.4)
2. Install [Rstudio](https://www.rstudio.com/products/rstudio/download/) (recommended)
3. Make sure you have Java installed in your machine. Check by typing `java -version` into Terminal or CommandPrompt. 
4. Generate an auth [token](https://github.com/settings/tokens) for the Analogizer repo, making sure to include all repo permissions. 
5. Use the following R commands.
```
install.packages('devtools')
library(devtools)
install_github('MacoskoLab/Analogizer', auth_token = '<token>')
```

### Troubleshooting (MacOS only -- recommended before step 5)
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

## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 