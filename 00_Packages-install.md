Installation de tous les packages
================
Vincent Noah
19 décembre 2020

  - [Mise à jour de la machine
    virtuelle.](#mise-à-jour-de-la-machine-virtuelle.)
  - [Installation de Dada2 et de Phyloseq avec
    BiocManager.](#installation-de-dada2-et-de-phyloseq-avec-biocmanager.)
      - [Installation de BiocManager.](#installation-de-biocmanager.)
      - [Installation de Dada2.](#installation-de-dada2.)
      - [Installation de Phyloseq.](#installation-de-phyloseq.)
  - [Installation de packages pour
    phyloseq.](#installation-de-packages-pour-phyloseq.)
      - [Installation de phangorn.](#installation-de-phangorn.)
      - [Installation de DECIPHER.](#installation-de-decipher.)
      - [Installation de gridExtra.](#installation-de-gridextra.)
  - [Installation de packages pour des analyses complémentaires sur
    phyloseq.](#installation-de-packages-pour-des-analyses-complémentaires-sur-phyloseq.)

# Mise à jour de la machine virtuelle.

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
```

# Installation de Dada2 et de Phyloseq avec BiocManager.

A partir des instructions
<https://benjjneb.github.io/dada2/dada-installation.html>

## Installation de BiocManager.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.11')
```

## Installation de Dada2.

``` r
BiocManager::install("dada2", version = "3.11")
```

## Installation de Phyloseq.

``` r
BiocManager::install("phyloseq", version = "3.11")
```

# Installation de packages pour phyloseq.

## Installation de phangorn.

``` r
BiocManager::install("phangorn")
```

## Installation de DECIPHER.

``` r
BiocManager::install("DECIPHER")
```

## Installation de gridExtra.

``` r
install.packages("gridExtra")
```

# Installation de packages pour des analyses complémentaires sur phyloseq.

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
```

``` r
install.packages(.cran_packages)
devtools::install_github(.github_packages)
BiocManager::install(.bioc_packages)                        
```
