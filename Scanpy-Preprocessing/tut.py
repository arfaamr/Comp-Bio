# -*- coding: utf-8 -*-

"""
SCANPY Preprocessing
Arfaa Rashid
February 19, 2023

Following tutorial to learn to preprocess data using Python

TUTORIAL: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
INPUT DATA: http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

"""


#----------------------------------------
#IMPORT PACKAGES

# installed scanpy with: conda install -c conda-forge scanpy python-igraph leidenalg [https://scanpy.readthedocs.io/en/stable/installation.html]

import numpy as np
import pandas as pd
import scanpy as sc

#Setup scanpy
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


#----------------------------------------
#READ INPUT

#Read in count matrix [input data file linked in header] into AnnData object
# obj x var: 2700 x 32738
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

adata.var_names_make_unique()               #appends number to duplicates

results_file = 'write/pbmc3k.h5ad'  # file that will store analysis results


#----------------------------------------
#PREPROCESSING

#Plots boxplot showing the 20 genes with highest fraction of counts in each single cell across all cells
sc.pl.highest_expr_genes(adata, n_top=20, )

#Basic filtering - removing outliers? ** in terms of frequency
# if cell does not express at least 200 genes, or a gene is not expressed in at least 3 cells, they are removed?
# [.pp is preprocessing]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)



#Use biological understanding of statistics indicating poor quality cells to perform quality control (QC) 

#Annotate mitochondrial cells as "mt" and calculate QC metrics avout those cells - inPlace, so adata itself is updated with the info the function returns
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#Plot violinplot, scatterplots of the QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

#Remove cells with too many mt genes expressed or too many total counts
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]



#Normalize to 10000 reads per cell so counts can be compared between cells
sc.pp.normalize_total(adata, target_sum=1e4)

#Logarithmize
sc.pp.log1p(adata)

#Identify highly-variable genes - wide range of gene expression
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)



#set the .raw to the normalized, logarithmized version of adata?
adata.raw = adata



#Actually filter adata. so what were we doing before? ** 
adata = adata[:, adata.var.highly_variable]
#Regress out effects of total counts per cell and percentage of mitochondrial genes expressed
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#Scale to unit variance, discard outliers
sc.pp.scale(adata, max_value=10)


#----------------------------------------
#PRINCIPAL COMPONENT ANALYSIS
# dimensionality reduction; reveals the main axes of variation and denoises data

#Perform PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')

#Inspect contribution of single PCs to variance to determine how many PCs to consider when computing neighborhood relations of cells ?
sc.pl.pca_variance_ratio(adata, log=True)

#Write results
adata.write(results_file)


#----------------------------------------
#NEIGHBOURHOOD GRAPH







#----------------------------------------
"""
how to "look at" a variable, like in RStudio? like actually whats in the dataframe, etc

what does it mean by "counts"? 

what does .raw do?  "freezes the state of the AnnData object." ? may not be necessary?
regress out? unit variance?

what is "neighbourhood"? - review how PCA actually works 

"""


