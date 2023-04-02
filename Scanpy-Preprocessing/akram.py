#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 12:19:47 2020

@author: Akram

---------------------------------------

Akram's preprocessing script from Compute Canada server, in Akram>SingleCellDataAnalysis>preprocessing.ipynb
Login info is in Drive

"""
import numpy as np
# data preprocessing
import os
#  IntegratedBatches/scRNAseq_Benchmark_datasets/Inter-dataset/PbmcBench/10Xv2
os.chdir("D:\_Drive E\_MyThesis\Data")
import pandas as pd
import scanpy as sc

pbmc10Xv2 = pd.read_csv("./10Xv2_pbmc1.csv",index_col=0)

data1 = pbmc10Xv2

def dataSetAnalysis(df):
    #view starting values of data set
    print("Dataset Head")
    print(df.head(3))
    print("=" * 30)
    
    # View features in data set
    print("Dataset Features")
    print(df.columns.values)
    print("=" * 30)
    
    # View How many samples and how many missing values for each feature
    print("Dataset Features Details")
    print(df.info())
    print("=" * 30)
    
    # view distribution of numerical features across the data set
    print("Dataset Numerical Features")
    print(df.describe())
    print("=" * 30)
    
    # view distribution of categorical features across the data set
    # print("Dataset Categorical Features")
    # print(df.describe(include=['O']))
    # print("=" * 30)

dataSetAnalysis(data1)

sc.settings.set_figure_params(dpi=80, facecolor='white')
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

results_file = './results/10Xv2_pbmc1.h5ad'# the file that will store the analysis results

cols = data1.columns
genelist = cols

numberOfGenes = len(genelist)
labels = pd.read_csv("./10Xv2_pbmc1Labels.csv")
dataSetAnalysis(labels)

# Join col(gene) names to X values
# rowNames = adata.obs_names
cellsIndex  = pd.DataFrame(data1.index)
cells = data1.iloc[:,0]
cells.shape
#type(cellsIndex)

labelWithIndex = pd.concat([cellsIndex,labels], axis=1)
labelWithIndex.to_csv("./results/labelWithIndex.csv", index = False)
labelWithCells = pd.read_csv("./results/labelWithIndex.csv",index_col=0)
concatenatedData = pd.concat([data1, labelWithCells], axis=1)
dataSetAnalysis(concatenatedData)
df = concatenatedData
dfNotencode = df
# Notice that 'assigned_cluster' contains string names to represent cell types. Let's encode them to integer
print("Before encoding: ")
len(np.unique(dfNotencode.iloc[:,-1]))
cellTypes = dfNotencode.iloc[:,-1]
cellTypes.to_csv("./results/cellTypesPBMC10Xv2.csv")
df.cellType = pd.factorize(df.cellType)[0]

np.unique(df.cellType)

print("/nAfter encoding: ")
print(df.iloc[0:2000,-1])#example number of rows

#########################
df.to_csv("./results/dataPBMC10Xv2Factorized.csv")

# Scanpy ********************************************************************
adata = sc.read_csv("./results/dataPBMC10Xv2Factorized.csv", first_column_names=True)
data =adata
type(data)
adata.obs_names 
adata.var_names
#adata.X
# adata.write_csvs("./anndata2.csv",skip_data=False, sep=',')
# Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
adata # n_obs × n_vars = 1937 × 20126

adata.obs_names_make_unique()
adata.var_names_make_unique()
# normalizing counts per cell
sc.pl.highest_expr_genes(adata, n_top=20, )
# # Basic filtering.
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# # With pp.calculate_qc_metrics, we can compute many metrics very efficiently.
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# # A violin plot of some of the computed quality measures:

# # the number of genes expressed in the count matrix
# # the total counts per cell
# # the percentage of counts in mitochondrial genes

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# # Remove cells that have too many mitochondrial genes expressed or too many total counts.
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# # Actually do the filtering by slicing the AnnData object.
adata = adata[adata.obs.n_genes_by_counts < 4000, :]
adata.var_names
adata = adata[adata.obs.pct_counts_mt < 10, :]


labelsBeforeNorm = adata.X[:,-1]
max(labelsBeforeNorm)
min(labelsBeforeNorm)
sc.pp.normalize_total(adata, target_sum=2e4)#['X'[:,1:]] #*********************If choosing target_sum=1e6, this is CPM normalization.
adata.var_names
max(adata.X[:,-1])
# Logarithmize the data.
sc.pp.log1p(adata)#**********************
adata.var_names
max(adata.X[:,-1])
adata.X[:,-1] = labelsBeforeNorm
max(adata.X[:,-1])
sc.pl.highest_expr_genes(adata, n_top=20, )
c = adata
# *********************************************************
# Identify highly-variable genes. (feature selection) 
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=.5)
sc.pl.highly_variable_genes(adata)
# sc.pl.highly_variable_genes(adata[:,0:-2])
max(adata.X[:,-1])
# Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene expression 
# for later use in differential testing and visualizations of gene expression. 
# This simply freezes the state of the AnnData object.
adata.var_names
adata.raw = adata

adataWithCluster = adata
adata.var.highly_variable.cellType = True

adata = adata[:, adata.var.highly_variable]


# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
# Regress out (mostly) unwanted sources of variation.
# a simple batch correction method is available via pp.regress_out()
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# n_obs × n_vars = 17654 × 21678
adata.var_names
max(adata.X[:,-1])
adata.X[:,-1] = labelsBeforeNorm
max(adata.X[:,-1])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)
adata.var_names
max(adata.X[:,-1])
adata.X[:,-1] = labelsBeforeNorm
max(adata.X[:,-1])

#sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca(adata, color='CST3')
# adata = adataWithCluster
adata.var_names
# Join col(gene) names to X values
rowNames = adata.obs_names

type(rowNames)
rows = pd.DataFrame(index=adata.obs_names)
type(rows)
np.savetxt("./results/barcodespbmc10xV2.csv", rows)
rows
hvgList = adata.var_names
n_hvg = len(adata.var_names)
hvg = pd.DataFrame(index=adata.var_names)
hvgString = ','.join(hvgList)
# Save adata as a csv file with row and col names
np.savetxt("./results/datapbmc10xV2NormWithColNames.csv", adata.X, delimiter=",", header = hvgString, comments='')

dataNorm = pd.read_csv('./results/datapbmc10xV2NormWithColNames.csv')
dataSetAnalysis(dataNorm)
dataNorm.set_index(rowNames, inplace=True)
dataSetAnalysis(dataNorm)


dataNorm.to_csv("./results/datapbmc10xV2NormWithColNamesRowNames-finalPreprocessedData.csv")
dataNormRows = pd.read_csv('./results/datapbmc10xV2NormWithColNamesRowNames-finalPreprocessedData.csv')
dataNormRows
dataSetAnalysis(dataNormRows)
np.unique(dataNormRows.cellType)
#########################
dataNormRows.isnull().values.any()


