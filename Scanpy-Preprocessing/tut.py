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
#DATASET ANALYSIS
# from Akram's script - view starting values, features of set

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


df = adata.to_df()
dataSetAnalysis(df)


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

#Annotate mitochondrial cells as "mt" and calculate QC metrics about those cells - inPlace, so adata itself is updated with the info the function returns
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#Plot violinplot/scatterplots of the QC metrics
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



#set the .raw to the normalized, logarithmized version of adata? to use later?
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
#Plot scatterplot of PCA results.?
sc.pl.pca(adata, color='CST3')

#Inspect contribution of single PCs to variance to determine how many PCs to consider when computing neighborhood relations of cells ?
sc.pl.pca_variance_ratio(adata, log=True)

#Write results
adata.write(results_file)


#----------------------------------------
#NEIGHBOURHOOD GRAPH

#Compute neighbourhood
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# use UMAP instead of tSNE. if there are still connectivity violoations, run following code
#  sc.tl.paga(adata)
#  sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#  sc.tl.umap(adata, init_pos='paga')

#Compute UMAP
sc.tl.umap(adata)

#Plot using raw [normalized & logmarithmized but not corrected] data and then with scaled & corrected data
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)


#----------------------------------------
#CLUSTERING NEIGHBOURHOOD
# use Leiden graph clustering method - clusters neighbourhood graph which we just produced

sc.tl.leiden(adata
             )
#Plot and write results
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
adata.write(results_file)


#----------------------------------------
#FIND MARKER GENES

#Compute ranking for differential genes in each cluster using t-test
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

adata.write(results_file)

#Alternatively, rank using multivariate logistic regression
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#Define marker genes.? are these the ones we found?
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']


adata = sc.read(results_file)
#Show top 10 ranked genes per cluster 0, 1, …, 7 in dataframe
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
#Get table with scores and groups
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

#Compare to a single cluster
adata.uns['log1p']["base"] = None       #***following line was giving keyerror: base, before this was added
sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)


#Reload with computed differential expression
adata = sc.read(results_file)
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

#To compare certain gene across groups:
sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')

#Mark cell types
new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
#adata.rename_categories('leiden', new_cluster_names)        #*** gives valueErr -> "leiden" and new_cluster_names are diff lengths? -> but looking at exported image, it has 8 categories, same as this?
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')



#Visualise marker genes
#??? Semicolons? **
sc.pl.dotplot(adata, marker_genes, groupby='leiden');
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90);



#----------------------------------------
#EXPORTING

adata.write(results_file, compression='gzip')  

adata.raw.to_adata().write('./write/pbmc3k_withoutX.h5ad')

#To export to csv
# Export single fields of the annotation of observations
# adata.obs[['n_counts', 'louvain_groups']].to_csv(
#     './write/pbmc3k_corrected_louvain_groups.csv')

# Export single columns of the multidimensional annotation
# adata.obsm.to_df()[['X_pca1', 'X_pca2']].to_csv(
#     './write/pbmc3k_corrected_X_pca.csv')

# Or export everything except the data using `.write_csvs`.
# Set `skip_data=False` if you also want to export the data.
# adata.write_csvs(results_file[:-5], )


#----------------------------------------
"""
what does it mean by "counts"? ln 72

what does .raw do?  "freezes the state of the AnnData object." ? may not be necessary?

regress out? unit variance?

what is "neighbourhood"? - review how PCA actually works 

numpy not used. why did we import it?

ERRs: ctrF ***
=========================

earlier part of Akram's looks diff, middle looks similar to tut i followed. mine has neighborhood graph/clustering


"""


