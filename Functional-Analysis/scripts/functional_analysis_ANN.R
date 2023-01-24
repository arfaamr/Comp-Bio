
#Arfaa Rashid
#Nov 30, 2022 
#functional_analysis.r
#https://hbctraining.github.io/Training-modules/DGE-functional-analysis/
  
# Program takes results of differential analysis on a gene set and performs functional analysis on them
# Uses over-representation analysis and functional class scoring methods. Visualizes results using various graphs

#**this version uses own annotations queried from EnsDb.Hsapiens.v75 to perform analysis as opposed to using given annotation file

#---------------------------
#Install packages

#install.packages(c("BiocManager", "devtools"))    
#BiocManager::install(c("clusterProfiler", "DOSE", "org.Hs.eg.db", "pathview", "AnnotationDbi", "EnsDb.Hsapiens.v75"))
#devtools::install_github("tidyverse/tidyverse")


#---------------------------
#Load libraries

library(dplyr)                  #tools for manipulating dataframes          [--not in workshop, but dplyr functions not found until this was loaded, prob due to some update]
library(tidyverse)              #collection of integrated packages for common data science functionalities, eg tibble

library(org.Hs.eg.db)           #tools to query orgDb -   gene feature info for particular organism, but only latest genomic build available 
library(EnsDb.Hsapiens.v75)     #tools to query EnsDb -   get transcript/gene info using gene ID; can specify version for correct genomic build
library(AnnotationDbi)          #annotation interface to connect various databases(?**)

library(clusterProfiler)        #tools for analysing/visualising data -   performs OR analysis w/ hypergeometric testing, performs GSEA
library(DOSE)                   #used with clusterProfiler to perform ORA(?**)
library(pathview)               #for visualisation -  integrate GSEA data into pathway images
library(ggrepel)                #for visualisation -  use text label geom in ggplot to repel labels & prevent overlaps


#---------------------------
#Read in given datafiles
#     - data is a list of genes
#     - link to data: https://hbctraining.github.io/Training-modules/DGE-functional-analysis/lessons/01_setting_up.html

# differential expression results
res_tableOE <- read.csv("data/Mov10oe_DE_results.csv", row.names = 1)   

# given annotation file [as opposed to annotations we queried ourselves later]
annotations_ahb <- read.csv("data/annotations_ahb.csv")

#Creates tibble (modified data frame) of DE results
res_tableOE_tb <- res_tableOE %>%
 rownames_to_column(var="gene") %>% 
 as_tibble()


#---------------------------
#Extract annotations from db
#   -> get info about each gene in res_tableOE_tb
#   - as opposed to using given annotation file -> to show we can query databases ourselves if annotations had not been given (?**)


#using org.Hs.eg.db

#return Ensembl IDs for a set of genes by getting gene names from res_tableOE and querying database for associated ENSEMBL, ENTREZ IDs
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db,
                                           keys = res_tableOE_tb$gene,                      #data used to query                         
                                           columns = c("ENSEMBL", "ENTREZID","GENENAME"),   #what we are looking for
                                           keytype = "SYMBOL")                             

#remove duplicates by keeping only first occurrences 
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

#some genes have NA - no annotation
#this is because we used a database for a different genomic build than our dataset was based on, so the above annotations are incorrect
#instead, use EnsDb.Hsapiens.v75 to get annotations for correct genomic build

#using EnsDb.Hsapiens.v75

annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                         keys = res_tableOE_tb$gene,
                                         columns = c("GENEID", "ENTREZID","GENEBIOTYPE"),   
                                         keytype = "SYMBOL")

# remove dups
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)
annotations_edb <- annotations_edb[non_duplicates_idx, ]

#no longer any NA entries, so it worked correctly
length(which(is.na(annotations_edb$GENEID)))


#---------------------------
#OR Analysis

#goal: determine which GO categories are overRep'd in a subset of under/over expressed genes

## Merge annotation df with initial results df
res_ids <- inner_join(res_tableOE_tb, annotations_edb, by=c("gene"="SYMBOL"))  #!ahb to edb to use own annotations queried from dbs rather than file given by workshop; gene_name to SYMBOL

## background set: all genes tested for significance
allOE_genes <- as.character(res_ids$GENEID)      #!gene_id to GENEID

## significant set: only genes with p-adjusted values<0.05 - ones with a significant difference [in ? geneexp?]
sigOE <- dplyr::filter(res_ids, padj < 0.05)
sigOE_genes <- as.character(sigOE$GENEID)        #!gene_id to GENEID

## Gene Ontology enrichment analysis 
#uses ensembl IDs to search the org.Hs.eg.db database, categorize by GO term, and perform over-rep analysis
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP",                    #specifies GO subcategory to compare by. BP - biological processes
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.02, 
                readable = TRUE)

# save results to rda and csv files
cluster_summary <- data.frame(ego) 
save(ego, file="results/ego_ANN.rda")
write.csv(cluster_summary,"results/clusterProfiler_Mov10oe_ANN.csv")


#---------------------------
#Visualization of ORA

pdf("results/plots_ANN.pdf")    #opens pdf file to write to 

# Creating Dotplot
dotplot(ego, showCategory=10, title = "OR Analysis")      #how can i find parameter names to manipulate graph..? args() returning ... no good documentation (?**)

# Creating Netplot - details the genes associated with one or more terms (by default gives the top 5 significant terms (by padj))
#extracting the log2 fold changes from results table by creating named vector, in order to colour graph by log2fold change
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene

options(ggrepel.max.overlaps = Inf)     #should prevent warning about too many unlabeled data points, but it's not working(?**). how to use ggrepel to repel them in the first place? 
#only shows some labels. has to do with warning? tried increasing max.overlaps, changing label size..

#create plot
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         color.params = list(foldChange=OE_foldchanges), 
         vertex.label.font=1,
         ggrepel.max.overlaps = Inf)

dev.off()                     #close pdf file


#---------------------------
#Gene Set Enrichment Analysis (FCS)

#goal: using log2fold vals to measure change, test differential expression over entire gene sets

## Remove NA values 
res_entrez <- dplyr::filter(res_ids, ENTREZID != "NA")      #!entrezid to ENTREZID
## Remove Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$ENTREZID) == F), ]    #!entrezid to ENTREZID
## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$ENTREZID                     #!entrezid to ENTREZID
## Sort fold changes
foldchanges <- sort(foldchanges, decreasing = TRUE)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges,     # ordered vector of fold changes named by entrez ID
                    organism = "hsa",           # supported organisms listed below [?] - search for homosapiens data in KEGG database?
                    nPerm = 1000,               # default number permutations [??]
                    minGSSize = 20,             # minimum gene set size - only look at categories with 20+ genes in them
                    pvalueCutoff = 0.05,        # padj cutoff value - only look at significantly differentially expressed genes
                    verbose = FALSE)

#extract and write GSEA results
gseaKEGG_results <- gseaKEGG@result
write.csv(gseaKEGG_results, "results/gseaOE_kegg.csv", quote=F)


## GSEA doesnt necessarily need to use KEGG pathways. could group by GO, etc as well
# - using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = foldchanges, 
                OrgDb = org.Hs.eg.db, 
                ont = 'BP', 
                nPerm = 1000, 
                minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE) 

gseaGO_results <- gseaGO@result


#---------------------------
#Visualizing GSEA

pdf("results/plots_GSEA_ANN.pdf")    #open pdf file to write to 

#kegg
# visualizing one (arbitrary?) pathway as GSEA plot
gseaplot(gseaKEGG, geneSetID = 'hsa03040')     #hsa03040 is id of a kegg pathway

# visualizing one pathway as pathway images
#detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts # err so skip...

pathview(gene.data = foldchanges,
         pathway.id = "hsa03040",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))
#^ ?** how to save images to results folder rather than in same folder? 
# how to rename so doesnt overwrite ones with old workshop annotiations?

#go
gseaplot(gseaGO, geneSetID = 'GO:0007423')

dev.off()



