
#Arfaa Rashid
#Nov 30, 2022 
#functional_analysis_GO.r
#https://hbctraining.github.io/Training-modules/DGE-functional-analysis/
  
#takes results of differential analysis on gene set and performs functional analysis on it using over-representation analysis and functional class scoring methods. Visualizes results of FA using various graphs

# *** this version performs ORA using different gene ontology terms - molecular func and cellular component rather than biological process as in the original


#---------------------------
#Install packages

#install.packages(c("BiocManager", "devtools"))    
#BiocManager::install(c("clusterProfiler", "DOSE", "org.Hs.eg.db", "pathview", "AnnotationDbi", "EnsDb.Hsapiens.v75"))
#devtools::install_github("tidyverse/tidyverse")


#---------------------------
#Load libraries

library(tidyverse)              #collection of integrated packages for common data science functionalities, eg tibble
library(org.Hs.eg.db)           #package to query annotation db - gene feature info for particular organism, but only latest genomic build available 
library(EnsDb.Hsapiens.v75)     #package to query annotation db - get transcript/gene info using gene ID; can specify version for correct genomic build
library(clusterProfiler)        #performs OR analysis w/ hypergeometric testing
library(ggrepel)                #for visualisation - text label geom for ggplot to repel labels & prevent overlaps
library(DOSE)
library(pathview)               #used for visualisation
library(AnnotationDbi)


#---------------------------
#Read in given datafiles -> data is a list of genes - link to data: https://hbctraining.github.io/Training-modules/DGE-functional-analysis/lessons/01_setting_up.html

res_tableOE <- read.csv("data/Mov10oe_DE_results.csv", row.names = 1)   
#differential expression results
                                                        
annotations_ahb <- read.csv("data/annotations_ahb.csv")

#Creates tibble (modified data frame)
res_tableOE_tb <- res_tableOE %>%
 rownames_to_column(var="gene") %>% 
 as_tibble()


#---------------------------
#Extract annotations from db -> get info about each gene in res_tableOE_tb


#using org.Hs.eg.db

#Return the Ensembl IDs for a set of genes - get gene names from res_tableOE and query database for their associated ENSEMBL IDs, ENTREZ IDs
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

# Return the Ensembl IDs for a set of genes
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                         keys = res_tableOE_tb$gene,
                                         columns = c("GENEID", "ENTREZID","GENEBIOTYPE"),   
                                         keytype = "SYMBOL")

# remove dups
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)
annotations_edb <- annotations_edb[non_duplicates_idx, ]

#no longer any NA entries, so it worked correctly
length(which(is.na(annotations_edb$GENEID)))

#***what is the point of all this ^ ? never used?
#*annotations_edb contains mostly same info as annotations_ahb, which workshop says to read in and use instead
#*seems to have more info
#---------------------------
#OR Analysis

#goal: determine which categories are overRep'd in a subset of under/over expressed genes*

## Merge annotation df with initial results df
res_ids <- inner_join(res_tableOE_tb, annotations_ahb, by=c("gene"="gene_name"))

## background set: all genes tested for significance in results           
allOE_genes <- as.character(res_ids$gene_id)    

## significant set: only genes with p-adjusted values<0.05 - ones with a significant difference [in ? geneexp?]
sigOE <- dplyr::filter(res_ids, padj < 0.05)
sigOE_genes <- as.character(sigOE$gene_id)

## Gene Ontology enrichment analysis 
#uses ensembl IDs to search the org.Hs.eg.db database, categorize by GO term, and perform over-rep analysis
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "CC",                    #specifies GO subcategory to compare by. BP - biological processes
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.02, 
                readable = TRUE)

cluster_summary <- data.frame(ego) 
save(ego, file="results/ego_CC.rda")
write.csv(cluster_summary,"results/clusterProfiler_Mov10oe_CC.csv")


#---------------------------
#Visualization of ORA

pdf("results/plots_CC.pdf")    #open pdf file to write to 

# Creating Dotplot
dotplot(ego, showCategory=10, title = "OR Analysis by CC")      #*** how can i find parameter names to manipulate graph..? args() returning ...


# Creating Netplot
## extracting the log2 fold changes from results table by creating named vector, in order to colour graph by log2fold change
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene

options(ggrepel.max.overlaps = Inf)     #***should prevent warning about too many unlabeled data points, but it's not working..

## cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         color.params = list(foldChange=OE_foldchanges), 
         vertex.label.font=1,
         ggrepel.max.overlaps = Inf)

#***only shows some labels. has to do with warning? tried increasing max.overlaps, changing label size..

dev.off()                     #close pdf file




