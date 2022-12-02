
#Arfaa Rashid
#Nov 30, 2022
#functional_analysis.r


#---------------------------
#Install packages
#install.packages(c("BiocManager", "devtools"))    #tidyverse already installed
#BiocManager::install(c("clusterProfiler", "DOSE", "org.Hs.eg.db", "pathview", "AnnotationDbi", "EnsDb.Hsapiens.v75"))


#---------------------------
#Load libraries
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(clusterProfiler)        #performs OR analysis w/ hypergeometric testing
library(DOSE)
library(pathview)  #used for visualization
library(AnnotationDbi)
library(org.Hs.eg.db)           #package links to the relevant dbs for human organisms
library(EnsDb.Hsapiens.v75)     #package with dbs for correct build - must install correct release of Ensembl to query?


#---------------------------
#Read in data - Bulk RNA-seq data
#https://hbctraining.github.io/Training-modules/DGE-functional-analysis/lessons/01_setting_up.html - link to data
res_tableOE <- read.csv("data/Mov10oe_DE_results.csv", row.names = 1)   
#differential expression results
                                                                        
#w/out row.names=1, first col is row names instead. parameter sets which col to take row names from
annotations_ahb <- read.csv("data/annotations_ahb.csv")

#Creates tibble (modified data frame)
res_tableOE_tb <- res_tableOE %>%                                       # %>% - fancy syntax. equivalent is as_tibble(rownames_to_column(res_tableOE, var="gene"))
 rownames_to_column(var="gene") %>% 
 as_tibble()


#---------------------------
#Extract annotations from db -> get info about each gene in res_tableOE_tb


#org.Hs.eg.db

#Return the Ensembl IDs for a set of genes  [returns list(table) of wanted info]
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = res_tableOE_tb$gene,  # data to use for retrieval                                 [want info about genes in gene col of our tibble]
                                           columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retrieve for given data  [specific info wanted]
                                           keytype = "SYMBOL") # type of data given in 'keys' argument                              [? not really sure. heading for keys col(what we gave) in returned table? use keytypes() to find..]

#sometimes duplicates are returned - will keep just first versions
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)    #indices for non dups
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

#some genes have NA - no annotation
# this is because we used a database for a different gene build than our dataset was based on, so the above annotations are incorrect.. 


#EnsDb.Hsapiens.v75

# Return the Ensembl IDs for a set of genes
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                         keys = res_tableOE_tb$gene,
                                         columns = c("GENEID", "ENTREZID","GENEBIOTYPE"),   
                                         keytype = "SYMBOL")

# remove dups
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)
annotations_edb <- annotations_edb[non_duplicates_idx, ]

#no longer any NA entries, so it worked right :>
length(which(is.na(annotations_edb$GENEID)))


#---------------------------
#OR Analysis

## Merge annotation df with initial results df
res_ids <- inner_join(res_tableOE_tb, annotations_ahb, by=c("gene"="gene_name"))      #list with cols from orig df and annotations[categories ?]

## background set: all genes tested for significance in results           
allOE_genes <- as.character(res_ids$gene_id)    

## significant set: only genes with p-adjusted values<0.05 -> ? we will determine which categories are overRep'd in a subset of under/over expressed genes?
sigOE <- dplyr::filter(res_ids, padj < 0.05)

sigOE_genes <- as.character(sigOE$gene_id)

## Gene Ontology enrichment analysis 
#gene is sig set, universe is background set; uses ensembl IDs to search the org.Hs.eg.db database. *does gene build no longer matter?
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP",                    #specifies GO subcategory to compare by. BP - biological processes
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

cluster_summary <- data.frame(ego) #what does this mean *? what are cols 5-7?
save(ego, file="results/ego.rda")
write.csv(cluster_summary,"results/clusterProfiler_Mov10oe.csv")


#---------------------------
#Visualization of ORA

pdf("results/plots.pdf")    #open pdf file to write to 

#dotplot
dotplot(ego, showCategory=50, title = "OR Analysis")      #* how can i find parameter names to manipulate graph..? args() returning ...


#netplot
## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
# (fold change: measurement for change b/w A, B. fold change of B WRT A = B/A)
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene

options(ggrepel.max.overlaps = Inf)     #*should prevent warning about too many unlabeled data points, but it's not working..

## cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         color.params = list(foldChange=OE_foldchanges), 
         vertex.label.font=1,
         ggrepel.max.overlaps = Inf)
#*only shows some labels. has to do with warning? tried increasing max.overlaps, changing label size

dev.off()                     #close pdf file


#---------------------------
#Gene Set Enrichment Analysis (FCS)

## Remove NA values 
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
## Remove Entrez duplicates [use Entrez IDs this time?]
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
## Sort fold changes in decreasing order  [why?*]
foldchanges <- sort(foldchanges, decreasing = TRUE)
