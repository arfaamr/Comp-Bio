
#Arfaa Rashid
#Intro to R
#individual scripts for each lesson have been compiled into one

#2 - Syntax and Data Structures
#Nov 22, 2022

#vectors
glengths <- c(4.6, 3000, 50000)
glengths

species <- c("ecoli", "human", "corn")
species

combined <- c(glengths, species)


#factors

expression <- c("low", "high", "medium", "high", "low", "medium", "high")
expression <- factor(expression)
#low, med, high are categories, and are assigned integer values alphabetically

samplegroup <- c("CTL", "CTL", "CTL", "KO", "KO", "KO", "OE", "OE", "OE")
samplegroup <- factor(samplegroup)
samplegroup

df <- data.frame(species, glengths)
df


#---------------------------------------------------
#2 - Libraries, Packages, Functions
#Nov 23, 2022

#install.packages("ggplot2") #run once
library(ggplot2)
library(tidyverse)
search()          #ggplot2 shows after library loaded
sessionInfo()

#install.packages("tidyverse")


#Nov 23, 2022
#4-data_wrangling1.r
#https://hbctraining.github.io/Intro-to-R/lessons/04_introR-data-wrangling.html


###READING DATA

metadata <- read.csv(file="./data/mouse_exp_design.csv", stringsAsFactors = TRUE)
#by default, converts cols with chars into factors
#supposedly, but according to str(), didnt? must include specification arg

metadata

str(metadata)     #structure function

##########

age <- c(3,5,6)
idx1 <- c(2,3)
idx2 <- c(3,5)

age[idx1]
age[idx2]
age[1:3]

age[age>4]

##EX
metadata[1]
WTSamplegroup = metadata[metadata[1]!="KO"]
WTSamplegroup
#^ Ex: Extract only those elements in samplegroup that are not KO
#i think this is right but idk?


### RELEVELING FACTORS
expression <- c("low", "high", "medium", "high", "low", "medium", "high")
expression <- factor(expression)
#low, med, high are categories, and are assigned integer values alphabetically

expression <- factor(expression, levels=c("low", "medium", "high"))


#Nov 24, 2022
#5-data_wrangling2.r
#https://hbctraining.github.io/Intro-to-R/lessons/04_introR-data-wrangling.html

metadata <- read.csv(file="./data/mouse_exp_design.csv", stringsAsFactors = TRUE)

metadata

metadata$genotype


#---------------------------------------------------
#Nov 27, 2022
#6: match_reorder.r


rpkm_data <- read.csv("./data/counts.rpkm.csv", stringsAsFactors = TRUE)
head(rpkm_data) #shows first few lines, to see what data looks like

#note that colnames in data match rownames in metadata. However, they are out of order

ncol(rpkm_data) == nrow(metadata) #is there corresponding data for each piece of metadata? TRUE -> yes

A = c(1,2,3,4)
B = c(4,3,2,1)

A == B
all(A == B)

x <- rownames(metadata)
y <- colnames(rpkm_data)

all(x %in% y) #TRUE
all(x == y) #FALSE
#all present, but in the wrong order. Must use a reordering method

#reordering---

important_genes <- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", "ENSMUSG00000030970")

#? wants to extract rows corresponding to above genes with %in%, but how? 
indices <- important_genes %in% rpkm_data
indices
head(rpkm_data)
x
head(rpkm_data[indices])
#rpkm_data[indices] 
typeof(rpkm_data)
rpkm_data[1,] #row1
rpkm_data[,]

#confusion ** RETURN

vect <- c('A','B','C','D')
vect[c(1,2)] #vector containing index 1, 2
vect
vect[c(2,1)] #returns element at 2, element at 1. essentially reversed
vect

#reorder vector

orig <- c("A", "B", "C", "D")
want <- c("D", "B", "A", "C")

reord_idx <- match(want,orig)
orig_reord <- orig[reord_idx]
orig_reord # orig vector (A B C D) has been ordered to wanted vector (D B A C)

#get a subset

first <- c('A', 'B', 'C', 'D')
second <- c('B', 'A', 'D')
idx <- match(first, second)

second[idx] #gives A B NA D, so reorders rest with placeholder for C

#reorder given dataset and metadata so they match

rownames(metadata)    #in order
colnames(rpkm_data)   #not in order

genomic_idx <- match(rownames(metadata), colnames(rpkm_data))
rpkm_ordered <- rpkm_data[,genomic_idx]
head(rpkm_data)
head(rpkm_ordered)


#---------------------------------------------------
#Nov 29, 2022
#7: data_visualization.r


head(rpkm_ordered$sample1)
mean(rpkm_ordered$sample1)
#we want a vector of size 12 containing means of all 12 samples, so it can be added as a col in the metadata
#use map()

library(purrr)

sample_means <- map_dbl(rpkm_ordered, mean)   #new vector for means (named index sample_means.)
#length(sample_means)


age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32) #new vector for age of mice

new_metadata <- data.frame(metadata, sample_means, age_in_days) #data.frame() adds cols
view(new_metadata)

## Open device for exporting plot to pdf
pdf("./figures/scatterplot.pdf")

ggplot(new_metadata) + geom_point(aes(x=age_in_days, y=sample_means,color=genotype, shape=celltype),size=2.25) + theme(axis.title = element_text(size=rel(1.5)),plot.title=element_text(hjust=0.5)) + xlab("Age (days)") + ylab("Mean expression") + ggtitle("Mice Gene Expression WRT Age")

ggplot(new_metadata) + geom_boxplot(aes(x=genotype, y=sample_means, shape=celltype, color=genotype)) + theme(axis.title = element_text(size=rel(1.5)),plot.title=element_text(hjust=0.5)) + xlab("Genotype") + ylab("Mean expression") + ggtitle("Mice Gene Expression WRT Genotype")

dev.off() #must close file to finish exporting

#? my scatterplot doesnt look exactly like theirs :( 
#Some y-values are off on the graph, but the metadata dataframe matches theirs

