
#Arfaa Rashid
#Nov 23, 2022
#4-data_wrangling.r
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

