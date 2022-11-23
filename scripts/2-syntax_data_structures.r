
#Intro to R
#Syntax and Data Structures
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