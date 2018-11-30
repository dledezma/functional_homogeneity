######################################################################################################################
#                                               Coexpresion_11.0                                                     #
######################################################################################################################

# Author: Luis Fernando Altamirano Pacheco
# Date: 16/04/2018

# Version 2
# Date: 18/04/2018
# This version incorporates to the analysis a plot as an output, in the x-axis of the plot is the regulon size and in 
# the y-axxis the correlation. It also takes the absolute value of each correlation before calculating the mean

# Version 3
# Date: 30/04/2018
# This version adds to the output:
#  - Two scatterplots of the median and the mean of the correlations vs the regulon size
#  - Two scatterplots for the transcription units, one with the median and another of the mean of the 
#    correlations, as in the previous plots for the regulons. 
#  - A table as a positive control with two columns. The first one for the name of every gene in the Colombos dataset, 
#    and the second for the self correlation
#
# Also a function with the "main code" is created, so it can be used for different data sets.

# Version 4
# Date: 10/05/2018
# The plots are generated again, all with bnumbers instead; the previous version computed the TU correlations with the 
# names.
# It is included in the analysis, as a positive control, a set of clusters of genes that are known to be coexpressed
# Histograms of the correlations are created, excluding the sets of genes with a size of one 

# Version 5
# Date: 13/02/2018
# A file with the comparision of the boxplots od the mean and the median of the TUS (bnumbers vs names) is added

# Version 6
# Date: 21/05/2018
# Modifications:
#   It Ads a summary table
#   The method to get the self correlations is modified
#   A summary table is created
#   The code can be run now at the terminal
# Usage:
# Rscript Coexpresion_6.0.R name/of/the/data/directory/ NameOfFile1 NameOfFile2 ... Name1 Name2 ...
# Where Name1 and Name2 are the names that are going to be assigned to the corresponding file in the analysis, those 
# are the names that are going to be in the outputs

# Version 7 
# Date: 25/05/2018
# Modifications:
#   The boxplots of the medians and the means are separated in two files, one for the means and other for the medians
# Usage:
# Rscript Coexpresion_7.0.R name/of/the/data/directory/ NameOfFile1 NameOfFile2 ... Name1 Name2 ...
# Where Name1 and Name2 are the names that are going to be assigned to the corresponding file in the analysis, those 
# are the names that are going to be in the outputs

# Version 8
# Date: 26/05/2018
# Modifications:
#   Wilcoxon tests are performed to all possible pairs of the distributions of the mean/median of the data sets
# Usage:
# Rscript Coexpresion_8.0.R name/of/the/data/directory/ NameOfFile1 NameOfFile2 ... Name1 Name2 ...
# Where Name1 and Name2 are the names that are going to be assigned to the corresponding file in the analysis, those 
# are the names that are going to be in the outputs

# Version 9
# Date: 26/06/2018, 02/07/2018
# Modifications:
#   Coments are added for better understanding of the code
#   The random controls are added
#   The boxplots pf the mean and median are generated also in a different format (with the base instead of ggplot2)
# Rscript Coexpresion_9.0.R name/of/the/data/directory/ NameOfFile1 NameOfFile2 ... Name1 Name2 ...

# Version 10
# Date: 04/07/2018
# Modifications
#    It uses 34 cores (maximum) instead of 1, only for some parts
#    Tiny modification (20/07/2018): The regulons with size = 1 are deleted from the control
# Rscript Coexpresion_10.0.R name/of/the/data/directory/ NameOfFile1 NameOfFile2 ... Name1 Name2 ...

# Version 11
# Date: 14/10/2018
# Modificcations
#    It is added twneo boxplot plot more, the saame mmean and median boxplots are plotted but without "GOs", 
#    "pathways" and their controls

# Important Notes
# The code is explecity multicore, to change the number of cores that the program should use, make sure to change the 
# number in every occurence of the parameter "mc.cores"

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(parallel)
#setwd("/home/feraltp/PGC/Coexpresion")

#################### Loading the data ####################

Colombos <- read.table("Data/Colomobos/colombos_ecoli_exprdata_20151029.txt", sep="\t", header=TRUE, row.names = 1)

ToList <-function(Genes){ # Genes: The path of a file with the next format:
	# SetName1	gene_1	gene_2	gene_3 ... gene_x
	# SetName2	gene_i	gene_i+1	...	gene_n
	List <- scan(Genes, what="", sep="\n") # scanning the file
	List <- strsplit(List, "[[:space:]]+") # stripping the file by space characters
	names(List) <- sapply(List, `[[`, 1) # Put the name of the set of genes (e.g. regulon, TU) as the name of each element in the list
	List <- lapply(List, `[`, -1) # Deleting the first element, i.e. the name of the list
}

args <- commandArgs(trailingOnly=TRUE)
l <- (length(args)-1)/2 # counting the number of gene sets
dir <- args[1] # directory where the gene sets can be found
sets <- args[2:(l+1)] # name of the files with the gene sets
sets <- lapply(paste(dir, sets, sep=""), ToList) # gene sets converted to list. This object is a list of lists, each list (a file of a gene set) has a set of vectors with the names of the genes for each group of genes. The name of each group is the name of each element of the list
names(sets) <- args[(l+2):(2*l+1)] # name given to each gene set

check <- function(X){  # This function cheks if all the genes of X are in the row names of Colomobos data set. The function returns those genes that are in Colombos
	X[X %in% rownames(Colombos)]
}

NoRepeats <-function(x){ # It deletes the repeated data from the symmetric matrix
	sort(x)[c(TRUE, FALSE)]
}
#################### Main function ####################

Correlation <- function(X, b=0){

	i = 1
	size <- vector()
	MeanCor <- vector()
	MedCor <- vector()
	names <- vector()
	df <- data.frame(Name = factor(), Cor = numeric(), Size=factor())
	for (regulon in names(X)){
		print(regulon)
		# Correlations are computed, the dataset is transposed because the correlations are made by column: 	
		if(   b==0 & length(check(X[[regulon]]))   |   b==1 & any(Colombos$Gene.name %in% X[[regulon]])   ){
			
			if(!b){
				RegCor <- abs(cor(t(Colombos[check(X[[regulon]]),-(1:2)]), use = "complete.obs", method = "spearman"))
			} else {
				RegCor <- abs(cor(t(Colombos[Colombos$Gene.name %in% X[[regulon]],-(1:2)]), use = "complete.obs", method = "spearman"))
			}
						
			diag(RegCor) <- NA #The diagonal is excluded...
			CorMean <- mean(RegCor, na.rm=TRUE) # ... so it can't affect the calculation of the mean correlations
			CorMed <- median(RegCor, na.rm=TRUE)			
			if(!b){
				len <- length(check(X[[regulon]])) # The length of the regulon is calculated
			} else {
				len <- sum(Colombos$Gene.name %in% X[[regulon]])
			}

			size[i] <-len
			MeanCor[i] <- CorMean
			MedCor[i] <- CorMed
			names[i] <- regulon
			i <- i+1

			if(length(as.vector(RegCor)) != 1){ # If the correlation matrix is not 1x1 ...
				df <- rbind(df, data.frame(Name = regulon, Cor = NoRepeats(na.omit(as.vector(RegCor))), Size=len )) # A data frame storing al the correlation coefficients for every posiible pair within the group of the set of genes, excluding the auto-correlations
			}
		}
	}
	CorStatistics <- data.frame(Names=names, RegulonSize=size, SpearmanCorrelationMean=MeanCor, SpearmanCorrelationMedian=MedCor)
	CorStatistics$SpearmanCorrelationMean[is.na(CorStatistics$SpearmanCorrelationMean)] <- 1 #Restoring the correlations that turned to NaN in previous steps
	CorStatistics$SpearmanCorrelationMedian[is.na(CorStatistics$SpearmanCorrelationMedian)] <- 1
	return(list(CorStatistics, df))
}

#################### Scatterplots and tables ####################

pdf("Results/Plots/PlotRegulonSizeVSCor.pdf") # The plots are saved in a pdf. The medians/means of the correlations are plotted against the size of the regulon, 
# i.e. the number of genes that compose it
# Each point represents the median or the mean of the correlations of  all possible pairs between two genes in each gene group. The self correlations 
# are excluded unless the size of the gene group is one.

correlations <- mclapply(sets, Correlation, mc.cores=20)

for (set in names(sets)){
	Reg <- correlations[[set]]
	CorData <- Reg[[1]] # A data frame with the median and mean of the correlations per gene group.
	BoxReg <- Reg[[2]] # A data frame with all correlations and the name of the gene group to which they correspond
	write.table(CorData, file =paste("Results/Tables/", set, "_", "correlations.txt", sep=""), sep = "\t", row.names = FALSE) # A file is created to save the correlations for each gene set
	print(ggplot(CorData, aes(RegulonSize,SpearmanCorrelationMean)) + geom_point() + ggtitle( paste("Mean", set, "Correlations") ) )
	print(ggplot(CorData, aes(RegulonSize,SpearmanCorrelationMedian)) + geom_point() + ggtitle(paste("Median", set, "Correlations")))

}

dev.off()

#----------------- Random Controls -------------------------

RandomPaths <- paste(dir, "random/random_", names(sets), sep="")
random <- lapply(RandomPaths, function(x) {
	paste(x, "/processed_ran_", 1:100, ".txt", sep="")
})  # list(c("/path1/processed_ran_1", "/path1/processed_ran_2"), c("/path2/processed_ran_1", "/path2/processed_ran_2"))

random <- lapply(random, function(x) {
	lapply(x, ToList)
	})
names(random) <- names(sets) # c("name1", "name2")

random <- lapply(names(random), function(x){ 
	all100 <- mclapply(random[[x]], Correlation, mc.cores=34) # A list of the results of the Correlation function for each file (gene set) of the 100 random files for just a type of gene set (e.g. TUs or Simple) 
	all100 <- lapply(all100, `[[`, 1) # Keeping just CorStatistics
	all100 <- lapply(all100, function(x){
		x[x$RegulonSize!=1,]
		})
	means <- lapply(all100, `[[`, "SpearmanCorrelationMean")
	medians <- lapply(all100, `[[`, "SpearmanCorrelationMedian")
	means <- do.call(c, means) # The vectors of the means are joined
	medians <- do.call(c, medians) # The vectors of the medians are joined
	return(data.frame(means=means, medians=medians, name=factor(x, levels = names(random)) ) )
	})

random <- do.call(rbind, random) # A dataframe with 3 columns: the name of the gene set, the medians and the means of the correlations

############################### Box Plots ###############################
quantdif <- function(x){ # Function that computes the "box size"
	quantile(x)["75%"]-quantile(x)["25%"]
}

i <-1
Mean <- vector("list", length = length(correlations))
Med <- vector("list", length = length(correlations))

for (set in names(sets)){
	Box <- correlations[[set]][[2]]
	Box$Name <- reorder(Box$Name, Box$Cor, quantdif)  # The data is ordered according to the "box size"
	wi <- 0.6*length(unique(Box$Name)) + 0.6*3
	pdf(paste("Results/Plots/Boxplots_",set,".pdf", sep=""), width=wi)
	print(ggplot(Box, aes(x = Name, y = Cor, fill=Size)) + scale_color_gradient(low="blue", high="red")+ geom_boxplot() + geom_point()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dev.off()

	CorData <- correlations[[set]][[1]]
	CorData <- CorData[CorData$RegulonSize !=1,] # Deleting the rows with regulon size equal to one: 

	Mean[[i]] <- CorData$SpearmanCorrelationMean
	Med[[i]] <- CorData$SpearmanCorrelationMedian

	i <- i+1
}

names(Mean) <- names(sets)
names(Med) <- names(sets)


# ---------- Wilcoxon test -----------------
wilcox <- function(x){ # Computes a symmetrical matrix of p-values between the correlations of every two gene sets
	WilMat <- matrix(nrow=length(x), ncol=length(x))
	colnames(WilMat) <- names(x)
	rownames(WilMat) <- names(x)
	for (set1 in names(x)){
		for(set2 in names(x)){
			WilMat[set1,set2] <- wilcox.test(x[[set1]],x[[set2]])$p.value
		}
	}
	return(WilMat)
}

wilcoxMean <- wilcox(Mean)
wilcoxMed <- wilcox(Med)

pdf("Results/Plots/WilcoxonTest.pdf")
heatmap.2(wilcoxMean, cellnote = round(wilcoxMean,2), notecol="black", density.info="none", trace="none", main="Mean", col=colorRampPalette(c("#00d0ff", "#9b0359"))(n = 300), margins = c(9, 9))
heatmap.2(wilcoxMed, cellnote = round(wilcoxMed,2), notecol="black", density.info="none", trace="none", main="Median", col=colorRampPalette(c("#00d0ff", "#9b0359"))(n = 300), margins = c(9, 9)) #00ff1d
dev.off()

# -----------------------------------------

Mean <- lapply(names(Mean), function(x){
	data.frame(category=factor(x, levels = names(Mean)), cor=Mean[[x]])
	})

Med <- lapply(names(Med), function(x){
	data.frame(category=factor(x, levels = names(Med)), cor=Med[[x]])
	})

Mean <- do.call(rbind, Mean)
Med <- do.call(rbind, Med)

# Sorting the factors to plot the boxplots as the order specified by the arguments
Mean$category = factor(Mean$category, args[(l+2):(2*l+1)])
Med$category = factor(Med$category, args[(l+2):(2*l+1)])
random$name = factor(random$name, args[(l+2):(2*l+1)])

# Mean
pdf("Results/Plots/Mean_Boxplots.pdf", width=14)
ggplot(Mean, aes(x = category, y = cor)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

labels <- character(length(levels(Mean$category))*2)
labels[c(T,F)] <- levels(Mean$category)
labels[c(F,T)] <- "Random"
boxplot(cor~category, data=Mean, frame=FALSE, at=seq(from=1,to=22,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "lightcoral", "darkgoldenrod2", "greenyellow", "paleturquoise4", "darkmagenta", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,22.7))
boxplot(means~name, data=random, frame=FALSE, at=seq(from=2,to=22,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(random$name))), col=rep("gainsboro", length(levels(random$name))), boxwex=0.4, notch=TRUE, add=TRUE)
axis(1, las=2, cex.axis=1, at=1:22, labels=labels, ,mgp=c(3,1,0))
axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
title(main="Means of Coexpressions (with Simplex)", font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)

MeanNoSimplex <- Mean[Mean$category!="Simplex",] # To delete the level "Simplex"
MeanNoSimplex$category <- factor(MeanNoSimplex$category) 
randomNoSimplex <- random[random$name!="Simplex",] 
randomNoSimplex$name <- factor(randomNoSimplex$name)
labelsNoSimplex <- character(length(levels(MeanNoSimplex$category))*2)
labelsNoSimplex[c(T,F)] <- levels(MeanNoSimplex$category)
labelsNoSimplex[c(F,T)] <- "Random"

boxplot(cor~category, data=MeanNoSimplex, frame=FALSE, at=seq(from=1,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "paleturquoise4", "darkmagenta", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,20.7))
boxplot(means~name, data=randomNoSimplex, frame=FALSE, at=seq(from=2,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoSimplex$name))), col=rep("gainsboro", length(levels(randomNoSimplex$name))), boxwex=0.4, notch=TRUE, add=TRUE)
axis(1, las=2, cex.axis=1, at=1:20, labels=labelsNoSimplex ,mgp=c(3,1,0))
axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)

	jpeg("Results/Plots/Mean_Boxplots_NoSimplex.jpeg") # Same last plot in the pdf file but in jpeg format
	boxplot(cor~category, data=MeanNoSimplex, frame=FALSE, at=seq(from=1,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "paleturquoise4", "darkmagenta", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,20.7))
	boxplot(means~name, data=randomNoSimplex, frame=FALSE, at=seq(from=2,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoSimplex$name))), col=rep("gainsboro", length(levels(randomNoSimplex$name))), boxwex=0.4, notch=TRUE, add=TRUE)
	axis(1, las=2, cex.axis=1, at=1:20, labels=labelsNoSimplex, mgp=c(3,1,0))
	axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
	title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)
	dev.off()

#Boxplot without GOs and pathways:
MeanNoPathGO <- MeanNoSimplex[MeanNoSimplex$category!="GOs" & MeanNoSimplex$category!="Pathways",] # To delete the levels "GOs" and "Pathways"
MeanNoPathGO$category <- factor(MeanNoPathGO$category) 
randomNoPathGO <- randomNoSimplex[randomNoSimplex$name!="GOs" & randomNoSimplex$name!="Pathways",] 
randomNoPathGO$name <- factor(randomNoPathGO$name)
labelsNoPathGO <- character(length(levels(MeanNoPathGO$category))*2)
labelsNoPathGO[c(T,F)] <- levels(MeanNoPathGO$category)
labelsNoPathGO[c(F,T)] <- "Random"

boxplot(cor~category, data=MeanNoPathGO, frame=FALSE, at=seq(from=1,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,16.7))
boxplot(means~name, data=randomNoPathGO, frame=FALSE, at=seq(from=2,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoPathGO$name))), col=rep("gainsboro", length(levels(randomNoPathGO$name))), boxwex=0.4, notch=TRUE, add=TRUE)
axis(1, las=2, cex.axis=1, at=1:16, labels=labelsNoPathGO ,mgp=c(3,1,0))
axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)

	jpeg("Results/Plots/Mean_Boxplots_NoSimplexGosPathways.jpeg") # Same last plot in the pdf file but in jpeg format
	boxplot(cor~category, data=MeanNoPathGO, frame=FALSE, at=seq(from=1,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,16.7))
	boxplot(means~name, data=randomNoPathGO, frame=FALSE, at=seq(from=2,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoPathGO$name))), col=rep("gainsboro", length(levels(randomNoPathGO$name))), boxwex=0.4, notch=TRUE, add=TRUE)
	axis(1, las=2, cex.axis=1, at=1:16, labels=labelsNoPathGO ,mgp=c(3,1,0))
	axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
	title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)
	dev.off()

dev.off()




#Median
pdf("Results/Plots/Median_Boxplots.pdf", width=14)
ggplot(Med, aes(x = category, y = cor)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

boxplot(cor~category, data=Med, frame=FALSE, at=seq(from=1,to=22,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "lightcoral", "darkgoldenrod2", "greenyellow", "paleturquoise4", "darkmagenta", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,22.7))
boxplot(medians~name, data=random, frame=FALSE, at=seq(from=2,to=22,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(random$name))), col=rep("gainsboro", length(levels(random$name))), boxwex=0.4, notch=TRUE, add=TRUE)
axis(1, las=2, cex.axis=1, at=1:22, labels=labels, mgp=c(3,1,0))
axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
title(main="Medians of Coexpressions (with Simplex)", font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)

MedNoSimplex <- Med[Med$category!="Simplex",] # To delete the level "Simplex"
MedNoSimplex$category <- factor(MedNoSimplex$category)

boxplot(cor~category, data=MedNoSimplex, frame=FALSE, at=seq(from=1,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "paleturquoise4", "darkmagenta", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,20.7))
boxplot(medians~name, data=randomNoSimplex, frame=FALSE, at=seq(from=2,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoSimplex$name))), col=rep("gainsboro", length(levels(randomNoSimplex$name))), boxwex=0.4, notch=TRUE, add=TRUE)
axis(1, las=2, cex.axis=1, at=1:20, labels=labelsNoSimplex, mgp=c(3,1,0))
axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)

	jpeg("Results/Plots/Median_Boxplots_NoSimplex.jpeg") # Same last plot in the previous pdf file but in jpeg format
	boxplot(cor~category, data=MedNoSimplex, frame=FALSE, at=seq(from=1,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "paleturquoise4", "darkmagenta", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,20.7))
	boxplot(medians~name, data=randomNoSimplex, frame=FALSE, at=seq(from=2,to=20,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoSimplex$name))), col=rep("gainsboro", length(levels(randomNoSimplex$name))), boxwex=0.4, notch=TRUE, add=TRUE)
	axis(1, las=2, cex.axis=1, at=1:20, labels=labelsNoSimplex, mgp=c(3,1,0))
	axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
	title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)
	dev.off()

#Boxplot without GOs and pathways:
MedNoPathGO <- MedNoSimplex[MedNoSimplex$category!="GOs" & MedNoSimplex$category!="Pathways",] # To delete the levels "GOs" and "Pathways"
MedNoPathGO$category <- factor(MedNoPathGO$category) 

boxplot(cor~category, data=MedNoPathGO, frame=FALSE, at=seq(from=1,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,16.7))
boxplot(medians~name, data=randomNoPathGO, frame=FALSE, at=seq(from=2,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoPathGO$name))), col=rep("gainsboro", length(levels(randomNoPathGO$name))), boxwex=0.4, notch=TRUE, add=TRUE)
axis(1, las=2, cex.axis=1, at=1:16, labels=labelsNoPathGO ,mgp=c(3,1,0))
axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)

	jpeg("Results/Plots/Median_Boxplots_NoSimplexGosPathways.jpeg") # Same last plot in the pdf file but in jpeg format
	boxplot(cor~category, data=MedNoPathGO, frame=FALSE, at=seq(from=1,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=TRUE, col=c("deepskyblue3", "tan3", "palevioletred3", "sienna1", "darkgoldenrod2", "greenyellow", "darkolivegreen4", "cadetblue3"), boxwex=0.4, notch=TRUE,  xlim=c(0.3,16.7))
	boxplot(medians~name, data=randomNoPathGO, frame=FALSE, at=seq(from=2,to=16,by=2), outpch=4, outcex=0.5, border="gray47", xaxt="n", yaxt="n", names=rep("Random", length(levels(randomNoPathGO$name))), col=rep("gainsboro", length(levels(randomNoPathGO$name))), boxwex=0.4, notch=TRUE, add=TRUE)
	axis(1, las=2, cex.axis=1, at=1:16, labels=labelsNoPathGO ,mgp=c(3,1,0))
	axis(2, las=2, cex.axis=1.3, at=seq(from=0, to=1, by=.1))
	title(font.lab=2, ylab="Spearman Correlation", cex.lab=1.3, cex.main=1.5)
	dev.off()

dev.off()




Corr <- lapply(correlations, function(x){
	x[[2]]$Cor
	})

#names(Corr) <- names(sets)

Corr <- lapply(names(Corr), function(x){
	data.frame(category=factor(x, levels = names(Corr)), cor=Corr[[x]])
	})

Corr <- do.call(rbind, Corr)

pdf("Results/Plots/all_pairs_boxplot.pdf")
ggplot(Corr, aes(x = category, y = cor)) + geom_boxplot()
dev.off()

############## Histograms ###############

pdf("Results/Plots/HistogramsCor.pdf")

ggplot(Corr, aes(x = cor, fill= category)) + geom_histogram(position="identity", alpha=0.5)

for (set in names(sets)){
	print(ggplot(correlations[[set]][[2]], aes(x = Cor)) + geom_histogram() + ggtitle(set))
}

dev.off()

pdf("Results/Plots/Reg_TU_Clus_HistogramCor.pdf")
ggplot(Corr[Corr$category %in% c("Regulon", "TU", "Cluster"),], aes(x = cor, fill= category)) + geom_histogram(position="identity", alpha=0.5)
dev.off()

############# Summary Table ######################  M -> mean   m -> median
MoM <- sapply(correlations, function(x){ # Mean of means
	mean(x[[1]]$SpearmanCorrelationMean[x[[1]]$RegulonSize!=1])
	})

moM <- sapply(correlations, function(x){ # Median of means
	median(x[[1]]$SpearmanCorrelationMean[x[[1]]$RegulonSize!=1])
	})

Mom <- sapply(correlations, function(x){ # Mean of medians 
	mean(x[[1]]$SpearmanCorrelationMedian[x[[1]]$RegulonSize!=1])
	})

mom <- sapply(correlations, function(x){ # Median of medians 
	median(x[[1]]$SpearmanCorrelationMedian[x[[1]]$RegulonSize!=1])
	})

SumTab <- data.frame(Name=names(sets), MeanOfMeans=MoM, MedianOfMeans=moM, MeanOfMedians=Mom, MedianOfMedians=mom)

write.table(SumTab, file = "Results/Tables/SummaryTable.txt", sep = "\t", row.names = FALSE)

############## Positive Control : All the genes are correlated, each with itself ##########

PosCon <- ToList(paste(dir, "processed_PositiveControl.txt", sep=""))

PosCon <- Correlation(PosCon)

#Table
write.table(PosCon[[1]], file =paste("Results/Tables/PositiveControl_Selfcorrelations.txt", sep=""), sep = "\t", row.names = FALSE)

# Boxplots
wi <- 0.6*length(unique(PosCon[[2]]$Name)) + 0.6*3
pdf("Results/Plots/Boxplots_PositiveControl.pdf",  width=wi)
print(ggplot(PosCon[[2]], aes(x = Name, y = Cor, fill=Size)) + scale_color_gradient(low="blue", high="red") + geom_boxplot() + geom_point())
dev.off()

PosCon <- data.frame(category=factor("Gene"), cor=PosCon[[2]]$Cor)
# Histogram
pdf("Results/Plots/Histogram_PositiveControl.pdf")
ggplot(PosCon, aes(x = cor, fill= category)) + geom_histogram() + ggtitle("Positive Control")
dev.off()