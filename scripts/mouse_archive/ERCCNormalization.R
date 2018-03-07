library(ggplot2)
library(scales)

# Read in the ERCC annotation, including the concentrations of each ERCC.
erccAnnotation <- read.table("Input/ERCC_Controls_Analysis.txt", sep="\t", header=TRUE)
erccConcentrations <- erccAnnotation[,c(2,4)]
colnames(erccConcentrations) <- c("ERCC", "Concentration")

for(f in c("All.FPKM", "All.FPKM.nolength", "All.FPKM.nolength2")) {

# Read in FPKM data
expressionData <- read.table(paste("Output/", f, sep=""), sep="\t", header=TRUE, row.name="tracking_id")

# Read the ERCC -> Transcript ID map
erccTranscripts <- read.table("Output/ERCC-Transcript.map", col.names=c("ERCC", "Transcript"), as.is=TRUE)

# Summarize split ERCC transcripts
# Determine which ERCCs have more than one transcript
duplicateERCCs <- which(duplicated(erccTranscripts$ERCC))

# From the bottom-up, fold the duplicated entry into the previous one.
# If there is more than one duplicate, then duplicate x will fold into duplicate x-1,
# which will fold into x-2, up until the original entry is reached.
statusColumns <- grep("status", colnames(expressionData))
for(duplicateERCC in sort(duplicateERCCs, decreasing=TRUE)) {
    originalTranscript <- erccTranscripts$Transcript[duplicateERCC-1]
    duplicateTranscript <- erccTranscripts$Transcript[duplicateERCC]
    numericColumns <- setdiff(9:ncol(expressionData), statusColumns)
    expressionData[originalTranscript, numericColumns] <- expressionData[originalTranscript, numericColumns] + expressionData[duplicateTranscript, numericColumns]
}

# Remove duplicated entries
if(length(duplicateERCCs) > 0) {
    expressionData <- expressionData[-which(rownames(expressionData) %in% erccTranscripts$Transcript[duplicateERCCs]),]
    erccTranscripts <- erccTranscripts[-duplicateERCCs,]
}

# Choose which measure (standard FPKM, low-confidence-interval FPKM) to use.
for( fpkm in c("FPKM_conf_lo", "FPKM.")) {
whichFPKM <- fpkm
#whichFPKM <- "FPKM_conf_lo"
statusData <- expressionData[,statusColumns]
hiDataRows <- apply(statusData=="HIDATA", 1, any)


#whichFPKM <- "FPKM."

# Subset only the FPKM data.
expressionFPKM <- expressionData[,grep(whichFPKM, colnames(expressionData), value=TRUE, fixed=TRUE)]

# Replace FPKM values in HIDATA rows by NA.
if(whichFPKM=="FPKM_conf_lo") {
    expressionFPKM[hiDataRows,] <- NA
}

# Build an ERCC-Concentration-FPKM table, as well as a non-ERCC expression table.
ERCCData <- cbind(ERCC=erccTranscripts$ERCC,
                  Concentration=erccConcentrations$Concentration[match(erccTranscripts$ERCC, erccConcentrations$ERCC)],
                  expressionFPKM[match(erccTranscripts$Transcript, rownames(expressionFPKM)),])
expressionNoERCC <- expressionFPKM[-which(rownames(expressionFPKM) %in% erccTranscripts$Transcript),]

# Get tissue/replicate information for each numeric column
sampleString <- gsub("coverage.", "", grep("coverage", colnames(expressionData), value=TRUE))
tissue <- gsub("(.*?)\\..", "\\1", sampleString, perl=TRUE)
repID <- gsub(".*?\\.(.)", "\\1", sampleString, perl=TRUE)

# Produce titration graphs for all replicates
outputFolder <- paste("Output/Analysis/", f, " - ", fpkm, "/", sep="");
dir.create(paste(outputFolder, "ERCC-Titration", sep=""), recursive=TRUE)
for(i in 3:ncol(ERCCData)) {
    erccDF <- data.frame(Concentration=log10(ERCCData$Concentration), FPKM=log10(ERCCData[,i]))
    graphTitle <- paste("ERCC Titration for ", tissue[i-2], "-", repID[i-2], sep="")
    g <- ggplot(erccDF, aes(x=Concentration, y=FPKM)) +
        geom_point() +
        labs(list(title=graphTitle, x="log10(ERCC Concentration)", y="log10(FPKM)")) +
        theme(axis.text.x=element_text(hjust=1),
              panel.background = element_rect(fill='#FFFFFF'),
              axis.line=element_line(colour="black", size=1, linetype="solid"),
              axis.text = element_text(colour="black", size=12),
              panel.grid.major.y = element_line(colour="#F0F0F0", size=1),
              panel.grid.minor.y = element_line(colour="#F4F4F4", size=1))          
    ggsave(paste(outputFolder, "ERCC-Titration/", graphTitle, ".png", sep=""), plot=g)
}

# Produce across-replicates cat-plots
library(matchBox)
source("Scripts/catPlotFunctions.R")

dir.create(paste(outputFolder, "Replicate-CAT", sep=""), recursive=TRUE)
for(t in unique(tissue)) {
    # Keep only samples from the same tissue.
    dataSubset <- cbind(1:nrow(expressionNoERCC), expressionNoERCC[,which(tissue %in% t)])
    
    # If these is more than one replicate:
    if(ncol(dataSubset) > 2) {
        # Rename columns according to replicate ID to make following steps simpler
        colIDs <- c("Id", "A", "B", "C")
        colnames(dataSubset) <- colIDs[1:ncol(dataSubset)]
        
        # Set arbitrary limit for the number of elements to compare.
        catLimit <- 5000

        # Compute CAT overlaps.
        catResults <- computeCat(dataSubset, idCol=1, method="equalRank", decreasing=TRUE, size=catLimit)
        
        # Build up a "pair" vector indicating which comparison each line of the data-frame represents.
        nComb <- choose(ncol(dataSubset)-1, 2)
        repPairs <- c("A-B", "A-C", "B-C")
        pairVec <- c()
        for(i in 1:nComb) {
            pairVec <- c(pairVec, rep(repPairs[i], catLimit))
        }
        
        # Concatenate all CAT results (if more than one)
        overlapVec <- c()
        for(i in 1:nComb) {
            overlapVec <- c(overlapVec, catResults[[i]]$cat[1:catLimit])
        }
        
        # Build data-frame and plot results.
        catDF <- data.frame(Overlap=overlapVec*100, Pair=pairVec, ListSize=1:catLimit)
        graphTitle <- paste("Replicate CAT for", t)
        ggplot(catDF, aes(x=ListSize, y=Overlap, group=Pair, color=Pair)) +
            geom_line() +
            ylim(c(0,100)) + 
            labs(list(title=graphTitle, x="List size", y="Overlap")) +
            theme(axis.text.x=element_text(hjust=1),
                  panel.background = element_rect(fill='#FFFFFF'),
                  axis.line=element_line(colour="black", size=1, linetype="solid"),
                  axis.text = element_text(colour="black", size=12),
                  panel.grid.major.y = element_line(colour="#F0F0F0", size=1),
                  panel.grid.minor.y = element_line(colour="#F4F4F4", size=1))   
        ggsave(paste(outputFolder, "Replicate-CAT/", graphTitle, ".png", sep=""))
    }
}

# Normalize FPKM values against ERCC values
expressionNorm <- t(apply(expressionNoERCC, 1, "/", apply(ERCCData[,-(1:2)], 2, sum)))
write.table(data.frame(tracking_id=rownames(expressionNorm), expressionNorm), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE,
            file=paste(outputFolder, "All.Norm.FPKM", sep=""))

# Tissue order
tissueOrder <- c()
if(ncol(expressionNoERCC) < 30) {
    tissueOrder <- c("UF", "GV", "1c", "2c", "4c", "8c", "Morula", "Early_blast", "Blast", "Expand_blast")
} else {
    tissueOrder <- c("GV", "MII", "1c", "2c", "4c", "8c_tot", "8c_tard", "Morula", "Jeune_blast", "Blast_exp")
}

errorBarData <- function(x) {
    results <- c(mean(x)-sd(x), mean(x)+sd(x))
    names(results) <- c("ymin", "ymax")
    return(results)
}

errorBarDataMinZero <- function(x) {
    results <- c(max(mean(x)-sd(x), 0), mean(x)+sd(x))
    names(results) <- c("ymin", "ymax")
    return(results)
}

plotNormalizedData <- function(indices, title) {
    tissueFactor <- factor(tissue, levels=tissueOrder)
    summedData <- apply(expressionNorm[indices,], 2, sum, na.rm=TRUE)
    summedFrame <- data.frame(Total=summedData, Tissue=tissueFactor, Rep=repID)
    
    # FPKM Point plot
    ggplot(summedFrame, aes(y=Total, x=Tissue)) + 
        geom_point(shape=18, size=3, colour="red") + 
        stat_summary (fun.y = mean, geom="line", linetype="dashed", mapping = aes (group = 1)) +
        stat_summary (fun.data="errorBarData", geom="errorbar", width=0.25, mapping = aes (group = 1)) +
        labs(title=paste("Expression profile for ", title, sep="")) + 
        scale_y_continuous(name="Transcript FPKM / Total ERCC FPKM") +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(paste(outputFolder, "Expression Profiles/Expression profile for ", title, " - points.png", sep=""))
    
    # FPKM Bar graph
    ggplot(summedFrame, aes(y=Total, x=Tissue)) + 
        stat_summary (fun.data="errorBarDataMinZero", geom="errorbar", width=0.25, size=1, mapping = aes (group = 1)) +
        stat_summary (fun.y = mean, geom="bar", mapping = aes (group = 1), colour="black", fill="black") +
        labs(title=paste("Expression profile for ", title, sep="")) + 
        scale_y_continuous(name="Transcript FPKM / Total ERCC FPKM", expand=c(0,0)) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              panel.background = element_rect(fill='#FFFFFF'),
              axis.line=element_line(colour="black", size=1, linetype="solid"),
              axis.text = element_text(colour="black", size=12),
              panel.grid.major.y = element_line(colour="#F0F0F0", size=1),
              panel.grid.minor.y = element_line(colour="#F4F4F4", size=1))
    ggsave(paste(outputFolder, "Expression Profiles/Expression profile for ", title, " - bars.png", sep=""))
    
    # Proportion graph (Point)
    proportionData <- apply(expressionNorm[indices,], 2, sum, na.rm=TRUE) / apply(expressionNorm, 2, sum, na.rm=TRUE)
    proportionFrame <- data.frame(Total=proportionData, Tissue=tissueFactor, Rep=repID)
    
    ggplot(proportionFrame, aes(y=Total, x=Tissue)) + 
        geom_point(shape=18, size=3, colour="red") + 
        stat_summary (fun.y = mean, geom="line", linetype="dashed", mapping = aes (group = 1)) +
        stat_summary (fun.data="errorBarData", geom="errorbar", width=0.25, mapping = aes (group = 1)) +
        labs(title=paste("Proportion of total expressed transcripts for ", title, sep="")) + 
        scale_y_continuous(name="Transcript FPKM / Total ERCC FPKM") +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(paste(outputFolder, "Expression Profiles/Proportion of total expressed transcripts for ", title, " - points.png", sep=""))
    
    # Proportion graph (Bar)
    ggplot(proportionFrame, aes(y=Total, x=Tissue)) + 
        geom_point(shape=18, size=3, colour="red") + 
        stat_summary (fun.y = mean, geom="line", linetype="dashed", mapping = aes (group = 1)) +
        stat_summary (fun.data="errorBarData", geom="errorbar", width=0.25, mapping = aes (group = 1)) +
        labs(title=paste("Proportion of total expressed transcripts for ", title, sep="")) + 
        scale_y_continuous(name="Transcript FPKM / Total ERCC FPKM") +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(paste(outputFolder, "Expression Profiles/Proportion of total expressed transcripts for ", title, " - bars.png", sep=""))    
}

dir.create(paste(outputFolder, "Expression Profiles", sep=""), recursive=TRUE)
plotNormalizedData(TRUE, "All")

}
}
