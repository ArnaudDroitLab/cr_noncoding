#####################################################################
###Function: getIntersection()
###Author: Luigi Marchionni
###The real function that computes the overlap between ordered lists of identifiers
###Note that the object x is a list containing to vectors of identifiers
###This function can compute the overlap between equal ranks or equally statistics

getIntersection <- function(x, size, method, decreasing) {

	##Extract the two vectors to be compared
	v1 <- x[[1]]
	v2 <- x[[2]]

	##check length
	if (length(v1)!=length(v2)) {
		stop("Vectors of identifiers to be compared must have equal length")
	} else {
		vOut <- length(v1)
	}

	##Using the increasing statistics methods
	if (method == "equalStat") {
		rStat <- c(v1, v2)

		##for decreasing ranking statistics
		if (decreasing) {
			rStat <- seq(from=max(rStat, na.rm=TRUE), to=min(rStat, na.rm=TRUE), length.out=vOut)

			##define the i and j indexes matrix
			rIndex <- sapply(rStat, function(rStat, v1, v2) {
				c(i=length(v1[abs(v1) >= rStat]), j=length(v2[abs(v2) >= rStat]))
			}, v1=v1, v2=v2)
		}
		##for increasing ranking statistics
		else {
			rStat <- seq(from=min(rStat, na.rm=TRUE), to=max(rStat, na.rm=TRUE), length.out=vOut)
			##define the i and j indexes matrix
			rIndex <- sapply(rStat, function(rStat, v1, v2) {
				c(i=length(v1[abs(v1) <= rStat]), j=length(v2[abs(v2) <= rStat]))
			}, v1=v1, v2=v2)
		}

		##compute cat
		cat <- vector("numeric", size)
		deg <- vector("numeric", size)
		for (z in 1:size) {
			ij <- rIndex[,z]
			int <- length(intersect(names(v1[1:ij["i"]]), names(v2[1:ij["j"]])))
			uni <- length(unique(union(names(v1[1:ij["i"]]), names(v2[1:ij["j"]]))))
			cat[z] <- int/uni
			deg[z] <- uni
		}
		out <- list(cat=cat, deg=deg)
	}

	##Using the equal Ranks method
	if (method == "equalRank") {
		cat <- vector("numeric", size)
		v1Names <- names(v1)
        v1Rank <- order(as.numeric(names(v1)))
        v2Names <- names(v2)
        v2Rank <- order(as.numeric(names(v2)))
        
        lastOverlap <- 0
        for (z in 1:size) {
#			cat[z] <- length(intersect(names(v1)[1:z], names(v2)[1:z]))/z
            if(v1Names[z]==v2Names[z]) {
                lastOverlap <- lastOverlap + 1
            } else {
                if(v2Rank[as.numeric(v1Names[z])] <= z) {
                    lastOverlap <- lastOverlap + 1
                }
                if(v1Rank[as.numeric(v2Names[z])] <= z) {
                    lastOverlap <- lastOverlap + 1
                }
            }
            cat[z] <- lastOverlap / z
            
		}
		out <- list(cat=cat)
	}

	##return
	return(out)
}


#####################################################################
#####################################################################
###Function: computeCat()
###Author: Luigi Marchionni
###This gets the the statistics and the identifiers out of dataframes or matrices
###and computes the  overlap after some validity checks

computeCat <- function (data,
			size=nrow(data),
			idCol=1, ref,
			method = c("equalRank", "equalStat"),
			decreasing = TRUE) {

	##check size
	if (size>nrow(data)) {
		size <- nrow(data)
	}

	##check data: colnames should be unique
	if (any(duplicated(colnames(data)))) {
		stop("Duplicated colnames(data) are not allowed")
	}

	##Check idCol: it  should be turned to a numeric if character
	if ( is.character(idCol) ) {
		if ( idCol %in% colnames(data) ) {
			idCol <- which( colnames(data) %in% idCol )
		}
	}

	##Check idCol: it should be a numeric corresponding to the first column if numeric
	if (is.numeric(idCol)) {
		if ( idCol != 1) {
			colOrder <- c(idCol, 1:ncol(data)[-idCol])
			data <- data[, colOrder]
		}
	} else {
		print("Invalid idCol argument! This should be a valid column name or index")
	}

	## prepare combinations combinations
	cData <- t(combn(colnames(data[,-idCol]),2))

	##filter if reference group is available
	if (!missing(ref)) {
		if (ref %in% colnames(data)) {
			refData <- apply(cData, 1, function(x,y) {any(x%in%y)}, y=ref)
			cData <- cData[refData,]
		} else {
			stop("Use a valid name for the reference group")
		}
	}

	##combine the identifiers ordered by the statistics
	tmp <- apply(cData, 1, function(x, y) {
		cnm <- colnames(y)
		##reorder identifiers and ranking statistics
		index1 <- which(cnm == x[1])
		index2 <- which(cnm == x[2])
		ord1 <- order(y[, index1], decreasing=decreasing)
		ord2 <- order(y[,  index2], decreasing=decreasing)
		X <- y[, index1][ord1]
		Y <- y[, index2][ord2]
		names(X) <- y[, idCol][ord1]
		names(Y) <- y[, idCol][ord2]
		##return
		ans <- list(X,Y)
	}, y=data)

	##create final names
	names(tmp) <- apply(cData, 1, function(x) {paste(x, collapse=".vs.")})

	##compute the overlap
	out <- lapply(tmp, getIntersection, size, method, decreasing)

	##return
	return(out)
}

#####################################################################
#####################################################################
###Function: calcHypPI()
###Author: Luigi Marchionni
###This function calculates the PI and probabilities using the hypergeometric distribution
###It is simulation based so it is very slow if the number of identifiers is very large
###This function will take more and more time when more and more genes are used:
###empirical evaluation showed that it takes 4 times more every doubling the genes

calcHypPI <- function(data,  expectedProp = 0.1, prob = c(0.999999,0.999,0.99,0.95) ) {

	## Reorder the percentiles and add the median
	prec <- sort(prob)
	otherprob <- 1-prob
	finalprob <- rep(NA, 2*length(prob))
	XX <- seq(1, length(finalprob), by=2)
	finalprob[XX] <- prob
	finalprob[XX+1] <- otherprob
	finalprob <- c(finalprob, 0.5)
	
	##How many features?
	nGenes <- nrow(data)

	###process expectedProp to define the number of white balls in the urn
	if (is.null(expectedProp)) {
		##if it is set to NULL by the user the number of white balls in the urn
		##will increase with the number of drawn balls (m == k)
		expected <- 1:nGenes
       } else {
	       ###it corresponds to the quantile
	       if (length(expectedProp) == 1) {
		       expected <- rep(round(nGenes*expectedProp), nGenes)
	       } else {
		       stop("Invalid 'expectedProp' argument: it should be a proportion or set equal to 'NULL'.")
	       }
       }

	##assemble numbers to be used with qhyper
	forQhyper <- cbind(mm=expected, nn=nGenes-expected, kk=1:nGenes)

	##Compute proportion for each quantile
	out <- t(apply(forQhyper, 1,
		      function(x, p) qhyper(p, x[1], x[2], x[3]), p=finalprob))/forQhyper[,1]
	colnames(out) <- format(finalprob)

	##Check if NA were generated
	if (any(is.na(out))) {
		stop("NaN were generated")
	} else {
		return(out)
	}
}


##################################################################
########################################################################
########################################################################
######Function: plotCat()
######Author: Luigi Marchionni
######This functions plots the computed CAT, generating the legend and the so on

plotCat <- function (catData,
		     whichToPlot = 1:length(catData),
		     preComputedPI,
		     size=500,
		     main="CAT-plot",
		     minYlim=0,
		     maxYlim=1,
		     col, pch, lty,
		     cex=1,
		     lwd=1,
		     spacePts=10,
		     cexPts=1,
		     legend=TRUE,
		     legendText,
		     where="center",
		     legCex=1,
		     plotLayout=layout(matrix(1:2, ncol = 2, byrow = TRUE), widths = c(0.7, 0.3)),
		     ...) {
	

	if ( any( ! whichToPlot %in% 1:length(catData) ) ) {
		stop("Argument 'WhichToPlot': incorrect indexes")
	} else {
		catData <- catData[whichToPlot]
	}

	##discover method used to generate catData
	if(length(catData[[1]]) == 1) {
		byRank=TRUE
		data <- sapply(catData, function(x) x[[1]] )
		if (missing(preComputedPI)) {
			ww <- matrix(NA, ncol=2, nrow=nrow(data))
			warning("No precomputed PI object was provided, hence PI will not be added to the plot")
		} else {
			ww <- preComputedPI
		}
	} else {
		byRank=FALSE
		data <- sapply(catData, function(x) x[[1]] )
		degData <- sapply(catData, function(x) x[[2]] )
		ww <- matrix(NA, ncol=2, nrow=nrow(data))
		if (!missing(preComputedPI)) {
			warning("Precomputed PI can be plotted only with CAT curves computed for equal ranks, not for equal statisitcs,\n hence no PI would be added to the plot")
		}
	}

	##set nr and nc
	nc <- ncol(data)
	nr <- nrow(data)

	##create color if needed
	if (missing(col)) {
		myCol <- sort(rainbow(nc))
	} else {
		myCol <- rep(col, length=nc)
		print("Recycling colors")
	}

	##generates points
	if (missing(pch)) {
		if (nc < 53) { myPch <- c(LETTERS, letters) }
		if (nc > 52 & nc < 121) { myPch <- c(1:25,33:127) }
		if (nc > 120) {
			myPch <- rep(c(1:25,33:127), length(nc))
			print("Recycling point types")
		}
	} else {
		myPch <- rep(pch, length=nc)
		print("Recycling point types")
	}

	##generates line types
	if (missing(lty)) {
		myLty <- 1:nc
	} else {
		myLty <- rep(lty, length=nc)
		print("Recycling line types")
	}

	##set maximal size to plot is more than number of row in data
	if (size > nr) {
		mySize <- nr
	} else {
		mySize <- size
	}

	##set the layout for the plot and the legend
	if (legend) {
		plotLayout
	}

	## set new par  parameters
	opar <- par(mar=c(4, 4, 3, 0))
	on.exit(par(opar))

	##start plotting
	plot(c(1:mySize, mySize:1), c(ww[1:mySize,1], ww[mySize:1,2]),
	     type="n",
	     ylim = c(minYlim, maxYlim), xlim = c(0, mySize),
 	     main = main, ylab = NA, xlab = NA,
	     cex.main = cex, cex.lab = cex, cex.axis = cex, axes = TRUE)

	##add polygon for cat curves by ranking
	if (byRank) {
		##prepare grey scale for PI if needed
		if ( ! all(is.na(ww)) ) {
			nGreys <- (ncol(ww)-1)/2
			myGreys <- grep("grey[23456789].",colors(),value=TRUE)
			greysSel <- ceiling(seq(from=1, to=length(myGreys), length.out=nGreys))
			ciCol <- rev(myGreys[greysSel])
		} else {
			ciCol <- NA
		}
			

		##create index to create polygons
		zz <- 1:(ncol(ww)-1)
		zz <- zz[matchBox:::is.odd(zz)]
		##print(zz)

		##plot polygons
		for (i in 1:length(zz) ) {
			polygon(c(1:mySize, mySize:1),
				c(ww[1:mySize, zz[i]], ww[mySize:1, zz[i]+1]),
				col=ciCol[i], border=NA)
		}
		##add 0.5 line
		points(ww[1:mySize, ncol(ww)], type = "l")
	}

	##plot and add points
	for (i in 1:nc) {
		points(data[1:mySize,i], type="l",
		       col=myCol[i], pch=myPch[i], lty=myLty[i], lwd=lwd)
		z <- seq(from=1+i, to=mySize, by=spacePts)
		points(z, data[z,i], col=myCol[i], pch=myPch[i], cex=cexPts)
	}

	##add text on axes
	mtext("Common Proportion", side=2, font=2, line=1+cex, cex=cex)
	mtext("List Size", side=1, font=2, line=1+cex, cex=cex)

	##add legend
	if (legend) {
		par(mar=c(3,0.5,3,0.5))
		plot.new()
		if (spacePts==0) {
			pch=NULL
		}
		##create legend text
		if (missing(legendText)) {
			legendText <- colnames(data)[1:nc]
		}
		legend(where, legend=legendText, title="Cat plot",
		       col=myCol, cex=legCex, lty=myLty, lwd=lwd, pch=myPch, ...)
	}
}


