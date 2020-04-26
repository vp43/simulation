logistic <- function(x) {return(1 / (1 + exp(-x)))}

generateCovariates = function(nSites,nVisits,nOccCovars=10,nDetCovars=10) {
  # The occupancy covariates are a matrix of dimension nSites x
  #  nOccCovars.  Since we are making up this example, we will use
  #  random normal numbers as occupancy covariates.
  oC = rnorm(nSites*nOccCovars)
  dim(oC) = c(nSites,nOccCovars)

  # The detection covariates are a matrix of dimension nSites*nVisits
  #  x nDetCovars.  This means there is a row for every observation.
  #  Again, we synthesize random normal numbers for this matrix.
  dC = rnorm(nSites*nVisits*nDetCovars, 0, 0.5)
  dim(dC) = c(nSites*nVisits,nDetCovars)

  index = array(0,c(nSites*nVisits,1))
  idx = 1
  for (i in 1:nSites) {	# for each site i
    for (t in 1:nVisits) {  # for each visit t to site i
      # set up the index value to indicate site i for this observation:
      index[idx] = i
      # increment idx and move to the next visit
      idx = idx + 1
    } # t
  } # i
  return(list(occCovars=oC, detCovars=dC, index=index))
}

generateSynData = function(occCovs,detCovs,index,random.flag,formula.type) {
	verbose = FALSE
	occCovs = cbind(INT=1,occCovs)
	detCovs = cbind(INT=0,detCovs)

# input occCovs and detCovs should include an intercept column in position 1
	nSites = dim(occCovs)[1]
	nOccCovs = dim(occCovs)[2]-1
	occCovNames = names(occCovs)
	nVisits = dim(detCovs)[1]
	nDetCovs = dim(detCovs)[2]-1
	detCovNames = names(detCovs)

# OCCUPANCY #
	# num. relevant covars
	if (random.flag) nRelOcc = sample(0:nOccCovs)[1]
	else nRelOcc = nOccCovs

	occFormula = "occ ~ INT"
	occCoeffs = c((rbeta(1,0.5,0.5)-0.5)*4)
	#occCoeffs = rnorm(1,0.5,1) 
	occVals = occCoeffs[1]*occCovs[,1]
	if (nRelOcc>0) {
		whichRelOcc = sort(sample(2:(nOccCovs+1))[1:nRelOcc])
		# covar relationships
		for (o in whichRelOcc) {
			thisCoeff = (rbeta(1,0.5,0.5)-0.5)*4
			#thisCoeff = rnorm(1,0,1)
			occCoeffs = c(occCoeffs,thisCoeff)
			if (is.numeric(occCovs[,o])) {
				if (random.flag) type = runif(1)
				else type = formula.type
				
				if (type<0.33) { # linear
					occFormula = paste(occFormula," + ",occCovNames[o],sep="")
					occVals = occVals + thisCoeff*occCovs[,o]
				} else if (type>=0.33 & type<0.66) { # quadratic
					occFormula = paste(occFormula," + ",occCovNames[o],"^2",sep="")	
					occVals = occVals + thisCoeff*occCovs[,o]^2
				} else { # exponential
					occFormula = paste(occFormula," + exp(-",occCovNames[o],")",sep="")
					occVals = occVals + thisCoeff*exp(-occCovs[,o])
				} # end if relationship types
			
			} else if (is.factor(occCovs[,o])) {
				vals = levels(occCovs[,o])
				nVals = length(vals)
				relVals = sample(1:(nVals-1))[1:sample(1:(nVals-1))[1]]
				occFormula = paste(occFormula," + ",occCovNames[o],paste(relVals,collapse=""),sep="")
				occVals = occVals + thisCoeff*(is.element(occCovs[,o],relVals))
			} else {
				print(paste("What kind of variable is occ cov. ",o,sep=""))
			} # end if variable types
		} # o

		if (random.flag) {
			if (nRelOcc>1) {
				# choose some interactions
				nPotIntsOcc = sample(0:(nRelOcc-1))[1]
				nIntsOcc = nPotIntsOcc
				if (nPotIntsOcc>0) {
					firstTerms = whichRelOcc
					for (i in 1:nPotIntsOcc) {
						thisFirst = sample(firstTerms)[1]
						firstTerms = firstTerms[-which(firstTerms==thisFirst)]
						thisSecond = sample(firstTerms)[1]
						if (!is.numeric(occCovs[,thisFirst]) | 
						!is.numeric(occCovs[,thisSecond])) {
						# don't deal with factors in interaction terms for now
						nIntsOcc = nIntsOcc-1
						next
					}
					thisCoeff = (rbeta(1,0.5,0.5)-0.5)*4
					occCoeffs = c(occCoeffs,thisCoeff)
					occFormula = paste(occFormula," + ",occCovNames[thisFirst],"*",occCovNames[thisSecond],sep="")
					occVals = occVals + thisCoeff*occCovs[,thisFirst]*occCovs[,thisSecond]
					} # i
				} # endif nPotIntsOcc>0
			} # endif nRelOcc>1	
		}	  
	} # endif nRelOcc>0
	
	if (verbose) {
		nRelOcc
		nIntsOcc
		occFormula
		layout(matrix(c(1,2),1,2))
		hist(occVals)
		hist(logistic(occVals))
		mean(logistic(occVals))
	}

# DETECTION #
	# num. relevant covars
	if (random.flag) nRelDet = sample(0:nDetCovs)[1]
	else nRelDet = nDetCovs
	
	detFormula = "det ~ INT"
	detCoeffs = c((rbeta(1,0.5,0.5)-0.5)*4)
	#detCoeffs = c(rnorm(1,0.5,1))
	detVals = detCoeffs[1]*occCovs[,1]
	if (nRelDet>0) {
	  whichRelDet = sort(sample(2:(nDetCovs+1))[1:nRelDet])
	  # covar relationships
	  for (d in whichRelDet) {
			thisCoeff = (rbeta(1,0.5,0.5)-0.5)*2
			#thisCoeff = rnorm(1,0.5,1)
			detCoeffs = c(detCoeffs,thisCoeff)
			if (is.numeric(detCovs[,d])) {
				if (random.flag) type = runif(1)
				else type = formula.type
				
				if (type<0.33) { # linear
					detFormula = paste(detFormula," + ",detCovNames[d],sep="")
					detVals = detVals + thisCoeff*detCovs[,d]
				} else if (type>=0.33 & type<0.66) { # quadratic
					detFormula = paste(detFormula," + ",detCovNames[d],"^2",sep="")	
					detVals = detVals + thisCoeff*detCovs[,d]^2
				} else { # exponential
					detFormula = paste(detFormula," + exp(-",detCovNames[d],")",sep="")
					detVals = detVals + thisCoeff*exp(-detCovs[,d])
				} # end if relationship types
			
			} else if (is.factor(detCovs[,d])) {
				vals = levels(detCovs[,d])
				nVals = length(vals)
				relVals = sample(1:(nVals-1))[1:sample(1:(nVals-1))[1]]
				detFormula = paste(detFormula," + ",detCovNames[d],paste(relVals,collapse=""),sep="")
				detVals = detVals + thisCoeff*(is.element(detCovs[,d],relVals))
			} else {
				print(paste("What kind of variable is det cov. ",d,sep=""))
			} # end if variable types
		} # d

		if (random.flag) {
			if (nRelDet>1) {
			 # choose some interactions
			 nPotIntsDet = sample(0:(nRelDet-1))[1]
			 nIntsDet = nPotIntsDet
			 if (nPotIntsDet>0) {
				firstTerms = whichRelDet
				for (i in 1:nPotIntsDet) {
					thisFirst = sample(firstTerms)[1]
				firstTerms = firstTerms[-which(firstTerms==thisFirst)]
				thisSecond = sample(firstTerms)[1]
				if (!is.numeric(detCovs[,thisFirst]) | 
					 !is.numeric(detCovs[,thisSecond])) {
					 # don't deal with factors in interaction terms for now
					 nIntsDet = nIntsDet-1
					 next
					 } 
						 thisCoeff = (rbeta(1,0.5,0.5)-0.5)*4
						 detCoeffs = c(detCoeffs,thisCoeff)
					 detFormula = paste(detFormula," + ",detCovNames[thisFirst],"*",detCovNames[thisSecond],sep="")
					 detVals = detVals + thisCoeff*detCovs[,thisFirst]*detCovs[,thisSecond]
				 } # i
			 } #endif nPotIntsDet>0
			} # endif nRelDet>1	
		}	  
	} # endif nRelDet>0
	
	if (verbose) {
		nRelDet
		nIntsDet
		detFormula
		layout(matrix(c(1,2),1,2))
		hist(detVals)
		hist(logistic(detVals))
		mean(logistic(detVals))
	}

# SIMULATE #
	occProbs = logistic(occVals)
	detProbs = logistic(detVals)
	Ztrue = as.numeric(runif(nSites)<occProbs)
	detHists = Ztrue[index]*(runif(nVisits)<detProbs)

# BUNDLE AND RETURN #
	simData = list(detHists=detHists,occCovars=occCovs[,2:(nOccCovs+1)],detCovars=detCovs[,2:(nOccCovs+1)],
			index=index,occProbs=occProbs,detProbs=detProbs,Ztrue=Ztrue,
			occForm=occFormula,occCoeffs=occCoeffs,detForm=detFormula,detCoeffs=detCoeffs)	
				
	return(simData)

} # end fn genRandDataset
