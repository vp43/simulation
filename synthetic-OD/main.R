#! /usr/bin/env Rscript
remove(list=ls())
source("genRandDataset.R")

########################################################################
# Data Generator for Occupancy-Detection modelss             
# Inputs: nSites					                                   
#					nVisits (to each site)                             
#					idx (for ID of multiple generations)               
#					formula.type ("linear", "quadratic", "exponential")
# Returns: Saved files of training/validation/test sets.
#					 Each dataset contains 
#					 - occCovars: occupancy covariates 
#          - detCovars: detection covariates
#          - occProbs: occupancy probabilitis
#          - detProbs: detection probabilitis
#          - trueOcc: sampled true occupancy status
#          - detHist: detection history based on trueOcc and detProbs
# How to Run: run(nSites,nVisits,idx,write.flag,formula.type) 
#             e.g. run(300,3,1,TRUE,"linear") 
#             e.g. run(600,5,2,FALSE,"quadratic")     		    	     	                       
########################################################################
run <- function(nSites, nVisits, idx, write.flag, formula.type="linear") {
	# Generate and load datasets
	myData <- generateData(nSites, nVisits, idx, formula.type)
	if (write.flag) saveData(myData)
}

generateData <- function(nSites, nVisits, idx, formula.type) {
	random.flag = FALSE 
	# If random.flalg = TRUE, randomly select relevant covariates and formula type.
	# If random.flalg = FALSE, the relevant covars is all covariates, 
	
	# Create directories
	dir_path  = "data/"
	if (!dir.exists(dir_path)) {dir.create(dir_path)}
	dir_path = paste(dir_path,nSites,"x",nVisits,"/", sep="")
	if (!dir.exists(dir_path)) {dir.create(dir_path)}
	dir_path = paste(dir_path,formula.type,"/", sep="")
	if (!dir.exists(dir_path)) {dir.create(dir_path)}
	save.path = paste(dir_path,idx,"/", sep="")
	print(save.path)
	if (!dir.exists(save.path)) {dir.create(save.path)}
	else stop("There is already data in it. Try it with a new version id.")
	
	# Generate covariates randomly
	randCovars = generateCovariates(nSites,nVisits)
	
	if (formula.type == "linear") {formula.id = 0.1}
	else if (formula.type == "quadratic") {formula.id = 0.5}
	else {formula.id = 0.9}
		 
	# Generate new synthetic data with given formula.
	synData = generateSynData(randCovars$occCovars, randCovars$detCovars, 
								            randCovars$index, random.flag, formula.id)
	print(synData$occForm)	
	
	# Check if the generated distributions of occupancy and detection are reasonable.
	threshold = 30 # prefered std on occupancy probabilitis
	Pz_balance = sd(hist(synData$occProbs, plot=FALSE)$counts)
	#Pd_balance = sd(hist(synData$detProbs, plot=FALSE)$counts)
	while ( Pz_balance > threshold & min(synData$detProbs) > 0.1) {
		synData = genRandDataset(demoData$occCovars, demoData$detCovars, demoData$index, 
							     random.flag, formula.id)
		Pz_balance = sd(hist(synData$occProbs, plot=FALSE)$counts)
		cat(c(Pz_balance, Pd_balance, "\n"))
		cat(c(min(synData$detProbs), "\n"))
	}
	#hist(synData$occProbs)
	#hist(synData$detProbs)
	
	# Data separation into training/validation/test sets
	trainData = list()
	validateData = list()
	testData = list()
	totVisits = nSites * nVisits
	
	# "detHists"
	trainData$detHists  = synData$detHists[1:(totVisits/3)]
	validateData$detHists  = synData$detHists[((totVisits/3)+1):(2*totVisits/3)]
	testData$detHists  = synData$detHists[((2*totVisits/3)+1):totVisits]

	# "occCovars"
	trainData$occCovars  = synData$occCovars[1:(nSites/3),]
	validateData$occCovars  = synData$occCovars[((nSites/3)+1):(2*nSites/3),]
	testData$occCovars  = synData$occCovars[((2*nSites/3)+1):nSites,]
	
	# "detCovars"
	trainData$detCovars  = synData$detCovars[1:(totVisits/3),]
	validateData$detCovars  = synData$detCovars[((totVisits/3)+1):(2*totVisits/3),]
	testData$detCovars  = synData$detCovars[((2*totVisits/3)+1):totVisits,]
	
	# "index"
	trainData$index  = synData$index[1:(totVisits/3)]
	validateData$index  = trainData$index
	testData$index  = trainData$index
	
	# "occProbs"
	trainData$occProbs  = synData$occProbs[1:(nSites/3)]
	validateData$occProbs  = synData$occProbs[((nSites/3)+1):(2*nSites/3)]
	testData$occProbs  = synData$occProbs[((2*nSites/3)+1):nSites]
	
	# "detProbs"
	trainData$detProbs  = synData$detProbs[1:(totVisits/3)]
	validateData$detProbs  = synData$detProbs[((totVisits/3)+1):(2*totVisits/3)]
	testData$detProbs  = synData$detProbs[((2*totVisits/3)+1):totVisits]
	
	# "Ztrue"
	trainData$Ztrue  = synData$Ztrue[1:(nSites/3)]
	validateData$Ztrue  = synData$Ztrue[((nSites/3)+1):(2*nSites/3)]
	testData$Ztrue  = synData$Ztrue[((2*nSites/3)+1):nSites]
	
	# "occForm"
	trainData$occForm  = synData$occForm
	validateData$occForm  = synData$occForm
	testData$occForm  = synData$occForm
	
	# "occCoeffs"
	trainData$occCoeffs  = synData$occCoeffs
	validateData$occCoeffs  = synData$occCoeffs
	testData$occCoeffs  = synData$occCoeffs
	
	# "detForm"
	trainData$detForm  = synData$detForm
	validateData$detForm  = synData$detForm
	testData$detForm  = synData$detForm
	
	# "detCoeffs"
	trainData$detCoeffs  = synData$detCoeffs
	validateData$detCoeffs  = synData$detCoeffs
	testData$detCoeffs  = synData$detCoeffs
	
	return(list(synData=synData, trainData=trainData, validateData=validateData, 
	            testData=testData, save.path=save.path))
}

saveData <- function(myData) {
	# Save distributions of generated data
	png(paste(myData$save.path,"occProb_dist.png",sep=""))
	hist(myData$testData$occProbs, xlab="True Occupancy Probabilities", ylab="Frequency")
	dev.off()
	
	png(paste(myData$save.path,"detProb_dist.png",sep=""))
	hist(myData$testData$detProbs, xlab="True Detection Probabilities", ylab="Frequency")
	dev.off()
	
	png(paste(myData$save.path,"sampling_dist.png",sep=""))
	plot(myData$testData$occProbs, myData$testData$Ztrue, 
	     xlab="True Occupancy Probabilities", ylab="True Occupancy (0/1)")
	dev.off()
	
	# Save training/validation/test sets
	trainData <- myData$trainData
	validateData <- myData$validateData
	testData <- myData$testData
	n.site <- nrow(trainData$occCovars)
	n.obs.total <- length(trainData$detHists)

	# Concatenate training and validation data for full training set:
	fullTrain = list()
	fullTrain$occCovars = rbind(trainData$occCovars,validateData$occCovars)
	fullTrain$detCovars = rbind(trainData$detCovars,validateData$detCovars)
	fullTrain$occProbs = append(trainData$occProbs,validateData$occProbs)
	fullTrain$detProbs = append(trainData$detProbs,validateData$detProbs)
	fullTrain$detHists = append(trainData$detHists,validateData$detHists)
	fullTrain$index = rbind(trainData$index,(validateData$index+n.site))
	n.site.full <- nrow(fullTrain$occCovars)
	n.obs.total.full <- length(fullTrain$detHists)

	# Save these new datsets as csv files for python codes.
	saveNewDataset(trainData, "train", myData$save.path)
	saveNewDataset(validateData, "valid", myData$save.path)
	saveNewDataset(testData, "test", myData$save.path)
	saveNewDataset(fullTrain, "fullTrain", myData$save.path)	
}

saveNewDataset <- function(dataList, data.type, save.path)  {
	# Save an object in a RData format
	save(dataList, file = paste(save.path,data.type,"Data.RData",sep="") )
	# Save an object to seperate files for Python programs
	write.table(dataList$occCovars, paste(save.path,data.type,"_occCovars.csv",sep=""), 
	            row.names=F, col.names=F, sep=",")
	write.table(dataList$detCovars, paste(save.path,data.type,"_detCovars.csv",sep=""), 
	            row.names=F, col.names=F, sep=",")
	write.table(dataList$occProbs, paste(save.path,data.type,"_occProbs.csv",sep=""), 
	            row.names=F, col.names=F, sep=",")
	write.table(dataList$Ztrue, paste(save.path,data.type,"_trueOcc.csv",sep=""),
             	row.names=F, col.names=F, sep=",")
	write.table(dataList$detProbs, paste(save.path,data.type,"_detProbs.csv",sep=""),
            	row.names=F, col.names=F, sep=",")
	write.table(dataList$detHists, paste(save.path,data.type,"_detHists.csv",sep=""),
	            row.names=FALSE, col.names=FALSE, sep=",")
}