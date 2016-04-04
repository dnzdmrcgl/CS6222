library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(VariantAnnotation)
library(ggplot2)

readBedFile <- function(bedFilePath) {
    df <- read.table(bedFilePath, header=F, stringsAsFactors=F)
    names(df) <- c('chr','start','end')
    gr <- with(df, GRanges(chr, IRanges(start, end)))
    return(gr)
}

loadTrainingData <- function(parentDirectory = '../data/training/', directory = 'syn1') {
  dataDirectory <- paste0(parentDirectory, directory, '/')
  tumorBam <- NULL
  normalBam <- NULL
  
  freebayesFilePath <- paste0(dataDirectory, directory, '_', 'freebayes.vcf')
  mutectFilePath <- paste0(dataDirectory, directory, '_', 'mutect.vcf')
  vardictFilePath <- paste0(dataDirectory, directory, '_', 'vardict.vcf')
  varscanFilePath <- paste0(dataDirectory, directory, '_', 'varscan.vcf')
  truthFilePath <- paste0(dataDirectory, directory, '_', 'truth.bed')
  
  freebayesVCF <- rowRanges(readVcf(freebayesFilePath, 'hg19'))
  mutectVCF <- rowRanges(readVcf(mutectFilePath, 'hg19'))
  vardictVCF <- rowRanges(readVcf(vardictFilePath, 'hg19'))
  varscanVCF <- rowRanges(readVcf(varscanFilePath, 'hg19'))
  
  truthSet <- readBedFile(truthFilePath)
  
  if ( substr(directory, 1, 3) == 'syn' ) {
    normalBamFilePath <- paste0(dataDirectory, directory, '_', 'normal.bam')
    tumorBamFilePath <- paste0(dataDirectory, directory, '_', 'tumor.bam')
    normalBam <- readGAlignments(BamFile(normalBamFilePath))
    tumorBam <- readGAlignments(BamFile(tumorBamFilePath))
  }
  
  return(list(freebayesVCF = freebayesVCF, mutectVCF = mutectVCF, vardictVCF = vardictVCF, 
              varscanVCF = varscanVCF, truthSet = truthSet, normalBam = normalBam, tumorBam = tumorBam))
}

getAllPredictions <- function(sampleData) {
  allPredictions <- unique(c(sampleData$freebayesVCF, sampleData$mutectVCF, sampleData$vardictVCF, sampleData$varscanVCF))
  
  allPredictions$freebayesSNP <- overlapsAny(allPredictions, sampleData$freebayesVCF, type = 'equal')
  allPredictions$mutectSNP <- overlapsAny(allPredictions, sampleData$mutectVCF, type = 'equal')
  allPredictions$vardictSNP <- overlapsAny(allPredictions, sampleData$vardictVCF, type = 'equal')
  allPredictions$varscanSNP <- overlapsAny(allPredictions, sampleData$varscanVCF, type = 'equal')
  
  allPredictions$truthTable <- overlapsAny(allPredictions, sampleData$truthSet, type = 'equal')
  
  return(allPredictions)
}

getRandomSample <- function(syn1.predictions, percentage) {
  index1 <- sample(1:length(syn1.predictions), as.integer(length(syn1.predictions) * percentage), replace = TRUE)
  index2 <- sample(1:length(syn1.predictions), as.integer(length(syn1.predictions) * percentage), replace = TRUE)
  index3 <- sample(1:length(syn1.predictions), as.integer(length(syn1.predictions) * percentage), replace = TRUE)
  index4 <- sample(1:length(syn1.predictions), as.integer(length(syn1.predictions) * percentage), replace = TRUE)
  
  return(list(index1 = index1, index2 = index2, index3 = index3, index4 = index4))
}

majorityVoting <- function(n, trueCounts, falseCounts, truthValues) {
  if(trueCounts[n] > falseCounts[n]) {
    return(1)
  } else if (trueCounts[n] == falseCounts[n]) {
    return(as.numeric(sample(truthValues[n,], 1)))
  } else {
    return (0)
  }
}

getPrediction <- function(query, sample.predictions, indexList) {
  overlapTable <- findOverlaps(query, sample.predictions, type = 'equal') 
  if (length(overlapTable) == 0) {
    return (0)
  }
  
  queryIndex <- subjectHits(overlapTable)
  
  table.index1 <- table(indexList$index1)
  table.index2 <- table(indexList$index2)
  table.index3 <- table(indexList$index3)
  table.index4 <- table(indexList$index4)
  
  times1 <- table.index1[match(queryIndex, names(table.index1))]
  times2 <- table.index2[match(queryIndex, names(table.index2))]
  times3 <- table.index3[match(queryIndex, names(table.index3))]
  times4 <- table.index4[match(queryIndex, names(table.index4))]
  
  times1[is.na(times1)] <- 0
  times2[is.na(times2)] <- 0
  times3[is.na(times3)] <- 0
  times4[is.na(times4)] <- 0
  
  times <- matrix(c(times1, times2, times3, times4), nrow = length(times1))
  truthValues <- as.matrix(mcols(sample.predictions[queryIndex])[,c('freebayesSNP', 'mutectSNP', 'vardictSNP', 'varscanSNP')])
  
  result <- times * truthValues
  trueCounts <- rowSums(result)
  falseCounts <- rowSums(times) - trueCounts
  
  pred <- sapply(1:length(trueCounts), majorityVoting, trueCounts, falseCounts, truthValues)
  return(pred)
}

getGroundTruth <- function(query, sample.predictions) {
  overlapTable <- findOverlaps(query, sample.predictions, type = 'equal')
  return(as.numeric(sample.predictions[subjectHits(overlapTable)]$truthTable))
}

evaluatePerformance <- function(predictions, labels) {
  true.positive <- sum((predictions + labels) == 2)
  test.outcome.positive <- sum(predictions)
  condition.positive <- sum(labels)
  
  precision <- true.positive / test.outcome.positive
  recall <- true.positive / condition.positive
  f1.score <- (2 * precision * recall) / (precision + recall)
  return(c(f1.score = f1.score, precision = precision, recall = recall))
}

multipleSamplingRun <- function(query, sample.predictions, numberOfSampling, samplingPercentage) {
  predictions <- rep(0, length(query))
  truePredictions <- getGroundTruth(query, sample.predictions)
  majorityVote <- as.integer(numberOfSampling / 2)
  for(i in 1:numberOfSampling) {
    sample.indexList <- getRandomSample(sample.predictions, samplingPercentage)
    predictions <- predictions + getPrediction(query, sample.predictions, sample.indexList)
  }
  predictions <- as.numeric(predictions > majorityVote)
  performance <- evaluatePerformance(predictions, truePredictions)
  return(performance)
}

getPerformance <- function(sample, sampleName, samplingRates, samplingNumber) {
  numberOfSampling <- samplingNumber[[1]]
  sampling.1 <- sapply(samplingRates, function(p) {multipleSamplingRun(sample, sample, numberOfSampling, samplingPercentage = p)})
  
  numberOfSampling <- samplingNumber[[2]]
  sampling.2 <- sapply(samplingRates, function(p) {multipleSamplingRun(sample, sample, numberOfSampling, samplingPercentage = p)})
  
  numberOfSampling <- samplingNumber[[3]]
  sampling.3 <- sapply(samplingRates, function(p) {multipleSamplingRun(sample, sample, numberOfSampling, samplingPercentage = p)})
  
  numberOfSampling <- samplingNumber[[4]]
  sampling.4 <- sapply(samplingRates, function(p) {multipleSamplingRun(sample, sample, numberOfSampling, samplingPercentage = p)})
  
  numberOfSampling <- samplingNumber[[5]]
  sampling.5 <- sapply(samplingRates, function(p) {multipleSamplingRun(sample, sample, numberOfSampling, samplingPercentage = p)})
  
  numberOfSampling <- samplingNumber[[6]]
  sampling.6 <- sapply(samplingRates, function(p) {multipleSamplingRun(sample, sample, numberOfSampling, samplingPercentage = p)})
  
  f1.scores <- c(sampling.1[1,], sampling.2[1,], sampling.3[1,], sampling.4[1,], sampling.5[1,], sampling.6[1,])
  precision <- c(sampling.1[2,], sampling.2[2,], sampling.3[2,], sampling.4[2,], sampling.5[2,], sampling.6[2,])
  recall <- c(sampling.1[3,], sampling.2[3,], sampling.3[3,], sampling.4[3,], sampling.5[3,], sampling.6[3,])
  samplingRatesAll <- rep(samplingRates, length(samplingNumber))
  samplingNumberAll <- rep(samplingNumber, rep(length(samplingRates), length(samplingNumber)))
  sampleNameAll <- rep(sampleName, length(samplingNumberAll))
  
  sample.performance <- data.frame(f1.scores, precision, recall, samplingRatesAll, samplingNumberAll, sampleNameAll)
  return(sample.performance)
}

plotPerformancePDF <- function(sample.performance, sampleName) {
  p <- ggplot(data = sample.performance, aes(y = f1.scores, x = samplingRatesAll, group = samplingNumberAll, color = samplingNumberAll))
  p <- p + geom_line() + geom_point() + scale_color_gradientn(colours = rainbow(6))
  pdf(paste0('./plots/', sampleName, '-f1score.pdf'), width = 15, height = 15)
  print(p)
  dev.off()
  
  p <- ggplot(data = sample.performance, aes(y = precision, x = samplingRatesAll, group = samplingNumberAll, color = samplingNumberAll))
  p <- p + geom_line() + geom_point() + scale_color_gradientn(colours = rainbow(6))
  pdf(paste0('./plots/', sampleName, '-precision.pdf'), width = 15, height = 15)
  print(p)
  dev.off()
  
  p <- ggplot(data = sample.performance, aes(y = recall, x = samplingRatesAll, group = samplingNumberAll, color = samplingNumberAll))
  p <- p + geom_line() + geom_point() + scale_color_gradientn(colours = rainbow(6))
  pdf(paste0('./plots/', sampleName, '-recall.pdf'), width = 15, height = 15)
  print(p)
  dev.off()
}

plotPerformance <- function(sample.performance, sampleName) {
  p <- ggplot(data = sample.performance, aes(y = f1.scores, x = samplingRatesAll, group = samplingNumberAll, color = samplingNumberAll))
  p <- p + geom_line() + geom_point() + scale_color_gradientn(colours = rainbow(6))
  jpeg(paste0('./plots/', sampleName, '-f1score.jpeg'), width = 450, height = 450)
  print(p)
  dev.off()
  
  p <- ggplot(data = sample.performance, aes(y = precision, x = samplingRatesAll, group = samplingNumberAll, color = samplingNumberAll))
  p <- p + geom_line() + geom_point() + scale_color_gradientn(colours = rainbow(6))
  jpeg(paste0('./plots/', sampleName, '-precision.jpeg'), width = 450, height = 450)
  print(p)
  dev.off()
  
  p <- ggplot(data = sample.performance, aes(y = recall, x = samplingRatesAll, group = samplingNumberAll, color = samplingNumberAll))
  p <- p + geom_line() + geom_point() + scale_color_gradientn(colours = rainbow(6))
  jpeg(paste0('./plots/', sampleName, '-recall.jpeg'), width = 450, height = 450)
  print(p)
  dev.off()
}

# Load the data from the specified input folder
# TODO: Since real data bam files will be added, update function so that it loads those bam files as well
syn1.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'syn1')
syn2.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'syn2')
syn3.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'syn3')
syn4.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'syn4')
syn5.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'syn5')
real1.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'real1')

# Create the table structure that contains the SNP calls for the union of all candidate SNPs
syn1.predictions <- getAllPredictions(syn1.data)
syn2.predictions <- getAllPredictions(syn2.data)
syn3.predictions <- getAllPredictions(syn3.data)
syn4.predictions <- getAllPredictions(syn4.data)
syn5.predictions <- getAllPredictions(syn5.data)
real1.predictions <- getAllPredictions(real1.data)

# Sample run with specified sampling number and rate for all samples
numberOfSampling <- 9
samplingPercentage <- 0.6

multipleSamplingRun(syn1.predictions, syn1.predictions, numberOfSampling, samplingPercentage)
multipleSamplingRun(syn2.predictions, syn2.predictions, numberOfSampling, samplingPercentage)
multipleSamplingRun(syn3.predictions, syn3.predictions, numberOfSampling, samplingPercentage)
multipleSamplingRun(syn4.predictions, syn4.predictions, numberOfSampling, samplingPercentage)
multipleSamplingRun(syn5.predictions, syn5.predictions, numberOfSampling, samplingPercentage)
multipleSamplingRun(real1.predictions, real1.predictions, numberOfSampling, samplingPercentage)

# Calculate the f1 score, precision and recall for each classifier on each sample data

evaluatePerformance(syn1.predictions$freebayesSNP, syn1.predictions$truthTable)
evaluatePerformance(syn1.predictions$mutectSNP, syn1.predictions$truthTable)
evaluatePerformance(syn1.predictions$vardictSNP, syn1.predictions$truthTable)
evaluatePerformance(syn1.predictions$varscanSNP, syn1.predictions$truthTable)

evaluatePerformance(syn2.predictions$freebayesSNP, syn2.predictions$truthTable)
evaluatePerformance(syn2.predictions$mutectSNP, syn2.predictions$truthTable)
evaluatePerformance(syn2.predictions$vardictSNP, syn2.predictions$truthTable)
evaluatePerformance(syn2.predictions$varscanSNP, syn2.predictions$truthTable)

evaluatePerformance(syn3.predictions$freebayesSNP, syn3.predictions$truthTable)
evaluatePerformance(syn3.predictions$mutectSNP, syn3.predictions$truthTable)
evaluatePerformance(syn3.predictions$vardictSNP, syn3.predictions$truthTable)
evaluatePerformance(syn3.predictions$varscanSNP, syn3.predictions$truthTable)

evaluatePerformance(syn4.predictions$freebayesSNP, syn4.predictions$truthTable)
evaluatePerformance(syn4.predictions$mutectSNP, syn4.predictions$truthTable)
evaluatePerformance(syn4.predictions$vardictSNP, syn4.predictions$truthTable)
evaluatePerformance(syn4.predictions$varscanSNP, syn4.predictions$truthTable)

evaluatePerformance(syn5.predictions$freebayesSNP, syn5.predictions$truthTable)
evaluatePerformance(syn5.predictions$mutectSNP, syn5.predictions$truthTable)
evaluatePerformance(syn5.predictions$vardictSNP, syn5.predictions$truthTable)
evaluatePerformance(syn5.predictions$varscanSNP, syn5.predictions$truthTable)

evaluatePerformance(real1.predictions$freebayesSNP, real1.predictions$truthTable)
evaluatePerformance(real1.predictions$mutectSNP, real1.predictions$truthTable)
evaluatePerformance(real1.predictions$vardictSNP, real1.predictions$truthTable)
evaluatePerformance(real1.predictions$varscanSNP, real1.predictions$truthTable)

#################################################################################
# Calculate the f1 scores, precision and recall for different samling rates and sampling times
# Store the calculated values in the performance folder in seralized manner
#################################################################################
samplingRates <- seq(from = 0.05, to = 0.95, by = 0.05)
samplingNumber <- seq(from = 1, to = 11, by = 2)

dir.create(file.path('./performance/'), showWarnings = FALSE)
sampleName <- 'syn1'
syn1.performance <- getPerformance(syn1.predictions, sampleName, samplingRates, samplingNumber)
saveRDS(syn1.performance, file = './performance/syn1Performace.rds')
sampleName <- 'syn2'
syn2.performance <- getPerformance(syn2.predictions, sampleName, samplingRates, samplingNumber)
saveRDS(syn2.performance, file = './performance/syn2Performace.rds')
sampleName <- 'syn3'
syn3.performance <- getPerformance(syn3.predictions, sampleName, samplingRates, samplingNumber)
saveRDS(syn3.performance, file = './performance/syn3Performace.rds')
sampleName <- 'syn4'
syn4.performance <- getPerformance(syn4.predictions, sampleName, samplingRates, samplingNumber)
saveRDS(syn4.performance, file = './performance/syn4Performace.rds')
sampleName <- 'syn5'
syn5.performance <- getPerformance(syn5.predictions, sampleName, samplingRates, samplingNumber)
saveRDS(syn5.performance, file = './performance/syn5Performace.rds')
sampleName <- 'real1'
real1.performance <- getPerformance(real1.predictions, sampleName, samplingRates, samplingNumber)
saveRDS(real1.performance, file = './performance/real1Performace.rds')


###########################
# The first part is to generate the predictions and caldulate the f1 scores, precision & recall
# Code below just loads the already calculated and stored scores and plots the graphs
###########################

syn1.performance <- readRDS('./performance/syn1Performace.rds')
syn2.performance <- readRDS('./performance/syn2Performace.rds')
syn3.performance <- readRDS('./performance/syn3Performace.rds')
syn4.performance <- readRDS('./performance/syn4Performace.rds')
syn5.performance <- readRDS('./performance/syn5Performace.rds')
real1.performance <- readRDS('./performance/real1Performace.rds')

dir.create(file.path('./plots'), showWarnings = FALSE)
plotPerformance(syn1.performance, 'syn1')
plotPerformance(syn2.performance, 'syn2')
plotPerformance(syn3.performance, 'syn3')
plotPerformance(syn4.performance, 'syn4')
plotPerformance(syn5.performance, 'syn5')
plotPerformance(real1.performance, 'real1')






