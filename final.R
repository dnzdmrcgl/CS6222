library(Rcpp)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(VariantAnnotation)
library(ggplot2)

################################################################
############### VCF DATA LOAD AND PROCESSING ###################
################################################################

# Converts the bed file into Granges object 
readBedFile <- function(bedFilePath) {
  df <- read.table(bedFilePath, header=F, stringsAsFactors=F)
  names(df) <- c('chr','start','end')
  gr <- with(df, GRanges(chr, IRanges(start, end)))
  return(gr)
}

# Loads the vcf files from the specified directory
# For training data, loads the truth set bed file
loadTrainingData <- function(parentDirectory = '../data/training/', directory = 'syn1') {
  dataDirectory <- paste0(parentDirectory, directory, '/')
  
  freebayesFilePath <- paste0(dataDirectory, directory, '_', 'freebayes.vcf')
  mutectFilePath <- paste0(dataDirectory, directory, '_', 'mutect.vcf')
  vardictFilePath <- paste0(dataDirectory, directory, '_', 'vardict.vcf')
  varscanFilePath <- paste0(dataDirectory, directory, '_', 'varscan.vcf')
  truthFilePath <- paste0(dataDirectory, directory, '_', 'truth.bed')
  
  freebayesVCF <- readVcf(freebayesFilePath, 'hg19')
  freebayes.ranges <- rowRanges(freebayesVCF)
  freebayes.info <- info(freebayesVCF)
  freebayes.geno <- geno(freebayesVCF)
  
  mutectVCF <- readVcf(mutectFilePath, 'hg19')
  mutect.ranges <- rowRanges(mutectVCF)
  mutect.info <- info(mutectVCF)
  mutect.geno <- geno(mutectVCF)
  
  vardictVCF <- readVcf(vardictFilePath, 'hg19')
  vardict.ranges <- rowRanges(vardictVCF)
  vardict.info <- info(vardictVCF)
  vardict.geno <- geno(vardictVCF)
  
  varscanVCF <- readVcf(varscanFilePath, 'hg19')
  varscan.ranges <- rowRanges(varscanVCF)
  varscan.info <- info(varscanVCF)
  varscan.geno <- geno(varscanVCF)
  
  if (directory == 'test') {
    truthSet <- readBedFile(truthFilePath) 
  } else {
    truthSet <- NULL
  }
  
  return(list(freebayes.ranges = freebayes.ranges, freebayes.info = freebayes.info, freebayes.geno = freebayes.geno,
              mutect.ranges = mutect.ranges, mutect.info = mutect.info, mutect.geno = mutect.geno,
              vardict.ranges = vardict.ranges, vardict.info = vardict.info, vardict.geno = vardict.geno,
              varscan.ranges = varscan.ranges, varscan.info = varscan.info, varscan.geno = varscan.geno,
              truthSet = truthSet))
}

# Merges the SNP callers results into a single Granges objects according to position
# For training data uses the truth set and labels the snp caller results accordingly
getAllPredictions <- function(sampleData, sampleName) {
  allPredictions <- unique(c(sampleData$freebayes.ranges, sampleData$mutect.ranges, sampleData$vardict.ranges, sampleData$varscan.ranges))
  allPredictions$freebayesSNP <- overlapsAny(allPredictions, sampleData$freebayes.ranges, type = 'equal')
  allPredictions$mutectSNP <- overlapsAny(allPredictions, sampleData$mutect.ranges, type = 'equal')
  allPredictions$vardictSNP <- overlapsAny(allPredictions, sampleData$vardict.ranges, type = 'equal')
  allPredictions$varscanSNP <- overlapsAny(allPredictions, sampleData$varscan.ranges, type = 'equal')
  
  if (is.null(sampleData$truthSet) == FALSE) {
    allPredictions$truthTable <- overlapsAny(allPredictions, sampleData$truthSet, type = 'equal')
  }
  
  dir.create(file.path('./predictions'), showWarnings = FALSE)
  saveRDS(allPredictions, paste0('./predictions/', sampleData, '_predictions.rds'))
  return(allPredictions)
}

# USAGE:
# syn1.data <- loadTrainingData(parentDirectory = '../data/training/', directory = 'syn1')
# syn1.predictions <- getAllPredictions(syn1.data)

################################################################
############## DBSNP ACCESS AND QUERY FUNCTION #################
################################################################

# Creates the mysql script that queries the UCSC dbSNP database version 146 for the called SNP position
createRsIdScript <- function(gr, snpFolder, sampleName) {
  sourceCpp('./UCSCBinFromRange.cpp')
  chromosomes <- seqnames(gr)
  start <- start(gr)
  end <- start + 1
  bins <- sapply(1:length(start), function(i) {binFromRangeStandard(start[i], end[i])})
  lines <- sapply(1:length(gr), function(i) {paste0('select chrom,chromStart,name from snp146 where bin = ', bins[i],' and chrom = \"chr', chromosomes[i], 
                                                    '\" and class = \"single\" and chromStart + 1 = ', start[i], ';')})
  write(lines, paste0(snpFolder, sampleName, '_snps.sql'))
}

# Checks whether the called SNPs are in the UCSC dbSNP list obtained and return a logical vector
isDbSNP <- function(gr, snpFolder, sampleName) {
  loadDbSNPs <- function(filename) {
    df <- read.table(filename, header=F, stringsAsFactors=F)
    names(df) <- c('chr','start','rsID')
    chr <- sapply(df$chr, substr, 4, 5)
    start <- df$start + 1
    rsID <- df$rsID
    gr <- GRanges(chr, IRanges(start, start))
    gr$rsID <- df$rsID
    return(gr)
  }
  
  snps <- loadDbSNPs(paste0(snpFolder, sampleName, '_snps.txt'))
  return(overlapsAny(gr, snps))
}

# Creates the sql script to query dbSNP
# Creates the shell script for the sql call
# Executes shell script via system call
# Saves the resulting vector in the dbSNP folder
dbSnpAnalysis <- function(sample.predictions, snpFolder = './dbSNP/', sampleName = 'syn1') {
  dir.create(file.path(snpFolder), showWarnings = FALSE)
  createRsIdScript(sample.predictions, snpFolder, sampleName)
  
  write('mysql -h genome-mysql.cse.ucsc.edu -u genome -A -D hg19 --skip-column-names < ', sampleName, '_snps.sql > ', sampleName, '_snps.txt', 
        file = paste0(snpFolder, sampleName, '_sql.sh'))

  system(paste0('cd ', snpFolder))
  system('sh ', sampleName, '_sql.sh')
  
  sample.dbSNP <- isDbSNP(sample.predictions, snpFolder, sampleName)
  saveRDS(sample.dbSNP, file = paste0(snpFolder, sampleName, '_dbSNP.rds'))
}

# USAGE:
# dbSnpAnalysis(sample.predictions = syn1.predictions, snpFolder = './dbSNP/', sampleName = 'syn1')

################################################################
############# MAPPING QUALITY ANALYSIS FROM BAM ################
################################################################

# Creates the shell script that calls samtools for the called SNP positions to determine the mean and rms mapping quality
createQualityScript <- function(qualityFolder = './quality/', filename, gr) {
  l <- length(gr)
  lines <- paste0(rep(paste0('samtools view ', filename, '.bam '), l), 
                  paste0(seqnames(gr), rep(':', l), start(gr), rep('-', l), end(gr)),  
                  rep(paste0(' | awk \'{sum+=$5} END { if (sum == 0)  print 0 ; else print sum/NR }\' >> ', filename, '_meanQuality.txt'), l)
  )
  write(lines, file = paste0('./scripts/', filename, '_meanQuality.sh'))
  
  lines <- paste0(rep(paste0('samtools view ', filename, '.bam '), l), 
                  paste0(seqnames(gr), rep(':', l), start(gr), rep('-', l), end(gr)),  
                  rep(paste0(' | awk \'{sum+=$5*$5} END { if (sum == 0)  print 0 ; else print sqrt(sum/NR) }\' >> ', filename, '_rmsQuality.txt'), l)
  )
  write(lines, file = paste0(qualityFolder, filename, '_rmsQuality.sh'))
}

# Reads the calculated MAPQ values for normal and tumor bam files into a data frame
getQualityValues <- function(qualityFolder = './quality/', sampleName) {
  normal.meanQuality <- read.table(file = paste0(qualityFolder, sampleName, '_normal_meanQuality.txt'))
  normal.rmsQuality <- read.table(file = paste0(qualityFolder, sampleName, '_normal_rmsQuality.txt'))
  tumor.meanQuality <- read.table(file = paste0(qualityFolder, sampleName, '_tumor_meanQuality.txt'))
  tumor.rmsQuality <- read.table(file = paste0(qualityFolder, sampleName, '_tumor_rmsQuality.txt'))
  dat <- cbind(normal.meanQuality, normal.rmsQuality, tumor.meanQuality, tumor.rmsQuality)
  names(dat) <- c('normal.meanQuality', 'normal.rmsQuality', 'tumor.meanQuality', 'tumor.rmsQuality')
  return(dat)
}

# Creates the shell scripts using samtools to calculate mapping quality at reported SNP locations
# Execute the shell scripts by system call
# Loads the quality values for the specified sample into a data frame and stores it in rds format
qualityAnalysis <- function(qualityFolder = './quality/', sampleName = 'syn1', sample.predictions) {
  dir.create(file.path(qualityFolder), showWarnings = FALSE)
  createQualityScript(qualityFolder, filename = paste0(sampleName, '_normal'), sample.predictions)
  createQualityScript(qualityFolder, filename = paste0(sampleName, '_tumor'), sample.predictions)
  
  system(paste0('cd ', qualityFolder))
  system(paste0('sh ', sampleName, '_normal_meanQuality.sh'))
  system(paste0('sh ', sampleName, '_normal_rmsQuality.sh'))
  system(paste0('sh ', sampleName, '_tumor_meanQuality.sh'))
  system(paste0('sh ', sampleName, '_tumor_rmsQuality.sh'))
  
  quality <- getQualityValues(qualityFolder, sampleName)
  saveRDS(quality, paste0(qualityFolder, sampleName, '_quality.rds'))
}

# USAGE:
# qualityAnalysis(qualityFolder = './quality/', sampleName = 'syn1', syn1.predictions)

################################################################
############# BAM FILE PROCESSING FOR COUNTING #################
################################################################

# Uses samtools to split the normal/tumor bam files into two bam files containing forward and reverse reads 
# Indices the original and split bam files for further calculations
createForwardReverseReadBams <- function(parentDirectory = '../data/training/', sampleName) {
  system(paste0('cd ', parentDirectory, sampleName))
  system(paste0('samtools index ', sampleName, '_normal.bam'))
  system(paste0('samtools view -b -F 16 ', sampleName, '_normal.bam > ', sampleName, '_normal_forward.bam'))
  system(paste0('samtools index ', sampleName, '_normal_forward.bam'))
  system(paste0('samtools view -b -f 16 ', sampleName, '_normal.bam > ', sampleName, '_normal_reverse.bam'))
  system(paste0('samtools index ', sampleName, '_normal_reverse.bam'))
}

# Calculates the frequency of each nucleotide at the positions specified in the input pileup result
# Stores the result in a data frame
pileupFrequency <- function(pileupResult) {
  nucleotides <- levels(pileupResult$nucleotide)
  res <- split(pileupResult, pileupResult$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    return(nuctab)
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupResult$nucleotide))
  res[2:ncol(res)] <- apply(res[2:ncol(res)], 2, as.numeric)
  return(res)
}

# Calculates the pileups for the called SNP locations on normal/tumor bam files of all/forward/reverse reads
# Gets the frequency at each location by using te pileups and stores the result in rds format in the frequencies folder
getPileupCounts <- function(gr, sampleName, parentDirectory = '../data/training/') {
  normalBam <- BamFile(paste0(parentDirectory, sampleName, '_normal.bam'))
  normalForwardBam <- BamFile(paste0(parentDirectory, sampleName, '_normal_forward.bam'))
  normalReverseBam <- BamFile(paste0(parentDirectory, sampleName, '_normal_reverse.bam'))
  
  tumorBam <- BamFile(paste0(parentDirectory, sampleName, '_tumor.bam'))
  tumorForwardBam <- BamFile(paste0(parentDirectory, sampleName, '_tumor_forward.bam'))
  tumorReverseBam <- BamFile(paste0(parentDirectory, sampleName, '_tumor_reverse.bam'))
  
  param <- ScanBamParam(which = gr)
  p_param <- PileupParam(max_depth=1000, distinguish_strand=FALSE)
  
  normal.pileup <- pileup(normalBam, scanBamParam=param, pileupParam=p_param)
  normalForward.pileup <- pileup(normalForwardBam, scanBamParam=param, pileupParam=p_param)
  normalReverse.pileup <- pileup(normalReverseBam, scanBamParam=param, pileupParam=p_param)
  
  tumor.pileup <- pileup(tumorBam, scanBamParam=param, pileupParam=p_param)
  tumorForward.pileup <- pileup(tumorForwardBam, scanBamParam=param, pileupParam=p_param)
  tumorReverse.pileup <- pileup(tumorReverseBam, scanBamParam=param, pileupParam=p_param)
  
  normal.freq <- pileupFrequency(normal.pileup)
  normalForward.freq <- pileupFrequency(normalForward.pileup)
  normalReverse.freq <- pileupFrequency(normalReverse.pileup)
  
  tumor.freq <- pileupFrequency(tumor.pileup)
  tumorForward.freq <- pileupFrequency(tumorForward.pileup)
  tumorReverse.freq <- pileupFrequency(tumorReverse.pileup)
  
  freq <- list(normal.freq = normal.freq, normalForward.freq = normalForward.freq, normalReverse.freq = normalReverse.freq, 
               tumor.freq = tumor.freq, tumorForward.freq = tumorForward.freq, tumorReverse.freq = tumorReverse.freq)
  
  dir.create(file.path('./frequency'), showWarnings = FALSE)
  saveRDS(freq, paste0('./frequency/', sampleName, '_frequency.rds'))
  return(freq)
}

# Uses frequency data to calculate the forward/reverse read counts supporting reference/variant in normal/tumor samples
# Calculates the read depth at called SNP positions in normal and tumor bam files
getReadCounts <- function(predictions, frequencies) {
  reference.nuc <- as.vector(predictions$REF)
  variant.nuc <- as.vector(unlist(predictions$ALT))
  
  tmp <- frequencies$normal.freq
  normal.gr <- GRanges(tmp$seqnames, IRanges(as.numeric(tmp$start), as.numeric(tmp$start)))
  tmp <- frequencies$normalForward.freq
  normalForward.gr <- GRanges(tmp$seqnames, IRanges(as.numeric(tmp$start), as.numeric(tmp$start)))
  tmp <- frequencies$normalReverse.freq
  normalReverse.gr <- GRanges(tmp$seqnames, IRanges(as.numeric(tmp$start), as.numeric(tmp$start)))
  
  tmp <- frequencies$tumor.freq
  tumor.gr <- GRanges(tmp$seqnames, IRanges(as.numeric(tmp$start), as.numeric(tmp$start)))
  tmp <- frequencies$tumorForward.freq
  tumorForward.gr <- GRanges(tmp$seqnames, IRanges(as.numeric(tmp$start), as.numeric(tmp$start)))
  tmp <- frequencies$tumorReverse.freq
  tumorReverse.gr <- GRanges(tmp$seqnames, IRanges(as.numeric(tmp$start), as.numeric(tmp$start)))
  
  normal.forward.reference <- sapply(1:length(normalForward.gr), function(i) {frequencies$normalForward.freq[i, reference.nuc[i]]})
  normal.reverse.reference <- sapply(1:length(normalReverse.gr), function(i) {frequencies$normalReverse.freq[i, reference.nuc[i]]})
  normal.forward.variant <- sapply(1:length(normalForward.gr), function(i) {frequencies$normalForward.freq[i, variant.nuc[i]]})
  normal.reverse.variant <- sapply(1:length(normalReverse.gr), function(i) {frequencies$normalReverse.freq[i, variant.nuc[i]]})
  
  tumor.forward.reference <- sapply(1:length(tumorForward.gr), function(i) {frequencies$tumorForward.freq[i, reference.nuc[i]]})
  tumor.reverse.reference <- sapply(1:length(tumorReverse.gr), function(i) {frequencies$tumorReverse.freq[i, reference.nuc[i]]})
  tumor.forward.variant <- sapply(1:length(tumorForward.gr), function(i) {frequencies$tumorForward.freq[i, variant.nuc[i]]})
  tumor.reverse.variant <- sapply(1:length(tumorReverse.gr), function(i) {frequencies$tumorReverse.freq[i, variant.nuc[i]]})
  
  normal.readDepth <- rowSums(frequencies$normal.freq[,c('A', 'C', 'G', 'T')])
  tumor.readDepth <- rowSums(frequencies$tumor.freq[,c('A', 'C', 'G', 'T')])
  
  getMissing <- function(predictions, gr, values) {
    tmp <- rep(0, length(predictions))
    tmp[queryHits(findOverlaps(predictions, gr))] <- values
    return(tmp)
  }
  
  normal.forward.reference <- getMissing(predictions, normalForward.gr, normal.forward.reference)
  normal.reverse.reference <- getMissing(predictions, normalReverse.gr, normal.reverse.reference)
  normal.forward.variant <- getMissing(predictions, normalForward.gr, normal.forward.variant)
  normal.reverse.variant <- getMissing(predictions, normalReverse.gr, normal.reverse.variant)
  
  tumor.forward.reference <- getMissing(predictions, tumorForward.gr, tumor.forward.reference)
  tumor.reverse.reference <- getMissing(predictions, tumorReverse.gr, tumor.reverse.reference)
  tumor.forward.variant <- getMissing(predictions, tumorForward.gr, tumor.forward.variant)
  tumor.reverse.variant <- getMissing(predictions, tumorReverse.gr, tumor.reverse.variant)
  
  normal.readDepth <- getMissing(predictions, normal.gr, normal.readDepth)
  tumor.readDepth <- getMissing(predictions, tumor.gr, tumor.readDepth)
  
  dat <- list(normal.forward.reference = normal.forward.reference, normal.reverse.reference = normal.reverse.reference,
              normal.forward.variant = normal.forward.variant, normal.reverse.variant = normal.reverse.variant,
              tumor.forward.reference = tumor.forward.reference, tumor.reverse.reference = tumor.reverse.reference,
              tumor.forward.variant = tumor.forward.variant, tumor.reverse.variant = tumor.reverse.variant,
              normal.readDepth = normal.readDepth, tumor.readDepth = tumor.readDepth)
  
  return(dat)
}

# USAGE:
# createForwardReverseReadBams(parentDirectory = '../data/training/', 'syn1')
# syn1.frequency <- getPileupCounts(syn1.predictions, 'syn1', parentDirectory = '../data/training/')
# syn1.readCounts <- getReadCounts(syn1.predictions, syn1.frequency)


################################################################
################## PREPROCESSING (COMPLETE) ####################
################################################################

# All preprocessing steps required for feature extraction for a single sample
preprocessing <- function(sampleName) {
  # Load vcf data from = different tool results
  sample.data <- loadTrainingData(parentDirectory = '../data/training/', directory = sampleName)
  # Merge different tool results into one Granges object with truth values
  sample.predictions <- getAllPredictions(sample.data)
  # For each SNP reported check whether it is in the UCSC dbSNP database
  dbSnpAnalysis(sample.predictions = sample.predictions, snpFolder = './dbSNP/', sampleName = sampleName)
  # Get mean and rms quality of the normal and tumor bam files for the reported SNP positions
  qualityAnalysis(qualityFolder = './quality/', sampleName = sampleName, sample.predictions)
  # Split the normal/tumor bam file into forward and reverse read bam files and index them
  createForwardReverseReadBams(parentDirectory = '../data/training/', sampleName = sampleName)
  # Calculate the pileups for reported SNP locations in original, reverse and forward read bam files
  sample.frequency <- getPileupCounts(sample.predictions, sampleName, parentDirectory = '../data/training/')
  # Find forward/reverse read counts supporting reference/variant in tumor/normal bam files & calculate read depth at SNP locations
  sample.readCounts <- getReadCounts(sample.predictions, sample.frequency)
}

################################################################
############### FEATURE EXTRACTION (COMPLETE) ##################
################################################################

# Extracting features by using the preprocessed data and storing it in the processedData folder as CSV files
extractFeatures <- function(sampleName) {
  sample.predictions <- readRDS(paste0('./predictions/', sampleName, '_predictions.rds'))
  sample.dbSNP <- readRDS(paste0('./dbSNP/', sampleName, '_dbSNP.rds'))
  sample.quality <- readRDS(paste0('./quality/', sampleName, '_quality.rds'))
  sample.frequencies <- readRDS(paste0('./frequency/', sampleName ,'_frequency.rds'))
  sample.readCounts <- getReadCounts(predictions = sample.predictions, frequencies = sample.frequencies)
  
  sample.features <- data.frame(freebayesSNP  = as.vector(sample.predictions$freebayesSNP), mutectSNP  = as.vector(sample.predictions$mutectSNP), 
                                varscanSNP  = as.vector(sample.predictions$varscanSNP), vardictSNP  = as.vector(sample.predictions$vardictSNP))
  
  sample.features$dbSNP <- sample.dbSNP
  sample.features$normalMAPQ <- sample.quality$normal.meanQuality
  sample.features$tumorMAPQ <- sample.quality$tumor.meanQuality
  
  sample.features$normal.forward.reference <- sample.readCounts$normal.forward.reference
  sample.features$normal.reverse.reference <- sample.readCounts$normal.reverse.reference
  sample.features$normal.forward.variant <- sample.readCounts$normal.forward.variant
  sample.features$normal.reverse.variant <- sample.readCounts$normal.reverse.variant
  
  sample.features$tumor.forward.reference <- sample.readCounts$tumor.forward.reference
  sample.features$tumor.reverse.reference <- sample.readCounts$tumor.reverse.reference
  sample.features$tumor.forward.variant <- sample.readCounts$tumor.forward.variant
  sample.features$tumor.reverse.variant <- sample.readCounts$tumor.reverse.variant
  
  sample.features$normal.readDepth <- sample.readCounts$normal.readDepth
  sample.features$tumor.readDepth <- sample.readCounts$tumor.readDepth
  
  sample.features$truthTable <- as.vector(sample.predictions$truthTable)

  dir.create(file.path('./processedData'), showWarnings = FALSE)
  write.csv(sample.features, file = paste0('./processedData/', sampleName ,'Features.csv'), row.names=FALSE, na="")
}

# USAGE: 
# extractFeatures('syn1')

################################################################
######################## MAIN SCRIPT ###########################
################################################################

#################################
######### PREPROCESSING #########
#################################

# Preprocesses the data and saves the intermediate results in the corresponding folders

# Training data
preprocessing(sampleName = 'syn1')
preprocessing(sampleName = 'syn2')
preprocessing(sampleName = 'syn3')
preprocessing(sampleName = 'syn4')
preprocessing(sampleName = 'syn5')
preprocessing(sampleName = 'real1')
preprocessing(sampleName = 'real2')

# Test data
preprocessing(sampleName = 'test')

#################################
###### FEATURE EXTRACTION #######
#################################

# Extracts the features and saves the final feature matrix as csv file into the processedData folder

# Training data
extractFeatures(sampleName = 'syn1')
extractFeatures(sampleName = 'syn2')
extractFeatures(sampleName = 'syn3')
extractFeatures(sampleName = 'syn4')
extractFeatures(sampleName = 'syn5')
extractFeatures(sampleName = 'real1')
extractFeatures(sampleName = 'real2')

# Test data
extractFeatures(sampleName = 'test')



