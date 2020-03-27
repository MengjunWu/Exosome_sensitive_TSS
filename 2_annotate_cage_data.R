#!/usr/local/bin/Rscript
###################################################
########### Sandelin Lab Cage Pipeline ############
######## Annotation of CAGE peaks (2 of 2) ########
######### Script for annotating Cage data #########
############ Kristoffer Vitting-Seerup ############
############### Evdoxia Karadoulama ###############
###################################################

# the idea is to load the Rdata file(s) prepared by the script "1_massage_offical_annotation.R" as well 
# as set of CAGE TCs. The cage TCs should then be compared to the annotation and get annotated
# in a ranked manner.
# The output should be identical to what the pipeline currently produce
# The input should be
#   - The working directory - allows me to load the file
#   - The organism
#   - The assembly

### Use argparse to make help file and parse input
if(TRUE) {
  suppressMessages( library('argparser') )
  
  # create parser
  argParse <- arg_parser(name = '2_annotate_cage_data.R', description = 'Description:\nPart 2 (of 2) of the Sandelin Lab Cage Pipeline Annotation scripts.\nThis script that annotates the cage clusters by comparing to the annotation of the GTF file prepare in part 1.')
  # add arguments
  argParse <- add_argument(parser = argParse, short = '-b', arg = '-bedFile',       help = 'The bed8 file produced by the pipeline file with the clusters to annotate. If NULL (default) the standard file outputted by the pipeline is used.')
  argParse <- add_argument(parser = argParse, short = '-g', arg = '-annotationGTF', help = 'The FULL path to the GTF (not gff) file that should be used for annoation (For example /.../refseq.gtf or /.../gencode.17.gtf).')
  argParse <- add_argument(parser = argParse, short = '-o', arg = '-outputFile',                  help = 'The FULL path and name to the file where the annoation should be outputted.')
  
  argParse <- add_argument(parser = argParse, arg = '-tssUpstreamWindow',          default = 100,  help = 'The maximum distance (in nucleotides) the peak of a CAGE cluster can be upstream of a annotated tss, and still be associated with that tss.')
  argParse <- add_argument(parser = argParse, arg = '-tssDownstreamWindow',        default = 100,  help = 'The maximum distance (in nucleotides) the peak of a CAGE cluster can be downstream of a annotated tss, and still be associated with that tss.')
  argParse <- add_argument(parser = argParse, arg = '-geneUpstreamPromoterWindow', default = 1000, help = 'The maximum distance (in nucleotides) the peak of a CAGE cluster can be upstream of a annotated gene, and still be associated as upstream of that gene.')
  
  # evaluate input
  inputedArguments <- parse_args(argParse, argv = commandArgs(trailingOnly = TRUE)) # creates a named list with arguments
  
  ### Test input
  inputTest <- sapply(inputedArguments[-c(2:3)], is.na) # first is always help and opts 
  if( any( inputTest ) ) {
    stop(paste('The following input arguments are missing: -', paste( names(inputTest)[which(inputTest)], collapse = ', -'), sep='') )
  }
  
  ### Note: Currently argparser ( v 0.1 ) gives a warning message for arguments starting with certain types of names
}

### Load dependencies
message('Loading dependencies...')
suppressMessages( library('rtracklayer')   )
suppressMessages( library('GenomicRanges') )
suppressMessages( library('data.table')    )
suppressMessages( library('plyr')          )

### Massage and test input
if(TRUE) {
  annotationPath <- paste( gsub('.gtf$', '', inputedArguments$annotationGTF, perl = T, ignore.case = T), '.parsed.Rdata', sep='')
  
  ### Test whether annotation exists
  if( ! file.exists(annotationPath) ) {
    
    ### If Noat and GTF file exist create the annoation
    if( ! file.exists(inputedArguments$annotationGTF) ) {
      stop('The supplied gtf file does not seem to exist')
    } else {
      message('The parsed version of the GTF file does not seem to exist - creating it now...')
      system(
        command = paste(
          '1_massage_official_annotation.R -gtf',
          inputedArguments$annotationGTF,
          sep=' '
        )
      )
      message('The passed version of the GTF file have now been created (and saved) - continuing with annotation...')
    }
  }
}

### Initialize helper functions
if(TRUE) {
  readBedFile <- function(bedPathAndFileName) {
    ### Read in the file
    aBedFile <- as.data.frame( fread(input=bedPathAndFileName, sep='\t', showProgress=F), stringsAsFactors = F)
    
    if( ncol(aBedFile) != 8) {
      stop('The bedfile supplied is not a bed8 file. Please revice.')
    } 
    
    ### Add names (and low for multiple different bed file formats) (could be more sophisticated using which types are in the different columns)
    setnames(aBedFile, colnames(aBedFile), new=c('chr','start','end','name', 'score','strand','peakStart','peakEnd'))
    
    ### Return result
    return(aBedFile)
  }
  invertGRangesStrand <- function(aGRange) {
    plusIndex  <- which(strand(aGRange) == '+')
    minusIndex <- which(strand(aGRange) == '-')
    
    # change
    strand(aGRange)[plusIndex]  <- '-'
    strand(aGRange)[minusIndex] <- '+'
    
    return(aGRange)
  }
  handleMultipleAssociations <- function(aDF) {
    isDuplicated <- duplicated( aDF$cageTSS) | duplicated( aDF$cageTSS, fromLast = T) 
    
    # devide into duplciated and non duplicated
    notDuplicatedData <- aDF[ which( ! isDuplicated ),]
    duplicatedData    <- aDF[ which(   isDuplicated ),]
    
    # loop over duplciated and deduplicate
    if(nrow(duplicatedData) != 0) {
      suppressMessages( 
        deDuplicatedData <- ddply(duplicatedData, .variables = 'cageTSS', function(localDF) {
          if(nrow(aDF) != 1) {
            localDF <- data.frame(
              cageTSS=localDF$cageTSS[1], 
              transcript_id=paste(localDF$transcript_id, collapse = ','), 
              gene_name    =paste(unique(localDF$gene_name), collapse = ','), 
              gene_id      =paste(unique(localDF$gene_id  ), collapse = ','),
              stringsAsFactors = F)
          }
          return(localDF)
        })
      )
      
      deDuplicatedData <- deDuplicatedData[, which(colnames(deDuplicatedData) %in% colnames(notDuplicatedData))]
      
      # combined deduplicated and nonduplcated data
      notDuplicatedData <- rbind(notDuplicatedData, deDuplicatedData)
      notDuplicatedData <- notDuplicatedData[sort.list(notDuplicatedData$cageTSS),]
      
    }
    
    return(notDuplicatedData)
  }
  
}

### Load and massage data
message('Step 1 of 3: Loading input data...')
if(TRUE) {
  ### Load required annotation
  load(file=annoationPath) 
  
  ### Load CAGE data
  cageClustersDf <- readBedFile(bedPathAndFileName = inputedArguments$bedFile)
  cageClustersDf[,c('tssAnnoation','transcript_id','gene_id','gene_name')] <- NA
  
  ### Convert to GRanges
  cageClustersGr <- GRanges(seqnames = cageClustersDf$chr, ranges = IRanges(cageClustersDf$peakStart, cageClustersDf$peakEnd), strand=cageClustersDf$strand, id=cageClustersDf$name)
  
  if( ! (strand(cageClustersGr) %in% c('+','-'))@values ) {
    warning('!!! There are cage clusters without strand annoation in the file !!!')
  }
}


### Prefilter all those that are not close to genes
message('Step 2 of 3: Classifying cage peaks according to annoation...')
if(TRUE) {
  ### Find overlap with genes - using the cage peaks
  extendedGeneEdges <- geneEdges
  start(extendedGeneEdges) <- start(extendedGeneEdges) - (as.integer(inputedArguments$geneUpstreamPromoterWindow) +1)
  end(extendedGeneEdges)   <- end(extendedGeneEdges)   + (as.integer(inputedArguments$geneUpstreamPromoterWindow) +1)
  start(extendedGeneEdges)[which(start(extendedGeneEdges) < 1)] <- 1
  
  annotationOverlapping <- suppressMessages( suppressWarnings( overlapsAny(cageClustersGr, extendedGeneEdges, ignore.strand = TRUE) ) )
  
  # remove those without overlap from GRanges
  toKeep <- which( annotationOverlapping )
  cageClustersGr   <- cageClustersGr  [toKeep,]
}

### Compare cage peaks to annoation
if(TRUE) {
  ####### Sense direction
  ### Primary TSS
  if(TRUE) {
    ### Extend cage peaks with the wanted window (since cage is extended upstream=downstream)
    cageTSS <- GenomicRanges::promoters(cageClustersGr, upstream = as.integer(inputedArguments$tssDownstreamWindow), downstream = as.integer(inputedArguments$tssUpstreamWindow))
    
    ### Find overlap - using the extended TSS
    primaryTSSoverlap <- suppressWarnings( as.data.frame( findOverlaps(cageTSS, primaryTSS) ))
    setnames(primaryTSSoverlap, c('cageTSS', 'transcript_id'))
    
    # Replace indexes with names
    primaryTSSoverlap$cageTSS       <- cageTSS$id              [ primaryTSSoverlap$cageTSS       ]
    primaryTSSoverlap$transcript_id <- primaryTSS$transcript_id[ primaryTSSoverlap$transcript_id ]
    primaryTSSoverlap$gene_name     <- primaryTSS$gene_name    [ match( primaryTSSoverlap$transcript_id , primaryTSS$transcript_id) ]
    primaryTSSoverlap$gene_id       <- primaryTSS$gene_id      [ match( primaryTSSoverlap$transcript_id , primaryTSS$transcript_id) ]
    
    # Handle that some might be associated with multiple transcripts
    primaryTSSoverlap <- handleMultipleAssociations(primaryTSSoverlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(primaryTSSoverlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- 'primary_tss'
    cageClustersDf$transcript_id [correspindingMatch] <- primaryTSSoverlap$transcript_id
    cageClustersDf$gene_name     [correspindingMatch] <- primaryTSSoverlap$gene_name
    cageClustersDf$gene_id       [correspindingMatch] <- primaryTSSoverlap$gene_id
    
    # remove identified from GRanges
    toKeep <- which( ! cageTSS$id %in% primaryTSSoverlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
    cageTSS          <- cageTSS         [toKeep,]
  }
  
  ### Alternative TSS
  if(TRUE) {
    ### Find overlap - using the extended TSS
    alternativeTSSoverlap <- suppressWarnings( as.data.frame( findOverlaps(cageTSS, alternativeTSS) ))
    setnames(alternativeTSSoverlap, c('cageTSS', 'transcript_id'))
    
    # Replace indexes with names
    alternativeTSSoverlap$cageTSS       <- cageTSS$id                  [alternativeTSSoverlap$cageTSS      ]
    alternativeTSSoverlap$transcript_id <- alternativeTSS$transcript_id[alternativeTSSoverlap$transcript_id]
    alternativeTSSoverlap$gene_name     <- alternativeTSS$gene_name    [ match( alternativeTSSoverlap$transcript_id , alternativeTSS$transcript_id) ]
    alternativeTSSoverlap$gene_id       <- alternativeTSS$gene_id      [ match( alternativeTSSoverlap$transcript_id , alternativeTSS$transcript_id) ]
    
    # Handle that some might be associated with multiple transcripts
    alternativeTSSoverlap <- handleMultipleAssociations(aDF = alternativeTSSoverlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(alternativeTSSoverlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- 'alternative_tss'
    cageClustersDf$transcript_id [correspindingMatch] <- alternativeTSSoverlap$transcript_id
    cageClustersDf$gene_name     [correspindingMatch] <- alternativeTSSoverlap$gene_name
    cageClustersDf$gene_id       [correspindingMatch] <- alternativeTSSoverlap$gene_id
    
    # remove identified from GRanges
    toKeep <- which( ! cageTSS$id %in% alternativeTSSoverlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
    cageTSS          <- cageTSS         [toKeep,]
  }
  
  ### 5UTR
  if(TRUE) {
    ### Find overlap - using just the cage peak
    utr5overlap <- suppressWarnings( as.data.frame( findOverlaps(cageClustersGr, utr5regGene) ))
    setnames(utr5overlap, c('cageTSS', 'gene_id'))
    
    # Replace indexes with names
    utr5overlap$cageTSS   <- cageClustersGr$id    [utr5overlap$cageTSS  ]
    utr5overlap$gene_id   <- utr5regGene$gene_id  [utr5overlap$gene_id  ]
    utr5overlap$gene_name <- utr5regGene$gene_name[match(utr5overlap$gene_id, utr5regGene$gene_id)]
    
    # Handle that some might be associated with multiple transcripts
    utr5overlap <- handleMultipleAssociations(utr5overlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(utr5overlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- '5utr'
    cageClustersDf$gene_name     [correspindingMatch] <- utr5overlap$gene_name
    cageClustersDf$gene_id       [correspindingMatch] <- utr5overlap$gene_id
    
    # remove identified from GRanges
    toKeep <- which( ! cageClustersGr$id %in% utr5overlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
  }
  
  ### Coding regions
  if(TRUE) {
    ### Find overlap - using just the cage peak
    codingOverlap <- suppressWarnings( as.data.frame( findOverlaps(cageClustersGr, codingRegions) ))
    setnames(codingOverlap, c('cageTSS', 'gene_id'))
    
    # Replace indexes with names
    codingOverlap$cageTSS   <- cageClustersGr$id      [codingOverlap$cageTSS  ]
    codingOverlap$gene_id   <- codingRegions$gene_id  [codingOverlap$gene_id  ]
    codingOverlap$gene_name <- codingRegions$gene_name[ match(codingOverlap$gene_id, codingRegions$gene_id) ]
    
    # Handle that some might be associated with multiple transcripts
    codingOverlap <- handleMultipleAssociations(codingOverlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(codingOverlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- 'coding'
    cageClustersDf$gene_id       [correspindingMatch] <- codingOverlap$gene_id
    cageClustersDf$gene_name     [correspindingMatch] <- codingOverlap$gene_name
    
    # remove identified from GRanges
    toKeep <- which( ! cageClustersGr$id %in% codingOverlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
  }
  
  ### 3UTR
  if(TRUE) {
    ### Find overlap - using just the cage peak
    utr3overlap <- suppressWarnings( as.data.frame( findOverlaps(cageClustersGr, utr3regGene) ))
    setnames(utr3overlap, c('cageTSS', 'gene_id'))
    
    # Replace indexes with names
    utr3overlap$cageTSS   <- cageClustersGr$id    [utr3overlap$cageTSS  ]
    utr3overlap$gene_id   <- utr3regGene$gene_id  [utr3overlap$gene_id  ]
    utr3overlap$gene_name <- utr3regGene$gene_name[match(utr3overlap$gene_id, utr3regGene$gene_id)]
    
    # Handle that some might be associated with multiple transcripts
    utr3overlap <- handleMultipleAssociations(utr3overlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(utr3overlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- '3utr'
    cageClustersDf$gene_id       [correspindingMatch] <- utr3overlap$gene_id
    cageClustersDf$gene_name     [correspindingMatch] <- utr3overlap$gene_name
    
    # remove identified from GRanges
    toKeep <- which( ! cageClustersGr$id %in% utr3overlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
  }
  
  ### Exonic
  if(TRUE) {
    ### Find overlap - using just the cage peak
    exonOverlap <- suppressWarnings( as.data.frame( findOverlaps(cageClustersGr, geneExons) ))
    setnames(exonOverlap, c('cageTSS', 'gene_id'))
    
    # Replace indexes with names
    exonOverlap$cageTSS   <- cageClustersGr$id  [exonOverlap$cageTSS]
    exonOverlap$gene_id   <- geneExons$gene_id  [exonOverlap$gene_id]
    exonOverlap$gene_name <- geneExons$gene_name[match(exonOverlap$gene_id, geneExons$gene_id)]
    
    # Handle that some might be associated with multiple transcripts
    exonOverlap <- handleMultipleAssociations(exonOverlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(exonOverlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- 'exonic'
    cageClustersDf$gene_id       [correspindingMatch] <- exonOverlap$gene_id
    cageClustersDf$gene_name     [correspindingMatch] <- exonOverlap$gene_name
    
    # remove identified from GRanges
    toKeep <- which( ! cageClustersGr$id %in% exonOverlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
  }
  
  ### Intron
  if(TRUE) {
    intronOverlap <- suppressWarnings( as.data.frame( findOverlaps(cageClustersGr, geneIntrons) ))
    setnames(intronOverlap, c('cageTSS', 'gene_id'))
    
    # Replace indexes with names
    intronOverlap$cageTSS   <- cageClustersGr$id    [intronOverlap$cageTSS  ]
    intronOverlap$gene_id   <- geneIntrons$gene_id  [intronOverlap$gene_id]
    intronOverlap$gene_name <- geneIntrons$gene_name[match(intronOverlap$gene_id, geneIntrons$gene_id)]
    
    # Handle that some might be associated with multiple transcripts
    intronOverlap <- handleMultipleAssociations(intronOverlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(intronOverlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- 'intronic'
    cageClustersDf$gene_id       [correspindingMatch] <- intronOverlap$gene_id
    cageClustersDf$gene_name     [correspindingMatch] <- intronOverlap$gene_name
    
    # remove identified from GRanges
    toKeep <- which( ! cageClustersGr$id %in% intronOverlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
  }
  
  ### Upstream of primary TSS
  if(TRUE) {
    ### Extend cage peaks with the wanted window (since cage is extended upstream=downstream)
    cageGeneUpstream <- GenomicRanges::promoters(cageClustersGr, upstream = 0, downstream = as.integer(inputedArguments$geneUpstreamPromoterWindow))
    
    ### Find overlap
    upstreamTSSoverlap <- suppressWarnings( as.data.frame( findOverlaps(cageGeneUpstream, primaryTSS) ))
    setnames(upstreamTSSoverlap, c('cageTSS', 'gene_id'))
    
    # Replace indexes with names
    upstreamTSSoverlap$cageTSS   <- cageGeneUpstream$id [upstreamTSSoverlap$cageTSS  ]
    upstreamTSSoverlap$gene_id   <- primaryTSS$gene_id[upstreamTSSoverlap$gene_id]
    upstreamTSSoverlap$gene_name <- primaryTSS$gene_name[match(upstreamTSSoverlap$gene_id, primaryTSS$gene_id)]
    
    # Handle that some might be associated with multiple transcripts
    upstreamTSSoverlap <- handleMultipleAssociations(upstreamTSSoverlap)
    
    ### annotate in result data.frame
    correspindingMatch <- match(upstreamTSSoverlap$cageTSS, cageClustersDf$name)
    cageClustersDf$tssAnnoation  [correspindingMatch] <- 'upstream'
    cageClustersDf$gene_id       [correspindingMatch] <- upstreamTSSoverlap$gene_id
    cageClustersDf$gene_name     [correspindingMatch] <- upstreamTSSoverlap$gene_name
    
    # remove identified from GRanges
    toKeep <- which( ! cageGeneUpstream$id %in% upstreamTSSoverlap$cageTSS)
    cageClustersGr   <- cageClustersGr  [toKeep,]
  }
}

### Make output
message('Step 3 of 3: Preparing output...')
if(TRUE) {
  ### Notes
  # the output should be a tap seperated file, without header. It has 1 row pr cluster and it has 11 collumns
  # The first 8 collumns correspond to a bed12 format
  # The last 3 are: Classification, Associated transcript, Associated gene.
  
  ### Add comment for which annotation file was used
  cageClustersDfCopy <- cageClustersDf[1:2,]
  cageClustersDfCopy[1,1] <- paste('#', date(), ': The annotation file used was:', inputedArguments$annotationGTF, sep=' ')
  cageClustersDfCopy[2,1] <- paste('#', paste(colnames(cageClustersDf), collapse = '\t'))
  cageClustersDfCopy[1:2,2:ncol(cageClustersDfCopy)] <- ''
  
  cageClustersDf <- rbind(cageClustersDfCopy, cageClustersDf)
  
  ### Write to file
  write.table(cageClustersDf , file=inputedArguments$outputFile, sep = '\t', quote = FALSE, row.names = F, col.names = F)
}

message('The classification of CAGE peaks have now been done.')
