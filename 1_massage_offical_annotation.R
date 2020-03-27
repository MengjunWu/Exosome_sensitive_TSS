#!/usr/local/bin/Rscript
# NOTE: run script from command with -h option to get help page

###################################################
########### Sandelin Lab Cage Pipeline ############
######## Annotation of CAGE peaks (1 of 2) ########
## Script for preparing the official annotation ###
############ Kristoffer Vitting-Seerup ############
###################################################

### Idea 
# This script should take the path to a gtf file as input.
# the GTF should then be massaged into the GRanges needed for the annotation and 
# save them to an Rdata file which will be loaded by the script which annotates the
# cage data

### Use argparse to make help file and parse input
if(TRUE) {
    suppressMessages( library('argparser') )
    
    ### Use argparse to make help file
    # create parser
    argParse <- arg_parser(name = '1_massage_official_annotation.R', description = 'Description:\nPart 1 (of 2) of the Sandelin Lab 
                            Cage Pipeline Annoation scripts. This script that massage the offical gene annoation so it can be used with 
                            part 2. This script produces two files with the same name as the GTF file, except they are called \'.Rdata\' 
                            and \'.parsed.Rdata\' instead of \'.gtf\'. The second file is used for the CAGE cluster annotation by the 
                            second script. The first file have two functions: First it allows for very fasts loading of the GRange version 
                            of the GTF file. Secondly it allows the user to manually parse non-standard GTF files and save them to an Rdata 
                            file, which is then used instead of the GTF file. The details for the Rdata version of the GTF file are as 
                            follows. It must be a standard GRanges object with chr start end and strand. Furthermore it must have 4 metadata 
                            collumns: 1) A collum called \'type\' indicating whether it is a Exon, a CDS or a UTR regions (other are ignored). 
                            2) A collumn called \'gene_id\' giving the unique gene id (fx ENSG00001. 3) A column called \'gene_name\' giving 
                            the gene name (Fx Rac1) and 4) a collumn called \'transcript_id\' with the unique transcript id. Lastly note that 
                            the GRange object must be called orgGTF.'
                           )
    
    # add arguments
    argParse <- add_argument(parser = argParse, short = '-g', arg = '-gtf',                           help = 'The FULL path to the GTF containing the official annotation.')
    argParse <- add_argument(parser = argParse,               arg = '-forceGTFimport', default=FALSE, help = 'A logic (TRUE/FALSE) which incates whether to import the GTF file even though a \".Rdata\" file containing the corresponding GRanges exists.')
    
    # evaluate input
    inputedArguments <- parse_args(argParse, argv = commandArgs(trailingOnly = TRUE)) # creates a named list with arguments
    
    ### Note: Currently argparser ( v 0.1 ) gives a warning message for arguments starting with certain types of names
    
    ### Test input
    inputTest <- sapply(inputedArguments[-c(1:2)], is.na) # first is always help and opts 
    if( any( inputTest ) ) {
        stop(paste('The following input arguments are missing: -', paste( names(inputTest)[which(inputTest)], collapse = ', -'), sep='') )
    }
 
    
}

### For devel
if(FALSE) {
    ### Manually create list with arguments
    inputedArguments <- list(
        gtf='/Volumes/BINF-Sandelin/projects/CAGE/genome/annotations/ASM294v2_r26/Schizosaccharomyces_pombe.ASM294v2.26.gtf'
        )
    
    '/Volumes/BINF-Sandelin/projects/CAGE/genome/annotations/ASM294v2_r26/Schizosaccharomyces_pombe.ASM294v2.26.'
    
    load('/Volumes/BINF-Sandelin/people/rtl144/test/gencodev19hg19.Rdata')# orgGTF
}

### Load dependencies
message('Loading dependencies...')
suppressMessages( library('rtracklayer')   )
suppressMessages( library('GenomicRanges') )
suppressMessages( library('XVector')       )
suppressMessages( library('plyr')          )

### Read in gtf file
if(TRUE) {
    gtfAsGRanges <- paste( gsub('.gtf$', '', inputedArguments$gtf, perl = T, ignore.case = T), '.Rdata', sep='')
    
    if( file.exists(gtfAsGRanges) & !as.logical(inputedArguments$forceGTFimport)) {
        
        message('Step 1 of 2: Importing the GRange format of the GTF file from the .Rdata file')
        load(gtfAsGRanges)
        
    } else {
        message('Step 1 of 2: Importing GTF. This may take a while...')
        
        orgGTF <- import(con = inputedArguments$gtf, format = 'gtf')
        orgGTF <- sort(orgGTF)
        
        # save for (potential) later usage
        save(orgGTF, file=gtfAsGRanges)
    }
    
    ### Test input
    if( ! all(c('gene_name','transcript_id') %in% colnames(mcols(orgGTF))) ) {
        stop('The gene name and transcript name was not found - please massage manually and save the Rdata file yourself')
    }
    
    ### Massage
    # subset to exon cds and utrs and only extract meta collums i need
    orgGTF <- orgGTF[which( tolower(orgGTF$type) %in% c('exon','cds','utr')) , c('type','gene_id','gene_name','transcript_id')]
    
    ### The GRanges that is imported should have the following collums
    # seqnames (chromosome name), ranges (start stop), stand - corresponding to standard GTF file format
    # 4 meta data collums called:
    # 'type' - containing basic annoation, exon cds, UTR etc - standard for GTF file
    # 'gene_id' - a unique id for each gene
    # 'gene_name' - containing the official gene name (hgnc, mgi etc)
    # 'transcript_id' - containing unique name of the transcript.
    
}

### Modify the GTF file to get a set of GRanges shaped excatly like I need (to avoid havint to do this on the each time)
message('Step 2 of 3: Massaging GTF. Usually takes a couple of minuts...')
if(TRUE) {
    ### Base gene and transcript regions 
    if(TRUE) {
        ### Gene edges - takes 6 secs
        geneEdges <- suppressMessages( unlist( range( split(orgGTF[,0], f=orgGTF$gene_id) ) ) )
        geneEdges$gene_id <- names(geneEdges)
        names(geneEdges) <- NULL
        geneEdges$gene_name <- orgGTF$gene_name[match(geneEdges$gene_id , orgGTF$gene_id)]
        
        ### Transcript edges - takes 16 secs
        transcriptEdges <- unlist( range(split(orgGTF[,0], f=orgGTF$transcript_id)) )
        transcriptEdges$transcript_id <- names(transcriptEdges)
        names(transcriptEdges) <- NULL
        transcriptEdges$gene_id   <- orgGTF$gene_id  [match(transcriptEdges$transcript_id, orgGTF$transcript_id)]
        transcriptEdges$gene_name <- orgGTF$gene_name[match(transcriptEdges$transcript_id, orgGTF$transcript_id)]
    }
    
    ### TSS
    if(TRUE) {
        ### Use gene edges to get most upstream TSS - takes < 1 sec
        primaryTSS <- GenomicRanges::promoters( geneEdges, upstream =0 , downstream = 1 )
        
        ### Extract all TSS
        transciptTSS <- GenomicRanges::promoters(transcriptEdges , upstream =0 , downstream = 1 )
        
        ### Remove those overlapping with primary TSS
        alternativeTSS <- transciptTSS[ which( ! overlapsAny(query = transciptTSS, subject = primaryTSS) ), ]
        ### Get primary TSS with annoation
        primaryTSS     <- transciptTSS[ which(   overlapsAny(query = transciptTSS, subject = primaryTSS) ), ]
        
    }
    
    ### Coding regions
    if(TRUE) {
        # takes 3 secs
        codingRegions <- orgGTF[which(orgGTF$type == 'CDS'),'gene_id']
        codingRegions <- unlist(reduce(split( codingRegions, f=codingRegions$gene_id)))
        codingRegions$gene_id <- names(codingRegions)
        names(codingRegions) <- NULL
        codingRegions$gene_name <- orgGTF$gene_name[match(codingRegions$gene_id , orgGTF$gene_id)]
    }
    
    ### UTR regions - Deviding into 5' and 3' have to be done one transcript at the time.
    if(TRUE) {
        # Takes the longest time - 4 min

        ### devide into 5' and 3' - this is done with data.frames since applying over these are much faster
        utrRegionsDF <- as.data.frame( orgGTF[which(orgGTF$type %in% c('CDS','UTR')) , c('type','transcript_id')] )
        utrRegionsDF$type <- as.vector(utrRegionsDF$type) # unfactor
        suppressMessages( 
            utrRegionsDF <- ddply(utrRegionsDF, .variables = 'transcript_id', .progress = 'none', .fun = function(aDF) {
                isPlusStrand <- aDF$strand[1] == '+'
                
                # Use a Rle object to change the UTR to 5UTR and 3UTR respectively
                myRle <- Rle(aDF$type)
                
                if(myRle@values[1] == 'UTR') {
                    if(isPlusStrand) {
                        myRle@values[1] <- '5UTR'
                    } else {
                        myRle@values[1] <- '3UTR'
                    }
                }
                if( myRle@values[length(myRle@values)] == 'UTR') {
                    if(isPlusStrand) {
                        myRle@values[length(myRle@values)] <- '3UTR'
                    } else {
                        myRle@values[length(myRle@values)] <- '5UTR'
                    }
                }
                
                # overwrite in data.frame
                aDF$type <- as.vector(myRle)
                
                # subset to only UTR regions
                aDF <- aDF[which(aDF$type != 'CDS'),]
                
                return(aDF)
            }) 
        )
        
        ### convert back to GRanges and add gene info
        utrRegionsGr <- GRanges(utrRegionsDF$seqnames, IRanges(utrRegionsDF$start, utrRegionsDF$end), strand=utrRegionsDF$strand, type=utrRegionsDF$type, transcript_id=utrRegionsDF$transcript_id)
        utrRegionsGr$gene_id <- orgGTF$gene_id[match(utrRegionsGr$transcript_id , orgGTF$transcript_id)]
        
        ### devide into 5UTR and 3UTR
        utr5reg <- utrRegionsGr[which(utrRegionsGr$type == '5UTR'),]
        utr3reg <- utrRegionsGr[which(utrRegionsGr$type == '3UTR'),]
        
        ### reduce on gene level
        utr5regGene <- sort( unlist( reduce( split(utr5reg, f=utr5reg$gene_id) )) )
        utr5regGene$gene_id <- names(utr5regGene)
        utr5regGene$gene_name <- orgGTF$gene_name[match( utr5regGene$gene_id, orgGTF$gene_id)]
        names(utr5regGene) <- NULL
        
        utr3regGene <- sort( unlist( reduce( split(utr3reg, f=utr3reg$gene_id) )) )
        utr3regGene$gene_id <- names(utr3regGene)
        utr3regGene$gene_name <- orgGTF$gene_name[match( utr3regGene$gene_id, orgGTF$gene_id)]
        names(utr3regGene) <- NULL
        
        ### remove those parts overlapping with CDS
        #utr5regGene <- GenomicRanges::setdiff(utr5regGene, codingRegions)
        #utr3regGene <- GenomicRanges::setdiff(utr3regGene, codingRegions)
        
    }
    
    ### Exons
    if(TRUE) {
        # takes 8 sec
        reducedGeneList <- reduce(split(orgGTF[,0], f=orgGTF$gene_id))
        geneExons <- unlist(reducedGeneList)
        geneExons$gene_id <- names(geneExons)
        geneExons$gene_name <- orgGTF$gene_name[match( geneExons$gene_id, orgGTF$gene_id)]
        names(geneExons) <- NULL
    }
    
    ### Intron
    if(TRUE) {
        # takes < 1 sec
        geneIntronsIRanges <- unlist( gaps(ranges( reducedGeneList )) )
        
        # convert the IRanges back to GRanges
        matchingGRangesIndex <- match(names(geneIntronsIRanges), orgGTF$gene_id)
        geneIntrons <- GRanges(
            seqnames=seqnames(orgGTF)[matchingGRangesIndex],
            ranges=geneIntronsIRanges,
            strand=strand(orgGTF)[matchingGRangesIndex],
            gene_id  =orgGTF$gene_id[matchingGRangesIndex],
            gene_name=orgGTF$gene_name[matchingGRangesIndex]
            )
        
        names(geneIntrons) <- NULL
        geneIntrons <- sort(geneIntrons)
        
        ### remove those parts overlapping with exons
        #geneIntrons <- GenomicRanges::setdiff(geneIntrons, geneExons)
    }
    
}

message('Step 3 of 3: Preparing output...')
### Save all files in a Rdata object
if(TRUE) {
    ### Sanity check that it worked
    testList <- list(
        geneEdges=geneEdges,
        transcriptEdges=transcriptEdges,
        primaryTSS=primaryTSS,
        alternativeTSS=alternativeTSS,
        utr5regGene=utr5regGene,
        codingRegions=codingRegions,
        utr3regGene=utr3regGene,
        geneExons=geneExons,
        geneIntrons=geneIntrons
        )
    grLength <- sapply(testList, length)
    
    isOfLengthNull <- names(grLength)[which( grLength == 0)]
    
    if( length(isOfLengthNull) ) {
        warning('The parsing of the GTF migh not be sucessfull. The folloing subsets were empty:', paste(isOfLengthNull, collapse = ','))
    }
    
    
    ### Create path
    parsedGRangesPath <- paste( gsub('.gtf$', '', inputedArguments$gtf, perl = T, ignore.case = T), '.parsed.Rdata', sep='')
    save(
        geneEdges,
        transcriptEdges,
        primaryTSS,
        alternativeTSS,
        utr5regGene,
        codingRegions,
        utr3regGene,
        geneExons,
        geneIntrons,
        file=parsedGRangesPath
        )
    
}


#