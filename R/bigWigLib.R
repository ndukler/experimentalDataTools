#' @rdname splitBw-methods
#' @name splitBw
methods::setMethod("splitBw",signature=c(binnedBw="binnedBigWig",f="character",to.matrix="logical"),definition = function(binnedBw,f="gid",to.matrix=FALSE){
    interTab=getExperimentTable(binnedBw@expDes)
    data.table::setkeyv(interTab,cols=f)
    l=unique(interTab[[f]])
    return(lapply(binnedBw@counts, function(x){
        if(to.matrix){
            out=lapply(as.list(l), function(cat) as.matrix(x[,interTab[cat]$experiment.id,with=FALSE]))
        } else {
            out=lapply(as.list(l), function(cat) x[,interTab[cat]$experiment.id,with=FALSE])
        }
        names(out)=l
        return(out)
    }))
})

#' Chromosome information collection
#' 
#' Collects the chromosomal lengths used to build all bigwigs experimental design object
#' @param expDes experimentalDesign object
#' @param which.chrom character vector of which chromosomes information should be collected on. Default value of NULL collects info for all chromosomes. Cannot be specified at the same time as regex.chrom.
#' @param regex.chrom can supply regex expression matching which chromosomes information should be collected on. Cannot be specified at the same time as which.chrom.
#' @name getChromInfo
#' @import rtracklayer
#' @export
#'
getChromInfo <- function (expDes, which.chrom = NULL, regex.chrom = NULL) 
{
    ci = rtracklayer::seqinfo(BigWigFile(expDes@expDes$filepath[1]))
    chromInfo = data.table::data.table(chrom = ci@seqnames, length = ci@seqlengths)
    data.table::setkey(chromInfo, "chrom")
    if (!is.null(which.chrom)) {
        chromInfo = chromInfo[which.chrom]
    }
    if (!is.null(regex.chrom)) {
        chromInfo = chromInfo[grep(regex.chrom, chrom)]
    }
    if(any(is.na(chromInfo$length))){
        missing.chrom=chromInfo[is.na(length)]$chrom
        stop(paste("Bigwig file is missing chromosomes present in the query file:",paste(missing.chrom,collapse=",")))
    }
    if (length(getIds(expDes)) > 1) {
        for (i in 2:nrow(expDes@expDes)) {
            cur = rtracklayer::seqinfo(BigWigFile(expDes@expDes$filepath[i]))
            cur = data.table::data.table(chrom = cur@seqnames, 
                length = cur@seqlengths)
            data.table::setkey(cur, "chrom")
            if (!is.null(which.chrom)) {
                cur = cur[which.chrom]
            }
            if (!is.null(regex.chrom)) {
                cur = cur[grep(regex.chrom, cur$chrom)]
            }
            if (any(cur$length != chromInfo$length)) {
                stop(paste("Bigwig", getFilepaths(expDes)[i], 
                  "was built with a different chromInfo file than", 
                  getFilepaths(expDes)[i]))
            }
        }
    }
    return(methods::new("chromInfo", info = chromInfo))
}

#' Imports bigWig(s)
#' 
#' Import bigwig files for later data retrieval. 
#' @param expDes experimentalDesign object
#' @param gReg.gr either a GRangesList or GRanges object containing genomic coordinates to import
#' @param as.type how to return the imported data
#' @param nthreads number of cores to use to import data (unimplemented)
#' @name importBwSelection
#' @import rtracklayer
#' @import doParallel
#' @import foreach
#' @export
importBwSelection <- function(expDes,gReg.gr,as.type='RleList',nthreads=1){
    fExist=file.exists(getFilepaths(expDes))
    if(sum(!fExist)>0){
        stop(paste(paste(getFilepaths(expDes)[!fExist],collapse=","),"do not exist."))
    }
    bwList=list()
    if(is.list(gReg.gr)){
        bwSel=IRanges::reduce(Reduce("c",gReg.gr))
    } else {
        bwSel=IRanges::reduce(gReg.gr)
    }
    ## Now check validity of query
    cInfo=getChromInfo(expDes,which.chrom=levels(GenomicRanges::seqnames(bwSel)))
    bwSel=BigWigSelection(bwSel)
    ## cl <- makeCluster(nthreads)
    ## registerDoParallel(cl)
    ## Loop over each id
    bwList <- foreach(i=getIds(expDes),.final = function(x) setNames(x,getIds(expDes)),.export=c("subsetExperiments","getExperimentTable","getFilepaths")) %do% {
        suppressMessages(library(data.table,quietly=TRUE,verbose=FALSE))
        sub.e=subsetExperiments(expDes,filters=list(experiment.id=i))
        bw=list()
        ## If that id has multiple strands associated with it you can capture all of them
        for(j in 1:nrow(getExperimentTable(sub.e))){
            write(paste("Reading in",getFilepaths(sub.e)[j],"..."),stderr())
            ## Try to load bigWigs and catch and report error if fails
            bw[[getStrand(sub.e)[j]]]=tryCatch({
                rtracklayer::import.bw(con=getFilepaths(sub.e)[j],selection=bwSel,as='RleList')[cInfo@info$chrom]
            }, warning = function(w) {
                print(w)
            }, error = function(e) {
                print(e)
                print(paste0("Unable to load ",getFilepaths(sub.e)[j],". Check error log for details"))
            })            
            write(paste("Finished importing:",getFilepaths(sub.e)[j]),stderr())
        }
        return(bw)
    }
    return(bwList)
}

#'  Imports bigWig(s)
#' 
#' Get read counts in bw regions for given list of genomic ranges. If a strand is specified for bins, only bigwigs matching that strand will be queried. If "*" strand, both "+" and "-" strands will be queried
#' @param bins GRangesList object
#' @param expDes experimentalDesign object
#' @param nthreads number of cores to use to import data (not implemented)
#' @name sumBwOverGR
#' @import rtracklayer
#' @export
sumBwOverGR <- function(bins,expDes,nthreads=1){
    ## Check to see if tu.list is GRanges object, then split into stand based GRanges objects
    if(class(bins)[1]!="GRangesList"){
        stop("Bins must be in the form of a GRangesList object")
    }
    ## Get the number of rows in each element of the GRanges list
    list.lengths=unlist(lapply(bins,length))
    ## save the bin names
    txids=lapply(bins,names)
    st=lapply(bins,GenomicRanges::strand) 
    ## Now combine the GRangesList into a single GRanges for efficient querying
    bins=unlist(bins)
    ## Reset bin names to guarentee that data returned in same order as queried in bins
    names(bins)=1:length(bins)
    ## Create matrix to place counts in
    rep.count.matrix=data.table::data.table(txid=1:length(bins))
    data.table::setkey(rep.count.matrix,"txid")
    ## split bins by strand
    bins=split(bins,rtracklayer::strand(bins))
    ## Get reads in each region per sample
    for(id in getIds(expDes)) {
        write(paste("Reading data from bigwig",id,"into memory..."),stdout())
        bw.list=importBwSelection(subsetExperiments(expDes,filters=list(experiment.id=id)),gReg.gr=IRanges::reduce(Reduce("c",bins)))
        write("Summing reads...",stdout())
        rep.count.matrix[,foo:=0]
        data.table::setnames(rep.count.matrix,"foo",id)
        ## Loop over chromosomes on both strands
        if(!is.null(bins[["+"]]) && length(bins[["+"]])>0){
            ## If the TU is on the plus stand get only reads from the plus strand
            for(chr in unique(GenomicRanges::seqnames(bins[["+"]]))) {
                read.sum=sum(Views(bw.list[[id]][["+"]][[as.character(chr)]], IRanges::ranges(bins[["+"]][GenomicRanges::seqnames(bins[["+"]]) == chr])))
                rep.count.matrix[names(read.sum),id]=as.numeric(read.sum)
            }
        }
        ## If the TU is on the minus stand get only reads from the minus strand
        if(!is.null(bins[["-"]]) && length(bins[["-"]])>0){
            for(chr in unique(GenomicRanges::seqnames(bins[["-"]]))) {
                read.sum=sum(IRanges::Views(bw.list[[id]][["-"]][[as.character(chr)]], ranges(bins[["-"]][GenomicRanges::seqnames(bins[["-"]]) == chr])))
                rep.count.matrix[names(read.sum),id]=as.numeric(read.sum)
            }
        }
        ## If the TU is on the star stand get reads from all strands
        if(!is.null(bins[["*"]]) && length(bins[["*"]])>0){
            for(chr in unique(GenomicRanges::seqnames(bins[["*"]]))) {
                read.sum=numeric(length(bins[["*"]][GenomicRanges::seqnames(bins[["*"]]) == chr]))
                for(s in names(bw.list[[id]])){
                    read.sum=read.sum+abs(sum(IRanges::Views(bw.list[[id]][[s]][[as.character(chr)]], ranges(bins[["*"]][GenomicRanges::seqnames(bins[["*"]]) == chr]))))
                }
                rep.count.matrix[as.numeric(names(read.sum)),id]=as.numeric(read.sum)
            }
        }        
    }
    rep.count.matrix[,txid:=NULL]
    return(new("binnedBigWig",counts=split(rep.count.matrix,rep(1:length(list.lengths),list.lengths)),sample.ids= colnames(rep.count.matrix),
               bin.ids=txids,expDes=expDes,strand=st))
}

#' @rdname seperateBinnedBwElements-methods
#' @name seperateBinnedBwElements
setMethod("seperateBinnedBwElements",signature=c(binnedBw="binnedBigWig"),definition = function(binnedBw){
    ## Experimental design must be copied to all elements of the list
    expD=binnedBw@expDes
    ## Sample ids must be copied to all elements of the list 
    sid=binnedBw@sample.ids
    out <- foreach::foreach(i=1:length(binnedBw@counts)) %do% {
        new("binnedBigWig",counts=list(binnedBw@counts[[i]]),sample.ids=sid,bin.ids=list(binnedBw@bin.ids[[i]]),expDes=expD,strand=list(binnedBw@strand[[i]]))
    }
    return(out)
})

