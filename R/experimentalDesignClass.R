#' Create experimental design object
#'
#' This function creates an experimental design object that can be used to manage data
#' @param experiment.id an identifier for the sample
#' @param strand The strand the data comes from +, - , or * for unknown or irrelevant
#' @param filepath where the files are located
#' @param metadata a data.frame containing all other variables of interest
#' @param id.vars a character vector specifying which variables in the metadata should be interperted as experimental conditions
#' @name experimentalDesign
#' @importFrom data.table :=
#' @export 
experimentalDesign <- function(experiment.id,filepath,strand=NULL,metadata=NULL,id.vars=NULL){
    ## Some checks and default settings
    if(length(filepath)!=length(experiment.id)){
        stop("There must be the same number of ids and filepaths.")
    }
    if(is.null(strand)){
        strand=rep("*",length(experiment.id))
    }
    if(!is.null(metadata) & nrow(metadata) != length(experiment.id)){
        stop("If metadata exists it must be the same length as ids")
    }
    if(!is.null(id.vars) & sum(!id.vars %in% colnames(metadata))>0 ){
        stop("All variables specified in id.vars must be present in metadata.")
    }
    
    ## Now build the experimental design table
    d=data.table::data.table(experiment.id=experiment.id,strand=strand,filepath=filepath)
    cols=character(0)
    if(!is.null(metadata)){
        d=cbind(d,metadata)
        cols=c(id.vars)
    }
    cols=c(cols,"strand")
    d[,gid:=apply(.SD,1,paste0,collapse=":"), .SDcols=cols]
    methods::new("experimentalDesign",expDes=d,id.vars=id.vars)
}

#' @rdname getFilepaths-methods
#' @name getFilepaths
methods::setMethod("getFilepaths",signature=c(expDes="experimentalDesign"),definition=function(expDes){
    return(expDes@expDes$filepath)
})

#' @rdname getIds-methods
#' @name getIds
methods::setMethod("getIds",signature=c(expDes="experimentalDesign"),definition=function(expDes){
    return(expDes@expDes$experiment.id)
})

#' @rdname getStrand-methods
#' @name getStrand
methods::setMethod("getStrand",signature=c(expDes="experimentalDesign"),definition=function(expDes){
    return(expDes@expDes$strand)
})

#' @rdname getGroupIds-methods
#' @name getGroupIds
methods::setMethod("getGroupIds",signature=c(expDes="experimentalDesign"),definition=function(expDes){
    return(expDes@expDes$gid)
})

#' @rdname getTable-methods
#' @name getTable
methods::setMethod("getTable",signature=c(expDes="experimentalDesign"),definition=function(expDes){
    return(expDes@expDes)
}) 

#' @rdname subset-methods
#' @name subset
methods::setMethod("subset",signature=c(expDes="experimentalDesign",filters="list"),definition=function(expDes,filters){
    if(sum(!names(filters) %in% colnames(getTable(expDes)))>0){
        stop(paste("Invalid filters present:",names(filters)[!names(filters) %in% colnames(getTable(expDes))]))
    }
    sub=getTable(expDes)
    ind = !logical(nrow(sub))
    for(n in names(filters)){
        ind = ind & sub[[n]] %in% filters[[n]]
    }
    return(methods::new("experimentalDesign",expDes=sub[ind,],id.vars=expDes@id.vars))    
})
