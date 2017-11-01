#' Split Data
#' 
#' Splits data from binnedBigWig into list by a column in the experimental design table.
#' @param binnedBw binnedBigWig object
#' @param f value to split binnedBigWig on
#' @param to.matrix LOGICAL. Should the output be converted from data.table to matrix form.
#' @docType methods
#' @name splitBw
#' @rdname splitBw-methods
#' @export
methods::setGeneric("splitBw", function(binnedBw,f,to.matrix) {
    methods::standardGeneric("splitBw")
})

#'  Split and clone binnedBigWig
#' 
#' Creates seperate binnedBigWig objects for each element of @counts list in masted binnedBigWig object. 
#' @param binnedBw binnedBigWig object
#' @docType methods
#' @name seperateBinnedBwElements
#' @rdname seperateBinnedBwElements-methods
#' @export
methods::setGeneric("seperateBinnedBwElements", function(binnedBw) {
    methods::standardGeneric("seperateBinnedBwElements")
})

#' Accessor function for experimental design object
#' 
#' Return all file paths in experimentalDesign object
#' @docType methods
#' @name getFilepaths
#' @rdname getFilepaths-methods
#' @export
methods::setGeneric("getFilepaths",valueClass="character",function(expDes){ methods::standardGeneric("getFilepaths") })

#' Accessor function for experimental design object
#' 
#' Return all experiment.ids in experimentalDesign object
#' @docType methods
#' @name getIds
#' @rdname getIds-methods
#' @export
methods::setGeneric("getIds",valueClass="character",function(expDes){ methods::standardGeneric("getIds") })

#' Accessor function for experimental design object
#' 
#' Return all strands associated with experimental data
#' @docType methods
#' @name getGroupIds
#' @rdname getGroupIds-methods
#' @export
methods::setGeneric("getGroupIds",valueClass="character",function(expDes){ methods::standardGeneric("getGroupIds") })

#' Accessor function for experimental design object
#'
#' Return all group ids (constructed based in shared id.vars). All entries with the same group.id are replicates for the given set of experimental conditions
#' @docType methods
#' @name getStrand
#' @rdname getStrand-methods
#' @export
methods::setGeneric("getStrand",valueClass="character",function(expDes){ methods::standardGeneric("getStrand") })

#' Accessor function for experimental design object
#' 
#' Return full experimental design table.
#' @docType methods
#' @name getTable
#' @rdname getTable-methods
#' @export
methods::setGeneric("getTable",valueClass="data.frame",function(expDes){ methods::standardGeneric("getTable") })

#' Subset experiments
#' 
#' Method for producing filtered experimental design object
#' @docType methods
#' @name subset
#' @rdname subset-methods
#' @export
methods::setGeneric("subset",valueClass="experimentalDesign",function(expDes,filters){ methods::standardGeneric("subset") })


