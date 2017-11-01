## Splits data from binnedBigWig into list by a column in the experimental design table
methods::setGeneric("splitBw", function(binnedBw,f,to.matrix) {
    methods::standardGeneric("splitBw")
})

## Splits binnedBigWig into binnedBigWig List by chain so it can be easily iterated over in parallel
methods::setGeneric("splitBwByChain", function(binnedBw) {
    methods::standardGeneric("splitBwByChain")
})

## Generic methods
methods::setGeneric("getFilepaths",valueClass="character",function(expDes){ methods::standardGeneric("getFilepaths") })
methods::setGeneric("getIds",valueClass="character",function(expDes){ methods::standardGeneric("getIds") })
methods::setGeneric("getGroupIds",valueClass="character",function(expDes){ methods::standardGeneric("getGroupIds") })
methods::setGeneric("getStrand",valueClass="character",function(expDes){ methods::standardGeneric("getStrand") })
methods::setGeneric("getTable",valueClass="data.frame",function(expDes){ methods::standardGeneric("getTable") })
methods::setGeneric("subset",valueClass="experimentalDesign",function(expDes,filters){ methods::standardGeneric("subset") })


