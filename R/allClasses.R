## Experimental design class
methods::setClass("experimentalDesign",
                  slots=c(expDes='data.table', id.vars = "character"),
                  validity=function(object){
                      errors <- character()
                      num.entries <- nrow(getExperimentalTable(object))
                      if (num.entries == 0) {
                          msg <- paste("No valid entries in experimental design table.")
                          errors <- c(errors, msg)
                      }
                      if (length(errors) == 0) TRUE else errors                      
                  })

## Define a class that will hold all information about chromosome length
methods::setClass("chromInfo", slots=c(info="data.table"))

## Make binnedBigWig class
methods::setClass("binnedBigWig", slots=c(counts="list",sample.ids="character",bin.ids="list",expDes="experimentalDesign",strand="list"))
