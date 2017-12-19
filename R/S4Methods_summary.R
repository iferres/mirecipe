# setGeneric("summary")

# summary.foo <- function(object, ...) {
   ## implement summary.foo
#   "I'm a foo"
# }
# setMethod("summary", "foo", summary.foo)

#Just once!
setGeneric('summary')


#Set S3 Method
summary.prokka <- function(object){
  
  object@stats
  
}

#Set S4 Method
setMethod('summary', 'prokka', summary.prokka)