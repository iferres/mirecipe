# setGeneric("summary")

# summary.foo <- function(object, ...) {
   ## implement summary.foo
#   "I'm a foo"
# }
# setMethod("summary", "foo", summary.foo)

#Just once!
# setGeneric('summary')


#Set S3 Method
#' @export
summary.prokka <- function(object){
  
  slot(object, 'stats')
  
}

# #Set S4 Method
# setMethod('summary', 'prokka', summary.prokka)


#' @export
summary.roary <- function(object){
  
  slot(object, 'stats')
  
}

#' @export
summary.progressiveMauve <- function(object){
  
  slot(object, 'stats')
  
}


# print.summary.progressiveMauve <- function(object){
#   
#   cat('Percentage of shared genome between rows and columns.\n')
#   object
#   
# }



