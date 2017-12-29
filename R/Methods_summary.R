# setGeneric("summary")

# summary.foo <- function(object, ...) {
   ## implement summary.foo
#   "I'm a foo"
# }
# setMethod("summary", "foo", summary.foo)

#Just once!
if (!isGeneric('summary'))
  setGeneric('summary', function(object, ...)
    standardGeneric('summary'))



summary.prokka <- function(object, ...){
  
  slot(object, 'stats')
  
}
setMethod('summary', 'prokka', summary.prokka)


summary.roary <- function(object, ...){
  
  slot(object, 'stats')
  
}
setMethod('summary', 'roary', summary.roary)



summary.progressiveMauve <- function(object, ...){
  
  slot(object, 'stats')
  
}
setMethod('summary', 'progressiveMauve', summary.progressiveMauve)



