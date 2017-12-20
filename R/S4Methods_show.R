# setMethod("show",
#           "student",
#           function(object) {
#             cat(object@name, "\n")
#             cat(object@age, "years old\n")
#             cat("GPA:", object@GPA, "\n")
#           }
# )

# prokka

#' @export
setMethod('show',
          'prokka',
          function(object){
            cat('An object of class "prokka"\n')
            x <- object@in.files
            p <- paste('Annotation for', length(x), 'genomes:\n ')
            cat(p)
            if(length(x)>4){
              cat(x[1:2], sep = '\n ')
              cat(' ...\n ')
              cat(rev(x)[1])
            }else{
              cat(x, sep = '\n ')
            }
          }
)
