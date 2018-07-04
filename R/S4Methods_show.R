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
              cat(rev(x)[1], sep = '\n')
            }else{
              cat(x, sep = '\n ')
            }
          }
)



# roary

#' @export
setMethod('show',
          'roary',
          function(object){
            cat('An object of class "roary"\n')
            x <- object@in.files
            p <- paste('Pangenome of', length(x), 'genomes:\n ')
            cat(p)
            if(length(x)>4){
              cat(x[1:2], sep = '\n ')
              cat(' ...\n ')
              cat(rev(x)[1], sep = '\n')
            }else{
              cat(x, sep = '\n ')
            }
          }
)




# progressiveMauve

#' @export
setMethod('show',
          'progressiveMauve',
          function(object){
            cat('An object of class "progressiveMauve"\n')
            x <- object@in.files
            p <- paste('Alignment for', length(x), 'genomes:\n ')
            cat(p)
            if(length(x)>4){
              cat(x[1:2], sep = '\n ')
              cat(' ...\n ')
              cat(rev(x)[1], sep = '\n')
            }else{
              cat(x, sep = '\n ')
            }
          }
)


# Core genome

#' @export
setMethod('show',
          'coreGenome',
          function(object){
            cat('An object of class "coreGenome"\n')
            x <- object@in.files
            p <- paste('Core-genome extracted from the following alignment:\n ')
            cat(p)
            cat(x, sep = '\n')
          }
)




# Gubbins

#' @export
setMethod('show',
          'gubbins',
          function(object){
            cat('An object of class "gubbins"\n')
            x <- object@in.files
            p <- paste('Recombination detected from the following alignment:\n ')
            cat(p)
            cat(x, sep = '\n')
          }
)



# Core clusters

#' @export
setMethod('show',
          'coreClusters', 
          function(object){
            cat('An object of class "coreClusters"\n')
            x <- names(object@in.files$prokka)
            p <- 'Core genome clusters of the following genomes:\n '
            cat(p)
            if(length(x)>4){
              cat(x[1:2], sep = '\n ')
              cat(' ...\n ')
              cat(rev(x)[1], sep = '\n')
            }else{
              cat(x, sep = '\n ')
            }
          })



# Aligned core clusters

#' @export
setMethod('show',
          'alignedCoreClusters', 
          function(object){
            cat('An object of class "alignedCoreClusters"\n')
            x <- object@in.files$genes
            p <- 'Aligned core clusters:\n '
            cat(p)
            if(length(x)>4){
              cat(x[1:2], sep = '\n ')
              cat(' ...\n ')
              cat(rev(x)[1])
            }else{
              cat(x, sep = '\n ')
            }
          })











