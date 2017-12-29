# s4 Classes
#' For developers.-
#' 
#' All class constructors should start with a '.', and with capital letter.
#' All class constructors must be exported and documented.
#' 
#' All classes should have at least the following slots:
#'  in.files - A vector of characters.
#'  out.files - A list of character vectors, one per each input file if
#'   necessary (i.e. if each input file returns a set of output files). If the
#'   function return one output for the whole run, all the resulting files
#'   should be listed in a list of length 1.
#'  stats - A matrix, data.frame, or table if possible. Something that can be
#'   returned by calling summary() (the generic function). Some computation is
#'   allowed, but nothing too complex. Not yet standarized.
#'  call - A character vector. The call(s) passed to system() or the call 
#'   captured by capture.output(match.call()) at the parent frame of the main
#'   function.
#'  prefix - Character. If missing, it should inherit the prefix of the
#'   previous steps. If there is no previous step and is missing, the time as
#'   implemented in the prokka function.
#'  Other slots are allowed but should be carefully considered. One example is
#'  the panmatrix returned by roary() function. It is not suitable for a 
#'  summary dispatch but is useful for plotting and extracting rich 
#'  information. 
#'  
#' The validation (validity) should evaluate two things:
#' First, if the expected out.files are present in the slot out.files. Second,
#' if the files present in the out.files slot exists in their respective paths.
#' The first validation is to check if the process went all fine. It's 
#' important when creating the final objects of each main functions. The second
#' one is to check if all the out.files of the previous process that may serve
#' as input for the current process still exist at the specified location. Any
#' function that takes as input another 'mirecipe' object should check the
#' validity of this one at the begining of the function using the validObject()
#' function (methods package).


#' @export
.Prokka <- setClass('prokka', 
                    slots = list(in.files = 'character',
                                 out.files = 'list',
                                 stats = 'data.frame',
                                 call = 'character', 
                                 prefix = 'character'), 
                    
                    validity = function(object){
                      
                      exts <- c("err", "faa", "ffn", "fna",
                                "fsa", "gbk", "gff", "log",
                                "sqn", "tbl", "txt")
                      
                      #Check if all expected file extensions exists in object
                      obj <- sapply(object@out.files, function(x){
                        ss <- sapply(strsplit(x, '[.]'), rev, simplify = FALSE)
                        all(sapply(ss, '[', 1) %in% exts)
                      })
                      
                      #Check if all object files exists in specified path
                      fex <- vapply(unlist(object@out.files),
                                    file.exists, 
                                    FUN.VALUE = NA)
                      
                      #Return validity
                      all(obj) & all(fex)
                      
                    }
)



#' @export
.Roary <- setClass('roary', 
                   slots = list(in.files = 'character',
                                out.files = 'list',
                                stats = 'data.frame',
                                call = 'character',
                                prefix = 'character',
                                panmatrix = 'matrix'),
                   
                   validity = function(object){
                     
                     val <- sapply(slot(object, 'out.files'),
                                   sapply,
                                   file.exists, 
                                   simplify = F)
                     
                     all(sapply(val, all))
                     
                   }
)



#' @export
.ProgressiveMauve <- setClass('progressiveMauve', 
                              slots = list(in.files = 'character',
                                           out.files = 'list',
                                           stats = 'matrix',
                                           call = 'character', 
                                           prefix = 'character'),
                              validity = function(object){
                                
                                val <- sapply(slot(object, 'out.files'),
                                              sapply,
                                              file.exists,
                                              simplify = F)
                                
                                all(sapply(val, all))
                                
                              }
)



#' @export
.CoreGenome <- setClass('coreGenome', 
                        slots = list(in.files = 'character',
                                     out.files = 'list',
                                     stats = 'matrix',
                                     call = 'character', 
                                     prefix = 'character'), 
                        validity = function(object){
                          val <- sapply(slot(object, 'out.files'), 
                                        sapply, 
                                        file.exists, 
                                        simplify = F)
                          all(sapply(val, all))
                        }
)