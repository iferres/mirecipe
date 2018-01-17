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
                      
                      esp <- c("err", "faa", "ffn", "fna",
                               "fsa", "gbk", "gff", "log",
                               "sqn", "tbl", "txt")
                      
                      #Check if all expected file extensions exists in object
                      obj <- sapply(object@out.files, function(x){
                        ss <- sapply(strsplit(x, '[.]'), rev, simplify = FALSE)
                        all(sapply(ss, '[', 1) %in% esp)
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
                     
                     esp <- c("accessory_binary_genes.fa", 
                              "accessory_binary_genes.fa.newick", 
                              "accessory_graph.dot", 
                              "accessory.header.embl", 
                              "accessory.tab", 
                              "blast_identity_frequency.Rtab", 
                              "clustered_proteins", 
                              "core_accessory_graph.dot", 
                              "core_accessory.header.embl",
                              "core_accessory.tab",
                              "gene_presence_absence.csv", 
                              "gene_presence_absence.Rtab", 
                              "number_of_conserved_genes.Rtab", 
                              "number_of_genes_in_pan_genome.Rtab", 
                              "number_of_new_genes.Rtab", 
                              "number_of_unique_genes.Rtab", 
                              "summary_statistics.txt")
                     
                     #Check if all expected files exists in object
                     obj <- sapply(object@out.files[[1]], function(x){
                       ss <- sapply(strsplit(x, '/'), rev, simplify = FALSE)
                       sapply(ss, '[', 1) %in% esp
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
.ProgressiveMauve <- setClass('progressiveMauve', 
                              slots = list(in.files = 'character',
                                           out.files = 'list',
                                           stats = 'matrix',
                                           call = 'character', 
                                           prefix = 'character'),
                              validity = function(object){
                                
                                esp <- c("alignment.xmfa", 
                                         "alignment.xmfa.backbone", 
                                         "alignment.xmfa.bbcols")
                                
                                #Check if all expected files exists in object
                                obj <- sapply(object@out.files[[1]], function(x){
                                  ss <- sapply(strsplit(x, '/'), rev, simplify = FALSE)
                                  sapply(ss, '[', 1) %in% esp
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
.CoreGenome <- setClass('coreGenome', 
                        slots = list(in.files = 'character',
                                     out.files = 'list',
                                     stats = 'matrix',
                                     call = 'character', 
                                     prefix = 'character'), 
                        validity = function(object){
                          
                          esp <- c("coreAliCoord.txt",
                                   "coreGenomeAlignment.fasta")
                          
                          #Check if all expected files exists in object
                          obj <- sapply(object@out.files[[1]], function(x){
                            ss <- sapply(strsplit(x, '/'), rev, simplify = FALSE)
                            sapply(ss, '[', 1) %in% esp
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
.Gubbins <- setClass('gubbins', 
                      slots = list(in.files = 'character',
                                   out.files = 'list',
                                   stats = 'matrix',
                                   call = 'character', 
                                   prefix = 'character'), 
                      validity = function(object){
                        
                        esp <- c(".recombination_predictions.embl",
                                 ".recombination_predictions.gff",
                                 ".branch_base_reconstruction.embl",
                                 ".summary_of_snp_distribution.vcf",
                                 ".per_branch_statistics.csv",
                                 ".filtered_polymorphic_sites.fasta",
                                 ".filtered_polymorphic_sites.phylip",
                                 ".final_tree.tre",
                                 ".node_labelled.final_tree.tre")
                          
                        #Check if all expected files exists in object
                        obj <- sapply(object@out.files[[1]], function(x){
                          ss <- sapply(strsplit(x, '/'), rev, simplify = FALSE)
                          paste0(prefix, sapply(ss, '[', 1)) %in% esp
                        })
                          
                        #Check if all object files exists in specified path
                        fex <- vapply(unlist(object@out.files),
                                      file.exists, 
                                      FUN.VALUE = NA)
                          
                        #Return validity
                        all(obj) & all(fex)
                      }
)











