#' @export
setClass('roary', 
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


#' @importFrom methods new
#' @export
roary <- function(prokka, 
                  out.dir = '.', 
                  prefix, 
                  cpus = 1L, 
                  cd = 99, 
                  i = 95){
  
  if (Sys.which('roary')==""){
    
    stop('roary is not in your $PATH.')
    
  }
  
  cl <- class(prokka)
  if (cl != 'prokka' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "prokka" is not of class "prokka" (mirecipe package).')
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  in.files <- vapply(slot(prokka, 'out.files'), function(x){
    grep('[.]gff$', x, value = TRUE)
  }, FUN.VALUE = '')
  
  if (!all(file.exists(in.files))){ 
    
    stop("cannot open the connection. One or more files doesn't exist.")
    
  }
  
  in.files <- normalizePath(in.files)
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- slot(prokka, 'prefix')
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_roary')
  
  # if(!dir.exists(dout)){
  #   dir.create(dout)
  # }
  
  out.files <- vector('list', 1L)
  
  run <- paste0('roary -cd ', cd, 
                ' -i ', i, 
                ' -p ', cpus,
                ' -f ', dout, 
                ' ', paste(in.files, collapse = ' '))

  system(run)
  
  out.files[[1L]] <- normalizePath(list.files(path = dout, 
                                        full.names = TRUE))
  #Stats
  gp <- grep('summary_statistics.txt', 
             out.files[[1L]], 
             fixed = TRUE, 
             value = TRUE)
  stats <- read.csv(gp, sep = '\t', header = FALSE)[, c(1,3)]
  colnames(stats) <- c('Description', 'Number')
  
  #Panmatrix
  gp <- grep('gene_presence_absence.csv', 
             out.files[[1L]], 
             fixed = TRUE, 
             value = TRUE)
  
  pm <- read.csv(gp)
  rownames(pm) <- pm$Gene
  pm <- pm[, c(15:dim(pm)[2])]
  
  panmatrix <- t(apply(pm, 2, function(j){
    
    vapply(j, function(k){ length(strsplit(k, '\t')[[1]]) }, 1L)
    
  }))
  
  res <- new('roary', 
             in.files = in.files, 
             out.files = out.files, 
             stats = stats, 
             call = run,
             prefix = prefix, 
             panmatrix = panmatrix)
  
}