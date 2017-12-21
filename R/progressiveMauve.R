#' @export
setClass('progressiveMauve', 
         slots = list(in.files = 'character',
                      out.files = 'list',
                      stats = 'data.frame',
                      call = 'character', 
                      prefix = 'character',
                      shared.percent = 'matrix'), 
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

progressiveMauve <- function(prokka, 
                             out.dir = '.', 
                             prefix, 
                             cpus = 1L){
  
  if (Sys.which('progressiveMauve')==""){
    
    stop('progressiveMauve is not in your $PATH.')
    
  }
  
  cl <- class(prokka)
  if (cl != 'prokka' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "prokka" is not of class "prokka" (mirecipe package).')
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  in.files <- vapply(slot(prokka, 'out.files'), function(x){
    grep('[.]fna$', x, value = TRUE)
  }, FUN.VALUE = '')
  
  if (!all(file.exists(in.files))){ 
    
    stop("cannot open the connection. One or more files doesn't exist.")
    
  }
  
  in.files <- normalizePath(in.files)
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- slot(prokka, 'prefix')
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_progressiveMauve')
  
  if(!dir.exists(dout)){
    dir.create(dout)
  }
  
  
  oal <- paste0(dout, '/alignment.xmfa')
  
  run <- paste0('progressiveMauve ',
                '--output=', oal, ' ', 
                paste0(in.files, collapse = ' '))
  
  system(run)
  
  out.files <- vector('list', 1L)
  
  out.files[[1L]] <- normalizePath(list.files(path = dout, 
                                              full.names = TRUE))
  #Stats
  gp <- grep('alignment.xmfa.backbone', 
             out.files[[1L]], 
             fixed = TRUE, 
             value = TRUE)
  bk <- read.csv(gp, sep = '\t', header = TRUE)
  
  m <- pcSharedGenome(prokka = prokka, backbone = bk)
  
  res <- new('progressiveMauve', 
             in.files = in.files, 
             out.files = out.files, 
             stats = stats, 
             call = run,
             prefix = prefix, 
             shared.percent = m)
  
}




#percentage of "shared" (homologous) genome between 2 genomes.
pcSharedGenome <- function(prokka, backbone){
  
  in.prokka <- slot(prokka, 'in.files')

  gnam <- sapply(strsplit(in.prokka, '/'), function(x) x[length(x)] )
  len <- slot(prokka, 'stats')$bases
  names(len) <- gnam
  
  sb <- backbone[, 2*(1:length(in.prokka))]
  
  cb <- combn(length(in.prokka), 2L)
  
  sh <- apply(cb, 2, function(x){
    which(sb[,x[1]]>0 & sb[,x[2]]>0)
    })
  
  su <- abs(backbone[, c(T,F)] - backbone[, c(F,T)])
  
  sha <- lapply(1:length(in.prokka), function(x){
    a <- (sum(su[sh[[x]], cb[1L, x]])/len[ cb[1L, x] ]) * 100 
    b <- (sum(su[sh[[x]], cb[2L, x]])/len[ cb[2L, x] ]) * 100
    c(a, b)
  })
  
  m <- matrix(100, 
              nrow = length(in.prokka), 
              ncol = length(in.prokka))
  colnames(m) <- gnam -> rownames(m)
  
  for (i in 1:length(sha)){
    m[names(sha[[i]])[1], names(sha[[i]])[2]] <- sha[[i]][1]
    m[names(sha[[i]])[2], names(sha[[i]])[1]] <- sha[[i]][2]
    }
  
  return(m)
}




