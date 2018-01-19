
#' @export
gubbins <- function(coreGenome, 
                    out.dir = '.',
                    prefix,
                    cpus = 1L){
  
  if (Sys.which('run_gubbins.py')==""){
    
    stop('run_gubbins.py is not in your $PATH.')
    
  }
  
  cl <- class(coreGenome)
  if (cl != 'coreGenome' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "coreGenome" is not of class "coreGenome" (mirecipe package).')
  }
  
  if(!validObject(coreGenome)){
    stop('Object "coreGenome" is not valid.')
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  in.files <- vapply(slot(coreGenome, 'out.files'), function(x){
    grep('[.]fasta$', x, value = TRUE)
  }, FUN.VALUE = '')
  
  in.files <- normalizePath(in.files)
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- slot(coreGenome, 'prefix')
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_gubbins')
  
  if(!dir.exists(dout)){
    dir.create(dout)
  }else{
    unlink(dout, recursive = TRUE)
    dir.create(dout)
  }
  
  
  # px <- paste0(dout, '/', prefix)
  
  #Workaround for putting out files in specified folder..
  wd <- getwd()
  on.exit(setwd(wd))
  
  setwd(dout)
  
  run <- paste0('run_gubbins.py -t fastree ', 
                '--threads ', cpus, ' ',
                '--prefix ', prefix, ' ',
                in.files)
  
  system(run)
  
  out.files <- vector('list', 1L)
  
  out.files[[1L]] <- normalizePath(list.files(path = dout, 
                                              full.names = TRUE))
  
  
  
  
  
  .Gubbins(in.files = in.files, 
           out.files = out.files, 
           # stats = stats, 
           call = run, 
           prefix = prefix)
  
  
  
}