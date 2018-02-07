
#' @importFrom parallel mclapply
#' @export
alignCoreClusters <- function(coreClusters, aligner = 'mafft', prefix, cpus = 1L){
  
  
  aligner <- match.arg(aligner, c('mafft'))
  if (aligner=='mafft' & Sys.which('mafft')==''){
    stop('mafft is not in your $PATH.')
  }
  
  
  cl <- class(coreClusters)
  if (cl != 'coreClusters' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "coreClusters" is not of class "coreClusters" (mirecipe package).')
  }
  
  if(!validObject(coreClusters)){
    stop('Object "coreClusters" is not valid.')
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  in.files <- slot(coreClusters, 'out.files')
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- slot(coreClusters, 'prefix')
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_alignedCoreClusters')
  
  if(!dir.exists(dout)){
    dir.create(dout)
    dir.create(paste0(dout, '/proteins'))
    dir.create(paste0(dout, '/genes'))
  }else{
    unlink(dout, recursive = TRUE)
    dir.create(dout)
    dir.create(paste0(dout, '/proteins'))
    dir.create(paste0(dout, '/genes'))
  }
  
  dirout <- paste0(dout, c('/proteins', '/genes'))
  out.files <- lapply(1:2, function(x){
    
    unlist(mclapply(in.files[[x]], function(y){
      align(y, aligner=aligner, extension = '.ali', dirout = dirout[x])
    }, mc.cores = cpus, mc.preschedule = FALSE))
    
  })
  names(out.files) <- c('proteins', 'genes')
  
  call <- capture.output(match.call())
  
  .AlignedCoreClusters(in.files = in.files, 
                       out.files = out.files, 
                       # stats = stats, 
                       call = call, 
                       prefix = prefix)
  
}


align <- function(fasta, aligner = 'mafft', extension = '.ali', dirout){
  
  
  outfile <- paste0(basename(fasta), extension)
  outfile <- paste0(dirout, '/', outfile)
  
  if(aligner=='mafft'){
    al <- mafft(fasta, outfile)
  }else{
    stop('Not implemented.')
  }
  
  al
}

mafft <- function(fasta, outfile){
  run <- paste0('ginsi --quiet ', fasta, ' > ', outfile)
  system(run)
  outfile
}

