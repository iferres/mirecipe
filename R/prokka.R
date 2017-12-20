setClass('prokka', 
         slots = list(in.files = 'character',
                      out.files = 'list',
                      stats = 'data.frame',
                      call = 'character', 
                      prefix = 'character'))





#' @export
prokka <- function(in.files, prefix, cpus = 1L){
  
  if (Sys.which('prokka')==""){
    
    stop('prokka is not in your $PATH.')
    
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  if (!all(file.exists(in.files))){ 
    
    stop("cannot open the connection. One or more files doesn't exist.")
    
  }
  
  if (missing(prefix)){
    
    prefix <- format(Sys.time(), "%b%d%H%M%S")
  
  }
  
  dout <- paste0('prokka_', prefix, '/')
  
  if(!dir.exists(dout)){
    dir.create(dout)
  }
  
  ln <- length(in.files)
  
  out.files <- vector('list', ln)
  names(out.files) <- in.files
  
  run <- vector('character', ln)
  names(run) <- in.files
  
  stats <- vector('list', ln)
  names(stats) <- in.files
  
  pb <- txtProgressBar(min = 0, max = length(in.files), style = 3)
  for (i in seq_along(in.files)){
    
    outdir <- paste0(dout, sub('[.]\\w+$', '', in.files[i]))
    
    px <- rev(strsplit(normalizePath(in.files[i]), '/')[[1]])[1]
    
    px <- sub('[.]\\w+$', '', px)
    
    run[i] <- paste0('prokka --quiet --outdir ',outdir, 
                     ' --prefix ',px,
                     ' --locustag ', px,
                     ' --cpus ',cpus,' ',
                     in.files[i])
    
    system(run[i])
    
    out.files[[i]] <- list.files(path = outdir, 
                                 full.names = TRUE)
    
    #Stats
    gp <- grep('[.]txt$', out.files[[i]], value = TRUE)
    rl <- readLines(gp)[-1]
    sp <- sapply(rl, strsplit, ': ')
    stats[[i]] <- as.integer(sapply(sp, '[', 2L))
    names(stats[[i]]) <- sapply(sp, '[', 1L)

    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  # Stats
  stats <- as.data.frame(do.call(rbind, stats))
  
  res <- new('prokka', 
             in.files = in.files,
             out.files = out.files,
             stats = stats,
             call = run, 
             prefix = prefix)
  
  return(res)
  
}