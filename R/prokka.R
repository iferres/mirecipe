

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
prokka <- function(in.files, 
                   out.dir = '.', 
                   prefix, 
                   cpus = 1L){
  
  if (Sys.which('prokka')==""){
    
    stop('prokka is not in your $PATH.')
    
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  if (!all(file.exists(in.files))){ 
    
    stop("cannot open the connection. One or more files doesn't exist.")
    
  }
  
  in.files <- normalizePath(in.files)
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- format(Sys.time(), "%b%d%H%M%S")
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_prokka/')
  
  if(!dir.exists(dout)){
    dir.create(dout)
  }
  
  ln <- length(in.files)
  
  out.files <- vector('list', ln)
  names(out.files) <- in.files
  
  run <- vector('character', ln)
  names(run) <- in.files
  
  stats <- vector('list', ln)
  names(stats) <- basename(in.files)
  vars <- c("contigs", "bases", "CDS", "tRNA", "tmRNA", "repeat_region")
  
  
  pb <- txtProgressBar(min = 0, max = length(in.files), style = 3)
  for (i in seq_along(in.files)){
    
    fn <- rev(strsplit(in.files[i], '/')[[1]])[1]
    
    px <- sub('[.]\\w+$', '', fn)
    
    outdir <- paste0(dout, '/', px)
    
    run[i] <- paste0('prokka --quiet --outdir ',outdir, 
                     ' --prefix ',px,
                     ' --locustag ', px,
                     ' --cpus ',cpus,' ',
                     in.files[i])
    
    system(run[i])
    
    out.files[[i]] <- list.files(path = outdir, 
                                 full.names = TRUE)
    
    out.files[[i]] <- normalizePath(out.files[[i]])
    
    #Stats
    gp <- grep('[.]txt$', out.files[[i]], value = TRUE)
    rl <- readLines(gp)[-1]
    sp <- strsplit(rl, ': ')
    stats[[i]] <- sapply(vars, function(x){
      gp <- grep(x, sp, fixed = T)
      ifelse(length(gp)!=0, as.integer(sp[[gp]][2]), NA_integer_)
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  # Stats
  stats <- as.data.frame(do.call(rbind, stats))
  
  
  #Return
  .Prokka(in.files = in.files,
          out.files = out.files,
          stats = stats,
          call = run, 
          prefix = prefix)
  
}