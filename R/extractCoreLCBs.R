

#' @export
extractCoreLCBs <- function(progressiveMauve,
                            out.dir = '.',
                            prefix,
                            msi=500L,
                            pco=0.1, #Less is less gaps
                            pal=0.9, #More is less gaps
                            nco=length(progressiveMauve@in.files)){
  
  
  cl <- class(progressiveMauve)
  if (cl != 'progressiveMauve' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "progressiveMauve" is not of class "progressiveMauve" (mirecipe package).')
  }
  # 
  # if( suppressWarnings(is.na(as.integer(cpus))) ){ 
  #   
  #   stop('First argument (cpus) is not an integer.\n')
  #   
  # }
  
  in.files <- vapply(slot(progressiveMauve, 'out.files'), function(x){
    grep('[.]xmfa$', x, value = TRUE)
  }, FUN.VALUE = '')
  
  if (!all(file.exists(in.files))){ 
    
    stop("cannot open the connection. One or more files doesn't exist.")
    
  }
  
  in.files <- normalizePath(in.files)
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- slot(progressiveMauve, 'prefix')
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_coreFasta')
  
  if(!dir.exists(dout)){
    dir.create(dout)
  }
  
  # 
  # oal <- paste0(dout, '/alignment.xmfa')
  # 
  # 
  # xmfa <- grep('xmfa$', slot(progressiveMauve, 'out.files')[[1]], value = TRUE)
  # 
  ncom <- (nco * 2L) + 2L
  rl <- readLines(in.files)
  rl <- rl[-(1L:ncom)]
  
  #index of chunks (blocks) of sequences in xmfa
  eq <- grep('^=', rl)
  vp <- vapply(1:length(eq), function(x){
    
    ini <- ifelse(eq[x]==eq[1], 1L, eq[x-1L]+1L)
    end <- eq[x]-1L
    c(ini, end)
    
  }, FUN.VALUE = c(1L,1L))
  
  # Stats about chunks
  ap <- t(apply(vp, 2, function(x){
    
    rr <- rl[x[1]:x[2]]
    gp <- grep('^>', rr)
    ln <- length(gp)
    nch <- nchar(paste0(rr[2L:(ifelse(ln>1L, gp[2L]-1, length(rr)-1L))],
                        collapse = ''))
    
    vs <- vapply(1:length(gp), function(y){
      
      ini <- gp[y] + 1L
      end <- ifelse(gp[y]!=rev(gp)[1], gp[y+1L] - 1L, length(rr))
      c(ini, end)
      
    }, FUN.VALUE = c(1L,1L))
    
    cb <- strsplit(apply(vs, 2, function(y){
      paste0(rr[y[1]:y[2]], collapse = '')
    }), '')
    
    cb <- do.call(rbind, cb)
    
    #Number of columns with a proportion of gaps lesser than pco.
    prco <- length(which(apply(cb, 2, function(z){
      #proportion of gaps in each column.
      length(which(z=='-'))/ln
    })<=pco))
    
    #Proportion of columns with less than pco gaps.
    pral <- prco/nch
    
    c(ln, nch, pral)
    
  }))
  
  #Save blocks which have a proportion of columns with less than pco gaps more
  #than pal.
  wh <- which(ap[,1]==nco & ap[,2]>=msi & ap[,3]>=pal)
  
  if (length(wh)==0L){
    stop('No LCBs pass filters.')
  }
  
  #Extract selected blocks for each genome.
  ck <- t(vp[, wh])
  ff <- sapply(1:nrow(ck), function(x){
    
    rr <- rl[ck[x,1]:ck[x,2]]
    gp <- grep('^>', rr)
    
    # chu <- sub('[.]xmfa$','.lcbs',xmfa)
    chu <- paste0(dout, '/coreAliCoord.txt')
    
    nch <- ap[wh[x], 2]
    fi <- paste(rr[gp], ': nch =',nch)
    cat(fi, sep = '\n', file = chu, append = TRUE)
    cat('=',sep = '\n', file = chu, append = TRUE)
    
    idx <- vapply(strsplit(sub('> ','',rr[gp]),':'), '[', 1, FUN.VALUE = '')
    
    ini <- gp+1L
    end <- vapply(1:length(gp), function(y){
      ifelse(gp[y]!=rev(gp)[1], gp[y+1L]-1L, length(rr))
    }, FUN.VALUE = 1L)
    
    m <- cbind(ini, end)
    fnam <- paste0(dout, '/', idx, '.core.fasta')
    
    fis <- vapply(1:nrow(m), function(y){
      cat(rr[m[y,1]:m[y,2]], file = fnam[y], sep = '', append = TRUE)
      fnam[y]
    }, FUN.VALUE = '')
    
    fis
    
  })
  
  tmps <- ff[,1]
  rm(ff)
  
  #Create outfile
  out <- paste0(dout, '/coreGenomeAlignment.fasta')
  file.create(out)
  
  #Concatenate blocks
  
  spl <- strsplit(progressiveMauve@in.files, '/')
  he <- sub('[.]\\w+$','',sapply(spl, function(y){rev(y)[1]}))
  vapply(1:length(tmps), function(y){
    cat(paste0('>',he[y], '\n'), file = out, append = TRUE)
    file.append(out, tmps[y])
    cat('\n', file = out, append = TRUE)
    TRUE
  }, FUN.VALUE = TRUE)
  
  file.remove(tmps)
  
  out.files <- vector('list', 1L)
  
  out.files[[1L]] <- normalizePath(list.files(path = dout, 
                                              full.names = TRUE))
  
  call <- capture.output(match.call())
  
  
  new('coreGenome',
      in.files = in.files, 
      out.files = out.files, 
      stats = stats,
      call = call, 
      prefix = prefix)
}








