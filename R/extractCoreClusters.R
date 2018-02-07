
#' @importFrom parallel mclapply
#' @importFrom seqinr read.fasta write.fasta
#' @export
extractCoreClusters <- function(prokka, 
                                coreGenome, 
                                roary, 
                                out.dir = '.', 
                                prefix, 
                                cpus = 1L){
  
  
  
  cl <- class(prokka)
  if (cl != 'prokka' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "prokka" is not of class "prokka" (mirecipe package).')
  }
  
  if(!validObject(prokka)){
    stop('Object "prokka" is not valid.')
  }
  
  cl <- class(coreGenome)
  if (cl != 'coreGenome' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "coreGenome" is not of class "coreGenome" (mirecipe package).')
  }
  
  if(!validObject(coreGenome)){
    stop('Object "coreGenome" is not valid.')
  }
  
  cl <- class(roary)
  if (cl != 'roary' | attr(cl, 'package') != 'mirecipe'){
    stop('Object "roary" is not of class "roary" (mirecipe package).')
  }
  
  if(!validObject(roary)){
    stop('Object "roary" is not valid.')
  }
  
  if( suppressWarnings(is.na(as.integer(cpus))) ){ 
    
    stop('First argument (cpus) is not an integer.\n')
    
  }
  
  #Check if all prefixes are the same. Stop if not
  pfx <- vapply(c(prokka, coreGenome, roary), function(x){
    slot(x, 'prefix')
    }, FUN.VALUE = NA_character_)
  if(! pfx[1]==pfx[2] & pfx[2]==pfx[3]){
    stop("Prefixes doesn't match.")
  }
  
  
  in.files <- lapply(c(prokka, coreGenome, roary), function(x){
    slot(x, 'out.files')
  })
  names(in.files) <- c('prokka','coreGenome','roary')
  # in.files <- normalizePath(in.files)
  
  out.dir <- normalizePath(out.dir)
  
  if (missing(prefix)){
    
    prefix <- pfx[1]
    
  }
  
  dout <- paste0(out.dir, '/', prefix, '_coreClusters')
  
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
  

  
  lcbs <- grep('coreAliCoord.txt', in.files$coreGenome[[1]], value = TRUE)
  gffs <- sapply(in.files$prokka, function(x){grep('gff$', x, value = TRUE)})
  roary_clusters <- grep('clustered_proteins$', in.files$roary[[1]], value = TRUE)
  
  clu <- getCoreClusters(lcbs = lcbs, 
                         gffs = gffs, 
                         roary_clusters = roary_clusters, 
                         cpus = cpus)
  
  ffa <- lapply(in.files$prokka, function(x){
    grep('faa$|ffn$', x, value = TRUE)
    })
  
  
  faas <- unlist(mclapply(sapply(ffa, '[', 1), 
                          read.fasta, 
                          'AA', 
                          TRUE, 
                          mc.cores=cpus), 
                 recursive = FALSE)
  names(faas) <- vapply(faas, 
                        attr, 
                        'name', 
                        FUN.VALUE = NA_character_)
  
  ffns <- unlist(mclapply(sapply(ffa, '[', 2), 
                          read.fasta, 
                          'DNA', 
                          TRUE, 
                          mc.cores=cpus), 
                 recursive = FALSE)
  names(ffns) <- vapply(ffns, 
                        attr, 
                        'name', 
                        FUN.VALUE = NA_character_)
  
  out.files <- mclapply(names(clu), function(X){
    
    fo <- paste0(dout,'/proteins/',X, '.faa')
    fa <- paste0(dout,'/genes/',X, '.ffn')
    
    write.fasta(faas[clu[[X]]], names = clu[[X]], file.out = fo, as.string = TRUE)
    write.fasta(ffns[clu[[X]]], names = clu[[X]], file.out = fa, as.string = TRUE)
    
    c(fo, fa)
  }, mc.cores = cpus)
  
  out.files <- lapply(1:2, function(X){
    vapply(out.files, '[', X, FUN.VALUE = NA_character_)
    })
  names(out.files) <- c('proteins', 'genes')
  
  call <- capture.output(match.call())
  
  .CoreClusters(in.files = in.files,
                out.files = out.files,
                #stats = stats, 
                call = call, 
                prefix = prefix)
  
}














# Support functions
getCoreClusters <- function(lcbs,
                            gffs,
                            roary_clusters,
                            prefix,
                            cpus = 1L){
  
  #Read gffs as list with data.frame and gene sequences
  x <- parallel::mclapply(gffs, readGff, mc.cores = cpus)
  names(x) <- sub('[.]gff$','', basename(gffs))
  
  #Read LCBs
  lcb <- readLCBS(lcbs)
  
  #Read roary's clusters
  clusters <- readRoaryClusters(roary_clusters)
  
  #Get which genes appear in the progressiveMauve core-alignment.
  system.time(coreGenes <- parallel::mclapply(names(x), function(i){
    
    tt <- x[[i]]
    ll <- do.call(rbind, lapply(lcb, function(j){j[i, ]}))
    tt[which(apply(tt, 1, function(j){
      any(apply(ll, 1, function(k){
        overlap(c(j[7], j[8]), c(k[2], k[3]))
      }))
    })), 2]
    
  }))
  names(coreGenes) <- names(x)
  
  #Get which roary's clusters appear in the progressiveMauve core-alignment.
  ul <- unlist(coreGenes)
  ulinx <- parallel::mclapply(clusters, function(y){
    all(y%in%ul)
  }, mc.cores = cpus)
  uliny <- unlist(ulinx)
  
  clu <- clusters[names(which(uliny))]
  clup <- lapply(clu, function(x){sapply(strsplit(x, '_'), '[', 1)})
  
  #Just keep clusters with one gene per genome.
  clu <- clu[sapply(clup, function(i){all(names(x)%in%i)})]
  clu <- clu[which(sapply(clu, length)==length(x))]
  
 clu 
}

readGff <- function(gff){
  rl <- readLines(gff)
  x <- getGffTable(rl)
  x
}

getGffTable <- function(rl){
  
  # rl <- readLines(gff)
  
  w <-which(grepl('^\\#\\#',rl))
  
  re <- rl[w][-c(1, length(w))]
  re <- do.call(rbind,lapply(strsplit(re,' '), '[', c(2:4)))
  re <- as.data.frame(re,stringsAsFactors=FALSE)
  re[,2] <- as.integer(re[,2])
  re[,3] <- as.integer(re[,3])
  
  re$V5 <- cumsum(re$V3)
  re$V4 <- re$V5 - re$V3 + 1L
  re <- re[, c(1,2,3,5,4)]
  
  upto <- rev(w)[1] - 1
  from <- rev(w)[2] + 1
  o <- rl[from:upto]
  
  lst <- strsplit(o,'\t')
  
  contig <- sapply(lst, function(x){ x[1] })
  type <- sapply(lst, function(x){ x[3] })
  from <- as.integer(sapply(lst, function(x){ x[4] }))
  to <- as.integer(sapply(lst, function(x){ x[5] }))
  strand <- sapply(lst, function(x){ x[7] })
  phase <- sapply(lst, function(x){ x[8] })
  attrib <- sapply(lst, function(x){ x[9] })
  
  metadata <- strsplit(attrib,';')
  
  id <- sapply(metadata,function(x){
    gp<-grep('ID=',x,value = T)
    if(length(gp)>0){sub('ID=','',gp)}else{''}
  })
  locustag <- sapply(metadata,function(x){
    gp<-grep('locus_tag=',x,value = T)
    if(length(gp)>0){sub('locus_tag=','',gp)}else{''}
  })
  gene <- sapply(metadata,function(x){
    gp<-grep('gene=',x,value = T)
    if(length(gp)>0){sub('gene=','',gp)}else{''}
  })
  product <- sapply(metadata,function(x){
    gp<-grep('product=',x,value = T)
    if(length(gp)>0){sub('product=','',gp)}else{''}
  })
  
  out <- data.frame(Contig=contig,
                    ID=id,
                    LocusTag=locustag,
                    Gene=gene,
                    Product=product,
                    Type=type,
                    From=from,
                    To=to,
                    Strand=strand,
                    Phase=phase,
                    stringsAsFactors = F)
  
  out$From <- apply(out, 1, function(x){
    re[which(re[, 1]==x[1]), 4] + as.integer(x[7]) - 1L
  })
  
  out$To <- apply(out, 1, function(x){
    re[which(re[ ,1]==x[1]), 4] + as.integer(x[8]) - 1L
  })
  
  out
  
}


getGffSeqs <- function(x, rl){
  gp <- grep('##FASTA',rl, fixed = TRUE) + 1L
  fna <- rl[gp:length(rl)]
  
  
  gph <- grep('^>', fna)
  fna <- paste0(fna[-gph], collapse = '')
  fna <- strsplit(fna, '')[[1]]
  
  sqs <- apply(x, 1, function(y){
    gen <- fna[y[[7]]:y[[8]]]
    if(y[[9]]=='-'){
      gen <- rev(seqinr::comp(gen, forceToLower = FALSE))
    }
    
    se <- paste0(gen, collapse = '')
    seqinr::as.SeqFastadna(se)
  })
  
  sqs <- as.list(sqs)
  names(sqs) <- x$LocusTag
  
  return(sqs)
}

readLCBS <- function(lcbs){
  
  lcb <- readLines(lcbs)
  
  eq <- grep('^=', lcb)
  vp <- vapply(1:length(eq), function(x){
    
    ini <- ifelse(eq[x]==eq[1], 1L, eq[x-1L]+1L)
    end <- eq[x]-1L
    c(ini, end)
    
  }, FUN.VALUE = c(1L,1L))
  
  dfs <- apply(vp, 2, function(x){
    sp <- strsplit(lcb[x[1]:x[2]], ':')
    sp2 <- lapply(sp, function(y){strsplit(y,' ')[[2]]})
    
    df <- vector('list', 5)
    names(df) <- c('fasta', 'start', 'end', 'strand', 'length')
    df$fasta <- sub('.fna', '', basename(sapply(sp2, '[', 3)))
    df$start <- as.integer(sapply(sp2, function(y) {
      strsplit(y[1],'-')[[1]][1]
    } ))
    df$end <- as.integer(sapply(sp2, function(y) {
      strsplit(y[1],'-')[[1]][2]
    } ))
    df$strand <- sapply(sp2, '[', 2)
    df$length <- as.integer(sapply(sp, function(y){
      sub(' nch = ','', y[3])
    }))
    
    df <- as.data.frame(df)
    rownames(df) <- sub('[.]fasta$','', df$fasta)
    
    df
  })
  
  dfs
}


readRoaryClusters <- function(roary_clusters){
  
  cl <- readLines(roary_clusters)
  sp <- strsplit(cl, ': ')
  
  out <- lapply(sp, function(x){
    strsplit(x[2],'\t')[[1]]
  })
  
  names(out) <- sapply(sp, '[', 1)
  
  return(out)
}



overlap <- function(a=c(1L, 1L), b=c(1L, 1L)){
  if(a[1]>a[2]){
    a <- rev(a)
  }
  if(b[1]>b[2]){
    b <- rev(b)
  }
  a[1]<=b[2] & a[2]>=b[1]
}









