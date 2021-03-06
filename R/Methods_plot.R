# setGeneric("plot")

# plot.foo <- function(object, ...) {
   ## implement plot.foo
#   "I'm a foo"
# }
# setMethod("plot", "foo", plot.foo)

if (!isGeneric('plot'))
  setGeneric('plot', function(x, y, ...)
    standardGeneric('plot'))




#' @export
plot.roary <- function(x,
                       hclustfun = hclust,
                       distfun = dist,
                       edgePar = list(lwd = 1.5),
                       margins = c(3, 0, 4, 5),
                       col = c('white','#E64B35FF'),
                       raster = TRUE,
                       labRow = NULL,
                       cexRow = 0.5,
                       ...){
  
  # cl <- class(x)
  # if (cl != 'roary' | attr(cl, 'package') != 'mirecipe'){
  #   stop('Object "roary" is not of class "roary" (mirecipe package).')
  # }
  # 
  # cl <- class(x)
  # if (cl != 'roary'){
  #   stop('Object "roary" is not of class "roary" (mirecipe package).')
  # }
  
  
  y <- slot(roary, 'panmatrix')
  
  y[which(y>1L)] <- 1L
  
  ### Distance and clustering ###
  dd <- distfun(y)
  hc <- hclustfun(dd)
  
  mm <- y[, order(colSums(y), decreasing = T)]
  mm <- as.matrix(mm)
  mm <- mm[hc$order, ]
  
  dhc <- as.dendrogram(hc)
  
  nr <- nrow(mm)
  nc <- ncol(mm)
  
  mat <- rbind(c(2,1))

  ### Start device functions ###
  op <- par(no.readonly = T)
  on.exit(dev.flush())
  on.exit(par(op), add = TRUE)
  
  dev.hold()
  
  ### Set layout ###
  
  
  layout(mat = mat,
         widths = c(1, 3),
         heights = c(1, 1),
         respect = F)
  
  ### Image ###
  par(mar = c(margins[1],
              margins[2],
              margins[3],
              margins[4]))
  
  image(x = 1:nc,
        y = 1:nr,
        z = t(mm),
        xlim = 0.5 + c(0, nc),
        ylim = 0.5 + c(0, nr),
        frame.plot=F,
        col = col,
        yaxt = 'n',
        xaxs= 'r',
        ylab = '',
        xlab = '',
        useRaster = raster,
        xpd = T)
  
  if(is.null(labRow)){
    labRow <- rownames(mm)
  }
  
  axis(4,
       1:nr,
       labels = labRow,
       las = 2,
       line = -1,
       tick = 0,
       cex.axis = cexRow)

  
  ### Dendrogram ###
  par(mar = c(margins[1], 0, margins[3]+0.5, 0), yaxs= 'i')
  plot(dhc,
       horiz = TRUE,
       # frame.plot = F,
       leaflab = 'n',
       edgePar = edgePar,
       axes=F)
  
  
}
setMethod("plot", signature(x = "roary", y = "missing"), 
          function(x, y, ...) plot.roary(x, ...))



#' @importFrom gplots heatmap.2
plot.progressiveMauve <- function(x, ...){
  
  ss <- slot(x, 'stats')
  heatmap.2(ss, trace = 'n', ...)
  
}
setMethod("plot", signature(x = "roary", y = "missing"), 
          function(x, y, ...) plot.roary(x, ...))



