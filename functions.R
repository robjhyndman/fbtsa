# Create nice R figures with minimal margins
# in landscape format suitable for slides and papers
savepdf <- function(file, width=16, height=10)
{
  fname <<- paste("figs/",file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54, pointsize=10, bg='white')
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}
# Crop pdf to remove all white space
endpdf <- function()
{
  #dev.off()
  crop::dev.off.crop(fname)
}

# Histograms
gghist <- function(data, mapping, ...)
{
  x <- data[[as.character(mapping$x[2])]]
  bw <- 0.2*bw.nrd0(x) + 0.8*bw.SJ(x)
  p <- ggplot(data, mapping) +
    geom_density(col=NA, fill="#cc5900", bw=bw)
  return(p)
}


# Produce plot of data and quantiles for electricity example
# To save passing large objects as arguments, this treats DT and qdemand
# as global variables
# Arguments:
#    id - the smartmetre id to plot
#    showquantiles - which quantiles to plot. By default, it plots the deciles

qdemandplot <- function(id, showquantiles=seq(0.1,0.9,by=0.1))
{
  library(data.table)
  library(ggplot2)
  idlist <- unique(DT[,id])
  if(id <= 500)
    id <- idlist[id]
  prob <- sort(unique(qdemand[,prob]))

  # Subset of DT
  j <- (DT[,id]==id)
  z <- DT[j, ]
  z[, tod:=z[,period]/2]
  z$dow <- factor(z$dow,levels=1:7,
                  labels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"))

  p1 <- ggplot(aes(y=demand, x=tod), data=z) +
    geom_point(shape=".") + facet_grid(~dow) +
    ylab("Demand (kWh)") + xlab("") +
    ggtitle(paste("Demand for ID:",id)) +
    guides(fill=FALSE) +
    scale_x_continuous(breaks=c(0,6,12,18,24))

  # Subset of qdemand
  j <- (qdemand[,id]==id)
  z <- qdemand[j, ]
  z[, tod:= ((z[,tow]-1) %% 48)/2 +1]
  z[, dow:= trunc((z[,tow]-1)/48) + 1]
  z$dow <- factor(z$dow, levels=1:7,
                  labels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"))

  if(!missing(showquantiles))
  {
    j <- z[,prob] %in% showquantiles
    z <- z[j,]
  }

  p2 <- ggplot(aes(y=demand, x=tod, colour=prob, group=prob), data=z) +
    geom_line() + facet_grid(~dow) +
    xlab("Time of day") + ylab("Quantiles") +
    scale_colour_gradientn(colours = rainbow(8),
                           name="Probability",
                           breaks=seq(0.1,0.9,by=0.2)) +
    theme(legend.position="bottom", legend.direction="horizontal",
          legend.key.width=unit(1,"cm"),
          legend.key.height=unit(.3,"cm"),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    scale_x_continuous(breaks=c(0,6,12,18,24))

  return(gridExtra::grid.arrange(p1,p2,ncol=1, heights=c(0.5,0.55)))
}

# Compute Jensen-Shannon distance
# based on quantiles q and p at probabilities prob
JS <- function(prob,q,p)
{
  # Compute approximate densities
  x <- seq(min(q,p),max(q,p), l=201)
  qpmf <- pmf(x,prob,q)
  ppmf <- pmf(x,prob,p)
  m <- 0.5 * (ppmf + qpmf)
  JS <- suppressWarnings(0.5*(sum(na.omit(ppmf*log(ppmf/m))) +
                              sum(na.omit(qpmf*log(qpmf/m)))))
  return(JS)
}

# Compute approximate discretized density (like a probability mass function)
# at each x (equally spaced) given quantiles q with probabilities p
pmf <- function(x, p, q)
{
  qcdf <- approx(q,p,xout=x,yleft=0,yright=1)$y
  qpmf <- c(0,diff(qcdf)/ (x[2]-x[1]))
  return(qpmf / sum(qpmf))
}

jsd <- function(qdemand)
{
  idlist <- unique(qdemand[,id])
  probs <- sort(unique(qdemand[,prob]))
  nid <- length(idlist)
  dmat <- matrix(0, nrow=nid, ncol=nid)
  rownames(dmat) <- colnames(dmat) <- idlist
  x <- seq(0, max(qdemand[,demand]), l=51)

  for(i in 2:nid)
    for(j in 1:(i-1))
    {
      tmp <- qdemand[id==idlist[i],]
      tmp[, demand2:=qdemand[id==idlist[j],demand]]
      dmat[i,j] <- sum(tmp[, JS(prob,demand,demand2), by=.(tow)]$V1)
      # for(dow in 1:7)
      #   for(period in 1:48)
      #   {
      #     tmp2 <- subset(tmp, period==period, dow==dow)
      #     js <- JS(prob,)
      #     }
    }

  # Create object of class "dist"
  return(as.dist(dmat + t(dmat)))
}


# Compute similarity matrix based on pairwise distances
# 3 different kernels are possible
# The h parameter is selected to be very large if the argument is omitted

similarity <- function(distances, h,
              kernel=c("Gaussian","Epanechnikov","Triangle"))
{
  if(class(distances) != "dist")
    stop("distances should be of class 'dist'")
  kernel <- match.arg(kernel)

  if(missing(h))
    h <- 1000*max(distances)

  distances <- as.matrix(distances)
  if(kernel=="Gaussian")
    sims <- exp(-(distances/h)^2)
  else if(kernel=="Epanechnikov")
    sims <- pmax(1-(distances/h)^2, 0)
  else
    sims <- pmax(1-distances/h, 0)

  return(sims)
}

# Compute similarity matrix based on pairwise distances
# Only implements Epanechnikov kernel
# distances must be a dist class.
# returns sparse object
sparsesimilarity <- function(distances, h)
{
  if(!is.element("dist",class(distances)))
    stop("distances must be of class 'dist'")

  if(missing(h))
    h <- 1000*max(distances)

  n <- attr(distances,"Size")
  k <- distances < h
  col <- rep(1:(n-1),(n-1):1)[k]
  row <- matrix(rep(1:n,n), n,n)[lower.tri(matrix(0,n,n))]
  row <- row[k]
  v <- 1-(distances[k]/h)^2
  sims <- Rcsdp::simple_triplet_sym_matrix(row,col,v,n=n)

  return(sims)
}


# Main function for finding embedding based on distances
# Arguments:
#  distances - an object of class "dist" (essential the lower triangle of distances matrix)
#  m         - embedding dimension
#  method    - which dimension reduction method to use. Default is LaplacianMDS which
#              uses Laplacian but with very large h. This is equivalent to "MDSiso" but much
#              faster.
#  h         - bandwidth for computing the similarity matrix. Only used for Laplacian
#              methods, apart from LaplacianMDS where hs is set to a large h.
#  ...    - any other arguments are passed to the function implementing the embedding method.

embedding <- function(distances, m=2,
  method=c("LaplacianMDS","Laplacian","Lrw","Lsym","Lsym2",
            "MDS","MDSiso","monoMDS","DPM","Rtsne"),
  h = median(distances), ...)
{
  method <- match.arg(method)
  if(class(distances)!="dist")
    stop("distances should be of class dist")
  if(m>6)
    stop("Maximum embedding is 6 dimensions")

  if(method == "LaplacianMDS")
  {
    #Laplacian eigenmap with large h
    w <- similarity(distances, ...)
    n <- NROW(w)
    D <- diag(rowSums(w))
    ei <- geigen::geigen(D-w,D,symmetric=TRUE)$vectors[,2:(m+1),drop=FALSE]
  }
  else if(method == "Laplacian")
  {
    #Laplacian eigenmap with regular h
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    D <- diag(rowSums(w))
    ei <- geigen::geigen(D-w,D,symmetric=TRUE)
    ei <- ei$vectors[,2:(m+1),drop=FALSE]
  }
  else if(method=="Lrw")
  {
    #Laplacian eigenmap. normalized Lrw
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    D <- diag(1/rowSums(w))
    ei <- geigen::geigen(diag(n) - D %*% w,D,symmetric=TRUE)
    ei <- ei$vectors[,2:(m+1),drop=FALSE]
  }
  else if(method=="Lsym")
  {
    #Laplacian eigenmap. normalized Lsym
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    wden <- rowSums(w)
    D <- diag(1/wden)
    Dhalf <- diag(1/sqrt(wden))
    ei <- geigen::geigen(diag(n) - Dhalf %*% w %*% Dhalf,D,symmetric=TRUE)
    ei <- ei$vectors[,2:(m+1),drop=FALSE]
  }
  else if(method=="Lsym2")
  {
    #Laplacian eigenmap. normalized Lsym
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    wden <- rowSums(w)
    D <- diag(1/wden)
    Dhalf <- diag(1/sqrt(wden))
    ei <- eigen(Dhalf %*% w %*% Dhalf,symmetric=TRUE)
    ei <- ei$vectors[,1:m,drop=FALSE]
  }
  else if(method=="Rtsne")
  {
    ei <- Rtsne::Rtsne(distances, dims=m, perplexity=9)$Y
  }
  else if(method=="MDS")
    ei <- mds(distances, ndim=m)$conf
  else if(method=="MDSiso")
  {
    # Multidimensional scaling
    mds <- MASS::isoMDS(distances, k=m)
    ei <- mds$points
  }
  else if(method=="monoMDS")
  {
    # Multidimensional scaling
    mds <- vegan::monoMDS(distances, k=m, model="local")
    ei <- mds$points
  }
  else if(method=="DPM")
  {
    # Density preserving map
    ei <- dpm(distances, m=m, ...)
    colnames(ei) <- paste("Comp",1:m, sep="")
    rownames(ei) <- attr(distances, "Labels")
    return(structure(scale(ei),class="embedding"))
  }

  else
    stop("Not implemented")

  colnames(ei) <- paste("Comp",1:m, sep="")
  rownames(ei) <- attr(distances,"Labels")

  # Scale embedding
  ei <- scale(ei)
  # Then take signs so medians are positive
  # Only purpose of this is to avoid arbitrary changes in sign for a component
  med_ei <- apply(ei, 2, median)
  j <- med_ei < 0
  ei[,j] <- -ei[,j]

  return(structure(list(y=ei,method=method,distances=distances),class="embedding"))
}

# Find outliers in matrix x
# embedded indicates if we look for outliers in embedded space
# or original space

outliers <- function(x,
  embedded=FALSE,
  method=c("kde","HDoutliers"),
  noutliers=NULL,
  pcoutliers=1,
  bandwidth=1e6)
{
  method <- match.arg(method)
  if(class(x)!="embedding")
    stop("This function is for objects of class 'embedding'")
  if(embedded)
  {
    # Extract embedded points
    x <- x$y
    m <- NCOL(x)
    # Check inputs
    if(m>3 & method=='kde')
    {
      warning("kde can't find outliers in more than 3d space. Switching to HDoutliers")
      method <- "HDoutliers"
    }
  }
  else if(method=="HDoutliers")
  {
    warning("HDoutliers only works on embedded space. Switching to kde")
    method <- "kde"
  }


  if(method=="kde")
  {
    # Work in original space
    if(!embedded)
      fxy <- kdedist(x$distances, bandwidth)
    else
      fxy <- kdeobs(x)
    return(kdeoutliers(fxy, noutliers=noutliers, pcoutliers=pcoutliers))
  }
  else
  {
    outliers <- HDoutliers::HDoutliers(x)
    names(outliers) <- rownames(x)[outliers]
    return(outliers)
  }
}

kdeoutliers <- function(fxy,
                        noutliers=NULL,
                        pcoutliers=1)
{
  if(is.null(noutliers))
  {
    ql <- quantile(fxy, prob=pcoutliers/100)
    noutliers <- sum(fxy < ql)
  }
  if(noutliers > 0)
  {
    outliers <- order(fxy)[seq(noutliers)]
    names(outliers) <- names(fxy)[outliers]
    return(outliers)
  }
  else
    return(NULL)
}


# Scatterplots of Laplacian eigenmaps
# If embedded=TRUE, then show HDRs and outliers from embedded space
# Else show from original space

plot.embedding <- function(embed,
  embedded=TRUE,
  m=NCOL(embed$y),
  noutliers=NULL,
  pcoutliers=1,
  outliermethod=c("kde","HDoutliers"),
  levels=c(1,50,99),
  showhdr=TRUE,
  kde.package=c("ash","ks"),
  main=paste("Embedding:",embed$method),
  bandwidth=1e6,
  labels=c("metres","rank"),
  ...)
{
  outliermethod <- match.arg(outliermethod)
  kde.package <- match.arg(kde.package)
  labels <- match.arg(labels)

  if(!embedded & outliermethod=="HDoutliers")
  {
    outliermethod <- "kde"
    warning("Using kde outlier method on original space")
  }
  data <- embed$y[,1:m]

  region <- NULL
  if(!showhdr)
    levels <- NULL

  # Find outliers and HDRs
  if(!is.null(levels) | outliermethod=="kde")
  {
    m <- NCOL(data)
    if(m > 2)
      kde.package <- "ks"
    if(embedded)
      fxy <- kdeobs(data, use_ash=(kde.package=='ash'))
    else
      fxy <- kdedist(embed$distances, bandwidth)
    if(!is.null(levels))
    {
      levels <- sort(levels)
      if(max(levels) < 1)
        levels <- levels*100
      ql <- quantile(fxy, prob=1-levels/100)
      region <- numeric(NROW(data)) + 100
      for(i in rev(seq_along(levels)))
        region[fxy > ql[i]] <- levels[i]
    }
    if(outliermethod=="kde")
      outliers <- kdeoutliers(fxy, noutliers,pcoutliers)
  }
  if(outliermethod=="HDoutliers")
  {
    outliers <- HDoutliers::HDoutliers(data)
    names(outliers) <- rownames(x)[outliers]
  }

  data <- as.data.frame(data)
  varnames <- colnames(data)

  if(length(outliers) > 0)
  {
    if(labels=="metres")
      labs <- rownames(data)[outliers]
    else
      labs <- 1:length(outliers)
  }

  if(m==1)
  {
    p <- ggplot2::ggplot(data,ggplot2::aes_string(varnames[1])) +
      ggplot2::geom_density(bw="SJ", fill="salmon", col=FALSE) +
      ggplot2::geom_rug()
    if(!is.null(outliers))
    {
      p <- p + ggplot2::annotate("text", x = data[[1]][outliers],
                                 y=rep(-max(fxy)/50,length(outliers)),
                                 label=labs, col='blue', cex=2.5)
    }
    if(!is.null(main))
      p <- p + ggplot2::ggtitle(main)
  }
  else if(m==2)
  {
    p <- annotatedplot(data, ggplot2::aes_string(x=varnames[1],y=varnames[2]),
                       outliers=outliers, labels=labs, region=region)
    if(!is.null(main))
      p <- p + ggplot2::ggtitle(main)
  }
  else
    p <- GGally::ggpairs(data,
                         title=main,
                         lower=list(continuous=GGally::wrap(annotatedplot,
                                                            outliers=outliers,
                                                            region=region, labels=labs, textsize=2.5)),
                         diag=list(continuous=mydensitydiag),
                         upper=list(continuous=GGally::wrap(annotatedplot,
                                                            outliers=outliers,
                                                            region=region, labels=labs, textsize=2.5))) +
        ggplot2::theme(text = ggplot2::element_text(size=10))
  return(p)
}


annotatedplot <- function(data, mapping, outliers=NULL, labels=NULL, region=NULL, textsize=2.5,...)
{
  xvar <- as.character(mapping$x)[2]
  yvar <- as.character(mapping$y)[2]
  xlim <- diff(range(data[[xvar]]))
  ylim <- diff(range(data[[yvar]]))

  if(!is.null(region))
  {
    # Construct region factor
    levels <- sort(unique(region[region < 100]), decreasing=TRUE)
    levels <- c(levels, 100)
    data$Region <- factor(region, levels=levels,
                          labels=c(paste(head(levels,-1)), ">99"))

    # Sort data so the larger regions go first (other than outliers)
    k <- region
    k[region==100] <- 0
    ord <- order(k, decreasing=TRUE)

    p <- ggplot2::ggplot(data[ord,], mapping) +
      ggplot2::geom_point(ggplot2::aes(col=data$Region[ord]))

    p <- p + ggplot2::scale_colour_manual(
      name="HDRs",
      breaks=c(paste(head(sort(levels),-1)), ">99"),
      values=c(RColorBrewer::brewer.pal(length(levels),"YlOrRd")[-1],"#000000"))
  }
  else
    p <- ggplot2::ggplot(data, mapping) + ggplot2::geom_point()
  if(!is.null(outliers))
  {
    if(is.null(labels))
      labels <- rownames(data)[outliers]
    p <- p + ggplot2::annotate("text", x = data[[xvar]][outliers]+xlim/50, y=data[[yvar]][outliers]+ylim/50,
                               label=labels, col='blue', cex=textsize)
  }
  return(p)
}


mydensitydiag <- function (data, mapping, ..., rescale = FALSE)
{
  p <- ggplot2::ggplot(data, mapping) + ggplot2::scale_y_continuous()
  if (identical(rescale, TRUE)) {
    p <- p + ggplot2::stat_density(aes(y = ..scaled.. * diff(range(x, na.rm = TRUE)) +
                                        min(x, na.rm = TRUE)),
                position = "identity", geom = "line", bw="SJ", , fill="salmon", col=FALSE, ...)
  }
  else {
    p <- p + ggplot2::geom_density(..., bw="SJ", fill="salmon", col=FALSE)
  }
  p$type <- "diag"
  p$subType <- "density"
  p
}

# Compute row means of sparse symmetric matrix
rowMeansSparse <- function(x)
{
  result <- numeric(x$n)
  for(i in seq(x$n))
  {
    k <- (x$i==i | x$j==i)
    result[i] <- sum(x$v[k])
  }
  return(result/x$n)
}

# Return number of non-zeros in each row of sparse symmetric matrix
rowNonZero <- function(x)
{
  result <- numeric(x$n)
  for(i in seq(x$n))
  {
    k <- (x$i==i | x$j==i)
    result[i] <- sum(k)
  }
  return(result)
}

# Add two sparse symmetric matrices of equal size
addSparse <- function(x,y)
{
  if(x$n != y$n)
    stop("Matrices not the same size")

  # Combine rows, colums and non-zero elements
  i <- c(x$i,y$i)
  j <- c(x$j,y$j)
  v <- c(x$v,y$v)
  # Find duplicates
  z <- duplicated(cbind(i,j))
  if(any(z))
  {
    #Split duplicates into separate vectors
    i2 <- i[z]
    j2 <- j[z]
    v2 <- v[z]
    i <- i[!z]
    j <- j[!z]
    v <- v[!z]
    # Add together any duplicate values
    for(k in seq_along(i2))
    {
      l <- which(i==i2[k] & j==j2[k])
      v[l] <- v[l] + v2[k]
    }
  }
  return(Rcsdp::simple_triplet_sym_matrix(i,j,v,n=x$n))
}


dij <- function(i,j,n)
{
  Rcsdp::simple_triplet_sym_matrix(
    i=c(i,j,i),
    j=c(i,j,j),
    v=c(1,1,-1),n=n)
}

dpm <- function(d, h, m)
{
  n <- attr(d,"Size")
  w <- sparsesimilarity(d, h=h)
  f <- rowMeansSparse(w) - rowNonZero(w)/n
  b <- c(f[f<0],0)

  A <- list()
  nA <- 0
  for(i in seq(n))
  {
    k <- (w$i==i | w$j==i)
    if(any(k))
    {
      idx <- which(k)
      i0 <- w$i[idx]
      j0 <- w$j[idx]
      tmp <- Rcsdp::.simple_triplet_zero_sym_matrix(w$n)
      if(length(idx)>0)
      {
        for(j in seq_along(idx))
          tmp <- addSparse(tmp, dij(i0[j],j0[j],n))
      }
      tmp$v <- -tmp$v/h^2/n
      nA <- nA+1
      A[[nA]] <- list(tmp)
    }
  }
  A[[nA+1]] <- list(matrix(1,n,n))
  C <- list(Rcsdp::.simple_triplet_diag_sym_matrix(1,n))
  tmp <- Rcsdp::csdp(C, A, b, list(type="s", size=n))

  if(tmp$status!=0)
    warning("Not converged")
  return(tmp$X[[1]][,1:m,drop=FALSE])
}


# Compute kernel density estimate at observations
kdeobs <- function(x, h=NULL, use_ash=TRUE)
{
  m <- NCOL(x)
  if(m==1)
  {
    if(is.null(h))
      h <- bw.SJ(x)
    den <- density(x, bw=h)
    fxy <- approx(den$x,den$y,xout=x)$y
  }
  else if(m==2)
  {
    kde.package <- ifelse(use_ash, "ash", "ks")
    den <- hdrcde::hdr.2d(x[,1], x[,2], prob=0.5, kde.package=kde.package, h=h)
    fxy <- den$fxy
  }
  else
  {
    if(is.null(h))
      h <- ks::Hpi.diag(x, binned = TRUE, nstage=1, optim.fun="optim")
    fxy <- ks::kde(x, eval.points=x,H=h)$estimate
  }

  names(fxy) <- rownames(x)
  return(fxy)
}

# Compute kernel density estimate from pairwise distances on original space
kdedist <- function(d, bandwidth)
{
  d <- as.matrix(d)
  fxy <- rowSums(exp(-d/bandwidth))
  return(fxy)
}
