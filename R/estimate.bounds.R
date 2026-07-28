# Helper to construct quantile-rule based continuous empirical CDFs (types 1-9)
.make_quantile_cdf <- function(x, qtype = 5) {
  x <- sort(as.numeric(x))
  n <- length(x)
  if (n == 0L) return(function(t) numeric(length(t)))
  p <- switch(as.character(qtype),
    "1" = (seq_len(n)) / n,
    "2" = (seq_len(n)) / n,
    "3" = (seq_len(n) - 0.5) / n,
    "4" = (seq_len(n)) / n,
    "5" = (seq_len(n) - 0.5) / n,
    "6" = (seq_len(n)) / (n + 1),
    "7" = if (n > 1L) (seq_len(n) - 1) / (n - 1) else 1,
    "8" = (seq_len(n) - 1/3) / (n + 1/3),
    "9" = (seq_len(n) - 3/8) / (n + 0.25),
    (seq_len(n) - 0.5) / n
  )
  ux <- unique(x)
  # Ties must be placed at the *midpoint* of the vertical jump they produce,
  # (s_{i-1} + s_i) / (2n), which is Appendix A of Ulrich, Miller & Schroeter
  # (2007) and is the mean -- not the maximum -- of the tied observations'
  # plotting positions.  Taking the maximum biased the estimated CDF upward by
  # (n_i - 1) / (2n) at every repeated value, which matters because RT is
  # recorded to the nearest millisecond.  Only reachable when `qtype` is
  # supplied; the default NULL path uses ecdf() and is unaffected.
  up <- vapply(ux, function(z) mean(p[x == z]), numeric(1))
  function(t) {
    out <- numeric(length(t))
    lo <- t < ux[1L]
    hi <- t > ux[length(ux)]
    mid <- !(lo | hi)
    out[hi] <- 1
    if (length(ux) == 1L) {
      out[mid] <- up
    } else if (any(mid)) {
      out[mid] <- stats::approx(ux, up, xout = t[mid], rule = 2)$y
    }
    out
  }
}


# Exact Gaussian kernel distribution estimator:
#   F_h(t) = n^-1 sum_i Phi((t - x_i) / h).
# Log-transformed Gaussian kernel distribution estimator:
#   F_h(t) = n^-1 sum_i Phi((log(t) - log(x_i)) / h) for t > 0, 0 for t <= 0.
# Bandwidth selector `density()` operates on log(x) to prevent boundary leakage.
.make_kernel_cdf <- function(x, bw = "nrd0") {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2L) {
    stop("Kernel CDF estimation needs at least two positive finite response times.",
         call. = FALSE)
  }
  tx <- log(x)
  h <- if (is.character(bw)) {
    if (length(bw) != 1L) stop("kernel_bw must have length one.", call. = FALSE)
    stats::density(tx, bw = bw, n = 2L)$bw
  } else {
    as.numeric(bw)
  }
  if (length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("kernel_bw must select or supply one positive finite bandwidth.",
         call. = FALSE)
  }
  out <- function(t) {
    t <- as.numeric(t)
    res <- numeric(length(t))
    pos <- t > 0
    if (any(pos)) {
      lt <- log(t[pos])
      res[pos] <- vapply(lt, function(z) mean(stats::pnorm((z - tx) / h)),
                         numeric(1))
    }
    res
  }
  attr(out, "bw") <- h
  attr(out, "range") <- range(x)
  out
}


estimate.bounds <- function (RT, CR = NULL, stopping.rule = c("OR","AND","STST"), assume.ID=FALSE,
                             numchannels=NULL, unified.space=FALSE, qtype = NULL,
                             cdf_method = NULL, kernel_bw = "nrd0")
{
  rule <- match.arg(stopping.rule, c("OR","AND","STST"))
  if (is.null(cdf_method)) {
    cdf_method <- if (is.null(qtype)) "ecdf" else "polygon"
  } else {
    cdf_method <- match.arg(cdf_method, c("ecdf", "polygon", "kernel"))
  }
  if (cdf_method == "polygon" && is.null(qtype)) qtype <- 5L
  
  nconds <- length(RT) #total number of input RTs
  n.channels <- numchannels
  #if numchannels for ID model not specified, what is number of channels in the UCIP model?
  if (is.null(numchannels) & !assume.ID) {
    n.channels <- nconds
  }
  
  #if n.channels <2, produce an error message that bounds cannot be 
  #estimated for less than 2 channels
  if (n.channels < 2){
    stop("Too few channels specified. Number of channels must be >=2 to estimate CDF bounds.")
  }
  #create correct/incorrect data if not provided
  if (is.null(CR) | length(CR) != nconds) {
    for (i in 1:nconds) {
      CR[[i]] <- rep(1, length(RT[[i]]))
    }
  }
  
  # Establish CDFs and survivor functions for each single-target distribution.
  keep <- lapply(seq_len(nconds), function(j) {
    RTx <- sort(RT[[j]], index.return = TRUE)
    RTj <- RTx$x
    CRj <- as.logical(CR[[j]])[RTx$ix]
    RTj[CRj & is.finite(RTj)]
  })
  G <- lapply(keep, function(x) {
    switch(cdf_method,
           ecdf = stats::ecdf(x),
           polygon = .make_quantile_cdf(x, qtype = qtype),
           kernel = .make_kernel_cdf(x, bw = kernel_bw))
  })
  times <- sort(unique(c(RT, recursive = TRUE)))
  times <- times[is.finite(times)]
  if (cdf_method == "kernel") {
    lower <- min(vapply(G, function(f) exp(log(attr(f, "range")[1L]) - 8 * attr(f, "bw")), numeric(1)))
    upper <- max(vapply(G, function(f) exp(log(attr(f, "range")[2L]) + 8 * attr(f, "bw")), numeric(1)))
    # A dense common grid preserves the smooth channel CDFs when the final
    # algebraic bounds are returned as approxfun objects.
    times <- sort(unique(c(times, seq(lower, upper, length.out = 4096L))))
  }
  nt <- length(times)
  Gmat <- vapply(G, function(f) f(times), numeric(nt))
  if (nconds == 1L) Gmat <- matrix(Gmat, ncol = 1L)
  Smat <- 1 - Gmat
  
  #given the stopping rule, estimate the CDF bounds
  if (rule == "OR"){  
    #two channel bounds
    if (n.channels == 2){
      if (assume.ID){
        if (unified.space){
          upper <- 2*Smat[,1]-1
          lower <- Smat[,1]
        }
        else {upper <- 2*Gmat[,1]
              lower <- Gmat[,1]}
      }
      #not ID, 2 channels
      else {
        if (unified.space){
          upper <- apply(Smat,1,sum)-1
          lower <- apply(Smat,1,min)
        }
        else{upper <- apply(Gmat, 1, sum)
             lower <- apply(Gmat, 1, max)} 
      }#end notID, 2 channels
    }#end (n.channels==2)
    #more than two channels
    if (n.channels > 2){
      if (assume.ID){
        #upper is 2Fn-1(t)-Fn-2(t)
        #lower is Fn-1(t)
        if (unified.space){
          upper <- 2*(Smat[,1]^(n.channels-1)) - Smat[,1]^(n.channels-2)
          lower <- Smat[,1]^(n.channels-1)
        }
        #note that min time Fn(t)=1-Gn(t) is 1-max time
        else {upper <- 2*(1-Gmat[,1]^(n.channels-1)) - (1-Gmat[,1]^(n.channels-2))
              lower <- 1-(Gmat[,1]^(n.channels-1)) }
      }  #end (assume.ID) 
      #not ID, >2 channels
      else {
        #find all n-1 distributions
        Fi.array <- matrix(rep(0,nt*nconds), nrow=nt, ncol=nconds)
        for (i.out in 1:nconds){
          Fi.array[,i.out] <- 1-apply(as.matrix(1-Gmat[,-i.out]), 1, prod)
        }
        i.out <- NULL
        
        num.cols <- (nconds*(nconds-1))/2
        Fij.array <- matrix(rep(0,nt*num.cols), nrow=nt, ncol=num.cols)
        Fij.array2 <- matrix(rep(0,nt*num.cols), nrow=nt, ncol=num.cols)
        column <- 1
        for (i.out in 1:(nconds-1)){
          for (j.out in (i.out+1):nconds){
            #CDF values for unified space
            Fij.array[,column] <- (1-apply(as.matrix(1-Gmat[,-i.out]), 1, prod)) + (1-apply(as.matrix(1-Gmat[,-j.out]), 1, prod)) - (1-apply(as.matrix(1-Gmat[,c(-i.out, -j.out)]), 1, prod))
            #survivor values for non-unified space
            Fij.array2[,column] <- apply(as.matrix(1-Gmat[,-i.out]), 1, prod) + apply(as.matrix(1-Gmat[,-j.out]), 1, prod) - apply(as.matrix(1-Gmat[,c(-i.out, -j.out)]), 1, prod)
            column <- column+1
          }
        if (unified.space){
          upper <- apply(Fij.array2, 1, max)
          lower <- apply(1-Fi.array, 1, min)
        }
        else{
          #min{i,j}[1-G^(i)_(n-1)(t) + 1-G^(j)_(n-1)(t) - 1-G^(i,j)_(n-2)(t)] for i!=j
          upper <- apply(Fij.array, 1, min)
          #max{i} F^(i)_(n-1)(t) = 1-G^(i)_(n-1)(t) -- max of all n-1 cdfs
          lower <- apply(Fi.array, 1, max)
          }
        }#end else [ie. not unified.space]
      }#end else [ie. not assume.ID]
     }#end (n.channels >2)
      if (unified.space){
        ucip <- apply(Smat, 1, prod)
        upper <- log(upper)/log(ucip)
        lower <- log(lower)/log(ucip)
      }
    }#end (rule=="OR")
    
    else if (rule == "AND"){
      if (n.channels == 2) {
        if (assume.ID){
          upper <- Gmat[,1]
          lower <- 2*Gmat[,1] - 1
        }
        else{upper <- apply(Gmat, 1, min)
             lower <- apply(Gmat, 1, sum) - 1
        }
      } #end (n.channels==2)
      if (n.channels > 2){
        if (assume.ID){
          #upper is G_(n-1)(t)
          #lower is 2G_(n-1)(t)-G_(n-2)(t)
          upper <- Gmat[,1]^(n.channels-1)
          lower <- (2*(Gmat[,1]^(n.channels-1))) - Gmat[,1]^(n.channels-2)
        }
        else {
          #min{i} G^(i)_(n-1)(t)  -- min of all n-1 cdfs
          Gi.array <- matrix(rep(0,nt*nconds), nrow=nt, ncol=nconds)
          for (i.out in 1:nconds){
            Gi.array[,i.out] <- apply(as.matrix(Gmat[,-i.out]), 1, prod)
          }
          upper <- apply(Gi.array, 1, min)
          i.out <- NULL
          
          #max{i,j}[G^(i)_(n-1)(t) + G^(j)_(n-1)(t) - G^(i,j)_(n-2)(t)] for i!=j
          num.cols <- (nconds*(nconds-1))/2
          Gij.array <- matrix(rep(0,nt*num.cols), nrow=nt, ncol=num.cols)
          column <- 1
          for (i.out in 1:(nconds-1)){
            for (j.out in (i.out+1):nconds){
              Gij.array[,column] <- apply(as.matrix(Gmat[,-i.out]), 1, prod) + apply(as.matrix(Gmat[,-j.out]), 1, prod) - apply(as.matrix(Gmat[,c(-i.out, -j.out)]), 1, prod)
              column <- column+1
            }
          }
          lower <- apply(Gij.array, 1, max)
        }
      } #end (n.channels>2)
      if (unified.space){
        ucip <- apply(Gmat, 1, prod)
        upper <- log(ucip)/log(upper)
        lower <- log(ucip)/log(lower)
      }
    }
    else if (rule == "STST"){
      if (assume.ID){
        upper <- n.channels * Gmat[,1]
        lower <- Gmat[,1]^(n.channels)
      }
      else {upper <- apply(Gmat, 1, sum)
            lower <- apply(Gmat, 1, prod)     
      }
      
      if (unified.space){
        ucip <- Gmat[,1]
        upper <- log(ucip)/log(upper)
        lower <- log(ucip)/log(lower)
      }
    }
    #create an interpolation function for each bound
    upper <- approxfun(times, upper)
    lower <- approxfun(times, lower)
    #return the bounds as a list
    return(list(Upper.Bound=upper, Lower.Bound=lower))
  }
