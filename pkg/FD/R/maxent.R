maxent <- function (c.means, c.mat, prior, tol = 1e-08, lambda = FALSE){
     
     # check input
    if (!is.numeric(c.means)) stop("c.means must be a numeric vector\n")
    if (!is.numeric(tol)) stop("tol must be numeric\n")
    if (!is.logical(lambda)) stop("lambda must be logical\n")
    if (!is.matrix(c.mat)){
       if (is.numeric(c.mat)){
          s.names <- names(c.mat)
          c.mat <- matrix(c.mat, nrow = 1, ncol = length(c.mat) )
          dimnames(c.mat) <- list("constraint", s.names)
         }
       else stop("if c.mat is not a matrix then it can only be a numeric vector\n")
     }
    s.names <- dimnames(c.mat)[[2]]
    c.names <- dimnames(c.mat)[[1]]
    dim.matrix <- dim(c.mat)
    if (dim.matrix[2] == 1 && dim.matrix[1] > 1){
        c.mat <- t(c.mat)
        dim.matrix <- dim(c.mat)
     }
    n.species <- dim.matrix[2]
    n.traits <- dim.matrix[1]
    if (missing(prior) ) prior <- rep(1 / n.species, n.species)
    if (length(prior) != n.species) stop("number of states in prior not equal to number in c.mat\n")
    if (any(is.na(c.means)) || any(is.na(c.mat)) || any(is.na(prior)) ) stop("no NA's allowed\n")
    n.constraints <- length(c.means)
    if (n.constraints != dim.matrix[1]) stop("number of constraint means not equal to number of constraints in c.mat\n")

    # FORTRAN loop  
    itscale <- .Fortran("itscale5", as.double(t(c.mat)), as.integer(n.species), as.integer(n.traits), as.double(c.means), as.double(prior), prob = double(n.species), entropy = double(1), niter = integer(1), as.double(tol), moments = double(n.traits), PACKAGE = "FD")

  # output
  res <- list()
  prob <- itscale$prob
  names(prob) <- s.names
  res$prob <- prob
  moments <- itscale$moments
  names(moments) <- c.names
  res$moments <- moments
  res$entropy <- itscale$entropy
  res$iter <- itscale$niter
  if (lambda){
    lambda <- coef(lm(log(prob) ~ t(c.mat) ) )
    names(lambda) <- c("intercept", c.names)
    res$lambda <- lambda
   }
  return(res)
}

