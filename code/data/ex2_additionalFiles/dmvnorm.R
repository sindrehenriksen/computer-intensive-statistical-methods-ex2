## Function for the density function of a multivariate normal distribution
## in canonical representation. 
##
## Arguments
## x : is the vector where we would like to evaluate the density
## b : b vector used in the canonical parameterisation of the normal distribution
## Q : precision matrix (with full rank)
## log : boolean indicating whether the log density should be returned
## memory : argument passed to the spam package (for the exercise this does not need to be changed).

dmvnorm.canonical <- function(x, b, Q, log=TRUE, memory=list(nnzcolindices=6467)){

    # some checks
    if (length(x) != NCOL(Q)) {
        stop("x and Q have non-conforming size")
    }
    if (length(b) != NROW(Q)) {
        stop("b and Q have non-conforming size")
    }
    # compute the log determinant
    logdet <- as.numeric(determinant(Q, logarithm=TRUE, memory=memory)$modulus)
    # get the mean
    mu <- solve.spam(Q, b, memory=memory)
    xmu <- (x-mu)
    # get the log-density
    logdens <- (- length(x) * log(2*pi) + logdet - t(xmu)%*%Q%*%xmu)/2

    if(log)
        return(logdens)
    exp(logdens)
}
