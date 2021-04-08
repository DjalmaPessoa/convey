#' J-divergence measure
#'
#' Estimate the J-divergence measure, an entropy-based measure of inequality
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param return.replicates Return the replicate estimates?
#' @param ... future expansion
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' This measure only allows for strictly positive variables.
#'
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Guilherme Jacob
#'
#' @note This function is experimental and is subject to change in later versions.
#'
#' @seealso \code{\link{svygei}}
#'
#' @references Nicholas Rohde (2016). J-divergence measurements of economic inequality.
#' J. R. Statist. Soc. A, v. 179, Part 3 (2016), pp. 847-870.
#' URL \url{http://onlinelibrary.wiley.com/doi/10.1111/rssa.12153/abstract}.
#'
#' Martin Biewen and Stephen Jenkins (2002). Estimation of Generalized Entropy
#' and Atkinson Inequality Indices from Complex Survey Data. \emph{DIW Discussion Papers},
#' No.345,
#' URL \url{https://www.diw.de/documents/publikationen/73/diw_01.c.40394.de/dp345.pdf}.
#'
#' @keywords survey
#'
#' @examples
#' library(survey)
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#'
#' # linearized design
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep(des_eusilc)
#'
#' svyjdiv( ~eqincome , design = subset( des_eusilc , eqincome > 0 ) )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' svyjdiv( ~eqincome , design = subset( des_eusilc_rep , eqincome > 0 ) )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyjdiv( ~py010n , design = subset( des_eusilc , py010n > 0 | is.na( py010n ) ) )
#' svyjdiv( ~py010n , design = subset( des_eusilc , py010n > 0 | is.na( py010n ) ), na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyjdiv( ~py010n , design = subset( des_eusilc_rep , py010n > 0 | is.na( py010n ) ) )
#' svyjdiv( ~py010n , design = subset( des_eusilc_rep , py010n > 0 | is.na( py010n ) ) , na.rm = TRUE )
#'
#' # database-backed design
#' library(RSQLite)
#' library(DBI)
#' dbfile <- tempfile()
#' conn <- dbConnect( RSQLite::SQLite() , dbfile )
#' dbWriteTable( conn , 'eusilc' , eusilc )
#'
#' dbd_eusilc <-
#' 	svydesign(
#' 		ids = ~rb030 ,
#' 		strata = ~db040 ,
#' 		weights = ~rb050 ,
#' 		data="eusilc",
#' 		dbname=dbfile,
#' 		dbtype="SQLite"
#' 	)
#'
#' dbd_eusilc <- convey_prep( dbd_eusilc )
#'
#' svyjdiv( ~eqincome , design = subset( dbd_eusilc , eqincome > 0 ) )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyjdiv <- function(formula, design, ...) {

  if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

  UseMethod("svyjdiv", design)

}

#' @rdname svyjdiv
#' @export
svyjdiv.survey.design <- function ( formula, design, na.rm = FALSE, deff=FALSE , linearized = FALSE , ... ) {

  # collect income data
  incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

  # fix missing values
  if (na.rm) {
    nas <- is.na(incvar)
    design <- design[nas == 0, ]
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
  }

  # collect sampling weights
  w <- 1/design$prob

  # test for positive income
  if ( any( incvar[ w > 0 ] <= 0 , na.rm = TRUE ) ) stop( "The J-divergence measure is defined for strictly positive variables only.  Negative and zero values not allowed." )

  # # method 1: compute components
  # U_0 <- list( value = sum( w ), lin = rep( 1, length( incvar ) ) )
  # U_1 <- list( value = sum( w * incvar ), lin = incvar )
  # T_0 <- list( value = sum( w * log( incvar ) ), lin = log( incvar ) )
  # T_1 <- list( value = sum( w * incvar * log( incvar ) ), lin = incvar * log( incvar ) )
  # list_all <- list(  U_0 = U_0, U_1 = U_1, T_0 = T_0, T_1 = T_1 )
  # estimate <- contrastinf( quote( ( T_1 / U_1 ) - ( T_0 / U_0 ) ) , list_all )
  # rval <- estimate$value
  # lin <- estimate$lin

  # method 2: compute point estimate
  estimate <- CalcJDiv( incvar , w )
  lin <- CalcJDiv_IF( incvar , w )
  names( lin ) <- rownames( design$variables )[ w > 0 ]

  # treat remaining missing
  if ( is.na( estimate ) ) {
    rval <- NA
    variance <- as.matrix(NA)
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr(rval, "statistic") <- "j-divergence"
    attr(rval, "var") <- variance
    return(rval)
  }

  # ensure length
  if ( length( lin ) != length( design$prob ) ) {
    tmplin <- rep( 0 , nrow( design$variables ) )
    tmplin[ w > 0 ] <- lin
    lin <- tmplin ; rm( tmplin )
    names( lin ) <- rownames( design$variables )
  }

  # compute variance
  variance <- survey::svyrecvar( lin/design$prob, design$cluster, design$strata, design$fpc, postStrata = design$postStrata )
  variance[ which( is.nan( variance ) ) ] <- NA
  colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

  # compute deff
  if ( is.character(deff) || deff) {
    nobs <- sum( weights( design ) != 0 )
    npop <- sum( weights( design ) )
    if (deff == "replace") vsrs <- survey::svyvar( lin , design, na.rm = na.rm) * npop^2/nobs
    else vsrs <- survey::svyvar( lin , design , na.rm = na.rm ) * npop^2 * (npop - nobs)/(npop * nobs)
    deff.estimate <- variance/vsrs
  }

  # keep necessary linearized functions
  lin <- lin[ 1/design$prob > 0 ]

  # coerce to matrix
  lin <- matrix( lin , nrow = length( lin ) , dimnames = list( names( lin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

  # build result object
  rval <- estimate
  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
  class(rval) <- c( "cvystat" , "svystat" )
  attr(rval, "var") <- variance
  attr(rval, "statistic") <- "j-divergence"
  if ( linearized ) attr(rval,"linearized") <- lin
  if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
  if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
  rval

}


#' @rdname svyjdiv
#' @export
svyjdiv.svyrep.design <- function ( formula, design, na.rm = FALSE, deff = FALSE , linearized = FALSE , return.replicates = FALSE , ... ) {

  # collect income data
  incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

  if(na.rm){
    nas<-is.na(incvar)
    design<-design[!nas,]
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
  }

  # collect sampling weights
  ws <- weights(design, "sampling")

  # check for positive incomes
  if ( any( incvar[ws != 0] <= 0, na.rm = TRUE ) ) stop( "The J-divergence measure is defined for strictly positive variables only.  Negative and zero values not allowed." )

  # compute point estimate
  estimate <- CalcJDiv( incvar, ws )

  # treat remaining missing
  if ( is.na( estimate ) ) {
    rval <- estimate
    variance <- as.matrix(NA)
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr(rval, "statistic") <- "j-divergence"
    attr(rval, "var") <- variance
    return(rval)
  }

  ### variance calculation

  # collect analysis weights
  wf <- weights(design, "analysis")

  # compute replicates
  qq <- apply( wf, 2 , function(wi) CalcJDiv( incvar , wi ) )

  # compute variance
  if ( any( is.na( qq ) ) ) variance <- as.matrix( NA ) else {
    variance <- survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = estimate )
    this.mean <- attr( variance , "means" )
    variance <- as.matrix( variance )
    attr( variance , "means" ) <- this.mean
  }
  colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

  # compute deff
  if ( is.character(deff) || deff || linearized ) {

    # compute linearized function
    lin <- CalcJDiv_IF( incvar , ws )

    # compute deff
    nobs <- length( design$pweights )
    npop <- sum( design$pweights )
    vsrs <- unclass( survey::svyvar( lin , design, na.rm = na.rm, return.replicates = FALSE, estimate.only = TRUE)) * npop^2/nobs
    if (deff != "replace") vsrs <- vsrs * (npop - nobs)/npop
    deff.estimate <- variance / vsrs

    # filter observation
    names( lin ) <- rownames( design$variables )

    # coerce to matrix
    lin <- matrix( lin , nrow = length( lin ) , dimnames = list( names( lin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

  }

  # build result object
  rval <- estimate
  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
  class(rval) <- c( "cvystat" , "svrepstat" )
  attr(rval, "var") <- variance
  attr(rval, "statistic") <- "j-divergence"
  if ( linearized ) attr( rval , "linearized" ) <- lin
  if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )

  # keep replicates
  if (return.replicates) {
    attr( qq , "scale") <- design$scale
    attr( qq , "rscales") <- design$rscales
    attr( qq , "mse") <- design$mse
    rval <- list( mean = rval , replicates = qq )
    class( rval ) <- c( "cvystat" , "svrepstat" )
  }

  # add design effect estimate
  if ( is.character(deff) || deff) attr( rval , "deff" ) <- deff.estimate

  # return object
  rval

}

#' @rdname svyjdiv
#' @export
svyjdiv.DBIsvydesign <-
  function (formula, design, ...) {

    design$variables <- getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)

    NextMethod("svyjdiv", design)
  }


# point estimate function
CalcJDiv <-  function( y , w ) {

  # filter observations
  y <- y[ w > 0]
  w <- w[ w > 0]

  # compute point esitmate
  N <- sum( w )
  mu <- sum( y * w ) / N
  jdiv <- ( ( y - mu ) / mu ) * log( y / mu )

  # compute point estimate
  sum( jdiv * w ) / N

}

# function to compute linearized function
CalcJDiv_IF <-  function( y , w ) {

  # filter observations
  y <- y[ w > 0]
  w <- w[ w > 0]

  # compute intermediate statistics
  N <- sum( w )
  Y <- sum( y * w )
  mu <- Y / N
  jdiv <- sum( w * ( ( y / mu ) - 1 ) * log( y / mu ) ) / N
  gei1 <- sum( w * (y/mu) * log(y/mu) ) / N

  # linearized function under fixed mean
  y.score <- ( ( ( y / mu ) - 1 ) * log( y / mu ) )
  lin.fixed <- ( y.score - jdiv ) / N

  # derivative wrt mean
  # djdiv.dmu <- 1/mu - sum( w * (y/mu) * ( log( y/mu ) + 1 ) ) / Y
  djdiv.dmu <- - gei1 / mu
  I.mu <- ( y - mu ) / N

  # compute final linearized function
  lin <- lin.fixed + djdiv.dmu * I.mu

  # return final linearized function
  return( lin )

}
