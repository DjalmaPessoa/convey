#' Gini coefficient
#'
#' Estimate the Gini coefficient, a measure of inequalty
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
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Djalma Pessoa and Anthony Damico
#'
#' @seealso \code{\link{svyarpr}}
#'
#' @references Guillaume Osier (2009). Variance estimation for complex indicators
#' of poverty and inequality. \emph{Journal of the European Survey Research
#' Association}, Vol.3, No.3, pp. 167-195,
#' ISSN 1864-3361, URL \url{http://ojs.ub.uni-konstanz.de/srm/article/view/369}.
#'
#' Jean-Claude Deville (1999). Variance estimation for complex statistics and estimators:
#' linearization and residual techniques. Survey Methodology, 25, 193-203,
#' URL \url{http://www5.statcan.gc.ca/bsolc/olc-cel/olc-cel?lang=eng&catno=12-001-X19990024882}.
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
#' svygini( ~eqincome , design = des_eusilc )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' svygini( ~eqincome , design = des_eusilc_rep )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svygini( ~ py010n , design = des_eusilc )
#' svygini( ~ py010n , design = des_eusilc , na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svygini( ~ py010n , design = des_eusilc_rep )
#' svygini( ~ py010n , design = des_eusilc_rep , na.rm = TRUE )
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
#' svygini( ~ eqincome , design = dbd_eusilc )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svygini <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svygini", design)

  }


#' @rdname svygini
#' @export
svygini.survey.design <-
  function( formula , design , na.rm=FALSE , deff=FALSE , linearized = FALSE , ... ) {

    # collect income data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if (na.rm) {
      nas <- is.na(incvar)
      design <- design[nas == 0, ]
      if (length(nas) > length(design$prob)) incvar <- incvar[nas == 0] else incvar[nas > 0] <- 0
    }

    # collect sampling weights
    w <- 1/design$prob

    # compute point estimate
    estimate <- CalcGini( incvar , w )

    # compute linearized function
    lin <- CalcGini_IF( incvar , w )

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
    attr(rval, "statistic") <- "gini"
    if ( linearized ) attr(rval,"linearized") <- lin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }

#' @rdname svygini
#' @export
svygini.svyrep.design <-
  function( formula , design , na.rm=FALSE , deff=FALSE , linearized=FALSE , return.replicates = FALSE , ...) {

    # collect data
    df <- model.frame(design)
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if( na.rm ) {
      nas <- is.na( incvar )
      design <- design[!nas,]
      df <- model.frame( design )
      incvar <- incvar[ !nas ]
    }

    # colelct sampling weights
    ws <- weights(design, "sampling")

    # compute point estimate
    estimate <- CalcGini(incvar, ws)

    # collect analysis weights
    ww <- weights(design, "analysis")

    # compute replicates
    qq <- apply(ww, 2, function(wi) CalcGini(incvar, wi))

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
      lin <- CalcGini_IF( incvar , ws )

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
    attr(rval, "statistic") <- "gini"
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


#' @rdname svygini
#' @export
svygini.DBIsvydesign <-
  function (formula, design, ...){

    design$variables <-
      getvars(
        formula,
        design$db$connection,
        design$db$tablename,
        updates = design$updates,
        subset = design$subset
      )

    NextMethod("svygini", design)
  }


# gini estimate function
CalcGini <-
  function(x, pw) {

    # filter observations
    x <- x[pw>0]
    pw <- pw[pw>0]

    # reorder
    pw <- pw[order(x)]
    x <- x[order(x)]

    # intermediate estimates
    N <- sum(pw)
    n <- length(x)
    big_t <- sum(x * pw)
    r <- cumsum(pw)
    Num <- sum((2 * r - 1) * x * pw)
    Den <- N * big_t

    # gini estimate
    (Num/Den) - 1

  }

# gini linearized function
CalcGini_IF <- function( x , pw ) {

  # filter observations
  x <- x[ pw > 0 ]
  pw <- pw[ pw > 0 ]

  # collect indices
  ind <- names( pw )

  # reorder observations
  ordx <- order(x)
  pw <- pw[ordx]
  x <- x[ordx]

  # population size
  N <- sum(pw)

  # total income
  Y <- sum(x * pw)

  # cumulative weight
  r <- cumsum(pw)

  # partial weighted function
  G <- cumsum(x * pw)
  T1 <- list( value = sum( r * x * pw ) , lin = ( Y - G + x * pw + r * x ) )
  T2 <- list( value = sum( x * pw ), lin = x)
  T3 <- list( value = sum( pw ) , lin = rep( 1 , length( x ) ) )

  # get T1
  list_all <- list(T1 = T1, T2 = T2, T3 = T3)
  GINI <- contrastinf( quote( ( 2 * T1 - T2 ) / ( T2 * T3 ) - 1 ) , list_all )
  lingini <- as.numeric( GINI$lin )

  # flip back to original order
  lingini <- lingini[order(ordx)]

  # add indices
  names(lingini) <- ind

  # return object
  return( lingini )

}
