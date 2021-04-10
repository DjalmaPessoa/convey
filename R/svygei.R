#' Generalized entropy index
#'
#' Estimate the generalized entropy index, a measure of inequality
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param epsilon a parameter that determines the sensivity towards inequality in the top of the distribution. Defaults to epsilon = 1.
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param influence Should a matrix of (weighted) influence functions be returned? (for compatibility with \code{\link[survey]{survey}})
#' @param return.replicates Return the replicate estimates?
#' @param ... future expansion
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' This measure only allows for strictly positive variables.
#'
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Guilherme Jacob, Djalma Pessoa and Anthony Damico
#'
#' @seealso \code{\link{svyatk}}
#'
#' @references Matti Langel (2012). Measuring inequality in finite population sampling.
#' PhD thesis: Universite de Neuchatel,
#' URL \url{https://doc.rero.ch/record/29204/files/00002252.pdf}.
#'
#' Martin Biewen and Stephen Jenkins (2002). Estimation of Generalized Entropy
#' and Atkinson Inequality Indices from Complex Survey Data. \emph{DIW Discussion Papers},
#' No.345,
#' URL \url{https://www.diw.de/documents/publikationen/73/diw_01.c.40394.de/dp345.pdf}.
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
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' # linearized design
#' svygei( ~eqincome , subset(des_eusilc, eqincome > 0), epsilon = 0 )
#' svygei( ~eqincome , subset(des_eusilc, eqincome > 0), epsilon = .5 )
#' svygei( ~eqincome , subset(des_eusilc, eqincome > 0), epsilon = 1 )
#' svygei( ~eqincome , subset(des_eusilc, eqincome > 0), epsilon = 2 )
#'
#' # replicate-weighted design
#' svygei( ~eqincome , subset(des_eusilc_rep, eqincome > 0), epsilon = 0 )
#' svygei( ~eqincome , subset(des_eusilc_rep, eqincome > 0), epsilon = .5 )
#' svygei( ~eqincome , subset(des_eusilc_rep, eqincome > 0), epsilon = 1 )
#' svygei( ~eqincome , subset(des_eusilc_rep, eqincome > 0), epsilon = 2 )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 0 )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 0, na.rm = TRUE )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = .5 )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = .5, na.rm = TRUE )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 1 )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 1, na.rm = TRUE )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 2 )
#' svygei( ~py010n , subset(des_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 2, na.rm = TRUE )
#'
#' # replicate-weighted design using a variable with missings
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = 0 )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = 0, na.rm = TRUE )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = .5 )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = .5, na.rm = TRUE )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = 1 )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = 1, na.rm = TRUE )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = 2 )
#' svygei( ~py010n , subset(des_eusilc_rep, py010n > 0 | is.na(py010n) ), epsilon = 2, na.rm = TRUE )
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
#' # database-backed linearized design
#' svygei( ~eqincome , subset(dbd_eusilc, eqincome > 0), epsilon = 0 )
#' svygei( ~eqincome , dbd_eusilc, epsilon = .5 )
#' svygei( ~eqincome , subset(dbd_eusilc, eqincome > 0), epsilon = 1 )
#' svygei( ~eqincome , dbd_eusilc, epsilon = 2 )
#'
#' # database-backed linearized design using a variable with missings
#' svygei( ~py010n , subset(dbd_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 0 )
#' svygei( ~py010n , subset(dbd_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 0, na.rm = TRUE )
#' svygei( ~py010n , dbd_eusilc, epsilon = .5 )
#' svygei( ~py010n , dbd_eusilc, epsilon = .5, na.rm = TRUE )
#' svygei( ~py010n , subset(dbd_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 1 )
#' svygei( ~py010n , subset(dbd_eusilc, py010n > 0 | is.na(py010n) ), epsilon = 1, na.rm = TRUE )
#' svygei( ~py010n , dbd_eusilc, epsilon = 2 )
#' svygei( ~py010n , dbd_eusilc, epsilon = 2, na.rm = TRUE )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svygei <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    #if( 'epsilon' %in% names( list(...) ) && list(...)[["epsilon"]] < 0 ) stop( "epsilon= cannot be negative." )

    UseMethod("svygei", design)

  }


#' @rdname svygei
#' @export
svygei.survey.design <-
  function ( formula, design, epsilon = 1, na.rm = FALSE, deff= FALSE , linearized = FALSE , influence = FALSE , ... ) {

    # collect income data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if (na.rm) {
      nas <- is.na(incvar)
      design <- design[!nas, ]
      if (length(nas) > length(design$prob) ) incvar <- incvar[nas == 0] else incvar[nas > 0] <- 0
    }

    # collect weights
    w <- 1/design$prob

    # check for strictly positive incomes
    if ( any( incvar[w > 0] <= 0 , na.rm = TRUE ) ) stop( "The GEI indices are defined for strictly positive variables only.\nNegative and zero values not allowed." )

    # compute value
    estimate <- CalcGEI( x = incvar, weights = w, epsilon = epsilon )

    # compute linearized functions
    lin <- CalcGEI_IF( x = incvar, weights = w, epsilon = epsilon )

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
    if ( is.character( deff ) || deff ) {
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
    attr(rval, "statistic") <- "gei"
    attr(rval,"epsilon")<- epsilon
    if ( linearized ) attr(rval,"linearized") <- lin
    if ( influence )  attr( rval , "influence" )  <- sweep( lin , 1 , design$prob[ is.finite( design$prob ) ] , "/" )
    if ( linearized | influence ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }


#' @rdname svygei
#' @export
svygei.svyrep.design <-
  function( formula , design , epsilon = 1 , na.rm=FALSE , deff= FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

    # collect income variable
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missings
    if(na.rm){
      nas<-is.na(incvar)
      design<-design[!nas,]
      df <- model.frame(design)
      incvar <- incvar[!nas]
    }

    # collect sampling weights
    ws <- weights(design, "sampling")

    # check for strictly positive incomes
    if ( any( incvar[ ws != 0 ] == 0, na.rm = TRUE) ) stop( "The GEI indices are defined for strictly positive variables only.\nNegative and zero values not allowed." )

    # compute point estimate
    estimate <- CalcGEI( x = incvar, weights = ws, epsilon = epsilon)

    # collect analysis weights
    ww <- weights(design, "analysis")

    # compute replicates
    qq <- apply( ww, 2 , function(wi) CalcGEI( incvar , wi , epsilon = epsilon ) )

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
      lin <- CalcGEI_IF( incvar , ws , epsilon )

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
    attr(rval, "statistic") <- "gei"
    attr(rval,"epsilon") <- epsilon
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


#' @rdname svygei
#' @export
svygei.DBIsvydesign <-
  function (formula, design, ...) {

    design$variables <- getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)

    NextMethod("svygei", design)
  }

# function for point estimates
CalcGEI <-
  function( x, weights, epsilon ) {

    # filter observations
    x <- x[weights > 0 ]
    weights <- weights[weights > 0 ]

    # intermediate stats
    N <- sum( weights )
    Y <- sum( weights*x )
    mu <- Y/N

    # branch on epsilon
    if ( epsilon == 0 ) {
      estimate <- - sum( weights * log( x/mu ) ) / N
    } else if ( epsilon == 1 ) {
      estimate <- sum( weights * (x/mu) * log( x/mu ) ) / N
    } else {
      estimate <- sum( weights * ( (x/mu)^epsilon - 1 ) ) / ( N * ( epsilon^2 - epsilon ) )
    }

    # return estimate
    return( estimate )

  }

# function for linearized functions
CalcGEI_IF <-
  function( x, weights, epsilon ) {

    # filter observations
    x <- x[ weights > 0 ]
    weights <- weights[ weights > 0 ]

    # compute intermediate statistics
    N <- sum( weights )
    X <- sum( weights * x )
    mu <- X/N

    # compute transformed income
    if ( epsilon == 0 ) {
      y <- - log( x / mu )
    } else if ( epsilon == 1 ) {
      y <- ( x / mu ) * log( x / mu )
    } else {
      y <- ( ( x / mu )^epsilon - 1 ) / ( epsilon^2 - epsilon )
    }
    gei.val <- sum( y * weights ) / N

    # linearized function under fixed mean
    lin.fixed <- ( y - gei.val ) / N

    # compute derivative wrt mean
    if ( epsilon == 0 ) {
      dgei.dmu <- 1/mu
    } else if ( epsilon == 1 ) {
      dgei.dmu <- sum( weights * (-x/mu^2) * ( log( x / mu ) + 1 ) ) / N
      # dgei.dmu <- -( gei.val + 1 ) / mu
    } else {
      dgei.dmu <- - sum( weights  *  (x/mu)^epsilon ) / ( X * ( epsilon - 1 ) )
    }

    # linearized function of the mean
    I.mu <- ( x - mu ) / N

    # linearized function
    lin <- lin.fixed + dgei.dmu * I.mu

    # add indices
    names( lin ) <- names( weights )

    # return linearized function
    return( lin )

  }


