#' Linearization of a variable quantile
#'
#' Computes the linearized variable of a quantile of variable.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param alpha the order of the quantile
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param influence Should a matrix of (weighted) influence functions be returned? (for compatibility with \code{\link[survey]{svyby}})
#' @param return.replicates Return the replicate estimates?
#' @param ... future expansion. not used.`
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
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
#'
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#' library(survey)
#' # linearized design
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep(des_eusilc)
#'
#' svyiqalpha( ~eqincome , design = des_eusilc, .50 )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' svyiqalpha( ~eqincome , design = des_eusilc_rep, .50 )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyiqalpha( ~ py010n , design = des_eusilc, .50 )
#' svyiqalpha( ~ py010n , design = des_eusilc , .50, na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyiqalpha( ~ py010n , design = des_eusilc_rep, .50 )
#' svyiqalpha( ~ py010n , design = des_eusilc_rep ,.50, na.rm = TRUE )
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
#' svyiqalpha( ~ eqincome , design = dbd_eusilc, .50 )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyiqalpha <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svyiqalpha", design)

  }

#' @rdname svyiqalpha
#' @export
svyiqalpha.survey.design <-
  function(formula, design, alpha=.5, na.rm=FALSE, deff= FALSE , linearized = FALSE , influence = FALSE , ...) {

    # collect data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if (na.rm) {
      nas <- is.na(incvar)
      design <- design[!nas, ]
      if (length(nas) > length(design$prob))
        incvar <- incvar[!nas]
      else incvar[nas] <- 0
    }

    # collect weights
    w <- 1/design$prob

    # compute value
    estimate <- computeQuantiles( incvar , w , alpha )

    # compute linearized functions
    # h <- h_fun(incvar, w)
    # Fprime <- densfun( formula = formula, design = design , estimate , h=h, FUN = "F", na.rm=na.rm )
    lin <- CalcQuantile_IF( incvar , w , alpha )

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
    if ( is.character(deff) || deff || linearized ) {
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
    attr(rval, "statistic") <- "quantile"
    if ( linearized ) attr( rval , "linearized" ) <- lin
    if ( influence )  attr( rval , "influence" )  <- sweep( lin , 1 , design$prob[ is.finite( design$prob ) ] , "/" )
    if ( linearized | influence ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }

#' @rdname svyiqalpha
#' @export
#'
svyiqalpha.svyrep.design <-
  function(formula, design, alpha=.5, na.rm=FALSE, deff = FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

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

    # compute point estimate
    estimate <- computeQuantiles( incvar , ws , alpha )

    # collect analysis weights
    ww <- weights(design, "analysis")

    # compute replicates
    qq <- apply( ww, 2 , function(wi) computeQuantiles( incvar , wi , alpha ) )

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
      lin <- CalcQuantile_IF( incvar , ws , alpha )

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
    attr(rval, "statistic") <- "quantile"
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

    # retorna objeto
    rval

  }


#' @rdname svyiqalpha
#' @export
svyiqalpha.DBIsvydesign <-
  function (formula, design, ...){

    if (!( "logical" %in% class(attr(design, "full_design"))) ){

      full_design <- attr( design , "full_design" )

      full_design$variables <-
        getvars(
          formula,
          attr( design , "full_design" )$db$connection,
          attr( design , "full_design" )$db$tablename,
          updates = attr( design , "full_design" )$updates,
          subset = attr( design , "full_design" )$subset
        )

      attr( design , "full_design" ) <- full_design

      rm( full_design )

    }

    design$variables <-
      getvars(
        formula,
        design$db$connection,
        design$db$tablename,
        updates = design$updates,
        subset = design$subset
      )

    NextMethod("svyiqalpha", design)
  }


# function for linearized functions
CalcQuantile_IF <- function(xx, w, p) {

  # filter observations
  xx <- xx[ w > 0 ]
  w <- w[w>0 ]

  # population size
  N <- sum(w)

  # calculate quantiles
  q_alpha <- computeQuantiles( xx , w , p )
  q_alpha <- as.numeric( q_alpha )

  # compute h
  h <- h_fun(xx, w)

  # compute Fprime
  Fprime <- compute_densfun( xx , w , q_alpha , h , "F" )

  # linearized variable
  iq <- -(1/(N * Fprime)) * ((xx <= q_alpha) - p )

  # add indices
  names(iq) <- names( w )

  # return object
  iq

}

# auxilliary function
compute_densfun <- function(incvar, w, x, h = NULL, FUN = "F" ) {

  if( !( FUN %in% c( "F" , "big_s" ) ) ) stop( "valid choices for `FUN=` are 'F' and 'big_s'" )

  # filter observations
  incvar <- incvar[w>0]
  w <- w[w>0]

  # run calculations
  N <- sum(w)
  if(is.null(h)) h <- h_fun(incvar,w)
  u <- (x - incvar)/h
  vectf <- exp(-(u^2)/2)/sqrt(2 * pi)
  if (FUN == "F")
    res <- sum(vectf * w)/(N * h) else {
      v <- w * incvar
      res <- sum(vectf * v)/h
    }
  res

}
