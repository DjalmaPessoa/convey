#' Linearization of the total below a quantile
#'
#' Computes the linearized variable of the total in the lower tail of
#' the distribution of a variable.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param alpha the order of the quantile
#' @param upper return the total in the total in the upper tail. Defaults to \code{FALSE}.
#' @param quantile return the upper bound of the lower tail
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param return.replicates Return the replicate estimates?
#' @param ... arguments passed on to `survey::svyquantile`
#'
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
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
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#' library(survey)
#' des_eusilc <- svydesign(ids = ~rb030, strata =~db040,  weights = ~rb050, data = eusilc)
#' des_eusilc <- convey_prep(des_eusilc)
#' svyisq(~eqincome, design=des_eusilc,.20 , quantile = TRUE)
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' svyisq( ~eqincome , design = des_eusilc_rep, .20 , quantile = TRUE )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyisq( ~ py010n , design = des_eusilc, .20 )
#' svyisq( ~ py010n , design = des_eusilc , .20, na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyisq( ~ py010n , design = des_eusilc_rep, .20 )
#' svyisq( ~ py010n , design = des_eusilc_rep , .20,  na.rm = TRUE )
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
#' svyisq( ~ eqincome , design = dbd_eusilc, .20 )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyisq <-
	function(formula, design, ...) {

		if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

		UseMethod("svyisq", design)

	}

#' @rdname svyisq
#' @export
svyisq.survey.design <-
	function(formula, design, alpha, quantile = FALSE, upper = FALSE , na.rm = FALSE, deff = FALSE , linearized = FALSE , ...) {

	  # test for convey_prep
		if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

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

	  # store quantile
	  q_alpha <- survey::svyquantile(x = formula, design = design, quantiles = alpha, method = "constant", na.rm = na.rm,...)

	  # compute value
	  estimate <- CalcISQ( incvar , w , alpha )
	  if (upper) estimate <- sum( w * incvar ) - estimate

	  # compute linearized functions
	  h <- h_fun(incvar, w)
	  Fprime0 <- densfun( formula = formula, design = design, q_alpha[[1]] , FUN = "F",     na.rm = na.rm )
	  Fprime1 <- densfun( formula = formula, design = design, q_alpha[[1]] , FUN = "big_s", na.rm = na.rm )
	  lin <- CalcISQ_IF( incvar , w, alpha , Fprime0 , Fprime1 )
	  if (upper) lin <- incvar[ w > 0 ] - lin

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
	    nobs <- sum( weights( design , "sampling" ) > 0 )
	    npop <- sum( weights( design , "sampling" ) )
	    if ( deff == "replace" ) vsrs <- survey::svyvar( lin , design, na.rm = na.rm) * npop^2/nobs
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
	  attr(rval, "statistic") <- "isq"
		if(quantile) attr(rval, "quantile") <- q_alpha
	  if ( linearized ) attr(rval,"linearized") <- lin
	  if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
	  if ( is.character(deff) || deff ) attr(rval,"deff") <- deff.estimate
	  rval

	}

#' @rdname svyisq
#' @export
svyisq.svyrep.design <-
	function(formula, design, alpha,quantile = FALSE, upper = FALSE , na.rm = FALSE, deff = FALSE , linearized=FALSE , return.replicates = FALSE , ...){

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

	  # store quantile
	  if (quantile) q_alpha <- survey::svyquantile(x = formula, design = design, quantiles = alpha, method = "constant", na.rm = na.rm,...)

	  # compute point estimate
	  q_alpha <- computeQuantiles( incvar , ws, alpha )
	  estimate <- CalcISQ( incvar , ws, alpha )
	  if (upper) estimate <- sum( ws * incvar ) - estimate

	  # collect analysis weights
	  ww <- weights(design, "analysis")

	  # compute replicates
	  qq <- apply( ww, 2 , function(wi) {
	    if ( upper ) sum( wi * incvar ) - CalcISQ( incvar , wi , alpha ) else CalcISQ( incvar , wi , alpha )
	  } )

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

	    # compute linearized functions
	    h <- h_fun( incvar, ws )
	    Fprime0 <- densfun( formula = formula, design = design , q_alpha[[1]] , FUN = "F",     na.rm = na.rm )
	    Fprime1 <- densfun( formula = formula, design = design , q_alpha[[1]] , FUN = "big_s", na.rm = na.rm )
	    lin <- CalcISQ_IF( incvar , ws, alpha , Fprime0 , Fprime1 )
	    if (upper) lin <- incvar[ ws > 0 ] - lin

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
	  attr(rval, "statistic") <- "isq"
	  if(quantile) attr(rval, "quantile") <- q_alpha
	  if ( linearized ) attr(rval,"linearized") <- lin
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

#' @rdname svyisq
#' @export
svyisq.DBIsvydesign <-
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

    NextMethod("svyisq", design)
  }


# function por point estimates
CalcISQ <- function( x , pw , alpha ) {

  # filter observations
  x <- x [ pw > 0 ]
  pw <- pw [ pw > 0 ]

  # compute quantile
  q_alpha <- computeQuantiles( x , pw , alpha )

  # compute total below quantile
  sum( x * ( x <= q_alpha ) * pw )

}

# function for linearized functions
CalcISQ_IF <- function( x , pw , alpha , Fprime0 , Fprime1 ) {

  # filter observations
  x <- x [ pw > 0 ]
  pw <- pw [ pw > 0 ]

  # population size
  N <- sum( pw )

  # compute quantile
  q_alpha <- computeQuantiles( x , pw , alpha )

  # linearization
  h <- h_fun(x, pw)
  # Fprime0 <- CalcDensFun( x , pw , q_alpha , h=h , FUN = "F" )
  # Fprime1 <- CalcDensFun( x , pw , q_alpha , FUN = "big_s" )
  iq <- -( 1 / ( N * Fprime0 ) ) * ( ( x <= q_alpha ) - alpha )
  isqalpha1 <- x * (x <= q_alpha)
  isqalpha <- isqalpha1 + Fprime1 * iq

  # add indices
  names( isqalpha ) <- names( pw )

  # return estimate
  isqalpha

}

# function for density estimation
CalcDensFun <- function( x , pw , q_alpha , h = NULL , FUN = "F" ) {

  # filter observations
  x <- x[pw>0]
  pw <- pw[pw>0]

  # intermediate estimates
  N <- sum(pw)
  if(is.null(h)) h <- h_fun(x,pw)

  # calculation
  u <- (x - x)/h
  vectf <- exp(-(u^2)/2)/sqrt(2 * pi)
  if (FUN == "F") {
    res <- sum(vectf * pw)/(N * h)
  } else {
    v <- pw * x
    res <- sum(vectf * v)/h
  }

  # final estimate
  res

}
