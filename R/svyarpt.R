#' At-risk-of-poverty threshold
#'
#' The standard definition is to use 60\% of the median income.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param quantiles income quantile quantiles, usually .50 (median)
#' @param percent fraction of the quantile, usually .60
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param return.replicates Return the replicate estimates?
#' @param ... arguments passed on to `survey::svyquantile`
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
#'
#' library(survey)
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#'
#' # linearized design
#'
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#' svyarpt( ~eqincome , design = des_eusilc )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#' svyarpt( ~eqincome , design = des_eusilc_rep )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyarpt( ~ py010n , design = des_eusilc )
#' svyarpt( ~ py010n , design = des_eusilc , na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyarpt( ~ py010n , design = des_eusilc_rep )
#' svyarpt( ~ py010n , design = des_eusilc_rep , na.rm = TRUE )
#'
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
#' svyarpt( ~ eqincome , design = dbd_eusilc )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyarpt <-
	function(formula, design, ...) {

		if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

		UseMethod("svyarpt", design)

	}

#' @rdname svyarpt
#' @export
svyarpt.survey.design <-
	function(formula, design, quantiles = 0.5, percent = 0.6,  na.rm = FALSE, deff=FALSE , linearized = FALSE , ...) {

		if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

		# if the class of the full_design attribute is just a TRUE, then the design is
		# already the full design.  otherwise, pull the full_design from that attribute.
		if ("logical" %in% class(attr(design, "full_design")))  full_design <- design else full_design <- attr(design, "full_design")

		# collect income data
		incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

		# treat missing values
		if(na.rm){
			nas<-is.na(incvar)
			design<-design[!nas,]
			if (length(nas) > length(design$prob))
			incvar <- incvar[!nas]
			else incvar[nas] <- 0
		}

		# collect sampling weights
		w <- 1/design$prob

		# compute quantile
		q_alpha <- svyiqalpha( formula, design = design, alpha = quantiles, method = "constant", na.rm = na.rm , linearized = TRUE , ...)
		rval <- percent * coef( q_alpha )[[1]]
		lin <- percent * attr( q_alpha , "linearized" )[,1]

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
		  if (deff == "replace") vsrs <- survey::svyvar( lin , design, na.rm = na.rm) * npop^2/nobs else vsrs <- survey::svyvar( lin , design , na.rm = na.rm ) * npop^2 * (npop - nobs)/(npop * nobs)
		  deff.estimate <- variance/vsrs
		}

		# keep necessary linearized functions
		lin <- lin[ 1/design$prob > 0 ]

		# coerce to matrix
		lin <- matrix( lin , nrow = length( lin ) , dimnames = list( names( lin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

		# build result object
		colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
		class( rval ) <- c( "cvystat" , "svystat" )
		attr( rval , "var" ) <- variance
		attr( rval , "statistic" ) <- "arpt"
		if ( linearized ) attr( rval , "linearized" ) <- lin
		if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
		if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
		rval

	}

#' @rdname svyarpt
#' @export
svyarpt.svyrep.design <-
	function(formula, design, quantiles = 0.5, percent = 0.6, na.rm = FALSE, deff=FALSE , linearized=FALSE , return.replicates = FALSE , ...) {

	  # check for convey_prep
		if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

		# if the class of the full_design attribute is just a TRUE, then the design is
		# already the full design.  otherwise, pull the full_design from that attribute.
		if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")


		# collect data
		df <- model.frame(design)
		incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

		# treat missing values
		if(na.rm){
			nas<-is.na(incvar)
			design<-design[!nas,]
			df <- model.frame(design)
			incvar <- incvar[!nas]
		}

		# collect weights
		w <- weights(design, "sampling")

		# point estimate
		quant_val <- computeQuantiles(incvar, w, p = quantiles)
		quant_val <- as.vector(quant_val)
		estimate <- percent * quant_val

		# collect analysis weights
		ww <- weights(design, "analysis")

		# compute replicates
		qq <- apply(ww, 2, function(wi) 0.6 * computeQuantiles(incvar, wi, p = quantiles))

		# compute variance
		if ( any( is.na( qq ) ) ) variance <- as.matrix( NA ) else {
		  variance <- survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = estimate )
		  this.mean <- attr( variance , "means" )
		  variance <- as.matrix( variance )
		  attr( variance , "means" ) <- this.mean
		}
		colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

		# compute deff
		if ( is.character(deff) || deff | linearized ) {

		  # compute linearized function
		  lin <- percent * CalcQuantile_IF( incvar , w , quantiles )

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
		attr(rval, "statistic") <- "arpt"
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

#' @rdname svyarpt
#' @export
svyarpt.DBIsvydesign <-
	function (formula, design, ...)
	{

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

		NextMethod("svyarpt", design)
	}
