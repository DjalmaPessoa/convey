#' Relative median poverty gap
#'
#' Estimate the difference between the at-risk-of-poverty threshold (\code{arpt}) and the median of incomes less than the \code{arpt} relative to the \code{arpt}.
#'
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param quantiles income quantile, usually .5 (median)
#' @param percent fraction of the quantile, usually .60
#' @param na.rm Should cases with missing values be dropped?
#' @param thresh return the poverty poverty threshold
#' @param poor_median return the median income of poor people
#' @param ... future expansion
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Djalma Pessoa and Anthony Damico
#'
#' @seealso \code{\link{svyarpt}}
#'
#' @references Guillaume Osier (2009). Variance estimation for complex indicators
#' of poverty and inequality. \emph{Journal of the European Survey Research
#' Association}, Vol.3, No.3, pp. 167-195,
#' ISSN 1864-3361, URL \url{https://ojs.ub.uni-konstanz.de/srm/article/view/369}.
#'
#' Jean-Claude Deville (1999). Variance estimation for complex statistics and estimators:
#' linearization and residual techniques. Survey Methodology, 25, 193-203,
#' URL \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X19990024882}.
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
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#'
#' svyrmpg( ~eqincome , design = des_eusilc, thresh = TRUE )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' svyrmpg( ~eqincome , design = des_eusilc_rep, thresh = TRUE )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyrmpg( ~ py010n , design = des_eusilc )
#' svyrmpg( ~ py010n , design = des_eusilc , na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyrmpg( ~ py010n , design = des_eusilc_rep )
#' svyrmpg( ~ py010n , design = des_eusilc_rep , na.rm = TRUE )
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
#' svyrmpg( ~ eqincome , design = dbd_eusilc )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyrmpg <-
	function(formula, design, ...) {

		if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

		UseMethod("svyrmpg", design)

	}


#' @rdname svyrmpg
#' @export
svyrmpg.survey.design <-
	function(formula, design, quantiles = 0.5, percent = 0.6, na.rm=FALSE, thresh = FALSE, poor_median = FALSE,...) {

		if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")


		# if the class of the full_design attribute is just a TRUE, then the design is
		# already the full design.  otherwise, pull the full_design from that attribute.
		if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

		incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

		if(na.rm){
			nas<-is.na(incvar)
			design<-design[!nas,]
			if (length(nas) > length(design$prob)) incvar <- incvar[!nas] else incvar[nas] <- 0
		}

		incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

		if(na.rm){
			nas<-is.na(incvec)
			full_design<-full_design[!nas,]
			if (length(nas) > length(full_design$prob)) incvec <- incvec[!nas] else incvec[nas] <- 0
		}

		ARPT <- svyarpt(formula = formula, full_design, quantiles = quantiles, percent = percent, na.rm = na.rm )
		arpt <- coef(ARPT)
		linarpt <- attr(ARPT, "lin")

		POORMED <- svypoormed(formula = formula, design = design, quantiles = quantiles, percent = percent, na.rm = na.rm)
		medp <- coef(POORMED)
		linmedp <- attr(POORMED, "lin")

		MEDP <- list(value = medp, lin = linmedp)
		ARPT <- list(value = arpt, lin = linarpt)
		list_all<- list(ARPT=ARPT, MEDP=MEDP)

		# linearize RMPG
		RMPG <- contrastinf( quote( ( ARPT - MEDP ) / ARPT ) , list_all )
		rval <- RMPG$value
		infun <- unlist(RMPG$lin)

		variance <- survey::svyrecvar(infun/full_design$prob, full_design$cluster,full_design$strata, full_design$fpc,postStrata = full_design$postStrata)

		colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
		class(rval) <- c( "cvystat" , "svystat" )
		attr( rval , "var" ) <- variance
		attr(rval, "lin") <- infun
		attr( rval , "statistic" ) <- "rmpg"
		if(thresh) attr(rval, "thresh") <- arpt
		if(poor_median)attr(rval, "poor_median") <- medp

		rval
	}


#' @rdname svyrmpg
#' @export
svyrmpg.svyrep.design <-
	function(formula, design, quantiles = 0.5, percent = 0.6,na.rm=FALSE, thresh = FALSE, poor_median = FALSE, ...) {

		if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

		# if the class of the full_design attribute is just a TRUE, then the design is
		# already the full design.  otherwise, pull the full_design from that attribute.
		if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

		df <- model.frame(design)
		incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

		if(na.rm){
			nas<-is.na(incvar)
			design<-design[!nas,]
			df <- model.frame(design)
			incvar <- incvar[!nas]
		}

		df_full<- model.frame(full_design)
		incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

		if(na.rm){
			nas<-is.na(incvec)
			full_design<-full_design[!nas,]
			df_full <- model.frame(full_design)
			incvec <- incvec[!nas]
		}

		wsf <- weights(full_design,"sampling")
		names(incvec) <- names(wsf) <- row.names(df_full)
		ind<- row.names(df)

		ComputeRmpg <-
			function(xf, wf, ind, quantiles, percent) {
				tresh <- percent * computeQuantiles(xf, wf, p = quantiles)
				x<-xf[ind]
				w<- wf[ind]
				indpoor <- (x <= tresh)
				medp <- computeQuantiles(x[indpoor], w[indpoor], p = 0.5)
				c(tresh, medp, 1 - (medp/tresh) )
			}

		ws <- weights(design, "sampling")
		Rmpg_val <- ComputeRmpg(xf = incvec, wf=wsf, ind= ind, quantiles = quantiles, percent = percent)
		rval <- Rmpg_val[3]

		wwf <- weights(full_design, "analysis")

		qq <-
			apply(wwf, 2, function(wi){
				names(wi)<- row.names(df_full)
				ComputeRmpg(incvec, wi, ind=ind, quantiles = quantiles,percent = percent)[3]
				})
		if(anyNA(qq))variance <- NA
		else variance <- survey::svrVar(qq, design$scale, design$rscales, mse = design$mse, coef = rval)

		variance <- as.matrix( variance )

		colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
		class(rval) <- c( "cvystat" , "svrepstat" )
		attr( rval , "var" ) <- variance
		attr(rval, "lin") <- NA
		attr( rval , "statistic" ) <- "rmpg"
		if(thresh) attr(rval, "thresh") <- Rmpg_val[1]
		if(poor_median)attr(rval, "poor_median") <- Rmpg_val[2]

		rval
	}

#' @rdname svyrmpg
#' @export
svyrmpg.DBIsvydesign <-
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

		NextMethod("svyrmpg", design)

	}
