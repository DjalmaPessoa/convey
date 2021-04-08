#' At-risk-of-poverty rate
#'
#' Estimate the proportion of persons with income below the at-risk-of-poverty threshold.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param quantiles income quantile, usually .50 (median)
#' @param percent fraction of the quantile, usually .60
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param return.replicates Return the replicate estimates?
#' @param ... arguments passed on to `svyarpt`
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
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#'
#' svyarpr( ~eqincome , design = des_eusilc )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' svyarpr( ~eqincome , design = des_eusilc_rep )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyarpr( ~ py010n , design = des_eusilc )
#' svyarpr( ~ py010n , design = des_eusilc , na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyarpr( ~ py010n , design = des_eusilc_rep )
#' svyarpr( ~ py010n , design = des_eusilc_rep , na.rm = TRUE )
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
#' svyarpr( ~ eqincome , design = dbd_eusilc )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyarpr <- function(formula, design, ...) {

  if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

  UseMethod("svyarpr", design)

}

#' @rdname svyarpr
#' @export
svyarpr.survey.design <-
  function(formula, design, quantiles = 0.5, percent = 0.6, na.rm=FALSE, deff = FALSE , linearized = FALSE , ...) {

    # test for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # extract values
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing
    if( na.rm ) {
      nas <-is.na( incvar )
      design <-design[ !nas, ]
      if ( length( nas ) > length( design$prob ) ) incvar <- incvar[!nas] else incvar[nas] <- 0
    }

    # collect domain index
    if( is.null( rownames( design$variables ) ) ) ind <- as.character( seq( length( design$prob ) ) ) else ind <- rownames( design$variables )[ 1/design$prob > 0 ]

    # collect (domain) weights
    w <- 1/design$prob

    # domain population
    N <- sum( w )

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # collect incomes
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing
    if( na.rm ){
      nas <- is.na( incvec )
      full_design <- full_design[!nas,]
      if (length(nas) > length(full_design$prob)) incvec <- incvec[!nas] else incvec[nas] <- 0
    }

    # collect full sample index
    if( is.null( names( full_design$prob ) ) ) ncom <- as.character( seq( length( full_design$prob ) ) ) else ncom <- names(full_design$prob)

    # full sample weights
    wf <- 1/full_design$prob

    # calculate h value
    htot <- h_fun( incvec , wf )

    # estimate at risk of poverty threshold
    ARPT <- svyarpt( formula = formula , design = full_design , quantiles = quantiles , percent = percent , na.rm = na.rm , linearized = TRUE , ... )
    arptv <- coef( ARPT )
    arptlin <- attr( ARPT , "linearized" )[,1]

    # ensure length of linearization vector
    if ( length( arptlin ) != length( full_design$prob ) ) {
      names( arptlin ) <- rownames( full_design$variables )[ 1/full_design$prob > 0 ]
      arptlin <- arptlin[ pmatch( rownames( full_design$variables ) , names( arptlin ) ) ]
      names( arptlin ) <- rownames( full_design$variables )
      arptlin[ is.na( arptlin ) ] <- 0
    }

    # value of arpr
    poor <- incvar <= arptv
    rval <- sum( poor * w ) / N

    # first term of linearized variable
    ID <- 1 * ( ncom %in% ind )
    arpr1lin <- ID * ( ( incvec <= arptv ) - rval ) / N

    # calculate Fprime for the domain
    Fprime <- densfun( formula = formula , design = design , arptv[[1]] , h=NULL , FUN = "F" , na.rm = na.rm )

    # calculate linearized variable
    arprlin <- arpr1lin + Fprime * arptlin

    # ensure length
    if ( length( arprlin ) != length( full_design$prob ) ) {
      tmplin <- rep( 0 , nrow( full_design$variables ) )
      tmplin[ w > 0 ] <- arprlin
      arprlin <- tmplin ; rm( tmplin )
      names( arprlin ) <- rownames( full_design$variables )
    }

    # compute variance
    variance <- survey::svyrecvar( arprlin/full_design$prob, full_design$cluster, full_design$strata, full_design$fpc, postStrata = full_design$postStrata )
    variance[ which( is.nan( variance ) ) ] <- NA
    colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # compute deff
    if ( is.character(deff) || deff) {
      nobs <- sum( weights( full_design , "sampling" ) > 0 )
      npop <- sum( weights( full_design , "sampling" ) )
      if (deff == "replace") vsrs <- survey::svyvar( arprlin , full_design, na.rm = na.rm) * npop^2/nobs
      else vsrs <- survey::svyvar( arprlin , full_design , na.rm = na.rm ) * npop^2 * (npop - nobs)/(npop * nobs)
      deff.estimate <- variance/vsrs
    }

    # keep necessary linearized functions
    arprlin <- arprlin[ 1/full_design$prob > 0 ]

    # coerce to matrix
    arprlin <- matrix( arprlin , nrow = length( arprlin ) , dimnames = list( names( arprlin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    # build result object
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "arpr"
    if ( linearized ) attr(rval, "linearized") <- arprlin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( arprlin ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }



#' @rdname svyarpr
#' @export
svyarpr.svyrep.design <-
  function(formula, design, quantiles = 0.5, percent = 0.6,na.rm=FALSE, deff =FALSE, linearized = FALSE , return.replicates = FALSE , ... ) {

    # test for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    ### calculate threshold

    # collect full sample incomes
    incvec <- model.frame( formula , full_design$variables , na.action = na.pass )[[1]]

    # treat missing values
    if( na.rm ){
      nas <- is.na(incvec)
      full_design <- full_design[!nas,]
      incvec <- model.frame( formula , full_design$variables , na.action = na.pass )[[1]]
    }

    # collect full sample sampling weights
    wsf <- weights(full_design,"sampling")

    # create indices
    names( incvec ) <- names( wsf ) <- row.names( full_design$variables )

    # compute threshold
    thresh <- percent * computeQuantiles( incvec , wsf, p = quantiles )

    ### compute at risk of poverty rate

    # collect domain sample income
    incvar <- model.frame( formula, design$variables, na.action = na.pass )[[1]]

    # treat missing values
    if( na.rm ){
      nas <- is.na( incvar )
      design <-design[ !nas, ]
      incvar <- model.frame( formula, design$variables, na.action = na.pass )[[1]]
    }

    # collect domain sample sampling weights
    wsd <- weights( design , "sampling" )

    # compute arpr
    Nd <- sum( wsd )
    rval <- sum( ( incvar <= thresh ) * wsd )/ Nd

    ### variance calculation

    # collect full analysis weights
    wwf <- weights(full_design, "analysis")

    # add domain indices
    ind <- rownames( full_design$variables ) %in% rownames( design$variables )

    # compute replicates
    qq <- apply( wwf, 2, function( wi ) {
      thr <- percent * computeQuantiles( incvec , wi, p = quantiles )
      sum( ( incvec[ ind ] <= thr ) * wi[ ind ] ) / sum( wi[ ind ] )
    } )

    # compute variance
    if ( any( is.na( qq ) ) ) variance <- as.matrix( NA ) else {
      variance <- survey::svrVar( qq , full_design$scale , full_design$rscales , mse = full_design$mse , coef = rval )
      this.mean <- attr( variance , "means" )
      variance <- as.matrix( variance )
      attr( variance , "means" ) <- this.mean
    }
    colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # compute deff
    if ( is.character(deff) || deff || linearized ) {

      # compute linearized threshold
      arptlin <- percent * CalcQuantile_IF( incvec , wsf , quantiles )

      # compute linearized function under fixed threshold
      ID <- 1 * ( rownames( full_design$variables ) %in% rownames( design$variables ) )
      arprlin1 <- ID * ( ( incvec <= thresh ) - rval ) / Nd

      # combines both
      Fprime <- densfun( formula , design = design , thresh , h = NULL , FUN = "F", na.rm = na.rm )
      arprlin <- arprlin1 + arptlin * Fprime

      # compute deff
      nobs <- length( full_design$pweights )
      npop <- sum( full_design$pweights )
      vsrs <- unclass( survey::svyvar( arprlin , full_design, na.rm = na.rm, return.replicates = FALSE, estimate.only = TRUE)) * npop^2/nobs
      if (deff != "replace") vsrs <- vsrs * (npop - nobs)/npop
      deff.estimate <- variance / vsrs

      # filter observation
      names( arprlin ) <- rownames( full_design$variables )

      # coerce to matrix
      arprlin <- matrix( arprlin , nrow = length( arprlin ) , dimnames = list( names( arprlin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    }

    # set up result object
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svrepstat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "arpr"
    if ( linearized ) attr(rval,"linearized") <- arprlin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( arprlin ) )

    # keep replicates
    if (return.replicates) {
      attr( qq , "scale") <- full_design$scale
      attr( qq , "rscales") <- full_design$rscales
      attr( qq , "mse") <- full_design$mse
      rval <- list( mean = rval , replicates = qq )
      class( rval ) <- c( "cvystat" , "svrepstat" )
    }

    # add design effect estimate
    if ( is.character(deff) || deff) attr( rval , "deff" ) <- deff.estimate

    # return object
    rval

  }

#' @rdname svyarpr
#' @export
svyarpr.DBIsvydesign <-
  function (formula, design, ...){

    # full design variables
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

    # domain design variables
    design$variables <-
      getvars(
        formula,
        design$db$connection,
        design$db$tablename,
        updates = design$updates,
        subset = design$subset
      )

    # move to next function
    NextMethod("svyarpr", design)

  }
