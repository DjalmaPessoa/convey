#' Quintile Share Ratio
#'
#' Estimate ratio of the total income received by the highest earners to the total income received by lowest earners, defaulting to 20%.
#'
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param alpha1 order of the lower quintile
#' @param alpha2 order of the upper quintile
#' @param na.rm Should cases with missing values be dropped?
#' @param upper_quant return the lower bound of highest earners
#' @param lower_quant return the upper bound of lowest earners
#' @param upper_tot return the highest earners total
#' @param lower_tot return the lowest earners total
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
#' svyqsr( ~eqincome , design = des_eusilc, upper_tot = TRUE, lower_tot = TRUE )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' svyqsr( ~eqincome , design = des_eusilc_rep, upper_tot = TRUE, lower_tot = TRUE )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyqsr( ~ db090 , design = des_eusilc )
#' svyqsr( ~ db090 , design = des_eusilc , na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyqsr( ~ db090 , design = des_eusilc_rep )
#' svyqsr( ~ db090 , design = des_eusilc_rep , na.rm = TRUE )
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
#' svyqsr( ~ eqincome , design = dbd_eusilc )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyqsr <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    if( 'alpha' %in% names( list(...) ) && list(...)[["alpha"]] > 0.5 ) stop( "alpha= cannot be larger than 0.5 (50%)" )

    UseMethod("svyqsr", design)

  }

#' @rdname svyqsr
#' @export
svyqsr.survey.design <-
  function(formula, design, alpha1 = 0.2 , alpha2 = ( 1 - alpha1 ) , na.rm=FALSE, upper_quant = FALSE, lower_quant = FALSE, upper_tot = FALSE, lower_tot = FALSE, deff = FALSE , linearized = FALSE , ...) {

    # test for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # colelct income data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar)
      design<-design[!nas,]
      if (length(nas) > length(design$prob)) incvar <- incvar[!nas] else incvar[nas] <- 0
    }

    # collect sampling weights
    w <- 1/design$prob

    # collect domain indices
    ind <- names(design$prob)

    # Linearization of S20
    S20 <- svyisq( formula = formula, design = design, alpha1, na.rm=na.rm, quantile = TRUE , deff = FALSE , linearized = TRUE )
    qS20 <- attr( S20, "quantile")
    totS20 <- coef( S20 )
    attributes( totS20 ) <- NULL
    S20 <- list(value= totS20[[1]], lin=attr(S20,"linearized")[,1])

    # treat missing
    if ( is.na( totS20 ) ) {
      rval <- as.numeric( NA )
      variance <- as.matrix( NA )
      colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      class(rval) <- c( "cvystat" , "svystat" )
      attr(rval, "var") <- variance
      attr(rval, "statistic") <- "qsr"
      return( rval )
    }

    # test division by zero
    if( S20$value == 0 ) stop( paste0( "division by zero. the alpha1=" , alpha1 , " percentile cannot be zero or svyqsr would return Inf" ) )

    # Linearization of S80C
    S80C <- svyisq( formula = formula, design = design, alpha2 , na.rm=na.rm , quantile = TRUE , upper = TRUE , deff = FALSE , linearized = TRUE )
    qS80C <- attr(S80C, "quantile")
    totS80C <- coef(S80C)
    attributes(totS80C) <- NULL
    S80C <- list(value=totS80C[[1]], lin=attr(S80C,"linearized")[,1])

    # ensure consistent lengths
    if ( length( unique( sapply( lapply( list( S80C , S20 ) , `[[` , "lin" ) , length ) ) ) != 1 ) stop()

    # LINEARIZED VARIABLE OF THE SHARE RATIO
    list_all <- list( S20 = S20 , S80C = S80C )
    QSR <- contrastinf( quote(S80C/S20), list_all)
    lin <- as.numeric(QSR$lin)
    names(lin) <- names(w)[w > 0]

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
    names( lin ) <- rownames( design$variables )[ w > 0 ]

    # coerce to matrix
    lin <- matrix( lin , nrow = length( lin ) , dimnames = list( names( lin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    # build result object
    rval <- as.numeric( QSR$value )
    attributes( rval ) <- NULL
    names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "qsr"
    if(upper_quant)  attr(rval, "upper_quant") <- qS80C
    if(lower_quant)  attr(rval, "lower_quant") <- qS20
    if(upper_tot)  attr(rval, "upper_tot") <- totS80C
    if(lower_tot)  attr(rval, "lower_tot") <- totS20
    if ( is.character(deff) || deff ) attr(rval,"deff") <- deff.estimate
    if ( linearized ) attr(rval, "linearized") <- lin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
    rval

  }

#' @rdname svyqsr
#' @export
svyqsr.svyrep.design <-
  function(formula, design, alpha1 = 0.2 , alpha2 = ( 1 - alpha1 ) , na.rm=FALSE, upper_quant = FALSE, lower_quant = FALSE, upper_tot = FALSE, lower_tot = FALSE, deff = FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

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

    # computation function
    ComputeQsr <-
      function(x, w, alpha1, alpha2) {
        quant_inf <- computeQuantiles(x, w, p = alpha1)
        quant_sup <- computeQuantiles(x, w, p = alpha2)
        rich <- (x > quant_sup) * x
        S80 <- sum(rich * w)
        poor <- (x <= quant_inf) * x
        S20 <- sum(poor * w)
        c( quant_sup, quant_inf, S80, S20, S80/S20)
      }

    # collect sampling weights
    ws <- weights(design, "sampling")

    # compute point estimate
    Qsr_val <- ComputeQsr(incvar, ws, alpha1 = alpha1, alpha2= alpha2)

    # treat missing
    if ( is.na( Qsr_val[[4]] ) ) {
      rval <- as.numeric( NA )
      variance <- as.matrix( NA )
      colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      class(rval) <- c( "cvystat" , "svystat" )
      attr(rval, "var") <- variance
      attr(rval, "statistic") <- "qsr"
      # attr(rval, "linearized") <- lin
      return( rval )
    }

    # test for division by zero
    if( Qsr_val[4] == 0 ) stop( paste0( "division by zero. the alpha1=" , alpha1 , " percentile cannot be zero or svyqsr would return Inf" ) )

    ### variance calculation

    # collect analysis weights
    ww <- weights(design, "analysis")

    # compute replicates
    qq <- apply( ww , 2 , function( wi ) ComputeQsr( incvar , w = wi , alpha1 = alpha1 , alpha2 = alpha2 )[5] )

    # compute variance
    if( anyNA( qq ) ) variance <- NA else variance <- survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = Qsr_val[5] )
    variance <- as.matrix( variance )
    colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # compute deff
    if ( is.character(deff) || deff || linearized ) {

      # Linearization of S20
      S20 <- svyisq(formula = formula, design = design, alpha1, na.rm=na.rm, quantile = TRUE , deff = FALSE , linearized = TRUE )
      qS20 <- attr(S20, "quantile")
      totS20 <- coef(S20)
      attributes(totS20) <- NULL
      S20 <- list( value= totS20[[1]], lin=attr(S20,"linearized")[,1] )

      # Linearization of S80C
      S80C <- svyisq( formula = formula, design = design, alpha2 , na.rm=na.rm , quantile = TRUE , upper = TRUE , deff = FALSE , linearized = TRUE )
      qS80C <- attr(S80C, "quantile")
      totS80C <- coef(S80C)
      attributes(totS80C) <- NULL
      S80C <- list( value=totS80C[[1]], lin=attr(S80C,"linearized")[,1] )


      # linearizatiion of the ratio
      list_all <- list( S20 = S20 , S80C = S80C )
      QSR <- contrastinf( quote( S80C/S20 ) , list_all )
      lin <- as.numeric( QSR$lin[,1] )
      names( lin ) <- rownames( design$variables )

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
    rval <- Qsr_val[[5]]
    attributes( rval ) <- NULL
    names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "qsr"
    class(rval) <- c( "cvystat" , "svrepstat" )
    if(upper_quant)  attr(rval, "upper_quant") <- Qsr_val[1]
    if(lower_quant)  attr(rval, "lower_quant") <- Qsr_val[2]
    if(upper_tot)  attr(rval, "upper_tot") <- Qsr_val[3]
    if(lower_tot)  attr(rval, "lower_tot") <- Qsr_val[4]
    if ( linearized ) attr(rval,"linearized") <- lin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )

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

#' @rdname svyqsr
#' @export
svyqsr.DBIsvydesign <-
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

    NextMethod("svyqsr", design)
  }
