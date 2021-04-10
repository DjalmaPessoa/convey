#' FGT indices decomposition (EXPERIMENTAL)
#'
#' Estimate the Foster et al. (1984) poverty class and its components
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param type_thresh type of poverty threshold. If "abs" the threshold is fixed and given the value
#' of abs_thresh; if "relq" it is given by percent times the quantile; if "relm" it is percent times the mean.
#' @param abs_thresh poverty threshold value if type_thresh is "abs"
#' @param g If g=2 estimates the average squared normalised poverty gap. This function is defined for g >= 2 only,
#' @param percent the multiple of the the quantile or mean used in the poverty threshold definition
#' @param quantiles the quantile used used in the poverty threshold definition
#' @param thresh return the poverty threshold value
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param influence Should a matrix of (weighted) influence functions be returned? (for compatibility with \code{\link[survey]{svyby}})
#' @param return.replicates Return the replicate estimates?
#' @param ... additional arguments. Currently not used.
#'
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' @return Object of class "\code{cvydstat}", with estimates for the FGT(g), FGT(0), FGT(1), income gap ratio and GEI(income gaps; epsilon = g) with a "\code{var}" attribute giving the variance-covariance matrix.
#' A "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Guilherme Jacob, Djalma Pessoa and Anthony Damico
#'
#' @note This function is experimental and is subject to change in later versions.
#'
#' @seealso \code{\link{svyfgt},\link{svyfgt},\link{svyfgt}}
#'
#' @references Oihana Aristondo, Cassilda Lasso De La vega and Ana Urrutia (2010).
#' A new multiplicative decomposition for the Foster-Greer-Thorbecke poverty indices.
#' \emph{Bulletin of Economic Research}, Vol.62, No.3, pp. 259-267.
#' University of Wisconsin. URL \url{http://dx.doi.org/10.1111/j.1467-8586.2009.00320.x}.
#'
#' James Foster, Joel Greer and Erik Thorbecke (1984). A class of decomposable poverty measures.
#' \emph{Econometrica}, Vol.52, No.3, pp. 761-766.
#'
#' Guillaume Osier (2009). Variance estimation for complex indicators
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
#'
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' # absolute poverty threshold
#' svyfgtdec(~eqincome, des_eusilc, g=2, abs_thresh=10000)
#' # poverty threshold equal to arpt
#' svyfgtdec(~eqincome, des_eusilc, g=2, type_thresh= "relq" , thresh = TRUE)
#' # poverty threshold equal to 0.6 times the mean
#' svyfgtdec(~eqincome, des_eusilc, g=2, type_thresh= "relm" , thresh = TRUE)
#'
#' # using svrep.design:
#' # absolute poverty threshold
#' svyfgtdec(~eqincome, des_eusilc_rep, g=2, abs_thresh=10000)
#' # poverty threshold equal to arpt
#' svyfgtdec(~eqincome, des_eusilc_rep, g=2, type_thresh= "relq" , thresh = TRUE)
#' # poverty threshold equal to 0.6 times the mean
#' svyfgtdec(~eqincome, des_eusilc_rep, g=2, type_thresh= "relm" , thresh = TRUE)
#'
#' \dontrun{
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
#' # absolute poverty threshold
#' svyfgtdec(~eqincome, dbd_eusilc, g=2, abs_thresh=10000)
#' # poverty threshold equal to arpt
#' svyfgtdec(~eqincome, dbd_eusilc, g=2, type_thresh= "relq" , thresh = TRUE)
#' # poverty threshold equal to 0.6 times the mean
#' svyfgtdec(~eqincome, dbd_eusilc, g=2, type_thresh= "relm" , thresh = TRUE)
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyfgtdec <-
  function(formula, design, ...) {

    warning("The svyfgtdec function is experimental and is subject to changes in later versions.")

    # if( 'type_thresh' %in% names( list( ... ) ) && !( list(...)[["type_thresh"]] %in% c( 'abs' ) ) ) stop( 'type_thresh= must be "abs". See ?svyfgtdec for more detail.' )
    if( 'type_thresh' %in% names( list( ... ) ) && !( list(...)[["type_thresh"]] %in% c( 'relq' , 'abs' , 'relm' ) ) ) stop( 'type_thresh= must be "relq" "relm" or "abs".  see ?svyfgt for more detail.' )

    if( !( 'g' %in% names(list(...)) ) ) stop( "g= parameter must be specified" )
    if( !is.na( list(...)[["g"]] ) && !( list(...)[["g"]] >= 2 ) ) stop( "this decomposition is defined for g >= 2 only." )

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svyfgtdec", design)

  }

#' @rdname svyfgtdec
#' @export
svyfgtdec.survey.design <-
  function(formula, design, g, type_thresh="abs",  abs_thresh=NULL, percent = .60, quantiles = .50, na.rm = FALSE, thresh = FALSE, deff = FALSE , linearized = FALSE , influence = FALSE , ...){

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # check for threshold type
    if( type_thresh == "abs" & is.null( abs_thresh ) ) stop( "abs_thresh= must be specified when type_thresh='abs'" )

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # compute fgt estimates
    fgt0 <- svyfgt( formula = formula, design=design, g=0, type_thresh=type_thresh, percent=percent, quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh = thresh , linearized = TRUE )
    fgt1 <- svyfgt( formula = formula, design=design, g=1, type_thresh=type_thresh, percent=percent, quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh = thresh , linearized = TRUE )
    fgtg <- svyfgt( formula = formula, design=design, g=g, type_thresh=type_thresh, percent=percent, quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh = thresh , linearized = TRUE )
    if (thresh) th <- attr( fgt0 , "thresh" )

    # keep threshold
    if ( thresh ) thresh.value <- attr( fgt0 , "thresh" )

    # income gap ratio
    fgt0 <- list( value = fgt0[[1]], lin = attr( fgt0 , "linearized" ) )
    fgt1 <- list( value = fgt1[[1]], lin = attr( fgt1 , "linearized" ) )
    igr <- contrastinf( quote( fgt1 / fgt0 ) , list( fgt0 = fgt0 , fgt1 = fgt1) )

    # generalized entropy index of poverty gaps
    # by residual
    fgtg <- list( value = fgtg[[1]], lin = attr( fgtg , "linearized" ) )
    gei_poor <- contrastinf( quote( ( fgtg / ( fgt0 * igr^g ) - 1 ) / ( g^2 - g ) ) , list( fgtg =fgtg , fgt0 = fgt0 , fgt1 = fgt1 , igr = igr , g = list( value = g , lin = rep( 0 , length( igr$lin ) ) ) ) )

    # matrix of linearized variables
    lin.matrix <- cbind( fgtg$lin, fgt0$lin, fgt1$lin , igr$lin , gei_poor$lin)
    lin.matrix <- as.matrix( lin.matrix )
    colnames( lin.matrix ) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )

    # ensure length
    if ( nrow( lin.matrix ) != length( full_design$prob ) ) {
      w <- 1/full_design$prob
      tmplin <- matrix( 0 , nrow = nrow( full_design$variables ) , ncol = ncol( lin.matrix ) )
      tmplin[ w > 0 , ] <- lin.matrix
      lin.matrix <- tmplin ; rm( tmplin )
      colnames( lin.matrix ) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )
    }

    # format objects
    estimates <- c( fgtg$value, fgt0$value, fgt1$value , igr$value , gei_poor$value )
    variance <- survey::svyrecvar( lin.matrix/full_design$prob , full_design$cluster, full_design$strata, full_design$fpc, postStrata = full_design$postStrata )

    # compute deff
    if ( is.character(deff) || deff ) {
      nobs <- sum( weights( full_design ) != 0 )
      npop <- sum( weights( full_design ) )
      if (deff == "replace") vsrs <- survey::svyvar( lin.matrix , full_design, na.rm = na.rm , estimate.only = TRUE ) * npop^2/nobs
      else vsrs <- survey::svyvar( lin.matrix , full_design , na.rm = na.rm , estimate.only = TRUE ) * npop^2 * (npop - nobs)/(npop * nobs)
      deff.estimate <- variance/vsrs
    }

    # add names
    rownames( lin.matrix ) <- rownames( full_design$variables )

    # build result object
    rval <- c( estimates )
    names( rval ) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- paste0("fgt", g , " decomposition")
    if (thresh) attr(rval, "thresh") <- th
    if ( linearized ) attr(rval,"linearized") <- lin.matrix
    if ( influence )  attr( rval , "influence" )  <- sweep( lin.matrix , 1 , full_design$prob[ is.finite( full_design$prob ) ] , "/" )
    if ( linearized | influence ) attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )
    if ( is.character(deff) || deff ) attr( rval , "deff") <- deff.estimate
    class(rval) <- c( "cvystat" , "svrepstat" , "svystat" )
    rval

  }


#' @rdname svyfgtdec
#' @export
svyfgtdec.svyrep.design <-
  function(formula, design, g, type_thresh="abs",  abs_thresh=NULL, percent = .60, quantiles = .50, na.rm = FALSE, thresh = FALSE, deff = FALSE , linearized = FALSE , return.replicates = FALSE , ...){

    # svyrep design h function
    h <- function(y,thresh,g) ( ( ( thresh - y ) / thresh )^g ) * ( y <= thresh )

    # svyrep design ComputeFGT function
    ComputeFGT <-
      function(y, w, thresh, g){
        y <- y[w>0]
        w <- w[w>0]
        N <- sum(w)
        sum( w * h( y , thresh , g ) ) / N
      }

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

    # test for threshold type
    if( type_thresh == "abs" & is.null( abs_thresh ) ) stop( "abs_thresh= must be specified when type_thresh='abs'" )

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # collect domain sample data
    df <- model.frame(design)
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar)
      design<-design[!nas,]
      df <- model.frame(design)
      incvar <- incvar[!nas]
    }

    # colelct somain weights
    ws <- weights(design, "sampling")

    # collect full sample data
    df_full<- model.frame(full_design)
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvec)
      full_design<-full_design[!nas,]
      df_full <- model.frame(full_design)
      incvec <- incvec[!nas]
    }

    # colelct full sample weights
    wsf <- weights(full_design,"sampling")
    names(incvec) <- names(wsf) <- row.names(df_full)

    # create domain indices
    ind <- row.names( df )

    # poverty threshold
    th <- NULL
    if(type_thresh=='relq') th <- percent * computeQuantiles( incvec, wsf, p = quantiles)
    if(type_thresh=='relm') th <- percent*sum(incvec*wsf)/sum(wsf)
    if(type_thresh=='abs') th <- abs_thresh

    # # estimates
    # fgt0 <- ComputeFGT(incvar, ws, g = 0 , thresh = th )
    # fgt1 <- ComputeFGT(incvar, ws, g = 1 , thresh = th )
    # fgtg <- ComputeFGT(incvar, ws, g = g , thresh = th )
    # igr <- fgt1/fgt0
    # gei_poor <- ( ( fgtg / ( fgt0 * igr^g ) - 1 ) / ( g^2 - g ) )
    # estimates <- c( fgtg, fgt0, fgt1, igr, gei_poor )

    # compute fgt estimates
    fgt0 <- svyfgt( formula = formula, design=design, g=0, type_thresh=type_thresh, percent=percent, quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh = thresh , linearized = TRUE )
    fgt1 <- svyfgt( formula = formula, design=design, g=1, type_thresh=type_thresh, percent=percent, quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh = thresh , linearized = TRUE )
    fgtg <- svyfgt( formula = formula, design=design, g=g, type_thresh=type_thresh, percent=percent, quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh = thresh , linearized = TRUE )

    # keep threshold
    if ( thresh ) thresh.value <- attr( fgt0 , "thresh" )

    # income gap ratio
    fgt0 <- list( value = fgt0[[1]], lin = attr( fgt0 , "linearized" ) )
    fgt1 <- list( value = fgt1[[1]], lin = attr( fgt1 , "linearized" ) )
    igr <- contrastinf( quote( fgt1 / fgt0 ) , list( fgt0 = fgt0 , fgt1 = fgt1) )

    # generalized entropy index of poverty gaps
    # by residual
    fgtg <- list( value = fgtg[[1]], lin = attr( fgtg , "linearized" ) )
    gei_poor <- contrastinf( quote( ( fgtg / ( fgt0 * igr^g ) - 1 ) / ( g^2 - g ) ) , list( fgtg =fgtg , fgt0 = fgt0 , fgt1 = fgt1 , igr = igr , g = list( value = g , lin = rep( 0 , length( igr$lin ) ) ) ) )

    # combine values and linearization
    estimates <- c( fgtg$value, fgt0$value, fgt1$value , igr$value , gei_poor$value )

    # matrix of linearized variables
    lin.matrix <- cbind( fgtg$lin, fgt0$lin, fgt1$lin , igr$lin , gei_poor$lin)
    lin.matrix <- as.matrix( lin.matrix )
    colnames( lin.matrix ) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )
    if ( nrow( lin.matrix ) != length( full_design$pweights ) ) {
      rownames( lin.matrix ) <- rownames( full_design$variables )[ full_design$pweights > 0 ]
      lin.matrix <- lin.matrix[ pmatch( rownames( full_design$variables ) , rownames(lin.matrix ) ) , ]
      lin.matrix[ full_design$pweights <= 0 , ] <- 0
    }
    lin.matrix <- lin.matrix[ full_design$pweights > 0 , ]

    # colelct full sample analysis weights
    ww <- weights( full_design, "analysis" )

    # get replicates
    qq <- apply( ww, 2, function( wi ) {
      if(type_thresh=='relq') thr <- percent * computeQuantiles( incvec , wi , p = quantiles )
      if(type_thresh=='relm') thr <- percent*sum( incvec*wi )/sum( wi )
      if(type_thresh=='abs')  thr <- abs_thresh
      c( fgtg = ComputeFGT( incvec , ifelse( rownames( df_full ) %in% ind , wi , 0 ) , g = g , thresh = thr ) ,
         fgt0 = ComputeFGT( incvec , ifelse( rownames( df_full ) %in% ind , wi , 0 ) , g = 0 , thresh = thr ) ,
         fgt1 = ComputeFGT( incvec , ifelse( rownames( df_full ) %in% ind , wi , 0 ) , g = 1 , thresh = thr ) )
    } )
    qq <- t( qq )
    qq <- cbind( qq , "igr" = qq[ , "fgt1" ] / qq[ , "fgt0" ] )
    qq <- cbind( qq , "gei_poor" = ( ( qq[ , "fgtg" ] / ( qq[ , "fgt0" ] * qq[ , "igr" ]^g ) - 1 ) / ( g^2 - g ) ) )
    colnames(qq) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )

    # compute variance
    if ( anyNA( qq ) ) {
      variance <- matrix( NA , ncol = 5 , nrow = 5  )
    } else {
      variance <- survey::svrVar(qq, full_design$scale, full_design$rscales, mse = full_design$mse, coef = estimates )
    }
    rownames( variance ) <- colnames( variance ) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )

    # compute deff
    if ( is.character(deff) || deff ) {
      nobs <- length( full_design$pweights )
      npop <- sum( full_design$pweights )
      vsrs <- unclass( survey::svyvar( lin.matrix , full_design, na.rm = na.rm, return.replicates = FALSE, estimate.only = TRUE)) * npop^2/nobs
      if (deff != "replace") vsrs <- vsrs * (npop - nobs)/npop
      deff.estimate <- variance / vsrs
    }

    # build result object
    rval <- c( estimates )
    names( rval ) <- c( paste0("fgt",g), "fgt0", "fgt1" , "igr" , paste0( "gei(poor;epsilon=",g,")" ) )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- paste0("fgt", g , " decomposition")
    if (thresh) attr(rval, "thresh") <- th
    class(rval) <- c( "cvystat" , "svrepstat" , "svystat" )
    if ( linearized ) attr( rval , "linearized") <- lin.matrix
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )

    # keep replicates
    if (return.replicates) {
      attr( qq , "scale") <- design$scale
      attr( qq , "rscales") <- design$rscales
      attr( qq , "mse") <- design$mse
      rval <- list( mean = rval , replicates = qq )
      class(rval) <- c( "cvystat" , "svrepstat" , "svystat" )
    }

    # add design effect estimate
    if ( is.character(deff) || deff ) attr( rval , "deff") <- deff.estimate

    # retorna objeto
    rval

  }


#' @rdname svyfgtdec
#' @export
svyfgtdec.DBIsvydesign <-
  function (formula, design, ...) {

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

    NextMethod("svyfgtdec", design)
  }

