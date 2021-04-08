#' Watts poverty index decomposition (EXPERIMENTAL)
#'
#' Estimate the Watts (1968) poverty measure and its components
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param type_thresh type of poverty threshold. If "abs" the threshold is fixed and given the value
#' of abs_thresh; if "relq" it is given by percent times the quantile; if "relm" it is percent times the mean.
#' @param abs_thresh poverty threshold value if type_thresh is "abs"
#' @param percent the multiple of the the quantile or mean used in the poverty threshold definition
#' @param quantiles the quantile used used in the poverty threshold definition
#' @param thresh return the poverty threshold value
#' @param na.rm Should cases with missing values be dropped?
#' @param ... additional arguments. Currently not used.
#'
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' @return Object of class "\code{cvydstat}", with estimates for the Watts index, FGT(0), Watts Poverty Gap Ratio, and Theil(poor incomes) with a "\code{var}" attribute giving the variance-covariance matrix.
#' A "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Guilherme Jacob, Djalma Pessoa and Anthony Damico
#'
#' @note This function is experimental and is subject to change in later versions.
#'
#' @seealso \code{\link{svywatts},\link{svyfgt},\link{svyfgt}}
#'
#' @references McKinley L. Blackburn (1989). Poverty measurement: an index related to a Theil measure of inequality.
#' \emph{Journal of Business & Economic Statistics}, Vol.7, No.4, pp. 475-481,
#' URL \url{http://amstat.tandfonline.com/doi/abs/10.1080/07350015.1989.10509760}.
#'
#' Satya R. Chakravarty, Joseph Deutsch and Jacques Silber (2008).
#' On the Watts multidimensional poverty index and its decomposition.
#' \emph{World Development}, Vol.36, No.6, pp.1067-1077.
#'
#' Harold W. Watts (1968). An economic definition of poverty.
#' \emph{Institute For Research on Poverty Discussion Papers}, n.5.
#' University of Wisconsin. URL \url{https://www.irp.wisc.edu/publications/dps/pdfs/dp568.pdf}.
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
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' # filter positive incomes
#' des_eusilc <- subset( des_eusilc , eqincome > 0 )
#' des_eusilc_rep <- subset( des_eusilc_rep , eqincome > 0 )
#'
#' # absolute poverty threshold
#' svywattsdec(~eqincome, des_eusilc, abs_thresh=10000)
#' # poverty threshold equal to arpt
#' svywattsdec(~eqincome, des_eusilc, type_thresh= "relq" , thresh = TRUE)
#' # poverty threshold equal to 0.6 times the mean
#' svywattsdec(~eqincome, des_eusilc, type_thresh= "relm" , thresh = TRUE)
#'
#' #  using svrep.design:
#' # absolute poverty threshold
#' svywattsdec(~eqincome, des_eusilc_rep, abs_thresh=10000)
#' # poverty threshold equal to arpt
#' svywattsdec(~eqincome, des_eusilc_rep, type_thresh= "relq" , thresh = TRUE)
#' # poverty threshold equal to 0.6 times the mean
#' svywattsdec(~eqincome, des_eusilc_rep, type_thresh= "relm" , thresh = TRUE)
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
#' dbd_eusilc <- subset( dbd_eusilc , eqincome > 0 )
#'
#' # absolute poverty threshold
#' svywattsdec(~eqincome, dbd_eusilc, abs_thresh=10000)
#' # poverty threshold equal to arpt
#' svywattsdec(~eqincome, dbd_eusilc, type_thresh= "relq" , thresh = TRUE)
#' # poverty threshold equal to 0.6 times the mean
#' svywattsdec(~eqincome, dbd_eusilc, type_thresh= "relm" , thresh = TRUE)
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svywattsdec <-
  function(formula, design, ...) {

    # warning("The svywattsdec function is experimental and is subject to changes in later versions.")

    if( 'type_thresh' %in% names( list( ... ) ) && !( list(...)[["type_thresh"]] %in% c( 'relq' , 'abs' , 'relm' ) ) ) stop( 'type_thresh= must be "relq" "relm" or "abs". see ?svywattsdec for more detail.' )

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svywattsdec", design)

  }

#' @rdname svywattsdec
#' @export
svywattsdec.survey.design <-
  function(formula, design, type_thresh="abs",  abs_thresh=NULL, percent = .60, quantiles = .50, na.rm = FALSE, thresh = FALSE , ...){

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # check for threshold type
    if( type_thresh == "abs" & is.null( abs_thresh ) ) stop( "abs_thresh= must be specified when type_thresh='abs'" )

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # domain
    if( is.null( names( design$prob ) ) ) ind <- as.character( seq( length( design$prob ) ) ) else ind <- names(design$prob)

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # compute poverty measures
    watts <- svywatts( formula=formula , design=design , type_thresh=type_thresh , percent=percent , quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh=thresh , linearized = TRUE )
    fgt0 <- svyfgt( formula=formula , design=design , g = 0 , type_thresh=type_thresh , percent=percent , quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh=thresh , linearized = TRUE )
    fgt1 <- svyfgt( formula=formula , design=design , g = 1 , type_thresh=type_thresh , percent=percent , quantiles=quantiles , abs_thresh=abs_thresh , na.rm=na.rm , thresh=thresh , linearized = TRUE )

    # Watts Poverty Gap Ratio
    fgt0 <- list( value = coef( fgt0 )[[1]] , lin = attr( fgt0 , "linearized" )[,1] )
    fgt1 <- list( value = coef( fgt1 )[[1]] , lin = attr( fgt1 , "linearized" )[,1] )
    W_pgr <- convey::contrastinf( quote( log( fgt0/(fgt0 - fgt1 ) ) ) , list( fgt0 = fgt0 , fgt1 = fgt1 ) )

    # Theil inequality index of incomes among the poor
    # by residual
    watts <- list( value = coef( watts )[[1]] , lin = attr( watts , "linearized" )[,1] )
    L_poor <- convey::contrastinf( quote( watts/fgt0 - W_pgr ) , list( watts = watts , fgt0 = fgt0 , W_pgr = W_pgr ) )
    L_poor$lin <- as.numeric( L_poor$lin )
    names( L_poor$lin ) <- names( watts$lin )

    # create linearization matrix
    lin.matrix <- do.call( cbind , lapply( list( watts , fgt0 , W_pgr , L_poor ) , `[[` , "lin" ) )
    rownames( lin.matrix ) <- names( watts$lin )
    colnames( lin.matrix ) <- c( "watts", "fgt0", "watts pov. gap ratio" , "theil(poor)" )

    # compile
    estimates <- c( watts$value, fgt0$value, W_pgr$value , L_poor$value )
    names( estimates ) <- c( "watts", "fgt0", "watts pov. gap ratio" , "theil(poor)" )

    # collect full weights
    wf <- 1/full_design$prob
    wf[ !( names( wf ) %in% rownames( lin.matrix ) ) ] <- 0

    # ensure length
    if ( nrow( lin.matrix ) != length( full_design$prob ) ) {
      tmplin <- matrix( 0 , nrow = nrow( full_design$variables ) , ncol = ncol( lin.matrix ) )
      tmplin[ wf > 0 , ] <- lin.matrix
      lin.matrix <- tmplin ; rm( tmplin )
      colnames( lin.matrix ) <- c( "watts", "fgt0", "watts pov. gap ratio" , "theil(poor)" )
    }

    # compute variance
    variance <- survey::svyrecvar( sweep( lin.matrix , 1 , wf , "*" ) , full_design$cluster, full_design$strata, full_design$fpc, postStrata = full_design$postStrata )

    # fix nan
    estimates[ is.nan( estimates ) ] <- NA
    variance[ which( is.nan( variance ) ) ] <- NA

    # build result object
    rval <- c( estimates )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "watts decomposition"
    class( rval ) <- c( "cvystat" , "svystat" )
    # attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )[ wf > 0 ]
    rval

  }


#' @rdname svywattsdec
#' @export
svywattsdec.svyrep.design <-
  function(formula, design, type_thresh="abs",  abs_thresh=NULL, percent = .60, quantiles = .50, na.rm = FALSE, thresh = FALSE, ...){

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # check for threshold type
    if( type_thresh == "abs" & is.null( abs_thresh ) ) stop( "abs_thresh= must be specified when type_thresh='abs'" )

    # svyrep design ComputeIndex functions
    ComputeWatts <-
      function( y , w , thresh ) {
        y <- y[ w > 0 ]
        w <- w[ w > 0 ]
        N <- sum(w)
        h <- function( y , thresh ) ifelse( y != 0 , ifelse( y <= thresh , log( thresh / y ) , 0 ) , 0 )
        sum( w * h( y , thresh ) ) / N
      }
    ComputeFGT <-
      function( y , w , g , thresh ) {
        y <- y[ w > 0 ]
        w <- w[ w > 0 ]
        N <- sum(w)
        h <- function(y,thresh,g) ( ( ( thresh - y ) / thresh )^g ) * ( y <= thresh )
        sum( w * h( y , thresh , g ) ) / N
      }

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # domain
    if( is.null( names( design$prob ) ) ) ind <- as.character( seq( length( design$prob ) ) ) else ind <- names(design$prob)

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    ### threshold calculation

    # collect income from full sample
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if (na.rm) {
      nas<-is.na( incvec )
      full_design<-full_design[!nas,]
      incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]
    }

    # collect full sample weights
    wsf <- weights( full_design , "sampling" )

    # poverty threshold
    if ( type_thresh=='relq' ) th <- percent * computeQuantiles( incvec, wsf, p = quantiles)
    if ( type_thresh=='relm' ) th <- percent * sum( incvec*wsf ) / sum( wsf )
    if ( type_thresh=='abs'  ) th <- abs_thresh

    ### domain poverty estimate

    # collect domain data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar)
      design<-design[!nas,]
      incvar <- incvar[!nas]
    }

    # collect domain sampling weights
    ws <- weights(design, "sampling")
    names( ws ) <- rownames( design$variables )

    # point estimates
    watts <- ComputeWatts(incvar, ws, thresh = th )
    fgt0 <- ComputeFGT(incvar, ws, g = 0 , thresh = th )
    fgt1 <- ComputeFGT(incvar, ws, g = 1 , thresh = th )
    w_pgr <- log( fgt0 / ( fgt0 - fgt1 ) )
    L_poor <- CalcGEI( incvar, ifelse( incvar <= th , ws , 0 ) , epsilon = 0 )

    ### variance calculation

    # collect full sample analysis weights
    wwf <- weights( full_design , "analysis" )

    # get replicates
    qq <- apply( wwf, 2, function(wi){

      # compute threshold
      if(type_thresh=='relq') thr <- percent * computeQuantiles( incvec , wi , p = quantiles )
      if(type_thresh=='relm') thr <- percent*sum( incvec*wi )/sum( wi )
      if(type_thresh=='abs')  thr <- abs_thresh

      # mark domains
      wsi <- ifelse( rownames( full_design$variables )[ wsf > 0 ] %in% rownames( design$variables )[ ws > 0 ] , wi , 0 )

      # compute estimates
      watts.rep  <- ComputeWatts( incvec , wsi , thresh = thr )
      fgt0.rep   <- ComputeFGT( incvec , wsi , g = 0 , thresh = thr )
      fgt1.rep   <- ComputeFGT( incvec , wsi , g = 1 , thresh = thr )
      w_pgr.rep  <- log( fgt0.rep / ( fgt0.rep - fgt1.rep ) )
      L_poor.rep <- CalcGEI( incvec, ifelse( incvec <= thr , wsi , 0 ) , epsilon = 0 )

      # combine esitmates
      c( watts.rep, fgt0.rep , w_pgr.rep , L_poor.rep )

    } )
    qq <- t( qq )
    colnames( qq ) <- c( "watts" , "fgt0" , "watts pov. gap ratio" , "theil(poor)" )
    qq[ is.nan( qq ) ] <- NA

    # combine estimates
    estimates <- c( watts, fgt0, w_pgr , L_poor )
    estimates[ is.nan( estimates ) ] <- NA

    # compute variance
    if ( anyNA( qq ) ) {
      variance <- diag( estimates )
      variance[,] <- NA
    } else {
      variance <- survey::svrVar( qq , full_design$scale, full_design$rscales, mse = full_design$mse, coef = estimates )
    }

    # build result object
    rval <- estimates
    names( rval ) <- c( "watts" , "fgt0" , "watts pov. gap ratio" , "theil(poor)" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "watts decomposition"
    class(rval) <- c( "cvystat" , "svrepstat" , "svystat" )
    rval

  }


#' @rdname svywattsdec
#' @export
svywattsdec.DBIsvydesign <-
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

    NextMethod("svywattsdec", design)
  }
