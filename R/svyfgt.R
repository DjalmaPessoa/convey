#' FGT measure of poverty
#'
#' Estimate the FGT measure.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param type_thresh type of poverty threshold. If "abs" the threshold is fixed and given the value
#' of abs_thresh; if "relq" it is given by percent times the quantile; if "relm" it is percent times the mean.
#' @param abs_thresh poverty threshold value if type_thresh is "abs"
#' @param g If g=0 estimates the headcount ratio; If g=1 estimates the average normalised poverty gap, and if g=2 estimates the average squared normalised poverty gap
#' @param percent the multiple of the the quantile or mean used in the poverty threshold definition
#' @param quantiles the quantile used used in the poverty threshold definition
#' @param thresh return the poverty threshold value
#' @param na.rm Should cases with missing values be dropped?
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned?
#' @param return.replicates Return the replicate estimates?
#' @param ... passed to \code{svyarpr} and \code{svyarpt}
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' The FGT poverty measures have three special cases.
#' When \code{g = 0}, the FGT measure is the headcount poverty rate, assigning the same "poverty-weight" to all persons below the poverty line.
#' When \code{g = 1}, it becomes the poverty gap ratio, a measure which accounts for the intensity of income shortfall among the poor.
#' When \code{g = 2}. it becomes the squared poverty gap ratio, a measure that also accounts for inequality of poverty intesity across the poor.
#' The \code{g} is a poverty sensitivity parameter, adding more weight to people with greater income shortfalls as it increases.
#'
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Djalma Pessoa, Anthony Damico and Guilherme Jacob
#'
#' @seealso \code{\link{svyarpt}}
#'
#' @references James Foster, Joel Greer and Erik Thorbecke (1984). A class of decomposable poverty measures.
#' \emph{Econometrica}, Vol.52, No.3, pp. 761-766.
#'
#' Y.G. Berger and C. J. Skinner (2003), Variance estimation for a low income proportion.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, Vol. 52, No. 4, pp. 457-468.
#' URL \url{https://doi.org/10.1111/1467-9876.00417}
#'
#' Buhong Zheng (2001). Statistical inference for poverty measures with relative poverty lines.
#' \emph{Journal of Econometrics}, Vol. 101, pp. 337-356.
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
#' # headcount ratio, poverty threshold fixed
#' svyfgt(~eqincome, des_eusilc, g=0,  abs_thresh=10000)
#' # poverty gap index, poverty threshold fixed
#' svyfgt(~eqincome, des_eusilc, g=1,  abs_thresh=10000)
#' # headcount ratio, poverty threshold equal to arpt
#' svyfgt(~eqincome, des_eusilc, g=0, type_thresh= "relq" , thresh = TRUE)
#' # poverty gap index, poverty threshold equal to arpt
#' svyfgt(~eqincome, des_eusilc, g=1, type_thresh= "relq", thresh = TRUE)
#' # headcount ratio, poverty threshold equal to .6 times the mean
#' svyfgt(~eqincome, des_eusilc, g=0, type_thresh= "relm", thresh = TRUE)
#' # poverty gap index, poverty threshold equal to 0.6 times the mean
#' svyfgt(~eqincome, des_eusilc, g=1, type_thresh= "relm" , thresh = TRUE)
#'
#' #  using svrep.design:
#' # headcount ratio, poverty threshold fixed
#' svyfgt(~eqincome, des_eusilc_rep, g=0,  abs_thresh=10000)
#' # poverty gap index, poverty threshold fixed
#' svyfgt(~eqincome, des_eusilc, g=1,  abs_thresh=10000)
#' # headcount ratio, poverty threshold equal to arpt
#' svyfgt(~eqincome, des_eusilc_rep, g=0, type_thresh= "relq" , thresh = TRUE)
#' # poverty gap index, poverty threshold equal to arpt
#' svyfgt(~eqincome, des_eusilc, g=1, type_thresh= "relq", thresh = TRUE)
#' # headcount ratio, poverty threshold equal to .6 times the mean
#' svyfgt(~eqincome, des_eusilc_rep, g=0, type_thresh= "relm" , thresh = TRUE)
#' # poverty gap index, poverty threshold equal to 0.6 times the mean
#' svyfgt(~eqincome, des_eusilc_rep, g=1, type_thresh= "relm", thresh = TRUE)
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
#'
#' dbd_eusilc <- convey_prep( dbd_eusilc )
#'
#' # headcount ratio, poverty threshold fixed
#' svyfgt(~eqincome, dbd_eusilc, g=0, abs_thresh=10000)
#' # poverty gap index, poverty threshold fixed
#' svyfgt(~eqincome, dbd_eusilc, g=1, abs_thresh=10000)
#' # headcount ratio, poverty threshold equal to arpt
#' svyfgt(~eqincome, dbd_eusilc, g=0, type_thresh= "relq", thresh = TRUE)
#' # poverty gap index, poverty threshold equal to arpt
#' svyfgt(~eqincome, dbd_eusilc, g=1, type_thresh= "relq")
#' # headcount ratio, poverty threshold equal to .6 times the mean
#' svyfgt(~eqincome, dbd_eusilc, g=0, type_thresh= "relm")
#' # poverty gap index, poverty threshold equal to 0.6 times the mean
#' svyfgt(~eqincome, dbd_eusilc, g=1, type_thresh= "relm")
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyfgt <-
  function(formula, design, ... ) {

    if( !( 'g' %in% names(list(...)) ) ) stop( "g= parameter must be specified" )

    if( !is.na( list(...)[["g"]] ) && !( ( list(...)[["g"]] >= 0 ) & ( list(...)[["g"]] %% 1 == 0 ) ) ) stop( "g= must be an integer greater or equal to 0" )

    if( 'type_thresh' %in% names( list( ... ) ) && !( list(...)[["type_thresh"]] %in% c( 'relq' , 'abs' , 'relm' ) ) ) stop( 'type_thresh= must be "relq" "relm" or "abs".  see ?svyfgt for more detail.' )

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svyfgt", design)

  }

#' @rdname svyfgt
#' @export
svyfgt.survey.design <-
  function(formula, design, g, type_thresh="abs",  abs_thresh=NULL, percent = .60, quantiles = .50, na.rm = FALSE, thresh = FALSE, deff = FALSE , linearized = FALSE , ...) {

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # check for threshold type
    if( type_thresh == "abs" & is.null( abs_thresh ) ) stop( "abs_thresh= must be specified when type_thresh='abs'" )

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # h function
    h <- function( y , thresh , g ) ( ( ( thresh - y ) / thresh )^g ) * ( y <= thresh )

    # ht function
    ht <- function( y , thresh , g ) ( g * ( ( ( thresh - y ) / thresh )^( g - 1 ) ) * ( y / ( thresh^2 ) ) ) * ( y <= thresh )

    ### threshold calculation

    # collect income from full sample
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na( incvec )
      full_design<-full_design[!nas,]
      if ( length(nas) > length( full_design$prob ) ) incvec <- incvec[!nas] else incvec[nas] <- 0
    }

    # collect full sample weights
    wf <- 1/full_design$prob

    # branch on threshold type and its linearized function
    if( type_thresh == 'relq' ) {

      ARPT <- svyiqalpha(formula = formula, full_design, alpha=quantiles,  na.rm=na.rm, linearized = TRUE , ...)
      th <- percent * coef(ARPT)[[1]]
      arptlin <- percent * attr( ARPT , "linearized" )[,1]

    } else if( type_thresh == 'relm') {

      Yf <- sum( wf * incvec )
      Nf <- sum( wf )
      muf <- Yf/Nf
      th <- percent * muf
      arptlin <- percent * ( incvec - muf ) / Nf
      arptlin <- arptlin[ wf > 0 ]
      names( arptlin ) <- rownames( full_design$variables )[ wf > 0 ]

    } else if ( type_thresh == 'abs' ) {
      th <- abs_thresh
      arptlin <- rep( 0 , sum( wf > 0 ) )
    }
    names( arptlin ) <- rownames( full_design$variables )[ wf > 0 ]

    ### domain poverty measure estimate

    # collect income from domain sample
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if ( na.rm ){
      nas <- is.na( incvar )
      design <- design[!nas,]
      if ( length(nas) > length( design$prob ) ) incvar <- incvar[!nas] else incvar[nas] <- 0
    }

    # collect full sample weights
    w <- 1/design$prob

    # domain population size
    Nd <- sum( w )

    # compute value
    estimate <- sum( w * h( incvar , th , g ) ) / Nd

    ### linearization

    # add linearized threshold
    ID <- 1 *( rownames( full_design$variables ) %in% rownames( design$variables )[1/design$prob > 0 ] )
    fgtlin1 <- ID * ( h( incvec , th , g ) - estimate ) / Nd
    fgtlin1 <- fgtlin1[ wf > 0 ]
    Fprime <- densfun( formula , design = design , th , h = NULL , FUN = "F", na.rm = na.rm )
    ahat <- sum( w * ht( incvar , th , g ) ) / Nd
    ahF <- ahat + h( th , th , g ) * Fprime
    if ( type_thresh %in% c( "relq" , "relm" ) ) fgtlin2 <- ahF * arptlin else fgtlin2 <- 0
    fgtlin <- fgtlin1 + fgtlin2
    names( fgtlin ) <- rownames( full_design$variables )[ wf > 0 ]

    # ensure length
    if ( length( fgtlin ) != length( full_design$prob ) ) {
      tmplin <- rep( 0 , nrow( full_design$variables ) )
      tmplin[ wf > 0 ] <- fgtlin
      fgtlin <- tmplin ; rm( tmplin )
      names( fgtlin ) <- rownames( full_design$variables )
    }

    # compute variance
    variance <- survey::svyrecvar( fgtlin/full_design$prob, full_design$cluster, full_design$strata, full_design$fpc, postStrata = full_design$postStrata )
    variance[ which( is.nan( variance ) ) ] <- NA
    colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # compute deff
    if ( is.character(deff) || deff) {
      nobs <- sum( weights( full_design ) != 0 )
      npop <- sum( weights( full_design ) )
      if (deff == "replace") vsrs <- survey::svyvar( fgtlin , full_design, na.rm = na.rm) * npop^2/nobs
      else vsrs <- survey::svyvar( fgtlin , full_design , na.rm = na.rm ) * npop^2 * (npop - nobs)/(npop * nobs)
      deff.estimate <- variance/vsrs
    }

    # keep necessary linearized functions
    fgtlin <- fgtlin[ 1/full_design$prob > 0 ]

    # coerce to matrix
    fgtlin <- matrix( fgtlin , nrow = length( fgtlin ) , dimnames = list( names( fgtlin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    # setup result object
    rval <- estimate
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- paste0("fgt",g)
    if(thresh) attr(rval, "thresh") <- th
    if ( linearized ) attr( rval , "linearized") <- fgtlin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( fgtlin ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }



#' @rdname svyfgt
#' @export
svyfgt.svyrep.design <-
  function(formula, design, g, type_thresh="abs", abs_thresh=NULL, percent = .60, quantiles = .50, na.rm = FALSE, thresh = FALSE, deff=FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

    # checked for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

    # check for threshold type
    if( type_thresh == "abs" & is.null( abs_thresh ) ) stop( "abs_thresh= must be specified when type_thresh='abs'" )

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # h function
    h <- function( y , thresh , g ) ( ( ( thresh - y ) / thresh )^g ) * ( y <= thresh )

    # ht function
    ht <- function( y , thresh , g ) ( g * ( ( ( thresh - y ) / thresh )^( g - 1 ) ) * ( y / ( thresh^2 ) ) ) * ( y <= thresh )

    # svyrep design ComputeFGT function
    ComputeFGT <-
      function(y, w, thresh, g){
        y <- y[w>0]
        w <- w[w>0]
        N <- sum(w)
        sum( w * h( y , thresh , g ) ) / N
      }

    # collect data from full sample
    df_full<- model.frame(full_design)
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvec)
      full_design<-full_design[!nas,]
      df_full <- model.frame(full_design)
      incvec <- incvec[!nas]
    }

    # collect weights from full sample
    wsf <- weights(full_design,"sampling")
    names(incvec) <- names(wsf) <- row.names(df_full)
    ind <- row.names(df)

    # poverty threshold
    if(type_thresh=='relq') th <- percent * computeQuantiles( incvec, wsf, p = quantiles)
    if(type_thresh=='relm') th <- percent*sum(incvec*wsf)/sum(wsf)
    if(type_thresh=='abs') th <- abs_thresh

    ### domain poverty measure

    # collect domain data
    df <- model.frame(design)
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar)
      design<-design[!nas,]
      df <- model.frame(design)
      incvar <- incvar[!nas]
    }

    # collect domain sampling weights
    ws <- weights(design, "sampling")
    names(ws) <- rownames(design$variables)

    # compute point estimate
    rval <- ComputeFGT(incvar, ws, g = g, th)

    ### variance calculation

    # collect full sample analysis weights
    wwf <- weights(full_design, "analysis")

    # calculate replicates
    qq <- apply( wwf, 2, function(wi) {

      names( wi ) <- row.names( df_full )
      if(type_thresh=='relq') thr <- percent * computeQuantiles( incvec , wi , p = quantiles )
      if(type_thresh=='relm') thr <- percent*sum( incvec*wi )/sum( wi )
      if(type_thresh=='abs')  thr <- abs_thresh
      wsi <- ifelse( names( wi ) %in% names( ws ) , wi , 0 )

      ComputeFGT( incvec , wsi , g = g , thr )

    } )

    # treat missing
    if ( anyNA( qq ) ) {
      variance <-  as.matrix( NA )
      names( rval ) <- names( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      class(rval) <- c( "cvystat" , "svrepstat" )
      attr( rval , "var" ) <- variance
      attr(rval, "statistic") <- paste0("fgt",g)
      if(thresh) attr(rval, "thresh") <- th
      return( rval )
    }

    # calculate variance
    variance <- survey::svrVar( qq, full_design$scale , full_design$rscales, mse = full_design$mse, coef = rval )

    # compute deff
    if ( is.character(deff) || deff || linearized ) {

      # compute threshold linearized function on full sample
      # branch on threshold type and its linearized function
      if( type_thresh == 'relq' ) arptlin <- percent * CalcQuantile_IF( incvec , wsf , quantiles )
      if( type_thresh == 'relm') {
        Yf <- sum( wsf * incvec )
        Nf <- sum( wsf )
        muf <- Yf/Nf
        arptlin <- percent * ( incvec - muf ) / Nf
        arptlin <- arptlin[ wsf > 0 ]
      }
      if ( type_thresh == 'abs' ) {
        arptlin <- rep( 0 , sum( wsf > 0 ) )
      }
      names( arptlin ) <- rownames( full_design$variables )[ wsf > 0 ]

      # compute poverty measure linearized function
      # under fixed poverty threshold
      Nd <- sum( ws )
      ID <- 1 *( rownames( full_design$variables ) %in% rownames( design$variables )[ ws > 0 ] )
      fgtlin1 <- ID * ( h( incvec , th , g ) - rval ) / Nd
      fgtlin1 <- fgtlin1[ wsf > 0]

      # combine both linearization
      Fprime <- densfun( formula , design = design , th , h = NULL , FUN = "F", na.rm = na.rm )
      ahat <- sum( ws * ht( incvar , th , g ) ) / Nd
      ahF <- ahat + h( th , th , g ) * Fprime
      if ( type_thresh %in% c( "relq" , "relm" ) ) fgtlin2 <- ahF * arptlin else fgtlin2 <- 0
      fgtlin <- fgtlin1 + fgtlin2

      # compute deff
      nobs <- length( full_design$pweights )
      npop <- sum( full_design$pweights )
      vsrs <- unclass( survey::svyvar( fgtlin , full_design, na.rm = na.rm, return.replicates = FALSE, estimate.only = TRUE)) * npop^2/nobs
      if (deff != "replace") vsrs <- vsrs * (npop - nobs)/npop
      deff.estimate <- variance / vsrs

      # add indices
      names( fgtlin ) <- rownames( full_design$variables )[ wsf > 0 ]

      # coerce to matrix
      fgtlin <- matrix( fgtlin , nrow = length( fgtlin ) , dimnames = list( names( fgtlin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    }

    # setup result object
    names( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svrepstat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- paste0("fgt",g)
    if(thresh) attr(rval, "thresh") <- th
    if ( linearized ) attr( rval , "linearized" ) <- fgtlin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( fgtlin ) )

    # keep replicates
    if (return.replicates) {
      attr( qq , "scale") <- full_design$scale
      attr( qq , "rscales") <- full_design$rscales
      attr( qq , "mse") <- full_design$mse
      rval <- list( mean = rval , replicates = qq )
      class( rval ) <- c( "cvystat" , "svrepstat" )
    }

    # add design effect estimate
    if ( is.character(deff) || deff ) attr( rval , "deff" ) <- deff.estimate

    # return object
    rval

  }

#' @rdname svyfgt
#' @export
svyfgt.DBIsvydesign <-
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

    NextMethod("svyfgt", design)
  }

