#' Relative median poverty gap
#'
#' Estimate the median of incomes less than the at-risk-of-poverty threshold (\code{arpt}).
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param quantiles income quantile, usually .5 (median)
#' @param percent fraction of the quantile, usually .60
#' @param na.rm Should cases with missing values be dropped?
#' @param ... arguments passed on to `survey::svyquantile`
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
#' svypoormed( ~eqincome , design = des_eusilc )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' svypoormed( ~eqincome , design = des_eusilc_rep )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svypoormed( ~ py010n , design = des_eusilc )
#' svypoormed( ~ py010n , design = des_eusilc , na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svypoormed( ~ py010n , design = des_eusilc_rep )
#' svypoormed( ~ py010n , design = des_eusilc_rep , na.rm = TRUE )
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
#' svypoormed( ~ eqincome , design = dbd_eusilc )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svypoormed <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svypoormed", design)

  }

#' @rdname svypoormed
#' @export
svypoormed.survey.design <-
  function(formula, design, quantiles = 0.5, percent = 0.6, na.rm=FALSE, ...) {

    # test for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # if the class of the full_design attribute is just a TRUE, then the design is
    # already the full design.  otherwise, pull the full_design from that attribute.
    if ("logical" %in% class(attr(design, "full_design"))) full_design <- design else full_design <- attr(design, "full_design")

    # collect domain sample income data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar)
      design<-design[!nas,]
      if (length(nas) > length(design$prob)) incvar <- incvar[!nas] else incvar[nas] <- 0
    }

    # collect domain weights
    w <- 1 / design$prob

    # treat lengths
    if( is.null( names( design$prob ) ) ) ind <- as.character( seq( length( design$prob ) ) ) else ind <- names(design$prob)

    # domain population size
    N <- sum(w)

    # collect full sample income data
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvec)
      full_design<-full_design[!nas,]
      if (length(nas) > length(full_design$prob)) incvec <- incvec[!nas] else incvec[nas] <- 0
    }

    # collect full sample weights
    wf <- 1 / full_design$prob

    # collect full sample indices
    if( is.null( names( full_design$prob ) ) ) ncom <- as.character( seq( length( full_design$prob ) ) ) else ncom <- names(full_design$prob)

    # point estimation
    htot <- h_fun(incvar, w)
    ARPT <- svyarpt(formula = formula, full_design, quantiles = quantiles, percent = percent, na.rm = na.rm)
    arpt <- coef(ARPT)
    linarpt <- attr(ARPT, "linearized")

    # treat remaining missing
    if(is.na(arpt)){
      rval <- NA
      variance <- as.matrix( NA )
      colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      class(rval) <- c( "cvystat" , "svystat" )
      attr( rval , "var" ) <- variance
      attr( rval , "statistic" ) <- "poormed"
      return( rval )
    }

    # below threshold sample domain
    nome <- terms.formula(formula)[[2]]
    dsub <- eval( substitute ( subset(design, incvar <= arpt ) , list( incvar = nome , arpt = as.numeric( arpt ) ) ) )
    # if( nrow( dsub ) == 0 ) stop( paste("zero records in the set of poor people.  determine the poverty threshold by running svyarpt on ~",nome ) )
    if( nrow( dsub ) <= 0 ){
      warning( paste("zero records in the set of poor people.  determine the poverty threshold by running svyarpt on ~",nome ) )
      rval <- NA
      variance <- as.matrix( 0 )
      colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      class(rval) <- c( "cvystat" , "svystat" )
      attr( rval , "var" ) <- variance
      attr(rval, "linearized") <- matrix( NA , nrow = sum( w > 0 ) , ncol = 1 , dimnames = list( names( w )[ w > 0 ] , NULL ) )
      attr(rval, "index") <- as.numeric( names( w )[ w > 0 ] )
      attr( rval , "statistic" ) <- "poormed"
      return( rval )
    }

    # compute median income of the poor
    # medp <- survey::svyquantile(x = formula, dsub, 0.5, method = "constant", na.rm=na.rm,...) # old
    ww <- 1 / dsub$prob
    incvar.top <- model.frame(formula, dsub$variables, na.action = na.pass)[[1]]
    medp <- computeQuantiles( incvar.top , ww , .5 )
    medp <- as.vector(medp)

    # add indices
    if (sum(1/design$prob==0) > 0) ID <- 1*(1/design$prob!=0) else ID <- 1 * ( ncom %in% ind )

    # compute at risk of poverty rate
    ARPR <- svyarpr( formula=formula, design= design, quantiles, percent, na.rm = na.rm )
    arpr <- coef( ARPR )
    ifarpr <- attr( ARPR, "linearized" )[,1]

    # ensure length
    if ( length( ifarpr ) != length( full_design$prob ) ) {
      tmplin <- rep( 0 , nrow( full_design$variables ) )
      tmplin[ rownames( full_design$variables ) %in% names( ifarpr ) ] <- ifarpr
      ifarpr <- tmplin ; rm( tmplin )
      names( ifarpr ) <- rownames( full_design$variables )
    }

    # compute density
    Fprimemedp <- densfun(formula = formula, design = design, medp, h = htot, FUN = "F", na.rm = na.rm )

    # linearize cdf of medp
    ifmedp <- ( 1 / N ) * ID * ( ( incvec <= medp ) - 0.5 * arpr )

    # linearize median of poor
    linmedp <- ( 0.5 * ifarpr - ifmedp ) / Fprimemedp

    # ensure length
    if ( length( linmedp ) != length( full_design$prob ) ) {
      tmplin <- rep( 0 , nrow( full_design$variables ) )
      tmplin[ rownames( full_design$variables ) %in% names( linmedp ) ] <- linmedp
      linmedp <- tmplin ; rm( tmplin )
      names( linmedp ) <- rownames( full_design$variables )
    }

    # computte variance
    variance <- survey::svyrecvar(linmedp/full_design$prob, full_design$cluster, full_design$strata, full_design$fpc, postStrata = full_design$postStrata)

    # keep necessary linearized functions
    linmedp <- linmedp[ 1/full_design$prob > 0 ]

    # coerce to matrix
    linmedp <- matrix( linmedp , nrow = length( linmedp ) , dimnames = list( names( linmedp ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    # build result object
    rval <- medp
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    attr( rval , "var" ) <- variance
    attr( rval , "statistic" ) <- "poormed"
    class(rval) <- c( "cvystat" , "svystat" )
    attr( rval , "linearized" ) <- linmedp
    attr( rval , "index" ) <- as.numeric( rownames( linmedp ) )
    rval

  }


#' @rdname svypoormed
#' @export
svypoormed.svyrep.design <-
  function(formula, design, quantiles = 0.5, percent = 0.6,na.rm=FALSE, ...) {

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

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

    # collect domain sample weights
    ws <- weights(design, "sampling")

    # colelct full sample data
    df_full <- model.frame(full_design)
    incvec <- model.frame(formula, full_design$variables, na.action = na.pass)[[1]]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvec)
      full_design<-full_design[!nas,]
      df_full <- model.frame(full_design)
      incvec <- incvec[!nas]
    }

    # collect full sample weighst
    wsf <- weights(full_design,"sampling")

    # create full sample indices
    names(incvec) <- names(wsf) <- rownames( full_design$variables )
    ind <- rownames( full_design$variables ) %in% rownames( design$variables )

    # computation function
    ComputePoormed <-
      function( xf , wf , ind , quantiles , percent ){
        tresh <- percent * computeQuantiles(xf, wf, p = quantiles)
        x <-xf[ind]
        w <- wf[ind]
        indpoor <- (x <= tresh)
        medp <- computeQuantiles( x[indpoor], w[indpoor], p = 0.5 )
        medp
      }

    # compute point estimate
    rval <- ComputePoormed(xf = incvec, wf=wsf, ind= ind, quantiles = quantiles, percent = percent)

    # collect domain sample weights
    ws <- weights(design, "sampling")

    # collect full sample analysis weights
    wwf <- weights(full_design, "analysis")

    # compute replicates
    qq <- apply( wwf , 2 , function( wi ) {
      names( wi ) <- rownames( df_full$variables )
      ComputePoormed( incvec , wi , ind=ind , quantiles = quantiles , percent = percent )
    } )

    # compute variance
    variance <- if ( all( is.na( qq ) ) ) NA else survey::svrVar(qq, design$scale, design$rscales, mse = design$mse, coef = rval)
    variance <- as.matrix( variance )
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # build result object
    class(rval) <- c( "cvystat" , "svrepstat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "poormed"
    rval

  }

#' @rdname svypoormed
#' @export
svypoormed.DBIsvydesign <-
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

    NextMethod("svypoormed", design)
  }

