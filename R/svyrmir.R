#' Relative median income ratio
#'
#' Estimate the ratio between the median income of people with age above 65 and the median income of people with age below 65.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param age formula defining the variable age
#' @param agelim the age cutpoint, the default is 65
#' @param quantiles income quantile, usually .5 (median)
#' @param na.rm Should cases with missing values be dropped?
#' @param med_old return the median income of people older than agelim
#' @param med_young return the median income of people younger than agelim
#' @param linearized Should a matrix of linearized variables be returned
#' @param ... arguments passed on to `survey::svyquantile`
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
#' library(survey)
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#'
#' # missing completely at random, missingness rate = .20
#' ind_miss <- rbinom(nrow(eusilc), 1, .20 )
#' eusilc$eqincome_miss <- eusilc$eqincome
#' is.na(eusilc$eqincome_miss)<- ind_miss==1
#'
#' # linearized design
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep(des_eusilc)
#'
#' svyrmir( ~eqincome , design = des_eusilc , age = ~age, med_old = TRUE )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' svyrmir( ~eqincome , design = des_eusilc_rep, age= ~age, med_old = TRUE )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svyrmir( ~ eqincome_miss , design = des_eusilc,age= ~age)
#' svyrmir( ~ eqincome_miss , design = des_eusilc , age= ~age, na.rm = TRUE )
#' # replicate-weighted design using a variable with missings
#' svyrmir( ~ eqincome_miss , design = des_eusilc_rep,age= ~age )
#' svyrmir( ~ eqincome_miss , design = des_eusilc_rep ,age= ~age, na.rm = TRUE )
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
#' svyrmir( ~eqincome , design = dbd_eusilc , age = ~age )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyrmir <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svyrmir", design)

  }

#' @rdname svyrmir
#' @export
svyrmir.survey.design  <-
  function(formula, design, age, agelim = 65, quantiles=0.5, na.rm=FALSE, med_old = FALSE, med_young = FALSE, linearized = FALSE , ...){

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # collect domain data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
    agevar <- model.frame(age, design$variables, na.action = na.pass)[[1]]
    x <- cbind(incvar,agevar)

    # treat missing values
    if(na.rm){
      nas<-rowSums(is.na(x))
      design<-design[nas==0,]
      if (length(nas) > length(design$prob)){
        incvar <- incvar[nas == 0]
        agevar <- agevar[nas==0]
      } else{
        incvar[nas > 0] <- 0
        agevar[nas > 0] <- 0
      }
    }

    # create indices
    if( is.null( names( design$prob ) ) ) names(design$prob) <- as.character( seq( length( design$prob ) ) )

    # collect domain weights
    w <- 1/design$prob

    # compute population size
    N <- sum(w)
    h <- h_fun(incvar,w)
    age.name <- terms.formula(age)[[2]]

    ### < 65 yo median income
    dsub1 <- eval( substitute( within_function_subset( design , age < agelim ) , list( age = age.name, agelim = agelim ) ) )
    if( sum( 1/dsub1$prob > 0 ) == 0 ) stop( "zero records in the set of non-elderly people" )
    if( "DBIsvydesign" %in% class( dsub1 ) ) {
      ind1 <- rownames( design$variables ) %in% which(is.finite( dsub1$prob ) )
    } else{
      ind1 <- rownames( design$variables ) %in% rownames( dsub1$variables )
    }
    q_alpha1 <- svyiqalpha( formula , dsub1 , quantiles , linearized = TRUE , ... )

    ### >= 65 yo income
    dsub2 <- eval( substitute( within_function_subset( design , age >= agelim ) , list( age = age.name, agelim = agelim ) ) )
    if( sum( 1/dsub2$prob > 0 ) == 0 ) stop( "zero records in the set of elderly people" )
    if( "DBIsvydesign" %in% class( dsub2 ) ) {
      ind2 <- rownames(design$variables) %in% which( is.finite( dsub2$prob ) )
    } else{
      ind2 <- rownames( design$variables ) %in% rownames( dsub2$variables )
    }

    # compute quantiles
    q_alpha2 <- svyiqalpha( formula , dsub2 , quantiles , linearized = TRUE , ... )

    # create objects for contrastinf
    MED1 <- list(value = q_alpha1[[1]] , lin = rep( 0 , nrow( design$variables ) ) )
    MED2 <- list(value = q_alpha2[[1]] , lin = rep( 0 , nrow( design$variables ) ) )

    # add partial linearizations
    MED1$lin [ rownames( design$variables ) %in% rownames( attr( q_alpha1 , "linearized" ) ) ] <- attr( q_alpha1 , "linearized" )[,1]
    MED2$lin [ rownames( design$variables ) %in% rownames( attr( q_alpha2 , "linearized" ) ) ] <- attr( q_alpha2 , "linearized" )[,1]
    list_all<- list(MED1=MED1, MED2=MED2)

    # linearize ratio of medians
    RMED <- contrastinf( quote( MED2/MED1 ),list_all )
    lin <- RMED$lin[,1]
    names( lin ) <- rownames( design$variables )

    # compute variance
    variance <- survey::svyrecvar( lin/design$prob, design$cluster, design$strata, design$fpc, postStrata = design$postStrata )
    variance[ which( is.nan( variance ) ) ] <- NA
    colnames( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # keep necessary linearized functions
    lin <- lin[ 1/design$prob > 0 ]

    # coerce to matrix
    lin <- matrix( lin , nrow = length( lin ) , dimnames = list( names( lin ) , strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]] ) )

    # build result object
    rval <- as.numeric( RMED$value )
    names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr(rval, "var") <- variance
    if ( linearized ) attr(rval,"linearized") <- lin
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin ) )
    if (med_old) attr( rval, "med_old") <- q_alpha2
    if (med_young) attr( rval, "med_young") <- q_alpha1
    attr( rval , "statistic" ) <- "rmir"
    rval

  }


#' @rdname svyrmir
#' @export
svyrmir.svyrep.design <-
  function(formula, design, age, agelim = 65, quantiles = 0.5, na.rm=FALSE, med_old = FALSE, med_young = FALSE,...) {

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

    # collect data
    df <- model.frame(design)
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
    agevar <- model.frame(age, design$variables, na.action = na.pass)[[1]]
    x <- cbind(incvar,agevar)

    # treat missing values
    if(na.rm){
      nas<-rowSums(is.na(x))
      design<-design[nas==0,]
      df <- model.frame(design)
      incvar <- incvar[nas==0]
      agevar<- agevar[nas==0]
    }

    # computation function
    ComputeRmir <-
      function(x, w, quantiles, age, agelim) {
        indb <- age < agelim
        inda <-  age >= agelim
        quant_below <- computeQuantiles( x[indb], w[indb], p = quantiles )
        quant_above <- computeQuantiles( x[inda], w[inda], p = quantiles )
        c(quant_above, quant_below, quant_above/quant_below)
      }

    # collect sampling weights
    ws <- weights(design, "sampling")

    # compute pont estimates
    Rmir_val <- ComputeRmir( incvar , ws , quantiles, agevar,  agelim )


    # treat missing
    if ( is.na( Rmir_val[[3]] ) ) {
      rval <- NA
      variance <-  as.matrix( NA )
      class(rval) <- c( "cvystat" , "svrepstat" )
      names( rval ) <- names( variance ) <- rownames( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
      attr( rval , "var" ) <- variance
      attr( rval , "statistic" ) <- "rmir"
      if (med_old) attr( rval, "med_old") <- Rmir_val[1]
      if (med_young) attr( rval, "med_young") <- Rmir_val[2]
      return( rval )
    }
    # collect analysis weights
    ww <- weights(design, "analysis")

    # compute replicates
    qq <- t( apply( ww, 2, function(wi) ComputeRmir( incvar , wi , quantiles, agevar , agelim ) ) )

    # compute variance
    variance <- survey::svrVar( qq[ , 3 ], design$scale, design$rscales, mse = design$mse, coef = Rmir_val[3] )
    names( variance ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # build result object
    rval <- Rmir_val[[3]]
    class(rval) <- c( "cvystat" , "svrepstat" )
    names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    attr( rval , "var" ) <- variance
    attr( rval , "statistic" ) <- "rmir"
    if (med_old) attr( rval, "med_old") <- Rmir_val[1]
    if (med_young) attr( rval, "med_young") <- Rmir_val[2]
    rval

  }

#' @rdname svyrmir
#' @export
svyrmir.DBIsvydesign <-
  function (formula, design, age, ...){

    if (!( "logical" %in% class(attr(design, "full_design"))) ){

      full_design <- attr( design , "full_design" )

      full_design$variables <-
        cbind(
          getvars(formula, attr( design , "full_design" )$db$connection, attr( design , "full_design" )$db$tablename, updates = attr( design , "full_design" )$updates, subset = attr( design , "full_design" )$subset),
          getvars(age, attr( design , "full_design" )$db$connection, attr( design , "full_design" )$db$tablename,
                  updates = attr( design , "full_design" )$updates, subset = attr( design , "full_design" )$subset)
        )

      attr( design , "full_design" ) <- full_design

      rm( full_design )

    }

    design$variables <-
      cbind(
        getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset),
        getvars(age, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)
      )

    NextMethod("svyrmir", design)
  }
