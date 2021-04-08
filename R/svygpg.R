#' Linearization of the gender pay (wage) gap
#'
#' Estimate the difference between the average gross hourly earnings of men and women expressed as a percentage of the average gross hourly earnings of men.
#'
#'
#' @param formula a formula specifying the gross hourly earnings variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param sex formula with a factor with labels 'male' and 'female'
#' @param na.rm Should cases with missing values be dropped?
#' @param ... future expansion
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#'@return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
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
#' library(laeken)
#' library(survey)
#' data(ses)
#' names( ses ) <- gsub( "size" , "size_" , tolower( names( ses ) ) )
#' des_ses <- svydesign(id=~1, weights=~weights, data=ses)
#' des_ses <- convey_prep(des_ses)
#'
#' # linearized design
#' svygpg(~earningshour, des_ses, ~sex)
#' # replicate-weighted design
#' des_ses_rep <-  as.svrepdesign( des_ses , type = "bootstrap" )
#' des_ses_rep <- convey_prep(des_ses_rep)
#'
#' svygpg(~earningshour, des_ses_rep, ~sex)
#'
#' \dontrun{
#'
#' # database-backed design
#' library(RSQLite)
#' library(DBI)
#' dbfile <- tempfile()
#' conn <- dbConnect( RSQLite::SQLite() , dbfile )
#' dbWriteTable( conn , 'ses' , ses )
#'
#' dbd_ses <- svydesign(id=~1, weights=~weights, data="ses", dbname=dbfile, dbtype="SQLite")
#' dbd_ses <- convey_prep( dbd_ses )
#'
#' svygpg(formula=~earningshour, design=dbd_ses, sex= ~sex)
#'
#' dbRemoveTable( conn , 'ses' )
#'
#' }
#'
#' @export
svygpg <-
  function(formula, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    UseMethod("svygpg", design)

  }

#' @rdname svygpg
#' @export
svygpg.survey.design <-
  function(formula, design, sex,  na.rm=FALSE,...) {

    # test for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your linearized survey design object immediately after creating it with the svydesign() function.")

    # collect data
    wagevar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
    sex <- model.matrix(eval(bquote(~0 + .(sex))), design$variables , na.action = na.pass )

    # treat missing values
    if( na.rm ) {
      nas<-rowSums( is.na( cbind( wagevar , sex ) ) )
      design<-design[nas==0,]
      if (length(nas) > length(design$prob)){
        wagevar <- wagevar[nas == 0]
        sex <- sex[nas==0,]
      } else{
        wagevar[nas > 0] <- 0
        sex[nas > 0,] <- 0
      }
    }

    # collect weights
    w <- 1 / design$prob
    ind <- names(design$prob)

    # fix names
    col_female <- grep("female", colnames(sex))
    col_male <- setdiff(1:2, col_female)

    # create linearization objects of totals
    INDM <- list(value = sum(sex[, col_male]*w), lin=sex[, col_male])
    INDF <- list(value = sum(sex[, col_female]*w), lin=sex[, col_female])
    TM <- list(value = sum(wagevar*sex[, col_male]*w), lin=wagevar*sex[, col_male])
    TF <- list(value = sum(wagevar*sex[, col_female]*w), lin=wagevar*sex[, col_female])
    list_all_tot <- list(INDM=INDM,INDF=INDF,TM=TM,TF=TF)

    # compute estimate
    IGPG <- contrastinf( quote( ( TM / INDM - TF / INDF ) / ( TM / INDM ) ) , list_all_tot )
    rval <- IGPG$value

    # calculate variance
    infun <- as.numeric( IGPG$lin )
    variance <- survey::svyrecvar(infun/design$prob, design$cluster, design$strata, design$fpc, postStrata = design$postStrata)
    variance[ is.nan( variance ) ] <- NA

    # add indices to linearized function
    names( infun ) <- rownames( design$variables )[ !is.infinite( design$prob ) ]

    # build result object
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]
    class(rval) <- c( "cvystat" , "svystat" )
    attr( rval , "var" ) <- variance
    attr(rval, "linearized" ) <- infun
    attr( rval , "index" ) <- as.numeric( names( infun ) )
    attr( rval , "statistic" ) <- "gpg"
    rval

  }


#' @rdname svygpg
#' @export
svygpg.svyrep.design <-
  function(formula, design, sex,na.rm=FALSE, ...) {

    # check for convey_prep
    if (is.null(attr(design, "full_design"))) stop("you must run the ?convey_prep function on your replicate-weighted survey design object immediately after creating it with the svrepdesign() function.")

    # collect data
    wagevar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
    sex <- model.matrix(eval(bquote(~0 + .(sex))), design$variables , na.action = na.pass )

    # treat missing values
    if( na.rm ) {
      nas<-rowSums( is.na( cbind( wagevar , sex ) ) )
      design<-design[nas==0,]
      if (length(nas) > length(design$prob)){
        wagevar <- wagevar[nas == 0]
        sex <- sex[nas==0,]
      } else{
        wagevar[nas > 0] <- 0
        sex[nas > 0,] <- 0
      }
    }

    # collect smapling weights
    ws <- weights(design, "sampling")
    design <- update(design, one = rep(1, length(wagevar)))

    # computation function
    ComputeGpg <-
      function(earn_hour, w, sex) {
        # filter observations
        earn_hour <- earn_hour[ w > 0 ]
        sex <- sex[ w > 0 , ]
        w <- w[ w > 0 ]

        # run computation
        col_female <- grep("female", colnames(sex))
        col_male <- setdiff(1:2, col_female)
        ind_men <- sex[, col_male]
        ind_fem <- sex[, col_female]
        med_men <- sum(ind_men * earn_hour * w)/sum(ind_men * w)
        med_fem <- sum(ind_fem * earn_hour * w)/sum(ind_fem * w)
        gpg <- (med_men - med_fem)/med_men
        gpg

      }

    # compute point estimate
    rval <- ComputeGpg( earn_hour = wagevar , w = ws , sex = sex )

    # collect analysis weights
    ww <- weights(design, "analysis")

    # compute replicates
    qq <- apply( ww , 2 , function(wi) ComputeGpg( wagevar , wi , sex = sex ) )

    # ompute variance
    if (anyNA(qq)) variance <- NA else variance <- survey::svrVar(qq, design$scale, design$rscales, mse = design$mse, coef = rval)
    variance <- as.matrix( variance )
    colnames( variance ) <- rownames( variance ) <-  names( rval ) <- strsplit( as.character( formula )[[2]] , ' \\+ ' )[[1]]

    # build result object
    class(rval) <- c( "cvystat" , "svrepstat" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "gpg"
    rval

  }

#' @rdname svygpg
#' @export
svygpg.DBIsvydesign <-
  function (formula, design, sex, ...){

    if (!( "logical" %in% class(attr(design, "full_design"))) ){

      full_design <- attr( design , "full_design" )

      full_design$variables <-
        cbind(
          getvars(formula, attr( design , "full_design" )$db$connection, attr( design , "full_design" )$db$tablename,updates = attr( design , "full_design" )$updates, subset = attr( design , "full_design" )$subset),

          getvars(sex, attr( design , "full_design" )$db$connection, attr( design , "full_design" )$db$tablename,updates = attr( design , "full_design" )$updates, subset = attr( design , "full_design" )$subset)
        )

      attr( design , "full_design" ) <- full_design

      rm( full_design )

    }

    design$variables <-
      cbind(
        getvars(formula, design$db$connection,design$db$tablename, updates = design$updates, subset = design$subset),

        getvars(sex, design$db$connection, design$db$tablename,updates = design$updates, subset = design$subset)
      )

    NextMethod("svygpg", design)
  }
