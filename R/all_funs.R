#'Computes the bandwidth needed to compute the derivative of the cdf function
#'
#'Using the whole sample, computes the bandwith used to get the linearized variable
#'
#' @param incvar income variable used in the estimation of the indicators
#' @param w vector of design weights
#' @return value of the bandwidth
#' @author Djalma Pessoa and Anthony Damico
#' @keywords survey
#' @export
h_fun <- function(incvar, w) {
  N <- sum(w)
  sd_inc <- sqrt((sum(w * incvar * incvar) - sum(w * incvar) * sum(w * incvar)/N)/N)
  h <- sd_inc/exp(0.2 * log(sum(w)))
  h
}

#'Estimate the derivative of the cdf function using kernel estimator
#'
#'computes the derivative of a function in a point using kernel estimation
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} from the \code{survey} library.
#' @param x the point where the derivative is calculated
#' @param h value of the bandwidth based on the whole sample
#' @param FUN if \code{F} estimates the derivative of the cdf function; if \code{big_s} estimates the derivative of total in the tails of the distribution
#' @param na.rm Should cases with missing values be dropped?
#' @param ... future expansion
#'
#' @return the value of the derivative at \code{x}
#'
#' @author Djalma Pessoa and Anthony Damico
#'
#' @keywords survey
#' @examples
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#' library(survey)
#' des_eusilc <- svydesign(ids = ~rb030, strata =~db040,  weights = ~rb050, data = eusilc)
#' des_eusilc <- convey_prep( des_eusilc )
#' densfun (~eqincome, design=des_eusilc, 10000, FUN="F" )
#' # linearized design using a variable with missings
#' densfun ( ~ py010n , design = des_eusilc, 10000, FUN="F" )
#' densfun ( ~ py010n , design = des_eusilc , 10000,FUN="F", na.rm = TRUE )
#'
#' @export
densfun <- function(formula, design, x, h = NULL, FUN = "F" , na.rm=FALSE, ...) {

  if( !( FUN %in% c( "F" , "big_s" ) ) ) stop( "valid choices for `FUN=` are 'F' and 'big_s'" )

  incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
  if(na.rm){
    nas<-is.na(incvar)
    design<-design[!nas,]
    if (length(nas) > length(design$prob))
      incvar <- incvar[!nas]
    else incvar[nas] <- 0
  }
  w <- weights(design , "sampling")
  N <- sum(w)
  if(is.null(h)) h <- h_fun(incvar,w)
  u <- (x - incvar)/h
  vectf <- exp(-(u^2)/2)/sqrt(2 * pi)
  if (FUN == "F")
    res <- sum(vectf * w)/(N * h) else {
      v <- w * incvar
      res <- sum(vectf * v)/h
    }
  res
}

#' Linearization of the cumulative distribution function (cdf) of a variable
#'
#' Computes the linearized variable of the cdf function in a point.
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param x the point where the cdf is calculated
#' @param na.rm Should cases with missing values be dropped?
#' @param ... future expansion
#'
#' @return Object of class "\code{cvystat}", which are vectors with a "\code{var}" attribute giving the variance and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Djalma Pessoa and Anthony Damico
#'
#' @seealso \code{\link{svyarpr}}
#'
#' @references Guillaume Osier (2009). Variance estimation for complex indicators
#'of poverty and inequality. \emph{Journal of the European Survey Research
#' Association}, Vol.3, No.3, pp. 167-195,
#' ISSN 1864-3361, URL \url{http://ojs.ub.uni-konstanz.de/srm/article/view/369}.

#'Jean-Claude Deville (1999). Variance estimation for complex statistics and estimators:
#' linearization and residual techniques. Survey Methodology, 25, 193-203,
#' URL \url{http://www5.statcan.gc.ca/bsolc/olc-cel/olc-cel?lang=eng&catno=12-001-X19990024882}.
#'
#' @keywords survey
#' @examples
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#' library(survey)
#' des_eusilc <- svydesign(ids = ~rb030, strata =~db040,  weights = ~rb050, data = eusilc)
#' des_eusilc <- convey_prep( des_eusilc )
#' icdf(~eqincome, design=des_eusilc, 10000 )
#' # linearized design using a variable with missings
#' icdf( ~ py010n , design = des_eusilc, 10000 )
#' icdf( ~ py010n , design = des_eusilc , 10000, na.rm = TRUE )
#' @export
icdf <- function(formula, design, x, na.rm = FALSE, ...) {

  incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
  ncom<- names(design$prob)
  if(na.rm){
    nas<-is.na(incvar)
    design<-design[!nas,]
    if (length(nas) > length(design$prob))
      incvar <- incvar[!nas]
    else incvar[nas] <- 0
  }
  w <- 1/design$prob
  #ind<- names(w)
  N <- sum(w)
  poor <- (incvar <= x) * 1
  value<- sum(poor*w)/N
  lin<-((incvar<=x)-value)/N
  rval <- value
  variance <- survey::svyrecvar(lin/design$prob, design$cluster,
                                design$strata, design$fpc, postStrata = design$postStrata)
  class(rval) <- c( "cvystat" , "svystat" )
  attr(rval, "linearized") <- lin
  attr(rval, "var") <- variance
  attr(rval, "statistic") <- "cdf"
  rval
}


# Functions U and big_t from Jenkins & Biewen:
U_fn <-
  function( x, weights, gamma ) {
    x <- x[weights != 0]

    weights <- weights[weights != 0]

    sum( weights * x^gamma )
  }

T_fn <-
  function( x, weights, gamma ) {
    x <- x[weights != 0]

    weights <- weights[weights != 0]

    sum( weights * x^gamma * log( x ) )
  }



# cvystat print method
#' @method print cvystat
#' @export
print.cvystat <- function(x, ...) {

  if (is.list(x)) x <- x[[1]]
  vv <- attr(x, "var")

  if (is.matrix(vv)) {
    m <- cbind(x, sqrt(diag(vv)))
  } else {
    m <- cbind(x, sqrt(vv))
  }

  # nattr <- length(names(attributes(x)))
  # if (nattr>5) {
  #   for(i in 6:nattr)
  #   {m <- cbind(m, attr(x, names(attributes(x)[i])))}
  #   colnames(m) <- c(attr(x, "statistic"), "SE", names(attributes(x))[6:nattr])
  # }
  # else {
  #   colnames(m) <- c(attr(x, "statistic"), "SE")
  # }

  colnames(m) <- c(attr(x, "statistic"), "SE")

  printCoefmat(m)

}


# cvystat vcov method
#' @export
vcov.cvystat <- function (object, ...) {
  if (is.list(object)) object <- object[[1]]
  as.matrix(attr(object, "var"))
}

# cvystat coef method
#' @export
coef.cvystat <- function(object, ...) {
  if (is.list(object)) object <- object[[1]]
  attr(object, "statistic") <- NULL
  attr(object, "deff") <- NULL
  attr(object, "var") <- NULL
  attr(object, "linearized") <- NULL
  attr(object, "influence") <- NULL
  attr(object, "index") <- NULL
  attr(object, "quantile") <- NULL
  attr(object, "epsilon") <- NULL
  attr(object, "dimensions") <- NULL
  attr(object, "parameters") <- NULL
  attr(object, "extra") <- NULL
  attr(object, "components") <- NULL
  attr(object, "group") <- NULL
  attr(object, "thresh") <- NULL
  unclass(object)
}


#' prepare svydesign and svyrep.design objects for the convey package
#'
#' stores the full survey design (needed for convey functions that use a global poverty threshold) within the design.  this function must be run immediately after the full design object creation with \code{svydesign} or \code{svrepdesign}
#'
#' @param design a survey design object of the library survey.
#'
#' @return the same survey object with a \code{full_design} attribute as the storage space for the unsubsetted survey design
#'
#' @author Djalma Pessoa and Anthony Damico
#'
#' @details  functions in the convey package that use a global poverty threshold require the complete (pre-subsetted) design in order to calculate variances correctly.  this function stores the full design object as a separate attribute so that functions from the \code{survey} package such as \code{subset} and \code{svyby} do not disrupt the calculation of error terms.
#'
#' @keywords survey
#'
#' @examples
#'
#' library(survey)
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#'
#' # linearized design: convey_prep must be run as soon as the linearized design has been created
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#' # now this linearized design object is ready for analysis!
#'
#' # # # CORRECT usage example # # #
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep( des_eusilc )
#' sub_eusilc <- subset( des_eusilc , age > 20 )
#' # since convey_prep() was run immediately after creating the design
#' # this will calculate the variance accurately
#' SE( svyarpt( ~ eqincome , sub_eusilc ) )
#' # # # end of CORRECT usage example # # #
#'
#' # # # INCORRECT usage example # # #
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' sub_eusilc <- subset( des_eusilc , age > 20 )
#' sub_eusilc <- convey_prep( sub_eusilc )
#' # since convey_prep() was not run immediately after creating the design
#' # this will make the variance wrong
#' SE( svyarpt( ~ eqincome , sub_eusilc ) )
#' # # # end of INCORRECT usage example # # #
#'
#' @export
convey_prep <- function(design) {

  if (!is.null(attr(design, "full_design")))stop("convey_prep has already been run on this design")

  if( as.character( design$call )[1] == 'subset' ) warning("this function must be run on the full survey design object immediately after the svydesign() or svrepdesign() call.")

  # store the full design within one of the attributes of the design
  attr(design, "full_design") <- design

  # store the full_design's full_design attribute as TRUE
  attr(attr(design, "full_design"), "full_design") <- TRUE

  class( design ) <- c( "convey.design" , class( design ) )

  design
}

#' @importFrom survey svyby
#' @export
svyby.convey.design <-
  function (formula, by, design, ...){

    if ( ( "DBIsvydesign" %in% class(design) ) & !( "logical" %in% class(attr(design, "full_design"))) ){

      full_design <- attr( design , "full_design" )

      if( 'sex' %in% names( list( ... ) ) ){

        full_design$variables <-
          cbind(
            getvars(formula, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset),
            getvars(by, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset) ,
            getvars(list( ... )[["sex"]], full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset)
          )

        design$variables <-
          cbind(
            getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset),
            getvars(by, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset) ,
            getvars(list( ... )[["sex"]], design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)
          )

      } else if( 'age' %in% names( list( ... ) ) ){

        full_design$variables <-
          cbind(
            getvars(formula, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset),
            getvars(by, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset) ,
            getvars(list( ... )[["age"]], full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset)
          )

        design$variables <-
          cbind(
            getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset),
            getvars(by, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset) ,
            getvars(list( ... )[["age"]], design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)
          )


      } else if( 'subgroup' %in% names( list( ... ) ) ){

        full_design$variables <-
          cbind(
            getvars(formula, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset),
            getvars(by, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset) ,
            getvars(list( ... )[["subgroup"]], full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset)
          )

        design$variables <-
          cbind(
            getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset),
            getvars(by, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset) ,
            getvars(list( ... )[["subgroup"]], design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)
          )


      } else {

        full_design$variables <-
          cbind(
            getvars(formula, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset),
            getvars(by, full_design$db$connection, full_design$db$tablename, updates = full_design$updates, subset = full_design$subset)
          )


        design$variables <-
          cbind(
            getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset),
            getvars(by, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)
          )


      }

      attr( design , "full_design" ) <- full_design

      rm( full_design )

    }

    # remove the "convey.design" and "DBIsvydesign" classes from the current object
    class(design) <- setdiff(class(design), "convey.design")
    class(design) <- setdiff(class(design), "DBIsvydesign")

    # method is unchanged for svyrep.design
    if ( "svyrep.design" %in% class( design ) ) {
      return( survey::svyby( formula , by, design , ... ) )
    }

    # method is unchanged for survey.design2 and !covmat
    if ( "survey.design" %in% class( design ) && ( is.null( list(...)[["covmat"]] ) || isFALSE( list(...)[["covmat"]] ) ) ) {
      return( survey::svyby( formula , by, design , ... ) )
    }

    # method is changes for survey.design2 and covmat
    if ( "survey.design2" %in% class( design ) && ( isTRUE( list(...)[["covmat"]] ) ) ) {
      custom.svyby <- function(formula , by , design , ... ) {
        dots <- list(...)
        FUN <- dots[ sapply( dots , is.function ) ][[1]]
        dots <- dots[ !sapply( dots , is.function ) ]
        dots$covmat <- FALSE
        dots$influence <- TRUE
        # do.call( survey:::svyby.survey.design2 , c( list( formula = formula , by = by , design = design , FUN = FUN ) , dots ) ) # survey's svyby.survey.design2
        do.call( convey:::svyby.survey.design2 , c( list( formula = formula , by = by , design = design , FUN = FUN ) , dots ) ) # convey custom svyby.survey.design2
      }

      # evaluate without covmat matrix
      rval <- custom.svyby( formula , by, design , ... )

      # collect influence functions
      inflmat <- attr( rval , "influence" )
      inflmat.full <- matrix( 0 , nrow = nrow( attr( design , "full_design")$variables ) , ncol = ncol( inflmat ) )
      inflmat.full[ rownames( attr( design , "full_design")$variables ) %in% rownames( inflmat ) , ] <- inflmat

      # compute covariance matrix
      covmat.mat <- survey::svyrecvar( inflmat.full , attr( design , "full_design" )$cluster, attr( design , "full_design" )$strata, attr( design , "full_design" )$fpc, postStrata = attr( design , "full_design" )$postStrata)

      # add covariance matrix
      attr( rval, "var" ) <- covmat.mat

      # return result
      return( rval )

    }

  }
