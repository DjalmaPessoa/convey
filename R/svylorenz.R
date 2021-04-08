#' Lorenz curve
#'
#' Estimate the Lorenz curve, an inequality graph
#'
#' @param formula a formula specifying the income variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param quantiles a sequence of probabilities that defines the quantiles sum to be calculated
#' @param ci Should the confidence interval be plotted? Defaults to \code{TRUE}.
#' @param alpha a number that especifies de confidence level for the graph.
#' @param empirical Should an empirical Lorenz curve be estimated as well? Defaults to \code{FALSE}.
#' @param plot Should the Lorenz curve be plotted? Defaults to \code{TRUE}.
#' @param add Should a new curve be plotted on the current graph?
#' @param curve.col a string defining the color of the curve.
#' @param na.rm Should cases with missing values be dropped? Defaults to \code{FALSE}.
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param return.replicates Return the replicate estimates?
#' @param ... additional arguments passed to \code{plot} methods
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' Notice that the 'empirical' curve is observation-based and is the one actually used to calculate the Gini index.
#' On the other hand, the quantile-based curve is used to estimate the shares, SEs and confidence intervals.
#'
#' This way, as the number of quantiles of the quantile-based function increases, the quantile-based curve approacches the observation-based curve.
#'
#' @return Object of class "\code{svyquantile}", which are vectors with a "\code{quantiles}" attribute giving the proportion of income below that quantile,
#' and a "\code{SE}" attribute giving the standard errors of the estimates.
#'
#' @author Guilherme Jacob, Djalma Pessoa and Anthony Damico
#'
#' @seealso \code{\link{svyquantile}}
#'
#' @references Milorad Kovacevic and David Binder (1997). Variance Estimation for Measures of Income
#' Inequality and Polarization - The Estimating Equations Approach. \emph{Journal of Official Statistics},
#' Vol.13, No.1, 1997. pp. 41 58. URL \url{https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/variance-estimation-for-measures-of-income-inequality-and-polarization---the-estimating-equations-approach.pdf}.
#'
#' Shlomo Yitzhaki and Robert Lerman (1989). Improving the accuracy of estimates of Gini coefficients.
#' \emph{Journal of Econometrics}, Vol.42(1), pp. 43-47, September.
#'
#' Matti Langel (2012). \emph{Measuring inequality in finite population sampling}. PhD thesis. URL \url{http://doc.rero.ch/record/29204}.
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
#' svylorenz( ~eqincome , des_eusilc, seq(0,1,.05), alpha = .01 )
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep( des_eusilc_rep )
#'
#' svylorenz( ~eqincome , des_eusilc_rep, seq(0,1,.05), alpha = .01 )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' svylorenz( ~py010n , des_eusilc, seq(0,1,.05), alpha = .01 )
#' svylorenz( ~py010n , des_eusilc, seq(0,1,.05), alpha = .01, na.rm = TRUE )
#' # demonstration of `curve.col=` and `add=` parameters
#' svylorenz( ~eqincome , des_eusilc, seq(0,1,.05), alpha = .05 , add = TRUE , curve.col = 'green' )
#' # replicate-weighted design using a variable with missings
#' svylorenz( ~py010n , des_eusilc_rep, seq(0,1,.05), alpha = .01 )
#' svylorenz( ~py010n , des_eusilc_rep, seq(0,1,.05), alpha = .01, na.rm = TRUE )
#'
#'
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
#' svylorenz( ~eqincome , dbd_eusilc, seq(0,1,.05), alpha = .01 )
#'
#' # highlithing the difference between the quantile-based curve and the empirical version:
#' svylorenz( ~eqincome , dbd_eusilc, seq(0,1,.5), empirical = TRUE, ci = FALSE, curve.col = "green" )
#' svylorenz( ~eqincome , dbd_eusilc, seq(0,1,.5), alpha = .01, add = TRUE )
#' legend( "topleft", c("Quantile-based", "Empirical"), lwd = c(1,1), col = c("red", "green"))
#' # as the number of quantiles increases, the difference between the curves gets smaller
#' svylorenz( ~eqincome , dbd_eusilc, seq(0,1,.01), empirical = TRUE, ci = FALSE, curve.col = "green" )
#' svylorenz( ~eqincome , dbd_eusilc, seq(0,1,.01), alpha = .01, add = TRUE )
#' legend( "topleft", c("Quantile-based", "Empirical"), lwd = c(1,1), col = c("red", "green"))
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline lines plot points polygon
#' @importFrom utils tail
#'
#' @export
svylorenz <- function(formula, design, ...) {

  if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

  UseMethod("svylorenz", design)

}


# thanks to alex for these wrapper functions
# http://stackoverflow.com/a/37518103/1759499
svylorenzplot_wrap <-
  function( cex = 0.1 , xlab = "Cumulative Population Share" , ylab = "Total Income Share" , ... ){
    plot(
      NULL ,
      NULL ,
      xlim = c( 0 , 1 ) ,
      ylim = c( 0 , 1 ) ,
      cex = cex ,
      xlab = xlab ,
      ylab = ylab ,
      ...
    )
  }

svylorenzlines_wrap <-
  function( x = x , y = y , pch = 16 , cex = 0.1 , lwd = 1 , col = col , ... ){
    lines( x , y , xlim = c( 0 , 1 ) , ylim = c( 0 , 1 ) , pch = pch , cex = cex , lwd = lwd , col = col , ... )
  }

svylorenzpoints_wrap <-
  function( x = x , y = y , pch = 16 , cex = 0.1 , lwd = 1 , col = col , ... ){

    points( x, y , xlim = c( 0 , 1 ) , ylim = c( 0 , 1 ) , pch = pch , cex = cex * 4 , lwd = lwd , col = col , ... )

  }

svylorenzpolygon_wrap <-
  function( x = x , y = y , col = col , border = NA , ... ){
    polygon( x , y , col = col , border = border , ... )
  }


#' @rdname svylorenz
#' @export
svylorenz.survey.design <- function ( formula , design, quantiles = seq(0,1,.1), empirical = FALSE, plot = TRUE, add = FALSE, curve.col = "red", ci = TRUE, alpha = .05, na.rm = FALSE , deff = FALSE , linearized = FALSE , ... ) {

  # collect income data
  incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

  # treat missing values
  if (na.rm) {
    nas <- is.na(incvar)
    design <- design[nas == 0, ]
    if (length(nas) > length(design$prob))
      incvar <- incvar[nas == 0]
    else incvar[nas > 0] <- 0
  }

  # collect weights
  w <- 1/design$prob

  # treat remaining missing
  if ( anyNA( incvar [ w > 0 ] ) ) {
    variance <- as.matrix(NA)
    cis <- array( rbind(rep(NA, length(quantiles)),rep(NA, length(quantiles))), dim = c(2, length(quantiles)), dimnames = list( c( "(lower", "upper)" ), as.character(quantiles) ) )
    rval <- t( matrix( data = rep(NA, length(quantiles)), nrow = length(quantiles), dimnames = list( as.character( quantiles ), as.character(formula)[2] ) ) )
    rval <- list(quantiles = rval, CIs = cis)
    attr(rval, "SE") <- rep(NA, length(quantiles))
    class(rval) <- c( "cvyquantile" , "svyquantile" )
    return(rval)
  }

  # calculate estimates
  N <- sum( w )
  Y <- sum( ifelse( w >0 , w * incvar , 0 ) )
  mu <- Y/N
  L.p <- as.numeric( lapply( quantiles, function( z ) CalcLorenz( incvar, w, z ) ) )
  GL.p <- L.p * mu

  # compute linearized
  lin.matrix <- do.call( cbind , lapply( quantiles, function( z ) CalcLorenz_IF( incvar, w , z ) ) )

  # treat out of sample
  if ( nrow( lin.matrix ) != length( design$prob ) ) {
    rownames( lin.matrix ) <- rownames( design$variables )[ w > 0 ]
    lin.matrix <- lin.matrix[ pmatch( rownames( design$variables ) , rownames(lin.matrix ) ) , ]
    lin.matrix[ w <= 0 , ] <- 0
  }

  # compute variance
  variance <- survey::svyrecvar( sweep( lin.matrix , 1 , 1/design$prob , "*" ) , design$cluster, design$strata, design$fpc, postStrata = design$postStrata )
  variance[ which( is.nan( variance ) ) ] <- NA
  rownames( variance ) <- colnames( variance ) <- paste0( "L(" , quantiles , ")" )
  se <- sqrt( diag( variance ) )

  # compute deff
  if ( is.character( deff ) || deff ) {
    nobs <- sum( weights( design , "sampling" ) > 0 )
    npop <- sum( weights( design , "sampling" ) )
    if ( deff == "replace" ) vsrs <- survey::svyvar( lin.matrix , design, na.rm = na.rm) * npop^2/nobs
    else vsrs <- survey::svyvar( lin.matrix , design , na.rm = na.rm ) * npop^2 * (npop - nobs)/(npop * nobs)
    deff.estimate <- diag( variance )/ vsrs
  }

  # keep necessary linearized functions
  lin.matrix <- lin.matrix[ 1/design$prob > 0 , ]
  rownames( lin.matrix ) <- rownames( design$variables )[ w > 0 ]
  colnames( lin.matrix ) <- paste0( "L(" , quantiles , ")" )

  # compute CIs
  CI.L <- L.p - se * qnorm( alpha , mean = 0, sd = 1, lower.tail = FALSE )
  CI.U <- L.p + se * qnorm( alpha , mean = 0, sd = 1, lower.tail = FALSE )
  cis <- structure( rbind( CI.L,CI.U ), .Dim = c(2L, length(quantiles), 1L), .Dimnames = list(c("(lower", "upper)"), as.character(quantiles),  as.character(formula)[2]))

  # build result object
  rval <- c( L.p )
  names( rval ) <- paste0( "L(" , quantiles , ")" )
  attr( rval , "var") <- variance
  attr( rval , "statistic") <- "lorenz"
  class( rval ) <- c( "cvystat" , "svystat" )
  if ( is.character(deff) || deff ) attr(rval,"deff") <- deff.estimate
  if ( linearized ) attr(rval,"linearized") <- lin.matrix
  if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )

  # calculate empirical curve
  if (empirical) empirical.lorenz <- emp.interp( incvar , w )

  # plot data
  if ( plot ) {

    plot_dots <- list( ... )

    # remove `deff` argument sent by svyby
    if( 'deff' %in% names( plot_dots ) ) plot_dots$deff <- NULL

    # overplot plot
    if ( !add ) do.call( svylorenzplot_wrap , plot_dots )

    # add error for special options
    if( any( c( 'xlim' , 'ylim' , 'col' ) %in% names( list( ... ) ) ) ) stop( "xlim=, ylim=, and col= parameters are fixed within `svylorenz`.  use curve.col= to change the line color" )

    # add reference line
    abline( 0 , 1 , ylim = c( 0 , 1 ) , plot_dots )

    # add empirical curve
    if( empirical ) {
      lines_dots <- plot_dots
      lines_dots$x <- empirical.lorenz$lorenz_x
      lines_dots$y <- empirical.lorenz$lorenz_y
      lines_dots$col = curve.col
      do.call( svylorenzlines_wrap , lines_dots )
    }

    # add main points
    points_dots <- plot_dots
    points_dots$x <- quantiles
    points_dots$y <- L.p
    points_dots$col <- curve.col

    do.call( svylorenzpoints_wrap , points_dots )

    # add confidence intervals
    if (ci) {
      X.Vec <- as.numeric( c(quantiles, tail(quantiles, 1), rev(quantiles), quantiles[1]) )
      Y.Vec <- as.numeric( c( CI.L, tail(CI.U, 1), rev(CI.U), CI.L[1] ) )

      polygon_dots <- plot_dots
      polygon_dots$x <- X.Vec
      polygon_dots$y <- Y.Vec
      polygon_dots$col <- adjustcolor( curve.col, alpha.f = .2)
      polygon_dots$border <- NA

      do.call( svylorenzpolygon_wrap , polygon_dots )

    }

  }

  # return final object
  return( rval )

}

#' @rdname svylorenz
#' @export
svylorenz.svyrep.design <- function(formula , design, quantiles = seq(0,1,.1), empirical = FALSE, plot = TRUE, add = FALSE, curve.col = "red", ci = TRUE, alpha = .05, na.rm = FALSE , deff = FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

  # collect income data
  incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]

  # treat missing values
  if(na.rm){
    nas<-is.na(incvar)
    design<-design[!nas,]
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[[1]]
  }

  # collect sampling weights
  ws <- weights(design, "sampling")

  # treat remaining missing
  if ( any( is.na( incvar [ ws > 0 ] ) ) ) {
    variance <- as.matrix(NA)
    cis <- array( rbind(rep(NA, length(quantiles)),rep(NA, length(quantiles))), dim = c(2, length(quantiles)), dimnames = list( c( "(lower", "upper)" ), as.character(quantiles) ) )
    rval <- t( matrix( data = rep(NA, length(quantiles)), nrow = length(quantiles), dimnames = list( as.character( quantiles ), as.character(formula)[2] ) ) )
    rval <- list(quantiles = rval, CIs = cis)
    attr(rval, "SE") <- rep(NA, length(quantiles))
    class(rval) <- c( "cvyquantile" , "svyquantile" )
    return(rval)
  }

  # compute point estimates
  L.p <- as.numeric( lapply( quantiles, function( z ) CalcLorenz( incvar, ws , z ) ) )
  rval <- t( matrix( data = L.p, dimnames = list( as.character( quantiles ) ) ) )

  # collect analysis weights
  ww <- weights(design, "analysis")

  # compute replicates
  qq <- t( apply( ww, 2, function(wi) as.numeric( lapply( quantiles, function( z ) CalcLorenz( incvar, wi , z ) ) ) ) )

  # compute variance
  if ( any(is.na(qq))) {
    variance <- as.matrix(NA)
    cis <- array( rbind(rep(NA, length(quantiles)),rep(NA, length(quantiles))), dim = c(2, length(quantiles)), dimnames = list( c( "(lower", "upper)" ), as.character(quantiles) ) )
    rval <- t( matrix( data = rep(NA, length(quantiles)), nrow = length(quantiles), dimnames = list( as.character( quantiles ), as.character(formula)[2] ) ) )
    rval <- list(quantiles = rval, CIs = cis)
    attr(rval, "SE") <- rep(NA, length(quantiles))
    class(rval) <- c( "cvyquantile" , "svyquantile" )
    return(rval)
  }
  variance <- survey::svrVar( qq, design$scale, design$rscales, mse = design$mse, coef = rval )
  rownames( variance ) <- colnames( variance ) <- paste0( "L(" , quantiles , ")" )

  # set up additional estimates
  se <- sqrt(diag(variance))

  # compute deff
  if ( is.character(deff) || deff ) {

    # compute linearized
    lin.matrix<- do.call( cbind , lapply( quantiles, function( z ) CalcLorenz_IF( incvar, ws , z ) ) )

    # add names
    rownames( lin.matrix ) <- rownames( design$variables )
    colnames( lin.matrix ) <- paste0( "L(" , quantiles , ")" )

    # compute deff
    nobs <- length( design$pweights )
    npop <- sum( design$pweights )
    vsrs <- unclass( survey::svyvar( lin.matrix , design, na.rm = na.rm, return.replicates = FALSE, estimate.only = TRUE)) * npop^2/nobs
    if (deff != "replace") vsrs <- vsrs * (npop - nobs)/npop
    deff.estimate <- diag( variance ) / vsrs

  }

  # compute CIs
  CI.L <- as.numeric( L.p - se * qnorm( alpha, mean = 0, sd = 1, lower.tail = FALSE ) )
  CI.U <- as.numeric( L.p + se * qnorm( alpha, mean = 0, sd = 1, lower.tail = FALSE ) )
  cis <- structure(rbind(CI.L,CI.U), .Dim = c(2L, length(quantiles), 1L), .Dimnames = list(c("(lower", "upper)"), as.character(quantiles),  as.character(formula)[2]))

  # set up result object
  rval <- c( L.p )
  names( rval ) <- paste0( "L(" , quantiles , ")" )
  attr( rval , "var") <- variance
  attr( rval , "statistic") <- "lorenz"
  class( rval ) <- c( "cvystat" , "svrepstat" )
  if ( linearized ) attr(rval,"linearized") <- lin.matrix
  if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )

  # empirical curve
  if (empirical) empirical.lorenz <- emp.interp( incvar , ws )

  # plot curve
  if ( plot ) {

    # set up options
    plot_dots <- list( ... )

    # remove `deff` argument sent by svyby
    if( 'deff' %in% names( plot_dots ) ) plot_dots$deff <- NULL

    # star new plot
    if ( !add ) do.call( svylorenzplot_wrap , plot_dots )

    # close unused option
    if( any( c( 'xlim' , 'ylim' , 'col' ) %in% names( list( ... ) ) ) ) stop( "xlim=, ylim=, and col= parameters are fixed within `svylorenz`.  use curve.col= to change the line color" )

    # add reference line
    abline( 0 , 1 , ylim = c( 0 , 1 ) , plot_dots )

    # plot empirical curve
    if( empirical ) {
      lines_dots <- plot_dots
      lines_dots$x <- empirical.lorenz$lorenz_x
      lines_dots$y <- empirical.lorenz$lorenz_y
      lines_dots$col = curve.col
      do.call( svylorenzlines_wrap , lines_dots )
    }

    # add point estimates of ordinates
    points_dots <- plot_dots
    points_dots$x <- quantiles
    points_dots$y <- L.p
    points_dots$col <- curve.col
    do.call( svylorenzpoints_wrap , points_dots )

    # add confidence intervals
    if (ci) {
      X.Vec <- as.numeric( c(quantiles, tail(quantiles, 1), rev(quantiles), quantiles[1]) )
      Y.Vec <- as.numeric( c( CI.L, tail(CI.U, 1), rev(CI.U), CI.L[1] ) )

      polygon_dots <- plot_dots
      polygon_dots$x <- X.Vec
      polygon_dots$y <- Y.Vec
      polygon_dots$col <- adjustcolor( curve.col, alpha.f = .2)
      polygon_dots$border <- NA
      do.call( svylorenzpolygon_wrap , polygon_dots )

    }

  }

  # keep replicates
  if (return.replicates) {
    attr( qq , "scale") <- design$scale
    attr( qq , "rscales") <- design$rscales
    attr( qq , "mse") <- design$mse
    rval <- list( mean = rval , replicates = qq )
    class( rval ) <- c( "cvystat" , "svrepstat" )
  }

  # add design effect estimate
  if ( is.character(deff) || deff ) attr(rval,"deff") <- deff.estimate

  # return final result object
  return(rval)

}

#' @rdname svylorenz
#' @export
svylorenz.DBIsvydesign <- function (formula, design, ...) {

  design$variables <- getvars(formula, design$db$connection, design$db$tablename, updates = design$updates, subset = design$subset)
  NextMethod("svylorenz", design)

}


# function for point estimates
CalcLorenz <- function( x , w = rep( 1 , length( x ) ) , p ) {

  # filter observations
  x <- x[ w > 0 ]
  w <- w[ w > 0 ]

  # population totals
  # N <- sum( w )
  Y <- sum( x * w )

  # weighted partial sum
  Y.p <- wtd.psum( x, w , p )

  # lorenz ordinate estimate
  Y.p / Y

}

# function for linearized functions
CalcLorenz_IF <- function( x , w , pc ) {

  # filter observations
  x <- x[ w > 0 ]
  w <- w[ w > 0 ]

  # population size
  N <- sum( w )
  Y <- sum( x * w )
  mu <- Y/N

  # compute interpolated quantile
  quant <- smooth.quantile( x, w , pc )

  # compute lorenz curve ordinate
  L.p <- CalcLorenz( x , w , pc )

  # compute linearized variable
  # branch on special cases
  lin <- if ( pc %in% c(0,1) ) rep( 0 , length( w ) ) else ( ( x - quant ) * ( x <= quant ) + ( pc * quant ) - ( x * L.p ) ) / Y

  # add indices
  names( lin ) <- names( w )

  # retur linearized variable
  lin

}

# # partial sum (1st definition):
# wtd.psum <- function (x, q = .5, weights = NULL ) {
#   indices <- weights != 0
#   x <- x[indices]
#   weights <- weights[indices]
#
#   x_thres <- wtd.qtl(x = x, q = q, weights = weights )
#
#   return( sum( weights * x * 1 * (x <= x_thres) ) )
#
# }

# partial sum (2nd definition)
wtd.psum <- function (x, weights = rep( 1 , length( x ) ) , q = .5 ) {

  # filter observations
  x <- x[ weights > 0 ]
  weights <- weights[ weights > 0 ]

  # reorder
  ordx <- order( x )
  x <- x[ ordx ]
  weights <- weights[ ordx ]

  # population size
  N <- sum( weights )

  # intermediate partial sums
  wsum <- cumsum( weights )
  wsum_1 <- wsum - weights
  alpha_k <- wsum / N
  t_k <- ( (q * N) - wsum_1 ) / weights

  # H function
  H_fn <- function(x) ifelse( x < 0 , 0 , ifelse( (0 <= x) & (x < 1) , x , ifelse( x >= 1 , 1 , NA ) ) )

  # compute estimate
  sum( weights * x * H_fn( t_k ) )

}

# quantile function:
smooth.quantile <- function (x , weights = rep( 1 , length( x) ) , q = .5 ) {

  # filter observations
  x <- x[ weights > 0 ]
  weights <- weights[ weights > 0 ]

  # reorder
  ordx <- order( x )
  x <- x[ ordx ]
  weights <- weights[ ordx ]

  # population size
  N <- sum(weights)

  # partial sums
  x_1 <- c( min( x ) , x[ - length( x ) ] )
  wsum <- cumsum( weights )
  wsum_1 <- wsum - weights
  alpha_k <- wsum / N
  t_k <- ( (q * N) - wsum_1 ) / weights

  # interpolation
  smth.vec <- ( x_1 + ( x - x_1 ) * t_k )

  # find k
  k <- which( ( wsum_1 < ( q * N) ) & ( (q * N) <= wsum) )

  # return estimate
  return( smth.vec[k] )

}

# empirical interpolation
emp.interp <- function( x , w ) {

  # filter observations
  x <- x[ w > 0 ]
  w <- w[ w > 0 ]

  # reorder
  ordx <- order( x )
  x <- x[ ordx ]
  w <- w[ ordx ]

  # interpolate x and y
  data.frame(
    # lorenz_x = ( 2*cumsum(w) - w ) / ( 2*sum(w) ) , # where is this from?
    lorenz_x = cumsum( w ) / sum(w) ,
    lorenz_y = cumsum( w * x ) / sum( w * x ) ,
    stringsAsFactors = FALSE )

}
