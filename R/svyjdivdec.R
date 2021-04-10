#' J-Divergence Decomposition
#'
#' Estimates the group decomposition of the generalized entropy index
#'
#' @param formula a formula specifying the income variable
#' @param subgroup a formula specifying the group variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param na.rm Should cases with missing values be dropped? Observations containing missing values in income or group variables will be dropped.
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
#' @param influence Should a matrix of (weighted) influence functions be returned? (for compatibility with \code{\link[survey]{survey}})
#' @param return.replicates Return the replicate estimates?
#' @param ... future expansion
#'
#' @details you must run the \code{convey_prep} function on your survey design object immediately after creating it with the \code{svydesign} or \code{svrepdesign} function.
#'
#' This measure only allows for strictly positive variables.
#'
#' @return Object of class "\code{cvydstat}", which are vectors with a "\code{var}" attribute giving the variance-covariance matrix and a "\code{statistic}" attribute giving the name of the statistic.
#'
#' @author Guilherme Jacob, Djalma Pessoa and Anthony Damico
#'
#' @note This function is experimental and is subject to change in later versions.
#'
#' @seealso \code{\link{svyjdiv}}
#'
#' @references Anthony F. Shorrocks (1984). Inequality decomposition by population subgroups.
#' \emph{Econometrica}, v. 52, n. 6, 1984, pp. 1369-1385.
#' URL \url{http://www.jstor.org/stable/1913511}.
#'
#' Nicholas Rohde (2016). J-divergence measurements of economic inequality.
#' \emph{J. R. Statist. Soc. A}, v. 179, Part 3 (2016), pp. 847-870.
#' URL \url{http://onlinelibrary.wiley.com/doi/10.1111/rssa.12153/abstract}.
#'
#' Martin Biewen and Stephen Jenkins (2002). Estimation of Generalized Entropy
#' and Atkinson Inequality Indices from Complex Survey Data. \emph{DIW Discussion Papers},
#' No.345,
#' URL \url{https://www.diw.de/documents/publikationen/73/diw_01.c.40394.de/dp345.pdf}.
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
#' des_eusilc <- convey_prep(des_eusilc)
#'
#' # replicate-weighted design
#' des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )
#' des_eusilc_rep <- convey_prep(des_eusilc_rep)
#'
#' # linearized design
#' svyjdivdec( ~eqincome , ~rb090 , subset(des_eusilc, eqincome > 0) )
#'
#' # replicate-weighted design
#' svyjdivdec( ~eqincome , ~rb090 , subset(des_eusilc_rep, eqincome > 0) )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' sub_des_eusilc <- subset(des_eusilc, py010n > 0 | is.na(py010n) )
#' svyjdivdec( ~py010n , ~rb090 , sub_des_eusilc )
#' svyjdivdec( ~py010n , ~rb090 , sub_des_eusilc , na.rm = TRUE )
#'
#' # replicate-weighted design using a variable with missings
#' sub_des_eusilc_rep <- subset(des_eusilc_rep, py010n > 0 | is.na(py010n) )
#' svyjdivdec( ~py010n , ~rb090 , sub_des_eusilc_rep )
#' svyjdivdec( ~py010n , ~rb090 , sub_des_eusilc_rep , na.rm = TRUE )
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
#' # database-backed linearized design
#' svyjdivdec( ~eqincome , ~rb090 , subset(dbd_eusilc, eqincome > 0) )
#'
#' # database-backed linearized design using a variable with missings
#' sub_dbd_eusilc <- subset(dbd_eusilc, py010n > 0 | is.na(py010n) )
#' svyjdivdec( ~py010n , ~rb090 , sub_dbd_eusilc )
#' svyjdivdec( ~py010n , ~rb090 , sub_dbd_eusilc , na.rm = TRUE )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svyjdivdec <-
  function( formula, subgroup, design, ...) {

    if( length( attr( terms.formula( formula ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `formula=` argument" )

    if( length( attr( terms.formula( subgroup ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `subgroup=` argument" )

    UseMethod("svyjdivdec", design)

  }

#' @rdname svyjdivdec
#' @export
svyjdivdec.survey.design <-
  function ( formula, subgroup, design, na.rm = FALSE, deff = FALSE , linearized = FALSE , influence = FALSE , ... ) {

    # collect income and subgroup data
    incvar <- model.frame( formula , design$variables , na.action = na.pass )[,]
    grpvar <- model.frame( subgroup , design$variables , na.action = na.pass , drop.unused.levels = TRUE )[,]

    # check format
    if ( class(grpvar) == "labelled" ) {
      stop( "This function does not support 'labelled' variables. Try factor().")
    }

    # treat missing values
    if (na.rm) {
      nas <- ( is.na(incvar) | is.na(grpvar ) )
      design <- design[nas == 0, ]
      incvar <- model.frame( formula , design$variables , na.action = na.pass )[,]
      grpvar <- model.frame( subgroup , design$variables , na.action = na.pass , drop.unused.levels = TRUE )[,]
    }

    # collect weights
    w <- 1/design$prob
    names( w ) <- rownames( design$variables )[ w>0 ]

    # test for positive incommes
    if ( any( incvar[ w != 0 ] <= 0, na.rm = TRUE) ) stop( "The J-divergence index is defined for strictly positive incomes only." )

    # add interaction
    grpvar <- interaction( grpvar )

    # total
    total.jdiv <- CalcJDiv( incvar, w  )

    # compute linearized function
    total.jdivlin <- CalcJDiv_IF( incvar, w )

    # create matrix of group-specific weights
    wg <- sapply( levels(grpvar) , function(z) ifelse( grpvar == z , w , 0 ) )

    # calculate group-specific GEI_0 and linearized functions
    grp.gei0 <- lapply( colnames( wg )  , function( this.group ) {
      wi <- wg[ , this.group ]
      statobj <- list(
        value = CalcGEI( x = incvar, weights = wi , epsilon = 0 ) ,
        lin = CalcGEI_IF( x = incvar, weights = wi , epsilon = 0 ) )
      lin2 <- rep( 0 , length( wi ) )
      lin2[ wi > 0 ] <- statobj$lin
      statobj$lin <- lin2
      statobj
    } )
    names( grp.gei0 ) <- colnames( wg )

    # calculate group-specific GEI_1 and linearized functions
    grp.gei1 <- lapply( colnames( wg )  , function( this.group ) {
      wi <- wg[ , this.group ]
      statobj <- list(
        value = CalcGEI( x = incvar, weights = wi , epsilon = 1 ) ,
        lin = CalcGEI_IF( x = incvar, weights = wi , epsilon = 1 ) )
      lin2 <- rep( 0 , length( wi ) )
      lin2[ wi > 0 ] <- statobj$lin
      statobj$lin <- lin2
      statobj
    } )
    names( grp.gei1 ) <- colnames( wg )

    # calculate gei0 within component weight
    grp.gei0.wgt <- lapply( colnames( wg ) , function(i) {

      wi <- wg[,i]
      this.linformula <- quote( ( N.g / N ) )
      contrastinf(
        this.linformula ,
        list( Y.g = list( value = sum( incvar * wi , na.rm = TRUE ) , lin = incvar * ( wi > 0 ) ) ,
              Y = list( value = sum( incvar * w , na.rm = TRUE ) , lin = incvar * ( w > 0 ) ) ,
              N.g = list( value = sum( wi , na.rm = TRUE ) , lin = ( wi > 0 ) ) ,
              N = list( value = sum( w , na.rm = TRUE ) , lin = ( w > 0 ) ) ) )

    } )
    names( grp.gei0.wgt ) <- colnames( wg )

    # calculate gei1 within component weight
    grp.gei1.wgt <- lapply( colnames( wg ) , function(i) {

      wi <- wg[,i]
      this.linformula <- quote( ( Y.g / Y ) )
      contrastinf(
        this.linformula ,
        list( Y.g = list( value = sum( incvar * wi , na.rm = TRUE ) , lin = incvar * ( wi > 0 ) ) ,
              Y = list( value = sum( incvar * w , na.rm = TRUE ) , lin = incvar * ( w > 0 ) ) ,
              N.g = list( value = sum( wi , na.rm = TRUE ) , lin = ( wi > 0 ) ) ,
              N = list( value = sum( w , na.rm = TRUE ) , lin = ( w > 0 ) ) ) )

    } )
    names( grp.gei1.wgt ) <- colnames( wg )

    # compute combined within components
    gei0.within.components <-
      list(
        value = sapply( grp.gei0.wgt , `[[` , "value" ) * sapply( grp.gei0 , `[[` , "value" ) ,
        lin = sweep( sapply( grp.gei0 , `[[` , "lin" ) , 2 , sapply( grp.gei0.wgt , `[[` , "value" ) , "*" ) +
          sweep( sapply( grp.gei0.wgt , `[[` , "lin" ) , 2 , sapply( grp.gei0 , `[[` , "value" ) , "*" ) )
    gei1.within.components <-
      list(
        value = sapply( grp.gei1.wgt , `[[` , "value" ) * sapply( grp.gei1 , `[[` , "value" ) ,
        lin = sweep( sapply( grp.gei1 , `[[` , "lin" ) , 2 , sapply( grp.gei1.wgt , `[[` , "value" ) , "*" ) +
          sweep( sapply( grp.gei1.wgt , `[[` , "lin" ) , 2 , sapply( grp.gei1 , `[[` , "value" ) , "*" ) )
    within.jdiv <- sum( gei0.within.components$value + gei1.within.components$value )
    within.jdivlin <- rowSums( gei0.within.components$lin + gei1.within.components$lin )[ w > 0 ]

    # compute between
    between.jdiv <- total.jdiv - within.jdiv
    between.jdivlin <- total.jdivlin - within.jdivlin

    # create vector of estimates
    estimates <- c( total.jdiv , within.jdiv , between.jdiv )
    names( estimates ) <- c( "total", "within", "between" )

    # create linearized matrix
    lin.matrix <- matrix( data = c( total.jdivlin , within.jdivlin , between.jdivlin ), ncol = 3, dimnames = list( NULL, c( "total", "within", "between" ) ) )
    rm( total.jdivlin , within.jdivlin , between.jdivlin )

    # ensure length
    if ( nrow( lin.matrix ) != length( design$prob ) ) {
      tmplin <- matrix( 0 , nrow = nrow( design$variables ) , ncol = ncol( lin.matrix ) )
      tmplin[ w > 0 , ] <- lin.matrix
      lin.matrix <- tmplin ; rm( tmplin )
      colnames( lin.matrix ) <- c( "total", "within", "between" )
    }
    rownames( lin.matrix ) <- rownames( design$variables )

    # compute variance
    variance <- survey::svyrecvar( sweep( lin.matrix , 1 , 1/design$prob , "*" ) , design$cluster, design$strata, design$fpc, postStrata = design$postStrata )
    variance[ which( is.nan( variance ) ) ] <- NA

    # compute deff
    if ( is.character(deff) || deff || linearized ) {
      nobs <- sum( weights( design ) != 0 )
      npop <- sum( weights( design ) )
      if (deff == "replace") vsrs <- survey::svyvar( lin.matrix , design, na.rm = na.rm) * npop^2/nobs
      else vsrs <- survey::svyvar( lin.matrix , design , na.rm = na.rm ) * npop^2 * (npop - nobs)/(npop * nobs)
      deff.estimate <- variance/vsrs
    }

    # keep necessary linearized functions
    lin.matrix <- lin.matrix[ 1/design$prob > 0 , ]

    # build result object
    rval <- c( estimates )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "jdiv decomposition"
    attr(rval,"group")<- as.character( subgroup )[[2]]
    class(rval) <- c( "cvystat" , "svystat" )
    if ( linearized ) attr(rval,"linearized") <- lin.matrix
    if ( influence )  attr( rval , "influence" )  <- sweep( lin.matrix , 1 , design$prob[ is.finite( design$prob ) ] , "/" )
    if ( linearized | influence ) attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }


#' @rdname svyjdivdec
#' @export
svyjdivdec.svyrep.design <-
  function( formula, subgroup, design, na.rm=FALSE, deff = FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

    # collect income and group data
    incvar <- model.frame( formula, design$variables, na.action = na.pass )[,]
    grpvar <- model.frame( subgroup, design$variables, na.action = na.pass , drop.unused.levels = TRUE )[,]

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar) | is.na(grpvar)
      design<-design[!nas,]
      incvar <- model.frame( formula, design$variables, na.action = na.pass )[,]
      grpvar <- model.frame( subgroup, design$variables, na.action = na.pass , drop.unused.levels = TRUE )[,]
    }

    # collect sampling weights
    ws <- weights(design, "sampling")

    # check for positive values
    if ( any( incvar[ ws != 0 ] <= 0, na.rm = TRUE) ) stop( "The J-divergence index is defined for strictly positive incomes." )

    # apply interaction
    grpmat <- model.matrix( update( subgroup , ~.+0 ) , design$variables )

    # compute point estimates
    N <- sum( ws )
    Y <- sum( incvar * ws )
    total.jdiv <- CalcJDiv( incvar, ws )
    grp.gei0   <- apply( grpmat , 2 , function( ind ) CalcGEI( incvar, ws * ind , epsilon = 0 ) )
    grp.gei1   <- apply( grpmat , 2 , function( ind ) CalcGEI( incvar, ws * ind , epsilon = 1 ) )
    grp.pshare <- colSums( sweep( grpmat , 1 , ws , "*" ) ) / N
    grp.sshare <- colSums( sweep( grpmat , 1 , ws*incvar , "*" ) ) / Y
    within.jdiv <- sum( grp.gei0 * grp.pshare + grp.gei1 * grp.sshare )
    between.jdiv <- total.jdiv - within.jdiv
    estimates <- c( total.jdiv , within.jdiv , between.jdiv )

    ### variance calculation

    # collect resampling weights
    wf <- weights( design , "analysis" )

    # compute replicates
    qq <- apply( wf , 2 , function( wi ) {

      N <- sum( wi )
      Y <- sum( incvar * wi )
      total.jdiv <- CalcJDiv( incvar, wi )
      grp.gei0   <- apply( grpmat , 2 , function( ind ) CalcGEI( incvar, wi * ind , epsilon = 0 ) )
      grp.gei1   <- apply( grpmat , 2 , function( ind ) CalcGEI( incvar, wi * ind , epsilon = 1 ) )
      grp.pshare <- colSums( sweep( grpmat , 1 , wi , "*" ) ) / N
      grp.sshare <- colSums( sweep( grpmat , 1 , wi*incvar , "*" ) ) / Y
      within.jdiv <- sum( grp.gei0 * grp.pshare + grp.gei1 * grp.sshare )
      between.jdiv <- total.jdiv - within.jdiv
      estimates <- c( total.jdiv , within.jdiv , between.jdiv )

    } )
    qq <- t( qq )
    dimnames( qq ) <- list( NULL, c( "total", "within", "between" ) )

    # compute variance
    if ( anyNA( qq ) ) {
      variance <- diag( estimates )
      variance[,] <- NA
    } else {
      variance <- survey::svrVar( qq , design$scale, design$rscales, mse = design$mse, coef = estimates )
    }

    # deff calculation
    if ( is.character(deff) || deff || linearized ) {

      ### compute linearized function matrix

      # compute linearized function
      total.jdivlin <- CalcJDiv_IF( incvar, ws )

      # add interaction
      grpvar <- interaction( grpvar )

      # create matrix of group-specific weights
      wg <- sapply( levels(grpvar) , function(z) ifelse( grpvar == z , ws , 0 ) )

      # calculate group-specific GEI_0 and linearized functions
      grp.gei0 <- lapply( colnames( wg )  , function( this.group ) {
        wi <- wg[ , this.group ]
        statobj <- list(
          value = CalcGEI( x = incvar, weights = wi , epsilon = 0 ) ,
          lin = CalcGEI_IF( x = incvar, weights = wi , epsilon = 0 ) )
        lin2 <- rep( 0 , length( wi ) )
        lin2[ wi > 0 ] <- statobj$lin
        statobj$lin <- lin2
        statobj
      } )
      names( grp.gei0 ) <- colnames( wg )

      # calculate group-specific GEI_1 and linearized functions
      grp.gei1 <- lapply( colnames( wg )  , function( this.group ) {
        wi <- wg[ , this.group ]
        statobj <- list(
          value = CalcGEI( x = incvar, weights = wi , epsilon = 1 ) ,
          lin = CalcGEI_IF( x = incvar, weights = wi , epsilon = 1 ) )
        lin2 <- rep( 0 , length( wi ) )
        lin2[ wi > 0 ] <- statobj$lin
        statobj$lin <- lin2
        statobj
      } )
      names( grp.gei1 ) <- colnames( wg )

      # calculate gei0 within component weight
      grp.gei0.wgt <- lapply( colnames( wg ) , function(i) {

        wi <- wg[,i]
        this.linformula <- quote( ( N.g / N ) )
        contrastinf(
          this.linformula ,
          list( Y.g = list( value = sum( incvar * wi , na.rm = TRUE ) , lin = incvar * ( wi > 0 ) ) ,
                Y = list( value = sum( incvar * ws , na.rm = TRUE ) , lin = incvar * ( ws > 0 ) ) ,
                N.g = list( value = sum( wi , na.rm = TRUE ) , lin = ( wi > 0 ) ) ,
                N = list( value = sum( ws , na.rm = TRUE ) , lin = ( ws > 0 ) ) ) )

      } )
      names( grp.gei0.wgt ) <- colnames( wg )

      # calculate gei1 within component weight
      grp.gei1.wgt <- lapply( colnames( wg ) , function(i) {

        wi <- wg[,i]
        this.linformula <- quote( ( Y.g / Y ) )
        contrastinf(
          this.linformula ,
          list( Y.g = list( value = sum( incvar * wi , na.rm = TRUE ) , lin = incvar * ( wi > 0 ) ) ,
                Y = list( value = sum( incvar * ws , na.rm = TRUE ) , lin = incvar * ( ws > 0 ) ) ,
                N.g = list( value = sum( wi , na.rm = TRUE ) , lin = ( wi > 0 ) ) ,
                N = list( value = sum( ws , na.rm = TRUE ) , lin = ( ws > 0 ) ) ) )

      } )
      names( grp.gei1.wgt ) <- colnames( wg )

      # compute combined within components
      gei0.within.components <-
        list(
          value = sapply( grp.gei0.wgt , `[[` , "value" ) * sapply( grp.gei0 , `[[` , "value" ) ,
          lin = sweep( sapply( grp.gei0 , `[[` , "lin" ) , 2 , sapply( grp.gei0.wgt , `[[` , "value" ) , "*" ) +
            sweep( sapply( grp.gei0.wgt , `[[` , "lin" ) , 2 , sapply( grp.gei0 , `[[` , "value" ) , "*" ) )
      gei1.within.components <-
        list(
          value = sapply( grp.gei1.wgt , `[[` , "value" ) * sapply( grp.gei1 , `[[` , "value" ) ,
          lin = sweep( sapply( grp.gei1 , `[[` , "lin" ) , 2 , sapply( grp.gei1.wgt , `[[` , "value" ) , "*" ) +
            sweep( sapply( grp.gei1.wgt , `[[` , "lin" ) , 2 , sapply( grp.gei1 , `[[` , "value" ) , "*" ) )
      within.jdiv <- sum( gei0.within.components$value + gei1.within.components$value )
      within.jdivlin <- rowSums( gei0.within.components$lin + gei1.within.components$lin )[ ws > 0 ]

      # compute between
      between.jdiv <- total.jdiv - within.jdiv
      between.jdivlin <- total.jdivlin - within.jdivlin

      # create matrix
      lin.matrix <- matrix( data = c( total.jdivlin , within.jdivlin , between.jdivlin ), ncol = 3, dimnames = list( NULL, c( "total", "within", "between" ) ) )
      rownames( lin.matrix ) <- rownames( design$variables )

      ## deff calculation

      # compute deff
      nobs <- length( design$pweights )
      npop <- sum( design$pweights )
      vsrs <- unclass( survey::svyvar( lin.matrix , design, na.rm = na.rm, return.replicates = FALSE, estimate.only = TRUE)) * npop^2/nobs
      if (deff != "replace") vsrs <- vsrs * (npop - nobs)/npop
      deff.estimate <- variance / vsrs

    }

    # build result object
    rval <- estimates
    names( rval ) <- c( "total", "within", "between" )
    attr(rval, "var") <- variance
    attr(rval, "statistic") <- "j-divergence decomposition"
    attr(rval,"group")<- as.character( subgroup )[[2]]
    class(rval) <- c( "cvystat" , "svrepstat" , "svystat" )
    if ( linearized ) attr( rval , "linearized" ) <- lin.matrix
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
    if ( is.character(deff) || deff) attr( rval , "deff" ) <- deff.estimate

    # return value
    rval

  }


#' @rdname svyjdivdec
#' @export
svyjdivdec.DBIsvydesign <-
  function (formula, subgroup, design, ...) {

    design$variables <-
      cbind(
        getvars(formula, design$db$connection,design$db$tablename, updates = design$updates, subset = design$subset),
        getvars(subgroup, design$db$connection, design$db$tablename,updates = design$updates, subset = design$subset)
      )

    NextMethod("svyjdivdec", design)

  }
