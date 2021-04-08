#' Generalized Entropy Index Decomposition
#'
#' Estimates the group decomposition of the generalized entropy index
#'
#' @param formula a formula specifying the income variable
#' @param subgroup a formula specifying the group variable
#' @param design a design object of class \code{survey.design} or class \code{svyrep.design} from the \code{survey} library.
#' @param epsilon a parameter that determines the sensivity towards inequality in the top of the distribution. Defaults to epsilon = 1.
#' @param na.rm Should cases with missing values be dropped? Observations containing missing values in income or group variables will be dropped.
#' @param deff Return the design effect (see \code{survey::svymean})
#' @param linearized Should a matrix of linearized variables be returned
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
#' @seealso \code{\link{svygei}}
#'
#' @references Anthony F. Shorrocks (1984). Inequality decomposition groups population subgroups.
#' \emph{Econometrica}, v. 52, n. 6, 1984, pp. 1369-1385.
#' URL \url{http://www.jstor.org/stable/1913511}.
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
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc, eqincome > 0 ) , epsilon = 0 )
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc, eqincome > 0 ) , epsilon = .5 )
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc, eqincome > 0 ) , epsilon = 1 )
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc, eqincome > 0 ) , epsilon = 2 )
#'
#' # replicate-weighted design
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc_rep, eqincome > 0 ) , epsilon = 0 )
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc_rep, eqincome > 0 ) , epsilon = .5 )
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc_rep, eqincome > 0 ) , epsilon = 1 )
#' svygeidec( ~eqincome , ~rb090 , subset( des_eusilc_rep, eqincome > 0 ) , epsilon = 2 )
#'
#' \dontrun{
#'
#' # linearized design using a variable with missings
#' sub_des_eusilc <- subset(des_eusilc, py010n > 0 | is.na(py010n) )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc , epsilon = 0 )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc , epsilon = 0, na.rm = TRUE )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc , epsilon = 1 )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc , epsilon = 1, na.rm = TRUE )
#'
#' # replicate-weighted design using a variable with missings
#' sub_des_eusilc_rep <- subset(des_eusilc_rep, py010n > 0 | is.na(py010n) )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc_rep , epsilon = 0 )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc_rep , epsilon = 0, na.rm = TRUE )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc_rep , epsilon = 1 )
#' svygeidec( ~py010n , ~rb090 , sub_des_eusilc_rep , epsilon = 1, na.rm = TRUE )
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
#' svygeidec( ~eqincome , ~rb090 , subset(dbd_eusilc, eqincome > 0) , epsilon = 0 )
#' svygeidec( ~eqincome , ~rb090 , subset(dbd_eusilc, eqincome > 0) , epsilon = .5 )
#' svygeidec( ~eqincome , ~rb090 , subset(dbd_eusilc, eqincome > 0) , epsilon = 1 )
#' svygeidec( ~eqincome , ~rb090 , subset(dbd_eusilc, eqincome > 0) , epsilon = 2 )
#'
#' # database-backed linearized design using a variable with missings
#' sub_dbd_eusilc <- subset(dbd_eusilc, py010n > 0 | is.na(py010n) )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = 0 )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = 0, na.rm = TRUE )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = .5 )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = .5, na.rm = TRUE )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = 1 )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = 1, na.rm = TRUE )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = 2 )
#' svygeidec( ~py010n , ~rb090 , sub_dbd_eusilc , epsilon = 2, na.rm = TRUE )
#'
#' dbRemoveTable( conn , 'eusilc' )
#'
#' dbDisconnect( conn , shutdown = TRUE )
#'
#' }
#'
#' @export
svygeidec <-
  function( formula, subgroup, design,  ...) {

    if( length( attr( terms.formula( subgroup ) , "term.labels" ) ) > 1 ) stop( "convey package functions currently only support one variable in the `subgroup=` argument" )

    # if( 'epsilon' %in% names( list(...) ) && list(...)[["epsilon"]] < 0 ) stop( "epsilon= cannot be negative." )

    UseMethod("svygeidec", design)

  }

#' @rdname svygeidec
#' @export
svygeidec.survey.design <-
  function ( formula, subgroup, design, epsilon = 1, na.rm = FALSE, deff = FALSE , linearized = FALSE , ... ) {

    # collect data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[,]
    grpvar <- model.frame( subgroup, design$variables, na.action = na.pass , drop.unused.levels = TRUE)[,]

    # check types
    if ( class(grpvar) == "labelled" ) {
      stop( "This function does not support 'labelled' variables. Try factor().")
    }

    # treat missing values
    if (na.rm) {
      nas <- ( is.na(incvar) | is.na(grpvar ) )
      design <- design[!nas, ]
      incvar <- model.frame(formula, design$variables, na.action = na.pass)[,]
      grpvar <- model.frame( subgroup, design$variables, na.action = na.pass , drop.unused.levels = TRUE)[,]
    }

    # collect sampling weights
    w <- 1/design$prob

    # check for strictly positive incomes
    if ( any(incvar[w != 0] <= 0,na.rm = TRUE) ) stop( "The GEI indices are defined for strictly positive variables only.\nNegative and zero values not allowed." )

    # add interaction
    grpvar <- interaction( grpvar )

    # total
    ttl.gei <- CalcGEI( x = incvar, weights = w, epsilon = epsilon )

    # compute linearized function
    ttl.lin <- CalcGEI_IF( x = incvar, weights = w, epsilon = epsilon )

    # treat domain
    ttl.lin <- ttl.lin[ pmatch( names( w ) , names( ttl.lin ) ) ]
    names( ttl.lin ) <- names( w )
    ttl.lin[ w<=0 ] <- 0

    # create matrix of group-specific weights
    wg <- sapply( levels(grpvar) , function(z) ifelse( grpvar == z , w , 0 ) )

    # calculate group-specific GEI and linearized functions
    grp.gei <- lapply( colnames( wg )  , function( this.group ) {
      wi <- wg[ , this.group ]
      statobj <- list(
        value = CalcGEI( x = incvar, weights = wi , epsilon = epsilon ) ,
        lin = CalcGEI_IF( x = incvar, weights = wi , epsilon = epsilon ) )
      lin2 <- rep( 0 , length( wi ) )
      lin2[ wi > 0 ] <- statobj$lin
      statobj$lin <- lin2
      statobj
    } )
    names( grp.gei ) <- colnames( wg )

    # calculate within component weight
    grp.gei.wgt <- lapply( colnames( wg ) , function(i) {

      wi <- wg[,i]

      if ( epsilon == 0 ) {
        this.linformula <- quote( ( N.g / N ) )
      } else if ( epsilon == 1 ) {
        this.linformula <- quote( ( Y.g / Y ) )
      } else {
        this.linformula <- substitute( quote( ( ( Y.g / Y )^epsilon ) * ( ( N.g / N )^( 1 - epsilon ) ) ) , list( epsilon = epsilon ) )
        this.linformula <- eval( this.linformula )
      }

      contrastinf(
        this.linformula ,
        list( Y.g = list( value = sum( incvar * wi , na.rm = TRUE ) , lin = incvar * ( wi > 0 ) ) ,
              Y = list( value = sum( incvar * w , na.rm = TRUE ) , lin = incvar * ( w > 0 ) ) ,
              N.g = list( value = sum( wi , na.rm = TRUE ) , lin = ( wi > 0 ) ) ,
              N = list( value = sum( w , na.rm = TRUE ) , lin = ( w > 0 ) ) ) )

    } )
    names( grp.gei.wgt ) <- colnames( wg )

    # calculate within component weight
    gei.within.components <-
      list(
        value = sapply( grp.gei.wgt , `[[` , "value" ) * sapply( grp.gei , `[[` , "value" ) ,
        lin = sweep( sapply( grp.gei , `[[` , "lin" ) , 2 , sapply( grp.gei.wgt , `[[` , "value" ) , "*" ) +
          sweep( sapply( grp.gei.wgt , `[[` , "lin" ) , 2 , sapply( grp.gei , `[[` , "value" ) , "*" ) )

    # compute within component
    wtn.gei <- sum( gei.within.components$value )
    within.lin <- rowSums( gei.within.components$lin )

    # between (residual)
    btw.gei <- ttl.gei - wtn.gei
    between.lin <- ttl.lin - within.lin

    # create vector of estimates
    estimates <- c( ttl.gei, wtn.gei, btw.gei )
    names( estimates ) <- c( "total", "within", "between" )

    # create linearized matrix
    lin.matrix <- matrix( data = c(ttl.lin, within.lin, between.lin), ncol = 3, dimnames = list( NULL, c( "total", "within", "between" ) ) )
    rm(ttl.lin, within.lin, between.lin)

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
    if ( is.character(deff) || deff ) {
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
    attr(rval,"group") <- as.character( subgroup )[[2]]
    attr(rval,"epsilon") <- epsilon
    class(rval) <- c( "cvystat" , "svystat" )
    if ( linearized ) attr(rval,"linearized") <- lin.matrix
    if ( linearized ) attr( rval , "index" ) <- as.numeric( rownames( lin.matrix ) )
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate
    rval

  }


#' @rdname svygeidec
#' @export
svygeidec.svyrep.design <-
  function( formula, subgroup, design, epsilon = 1, na.rm=FALSE, deff = FALSE , linearized = FALSE , return.replicates = FALSE , ...) {

    # between inequality function
    fun.btw.gei <- function( y , w , grp , epsilon ) {

      y <- y[ w>0 ]
      grp <- grp[ w>0 ]
      w <- w[ w>0 ]

      N <- sum( w )
      Y <- sum( y * w )
      mu <- Y/N
      N.g <- tapply( w , grp , sum , na.rm = TRUE )
      Y.g <- tapply( w * y , grp , sum , na.rm = TRUE )
      mu.g <- Y.g / N.g

      if ( epsilon == 0 ) {
        estimate <- - sum( N.g * log( mu.g / mu ) ) / N
      } else if ( epsilon == 1) {
        estimate <- sum( ( Y.g / Y ) * log( mu.g / mu ) )
      } else {
        estimate <- sum( ( N.g / N ) * ( ( mu.g / mu )^epsilon - 1 ) ) / ( epsilon^2 - epsilon )
      }
      estimate
    }

    # within inequality function
    fun.wtn.gei <- function( y , w , grp , epsilon ) {

      y <- y[ w>0 ]
      grp <- grp[ w>0 ]
      w <- w[ w>0 ]

      N <- sum( w )
      Y <- sum( y * w )
      mu <- Y/N
      N.g <- tapply( w , grp , sum , na.rm = TRUE )
      Y.g <- tapply( w * y , grp , sum , na.rm = TRUE )
      s.g <- Y.g / Y
      p.g <- N.g / N
      gei.g <- sapply( levels( grp ) , function( grpv ) CalcGEI( y , ifelse( grp == grpv , w , 0 ) , epsilon = epsilon ) )

      if ( epsilon == 0 ) {
        estimate <- sum( p.g * gei.g )
      } else if ( epsilon == 1) {
        estimate <- sum( s.g * gei.g )
      } else {
        estimate <- sum( (s.g^epsilon) * (p.g^(1-epsilon)) * gei.g )
      }

      estimate

    }

    # collect data
    incvar <- model.frame(formula, design$variables, na.action = na.pass)[,]
    grpvar <- model.frame( subgroup, design$variables, na.action = na.pass , drop.unused.levels = TRUE)[,]

    # check types
    if ( class(grpvar) == "labelled" ) {
      stop( "This function does not support 'labelled' variables. Try factor().")
    }

    # treat missing values
    if(na.rm){
      nas<-is.na(incvar) | is.na(grpvar)
      design<-design[!nas,]
      df <- model.frame(design)
      incvar <- model.frame(formula, design$variables, na.action = na.pass)[,]
      grpvar <- model.frame( subgroup, design$variables, na.action = na.pass , drop.unused.levels = TRUE)[,]
    }

    # collect samling weights
    ws <- weights(design, "sampling")

    # check for strictly positive incomes
    if ( any( incvar[ws != 0] <= 0 , na.rm = TRUE ) ) stop( "The GEI indices are defined for strictly positive variables only.\nNegative and zero values not allowed." )

    # create interaction
    grpvar <- interaction( grpvar )

    # collect analysis weights
    ww <- weights(design, "analysis")
    qq.ttl.gei <- apply(ww, 2, function(wi) CalcGEI(incvar, wi, epsilon = epsilon) )

    ### point estimates

    # total inequality
    ttl.gei <- CalcGEI( x = incvar, weights = ws , epsilon = epsilon )
    btw.gei <- fun.btw.gei( incvar , ws , grpvar , epsilon )
    # wtn.gei <- fun.wtn.gei( incvar , ws , grpvar , epsilon )
    wtn.gei <- ttl.gei - btw.gei
    estimates <- c( ttl.gei, wtn.gei, btw.gei )
    stopifnot( all.equal ( fun.btw.gei( incvar , ws , grpvar , epsilon ) , btw.gei , tolerance = 1e-10 ) )

    ### variance estimation

    # create matrix of replicates
    qq <- apply( ww , 2 , function( wi ) {
      ttl.rep <- CalcGEI( incvar, wi , epsilon )
      btw.rep <- fun.btw.gei( incvar , wi , grpvar , epsilon )
      wtn.rep <- ttl.rep - btw.rep
      c( ttl.rep , wtn.rep , btw.rep )
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

    # compute deff
    if ( is.character(deff) || deff || linearized ) {

      ### compute linearization

      # compute linearized function
      ttl.lin <- CalcGEI_IF( x = incvar, weights = ws, epsilon = epsilon )

      # treat domain
      ttl.lin <- ttl.lin[ pmatch( names( ws ) , names( ttl.lin ) ) ]
      names( ttl.lin ) <- names( ws )
      ttl.lin[ ws<=0 ] <- 0

      # create matrix of group-specific weights
      wg <- sapply( levels(grpvar) , function(z) ifelse( grpvar == z , ws , 0 ) )

      # calculate group-specific GEI and linearized functions
      grp.gei <- lapply( colnames( wg )  , function( this.group ) {
        wi <- wg[ , this.group ]
        statobj <- list(
          value = CalcGEI( x = incvar, weights = wi , epsilon = epsilon ) ,
          lin = CalcGEI_IF( x = incvar, weights = wi , epsilon = epsilon ) )
        lin2 <- rep( 0 , length( wi ) )
        lin2[ wi > 0 ] <- statobj$lin
        statobj$lin <- lin2
        statobj
      } )
      names( grp.gei ) <- colnames( wg )

      # calculate within component weight
      grp.gei.wgt <- lapply( colnames( wg ) , function(i) {

        wi <- wg[,i]

        if ( epsilon == 0 ) {
          this.linformula <- quote( ( N.g / N ) )
        } else if ( epsilon == 1 ) {
          this.linformula <- quote( ( Y.g / Y ) )
        } else {
          this.linformula <- substitute( quote( ( ( Y.g / Y )^epsilon ) * ( ( N.g / N )^( 1 - epsilon ) ) ) , list( epsilon = epsilon ) )
          this.linformula <- eval( this.linformula )
        }

        contrastinf(
          this.linformula ,
          list( Y.g = list( value = sum( incvar * wi , na.rm = TRUE ) , lin = incvar * ( wi > 0 ) ) ,
                Y = list( value = sum( incvar * ws , na.rm = TRUE ) , lin = incvar * ( ws > 0 ) ) ,
                N.g = list( value = sum( wi , na.rm = TRUE ) , lin = ( wi > 0 ) ) ,
                N = list( value = sum( ws , na.rm = TRUE ) , lin = ( ws > 0 ) ) ) )

      } )
      names( grp.gei.wgt ) <- colnames( wg )

      # calculate within component weight
      gei.within.components <-
        list(
          value = sapply( grp.gei.wgt , `[[` , "value" ) * sapply( grp.gei , `[[` , "value" ) ,
          lin = sweep( sapply( grp.gei , `[[` , "lin" ) , 2 , sapply( grp.gei.wgt , `[[` , "value" ) , "*" ) +
            sweep( sapply( grp.gei.wgt , `[[` , "lin" ) , 2 , sapply( grp.gei , `[[` , "value" ) , "*" ) )

      # compute within component
      wtn.gei <- sum( gei.within.components$value )
      within.lin <- rowSums( gei.within.components$lin )

      # between (residual)
      btw.gei <- ttl.gei - wtn.gei
      between.lin <- ttl.lin - within.lin

      # create vector of estimates
      estimates <- c( ttl.gei, wtn.gei, btw.gei )
      names( estimates ) <- c( "total", "within", "between" )

      # treat out of sample

      # create linearized matrix
      lin.matrix <- matrix( data = c(ttl.lin, within.lin, between.lin), ncol = 3, dimnames = list( NULL, c( "total", "within", "between" ) ) )
      rownames( lin.matrix ) <- rownames( design$variables )
      rm(ttl.lin, within.lin, between.lin)

      ### compute deff
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
    attr(rval, "statistic") <- "gei decomposition"
    attr(rval,"epsilon")<- epsilon
    attr(rval,"group")<- as.character( subgroup )[[2]]
    class(rval) <- c( "cvystat" , "svrepstat" , "svystat" )
    if ( linearized ) attr(rval,"linearized") <- lin.matrix
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
    if ( is.character(deff) || deff) attr( rval , "deff") <- deff.estimate

    # retorna objeto
    rval

  }


#' @rdname svygeidec
#' @export
svygeidec.DBIsvydesign <-
  function (formula, subgroup, design, ...) {

    design$variables <-
      cbind(
        getvars( formula , design$db$connection , design$db$tablename , updates = design$updates , subset = design$subset ),
        getvars( subgroup , design$db$connection , design$db$tablename , updates = design$updates , subset = design$subset )
      )

    NextMethod("svygeidec", design)

  }

