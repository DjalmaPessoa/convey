# load libraries
library( survey )
library( convey )
library( laeken )
library( testthat )
# library( vardpoor )

# loop through g and threshold type
for ( this.threshold in c( "abs", "relq","relm") ) for ( this.g in c(1.5,2) )  {

  # return test context
  context( paste0( "svyrich-" , this.threshold ," FGTT2 g=", this.g , " output survey.design and svyrep.design" ) )

  ### test 1: test if funtion works on unweighted objects

  # load data
  data("api")

  # set up convey design
  expect_warning( dstrat1<-convey_prep(svydesign( id=~1,data=apistrat) ) )

  # perform tests
  test_that( paste0( "svyrich g=", this.g , " works on unweighted designs"), {
    expect_false( is.na ( coef( svyrich( ~api00, design=dstrat1 , g = this.g , abs_thresh = 600 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE ) ) ) )
    expect_false( is.na ( SE( svyrich( ~api00, design=dstrat1 , g = this.g , abs_thresh = 600 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE ) ) ) )
  } )

  ### test 2: income data from eusilc --- data.frame-backed design object

  # collect and format data
  data( eusilc )
  names( eusilc ) <- tolower( names( eusilc ) )

  # set up survey design objects
  des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 , weights = ~rb050 , data = eusilc )
  des_eusilc_rep <-as.svrepdesign( des_eusilc , type= "bootstrap" , replicates = 50 )

  # prepare for convey
  des_eusilc <- convey_prep( des_eusilc )
  des_eusilc_rep <- convey_prep( des_eusilc_rep )

  # calculate estimates
  a1 <- svyrich( ~eqincome , des_eusilc , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
  a2 <- svyby( ~eqincome , ~hsize, des_eusilc , svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )
  a2.nocov <- svyby( ~eqincome , ~hsize, des_eusilc , svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = FALSE )
  b1 <- svyrich( ~eqincome , des_eusilc_rep , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
  b2 <- svyby( ~eqincome , ~hsize, des_eusilc_rep , svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )
  b2.nocov <- svyby( ~eqincome , ~hsize, des_eusilc_rep , svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = FALSE )

  # calculate auxilliary tests statistics
  cv_diff1 <- abs( cv( a1 ) - cv( b1 ) )
  se_diff2 <- max( abs( SE( a2 ) - SE( b2 ) ) , na.rm = TRUE )

  # perform tests
  test_that( "output svyrich" , {
    expect_is( coef( a1 ) ,"numeric" )
    expect_is( coef( a2 ) , "numeric" )
    expect_is( coef( b1 ) ,"numeric" )
    expect_is( coef( b2 ) ,"numeric" )
    expect_equal( coef( a1 ) , coef( b1 ) )
    expect_equal( coef( a2 ) , coef( b2 ) )
    # expect_lte( cv_diff1 , coef(a1) * 0.20 )         # the difference between CVs should be less than 5% of the coefficient, otherwise manually set it
    expect_lte( se_diff2 , max( coef(a2) ) * 0.20 )  # the difference between CVs should be less than 10% of the maximum coefficient, otherwise manually set it
    expect_is( SE( a1 ) , "matrix" )
    expect_is( SE( a2 ) , "numeric" )
    expect_is( SE( b1 ) , "numeric" )
    expect_is( SE( b2 ) , "numeric" )
    expect_lte( confint( a1 )[1] ,  coef( a1 ) )
    expect_gte( confint( a1 )[2] , coef( a1 ) )
    expect_lte( confint( b1 )[,1] , coef( b1 ) )
    expect_gte( confint( b1 )[2] , coef( b1 ) )
    expect_equal( sum( confint( a2 )[,1] <= coef( a2 ) ) , length( coef( a2 ) ) )
    expect_equal( sum( confint( a2 )[,2] >= coef( a2 ) ) , length( coef( a2 ) ) )
    expect_equal( sum( confint( b2 )[,1] <= coef( b2 ) ) , length( coef( b2 ) ) )
    expect_equal( sum( confint( b2 )[,2] >= coef( b2 ) ) , length( coef( b2 ) ) )

    # check equality of linearized variables
    expect_equal( attr( a1 , "linearized" ) , attr( b1 , "linearized" ) )

    # check equality vcov diagonals
    expect_equal( diag( vcov( a2 ) ) , suppressWarnings( diag( vcov( a2.nocov ) ) ) )
    expect_equal( diag( vcov( b2 ) ) , suppressWarnings( diag( vcov( b2.nocov ) ) ) )

  } )

  ### test 2: income data from eusilc --- database-backed design object

  # perform tests
  test_that("database svyrich",{

    # skip test on cran
    skip_on_cran()

    # load libraries
    library(RSQLite)
    library(DBI)

    # set-up database
    dbfile <- tempfile()
    conn <- dbConnect( RSQLite::SQLite() , dbfile )
    dbWriteTable( conn , 'eusilc' , eusilc )

    # database-backed design
    dbd_eusilc <-
      svydesign(
        ids = ~rb030 ,
        strata = ~db040 ,
        weights = ~rb050 ,
        data="eusilc",
        dbname=dbfile,
        dbtype="SQLite"
      )

    # prepare for convey
    dbd_eusilc <- convey_prep( dbd_eusilc )

    # calculate estimates
    c1 <- svyrich( ~eqincome , dbd_eusilc , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
    c2 <- svyby( ~eqincome , ~hsize, dbd_eusilc , svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )

    # remove table and close connection to database
    dbRemoveTable( conn , 'eusilc' )
    dbDisconnect( conn )

    # peform tests
    expect_equal( coef( a1 ) , coef( c1 ) )
    expect_equal( coef( a2 ) , coef( c2 ) )
    expect_equal( SE( a1 ) , SE( c1 ) )
    expect_equal( SE( a2 ) , SE( c2 ) )
    expect_equal( deff( a1 ) , deff( c1 ) )
    expect_equal( deff( a2 ) , deff( c2 ) )
    expect_equal( vcov( a2 ) , vcov( c2 ) )

    # check equality of linearized variables
    expect_equal( attr( c1 , "linearized" ) , attr( a1 , "linearized" ) )
    expect_equal( attr( c1 , "linearized" ) , attr( b1 , "linearized" ) )

    # check equality of indices
    expect_equal( attr( a1 , "index" ) , attr( c1 , "index" ) )

    # check equality of linearized variables
    expect_equal( attr( a1 , "linearized" ) , attr( b1 , "linearized" ) )

  } )

  ### test 3: compare subsetted objects to svyby objects

  # calculate estimates
  sub_des <- svyrich( ~eqincome , design = subset( des_eusilc , hsize == 1) , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
  sby_des <- svyby( ~eqincome, by = ~hsize, design = des_eusilc, FUN = svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )
  sub_rep <- svyrich( ~eqincome , design = subset( des_eusilc_rep , hsize == 1) , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
  sby_rep <- svyby( ~eqincome, by = ~hsize, design = des_eusilc_rep, FUN = svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )

  # perform tests
  test_that("subsets equal svyby",{

    # domain vs svyby: coefficients must be equal
    expect_equal( as.numeric( coef( sub_des ) ) , as.numeric( coef( sby_des ) )[1] )
    expect_equal( as.numeric( coef( sub_rep ) ) , as.numeric( coef( sby_rep ) )[1] )

    # domain vs svyby: SEs must be equal
    expect_equal( as.numeric( SE( sub_des ) ) , as.numeric( SE( sby_des ) )[1] )
    expect_equal( as.numeric( SE( sub_rep ) ) , as.numeric( SE( sby_rep ) )[1] )

    # domain vs svyby: DEffs must be equal
    expect_equal( as.numeric( deff( sub_des ) ) , as.numeric( deff( sby_des ) )[1] )
    expect_equal( as.numeric( deff( sub_rep ) ) , as.numeric( deff( sby_rep ) )[1] )

    # domain vs svyby and svydesign vs svyrepdesign:
    # coefficients should match across svydesign
    expect_equal( as.numeric( coef( sub_des ) ) , as.numeric( coef( sby_rep ) )[1] )

    # domain vs svyby and svydesign vs svyrepdesign:
    # coefficients of variation should be within five percent
    cv_diff <- abs( cv( sub_des ) - cv( sby_rep )[1] )
    expect_lte( cv_diff , .5 )

    # check equality of linearized variables
    expect_equal( attr( sub_des , "linearized" ) , attr( sub_rep , "linearized" ) )

    # check equality of linearized variables
    expect_equal( vcov( sub_des )[1] , vcov( sby_des )[1,1] )
    expect_equal( vcov( sub_rep )[1] , vcov( sby_rep )[1,1] )

  } )

  ### test 4: compare subsetted objects to svyby objects

  # compare database-backed designs to non-database-backed designs
  test_that("dbi subsets equal non-dbi subsets",{

    # skip test on cran
    skip_on_cran()

    # load libraries
    library( RSQLite )
    library( DBI )

    # set up database
    dbfile <- tempfile()
    conn <- dbConnect( RSQLite::SQLite() , dbfile )
    dbWriteTable( conn , 'eusilc' , eusilc )

    # create database-backed design (with survey design information)
    dbd_eusilc <-
      svydesign(
        ids = ~rb030 ,
        strata = ~db040 ,
        weights = ~rb050 ,
        data="eusilc",
        dbname=dbfile,
        dbtype="SQLite"
      )

    # create a hacky database-backed svrepdesign object
    # mirroring des_eusilc_rep
    dbd_eusilc_rep <-
      svrepdesign(
        weights = ~ rb050,
        repweights = attr( des_eusilc_rep , "full_design" )$repweights ,
        scale = attr( des_eusilc_rep , "full_design" )$scale ,
        rscales = attr( des_eusilc_rep , "full_design" )$rscales ,
        type = "bootstrap" ,
        data = "eusilc" ,
        dbtype="SQLite" ,
        dbname = dbfile ,
        combined.weights = FALSE
      )

    # prepare for convey
    dbd_eusilc <- convey_prep( dbd_eusilc )
    dbd_eusilc_rep <- convey_prep( dbd_eusilc_rep )

    # calculate estimates
    sub_dbd <- svyrich( ~eqincome , design = subset( dbd_eusilc , hsize == 1) , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
    sby_dbd <- svyby( ~eqincome, by = ~hsize, design = dbd_eusilc, FUN = svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )
    sub_dbr <- svyrich( ~eqincome , design = subset( dbd_eusilc_rep , hsize == 1) , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , linearized = TRUE )
    sby_dbr <- svyby( ~eqincome, by = ~hsize, design = dbd_eusilc_rep, FUN = svyrich , g = this.g , abs_thresh = 20000 , type_thresh = this.threshold , type_measure = "FGTT2" , deff = TRUE , covmat = TRUE )

    # remove table and disconnect from database
    dbRemoveTable( conn , 'eusilc' )
    dbDisconnect( conn )

    # perform tests
    expect_equal( coef( sub_des ) , coef( sub_dbd ) )
    expect_equal( coef( sub_rep ) , coef( sub_dbr ) )
    expect_equal( SE( sub_des ) , SE( sub_dbd ) )
    expect_equal( SE( sub_rep ) , SE( sub_dbr ) )
    expect_equal( deff( sub_des ) , deff( sub_dbd ) )
    expect_equal( deff( sub_rep ) , deff( sub_dbr ) )
    expect_equal( vcov( sby_des ) , vcov( sby_dbd ) )
    expect_equal( vcov( sby_rep ) , vcov( sby_dbr ) )

    # compare database-backed subsetted objects to database-backed svyby objects
    # dbi subsets equal dbi svyby
    expect_equal( as.numeric( coef( sub_dbd ) ) , as.numeric( coef( sby_dbd ) )[1] )
    expect_equal( as.numeric( coef( sub_dbr ) ) , as.numeric( coef( sby_dbr ) )[1] )
    expect_equal( as.numeric( SE( sub_dbd ) ) , as.numeric( SE( sby_dbd ) )[1] )
    expect_equal( as.numeric( SE( sub_dbr ) ) , as.numeric( SE( sby_dbr ) )[1] )
    expect_equal( as.numeric( deff( sub_dbd ) ) , as.numeric( deff( sby_dbd ) )[1] )
    expect_equal( as.numeric( deff( sub_dbr ) ) , as.numeric( deff( sby_dbr ) )[1] )

    # check equality of linearized variables
    expect_equal( attr( sub_dbd , "linearized" ) , attr( sub_des , "linearized" ) )
    expect_equal( attr( sub_dbr , "linearized" ) , attr( sub_rep , "linearized" ) )

  } )

}
