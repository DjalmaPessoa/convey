# load libraries
library( survey )
library( convey )
library( laeken )
library( testthat )
# library( vardpoor )

# return test context
context( "svyfgtdec-abs  decomposition output survey.design and svyrep.design")

### test 2: income data from eusilc --- data.frame-backed design object

# collect and format data
data( eusilc )
names( eusilc ) <- tolower( names( eusilc ) )

# set up survey design objects
des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 , weights = ~rb050 , data = eusilc )
des_eusilc_rep <- as.svrepdesign( des_eusilc , type= "bootstrap" , replicates = 50 )

# prepare for convey
des_eusilc <- convey_prep( des_eusilc )
des_eusilc_rep <- convey_prep( des_eusilc_rep )

# calculate estimates
a1 <- svyfgtdec( ~eqincome , des_eusilc , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
a2 <- svyby( ~eqincome , ~rb090, des_eusilc , svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
b1 <- svyfgtdec( ~eqincome , des_eusilc_rep , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
b2 <- svyby( ~eqincome , ~rb090, des_eusilc_rep , svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
d1 <- svyfgt( ~eqincome , des_eusilc , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
d2 <- svyby( ~eqincome , ~rb090, des_eusilc , svyfgt , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
e1 <- svyfgt( ~eqincome , des_eusilc , g = 0 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
e2 <- svyby( ~eqincome , ~rb090, des_eusilc , svyfgt , g = 0 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
f1 <- svyfgt( ~eqincome , des_eusilc , g = 1 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
f2 <- svyby( ~eqincome , ~rb090, des_eusilc , svyfgt , g = 1 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )

# calculate auxilliary tests statistics
cv_diff1 <- max( abs( cv( a1 ) - cv( b1 ) ) )
se_diff2 <- max( abs( SE( a2 ) - SE( b2 ) ) , na.rm = TRUE )

# perform tests
test_that( "output svyfgtdec" , {
  expect_is( coef( a1 ) ,"numeric" )
  expect_is( coef( a2 ) , "numeric" )
  expect_is( coef( b1 ) ,"numeric" )
  expect_is( coef( b2 ) ,"numeric" )
  expect_equal( coef( a1 ) , coef( b1 ) )
  expect_equal( coef( a2 ) , coef( b2 ) )
  # expect_lte( cv_diff1 , (coef(a1)[[1]]) * 0.05 )         # the difference between CVs should be less than 5% of the coefficient, otherwise manually set it
  expect_lte( se_diff2 , max( coef(a2) ) * 0.05 ) # the difference between CVs should be less than 10% of the maximum coefficient, otherwise manually set it
  expect_equal( sum( confint( a2 )[,1] <= coef( a2 ) ) , length( coef( a2 ) ) )
  expect_equal( sum( confint( a2 )[,2] >= coef( a2 ) ) , length( coef( a2 ) ) )
  expect_equal( sum( confint( b2 )[,1] <= coef( b2 ) ) , length( coef( b2 ) ) )
  expect_equal( sum( confint( b2 )[,2] >= coef( b2 ) ) , length( coef( b2 ) ) )
  expect_equal( coef( a1 )[[1]] , coef( d1 )[[1]] )
  expect_equal( SE( a1 )[[1]] , SE( d1 )[[1]] )
  expect_equal( coef( a1 )[[2]] , coef( e1 )[[1]] )
  expect_equal( SE( a1 )[[2]] , SE( e1 )[[1]] )
  expect_equal( coef( a1 )[[3]] , coef( f1 )[[1]] )
  expect_equal( SE( a1 )[[3]] , SE( f1 )[[1]] )
} )

### test 2: income data from eusilc --- database-backed design object

# perform tests
test_that("database svyfgtdec",{

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
  c1 <- svyfgtdec( ~eqincome , dbd_eusilc , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
  c2 <- svyby( ~eqincome , ~rb090, dbd_eusilc , svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )

  # remove table and close connection to database
  dbRemoveTable( conn , 'eusilc' )
  dbDisconnect( conn )

  # peform tests
  expect_equal( coef( a1 ) , coef( c1 ) )
  expect_equal( coef( a2 ) , coef( c2 )[ match( names( coef( c2 ) ) , names( coef( a2 ) ) ) ] )
  expect_equal( SE( a1 ) , SE( c1 ) )
  expect_equal( SE( a2 ) , SE( c2 )[ 2:1 , ] )
  

} )

### test 3: compare subsetted objects to svyby objects

# calculate estimates
sub_des <- svyfgtdec( ~eqincome , design = subset( des_eusilc , rb090 == "male" ) , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
sby_des <- svyby( ~eqincome, by = ~rb090, design = des_eusilc, FUN = svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
sub_rep <- svyfgtdec( ~eqincome , design = subset( des_eusilc_rep , rb090 == "male" ) , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
sby_rep <- svyby( ~eqincome, by = ~rb090, design = des_eusilc_rep, FUN = svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )

# perform tests
test_that("subsets equal svyby",{

  # domain vs svyby: coefficients must be equal
  expect_equal( as.numeric( coef( sub_des ) ) , as.numeric( coef( sby_des[1,] ) ) )
  expect_equal( as.numeric( coef( sub_rep ) ) , as.numeric( coef( sby_rep[1,] ) ) )

  # domain vs svyby: SEs must be equal
  expect_equal( as.numeric( SE( sub_des ) ) , as.numeric( SE( sby_des[1,] ) ) )
  expect_equal( as.numeric( SE( sub_rep ) ) , as.numeric( SE( sby_rep[1,] ) ) )

  # domain vs svyby and svydesign vs svyrepdesign:
  # coefficients should match across svydesign
  expect_equal( as.numeric( coef( sub_des ) ) , as.numeric( coef( sby_rep[1,] ) ) )

  # domain vs svyby and svydesign vs svyrepdesign:
  # coefficients of variation should be within five percent
  cv_diff <- max( abs( cv( sub_des ) - cv( sby_rep )[1,] ) )
  expect_lte( cv_diff , .05 )

  

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
  sub_dbd <- svyfgtdec( ~eqincome , design = subset( dbd_eusilc , rb090 == "male" ) , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
  sby_dbd <- svyby( ~eqincome, by = ~rb090, design = dbd_eusilc, FUN = svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
  sub_dbr <- svyfgtdec( ~eqincome , design = subset( dbd_eusilc_rep , rb090 == "male" ) , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )
  sby_dbr <- svyby( ~eqincome, by = ~rb090, design = dbd_eusilc_rep, FUN = svyfgtdec , g = 2 , abs_thresh = 7000 , type_thresh = "abs" , deff = TRUE )

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
  expect_equal( vcov( sub_des ) , vcov( sub_dbd ) )
  expect_equal( vcov( sub_rep ) , vcov( sub_dbr ) )


  # compare database-backed subsetted objects to database-backed svyby objects
  # dbi subsets equal dbi svyby
  expect_equal( as.numeric( coef( sub_dbd ) ) , as.numeric( coef( sby_dbd[2,] ) ) )
  expect_equal( as.numeric( coef( sub_dbr ) ) , as.numeric( coef( sby_dbr[2,] ) ) )
  expect_equal( as.numeric( SE( sub_dbd ) ) , as.numeric( SE( sby_dbd[2,] ) ) )
  expect_equal( as.numeric( SE( sub_dbr ) ) , as.numeric( SE( sby_dbr[2,] ) ) )

} )
