# load libraries
library( survey )
library( convey )
library( laeken )
library( testthat )
# library( vardpoor )

# return test context
context( "svygpg output survey.design and svyrep.design")

### test 2: income data from ses --- data.frame-backed design object

# collect and format data
data( ses )
rownames( ses ) <- NULL
names( ses ) <- gsub( "size" , "size_" , tolower( names( ses ) ) )

# set up survey design objects
des_ses <- svydesign( id=~1, weights=~weights, data=ses)
des_ses_rep <- as.svrepdesign( des_ses , type = "bootstrap" )

# prepare for convey
des_ses <- convey_prep( des_ses )
des_ses_rep <- convey_prep( des_ses_rep )

# calculate estimates
a1 <- svygpg( ~earningshour , des_ses , sex = ~sex )
a2 <- svyby( ~earningshour , ~location, des_ses, svygpg , epsilon =.5 , sex = ~sex )
b1 <- svygpg( ~earningshour , des_ses_rep , sex = ~sex )
b2 <- svyby( ~earningshour , ~location, des_ses_rep , svygpg , sex = ~sex )

# calculate auxilliary tests statistics
cv_diff1 <- abs( cv( a1 ) - cv( b1 ) )
se_diff2 <- max( abs( SE( a2 ) - SE( b2 ) ) , na.rm = TRUE )

# perform tests
test_that( "output svygpg" , {
  expect_is( coef( a1 ) ,"numeric" )
  expect_is( coef( a2 ) , "numeric" )
  expect_is( coef( b1 ) ,"numeric" )
  expect_is( coef( b2 ) ,"numeric" )
  expect_equal( coef( a1 ) , coef( b1 ) )
  expect_equal( coef( a2 ) , coef( b2 ) )
  expect_lte( cv_diff1 , coef(a1) * 0.05 )         # the difference between CVs should be less than 5% of the coefficient, otherwise manually set it
  # expect_lte( se_diff2 , max( coef(a2) ) * 0.05 ) # the difference between CVs should be less than 10% of the maximum coefficient, otherwise manually set it
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
} )

### test 2: income data from ses --- database-backed design object

# perform tests
test_that("database svygpg",{

  # skip test on cran
  skip_on_cran()

  # load libraries
  library(RSQLite)
  library(DBI)

  # set-up database
  dbfile <- tempfile()
  conn <- dbConnect( RSQLite::SQLite() , dbfile )
  dbWriteTable( conn , 'ses' , ses )

  # database-backed design
  dbd_ses <-
    svydesign(
      ids = ~1 ,
      weights = ~weights ,
      data="ses",
      dbname=dbfile,
      dbtype="SQLite"
    )

  # prepare for convey
  dbd_ses <- convey_prep( dbd_ses )

  # calculate estimates
  c1 <- svygpg( ~ earningshour , dbd_ses , sex = ~sex )
  c2 <- svyby( ~ earningshour , ~location , dbd_ses , FUN = svygpg , sex = ~sex )

  # remove table and close connection to database
  dbRemoveTable( conn , 'ses' )
  dbDisconnect( conn )

  # peform tests
  expect_equal( coef( a1 ) , coef( c1 ) )
  expect_equal( coef( a2 ) , coef( c2 ) )
  expect_equal( SE( a1 ) , SE( c1 ) )
  expect_equal( SE( a2 ) , SE( c2 ) )

  

} )

### test 3: compare subsetted objects to svyby objects

# calculate estimates
sub_des <- svygpg( ~earningshour , design = subset( des_ses , location == "AT1" ) , sex = ~sex )
sby_des <- svyby( ~earningshour, by = ~location, design = des_ses, FUN = svygpg , sex = ~sex  )
sub_rep <- svygpg( ~earningshour , design = subset( des_ses_rep , location == "AT1" ) , sex = ~sex  )
sby_rep <- svyby( ~earningshour, by = ~location, design = des_ses_rep, FUN = svygpg , sex = ~sex  )

# perform tests
test_that("subsets equal svyby",{

  # domain vs svyby: coefficients must be equal
  expect_equal( as.numeric( coef( sub_des ) ) , as.numeric( coef( sby_des ) )[1] )
  expect_equal( as.numeric( coef( sub_rep ) ) , as.numeric( coef( sby_rep ) )[1] )

  # domain vs svyby: SEs must be equal
  expect_equal( as.numeric( SE( sub_des ) ) , as.numeric( SE( sby_des ) )[1] )
  expect_equal( as.numeric( SE( sub_rep ) ) , as.numeric( SE( sby_rep ) )[1] )

  # domain vs svyby and svydesign vs svyrepdesign:
  # coefficients should match across svydesign
  expect_equal( as.numeric( coef( sub_des ) ) , as.numeric( coef( sby_rep ) )[1] )

  # domain vs svyby and svydesign vs svyrepdesign:
  # coefficients of variation should be within five percent
  cv_diff <- abs( cv( sub_des ) - cv( sby_rep )[1] )
  expect_lte( cv_diff , .5 )

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
  dbWriteTable( conn , 'ses' , ses )

  # create database-backed design (with survey design information)
  dbd_ses <-
    svydesign(
      ids = ~1 ,
      weights = ~weights ,
      data="ses",
      dbname=dbfile,
      dbtype="SQLite"
    )

  # create a hacky database-backed svrepdesign object
  # mirroring des_ses_rep
  dbd_ses_rep <-
    svrepdesign(
      weights = ~ weights,
      repweights = attr( des_ses_rep , "full_design" )$repweights ,
      scale = attr( des_ses_rep , "full_design" )$scale ,
      rscales = attr( des_ses_rep , "full_design" )$rscales ,
      type = "bootstrap" ,
      data = "ses" ,
      dbtype="SQLite" ,
      dbname = dbfile ,
      combined.weights = FALSE )

  # prepare for convey
  dbd_ses <- convey_prep( dbd_ses )
  dbd_ses_rep <- convey_prep( dbd_ses_rep )

  # calculate estimates
  sub_dbd <- svygpg( ~earningshour , design = subset( des_ses , location == "AT1" ) , sex = ~sex  )
  sby_dbd <- svyby( ~earningshour, by = ~location, design = des_ses, FUN = svygpg , sex = ~sex  )
  sub_dbr <- svygpg( ~earningshour , design = subset( des_ses_rep , location == "AT1" ) , sex = ~sex  )
  sby_dbr <- svyby( ~earningshour, by = ~location, design = des_ses_rep, FUN = svygpg , sex = ~sex  )

  # remove table and disconnect from database
  dbRemoveTable( conn , 'ses' )
  dbDisconnect( conn )

  # perform tests
  expect_equal( coef( sub_des ) , coef( sub_dbd ) )
  expect_equal( coef( sub_rep ) , coef( sub_dbr ) )
  expect_equal( SE( sub_des ) , SE( sub_dbd ) )
  expect_equal( SE( sub_rep ) , SE( sub_dbr ) )

  # compare database-backed subsetted objects to database-backed svyby objects
  # dbi subsets equal dbi svyby
  expect_equal( as.numeric( coef( sub_dbd ) ) , as.numeric( coef( sby_dbd ) )[1] )
  expect_equal( as.numeric( coef( sub_dbr ) ) , as.numeric( coef( sby_dbr ) )[1] )
  expect_equal( as.numeric( SE( sub_dbd ) ) , as.numeric( SE( sby_dbd ) )[1] )
  expect_equal( as.numeric( SE( sub_dbr ) ) , as.numeric( SE( sby_dbr ) )[1] )

} )
