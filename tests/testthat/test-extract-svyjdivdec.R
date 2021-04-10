# load libraries
library( survey )
library( convey )
library( laeken )
library( testthat )
# library( vardpoor )

# return test context
context( "svyjdivdec decomposition output survey.design and svyrep.design" )

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

# filter postive incomes
des_eusilc <- subset( des_eusilc , eqincome >0 )
des_eusilc_rep <- subset( des_eusilc_rep , eqincome > 0 )

# calculate estimates
a1 <- svyjdivdec( ~eqincome , des_eusilc , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
a2 <- svyby( ~eqincome , ~rb090, des_eusilc , svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )
a2.nocov <- svyby( ~eqincome , ~rb090, des_eusilc , svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = FALSE )
b1 <- svyjdivdec( ~eqincome , des_eusilc_rep , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
b2 <- svyby( ~eqincome , ~rb090, des_eusilc_rep , svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )
b2.nocov <- svyby( ~eqincome , ~rb090, des_eusilc_rep , svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = FALSE )
d1 <- svyjdiv( ~eqincome , des_eusilc , deff = TRUE , linearized = TRUE )
d2 <- svyby( ~eqincome , ~rb090, des_eusilc , svyjdiv , deff = TRUE )
e1 <- svygeidec( ~eqincome , des_eusilc , subgroup = ~db040 , epsilon = 0 , deff = TRUE , linearized = TRUE )
f1 <- svygeidec( ~eqincome , des_eusilc , subgroup = ~db040 , epsilon = 1 , deff = TRUE , linearized = TRUE )

# calculate auxilliary tests statistics
cv_diff1 <- max( abs( cv( a1 ) - cv( b1 ) ) )
se_diff2 <- max( abs( SE( a2 ) - SE( b2 ) ) , na.rm = TRUE )

# perform tests
test_that( "output svyjdivdec" , {
  expect_is( coef( a1 ) , "numeric" )
  expect_is( coef( a2 ) , "numeric" )
  expect_is( coef( b1 ) , "numeric" )
  expect_is( coef( b2 ) , "numeric" )
  expect_equal( coef( a1 ) , coef( b1 ) )
  expect_equal( coef( a2 ) , coef( b2 ) )
  # expect_lte( cv_diff1 , (coef(a1)[[1]]) * 0.05 )         # the difference between CVs should be less than 5% of the coefficient, otherwise manually set it
  expect_lte( se_diff2 , max( coef(a2) ) * 0.05 ) # the difference between CVs should be less than 10% of the maximum coefficient, otherwise manually set it
  expect_equal( sum( confint( a2 )[,1] <= coef( a2 ) ) , length( coef( a2 ) ) )
  expect_equal( sum( confint( a2 )[,2] >= coef( a2 ) ) , length( coef( a2 ) ) )
  expect_equal( sum( confint( b2 )[,1] <= coef( b2 ) ) , length( coef( b2 ) ) )
  expect_equal( sum( confint( b2 )[,2] >= coef( b2 ) ) , length( coef( b2 ) ) )
  expect_equal( coef( a1 )[[1]] , coef( d1 )[[1]] )
  expect_equal( as.numeric( coef( a2 )[1:2] ) , as.numeric( coef( d2 ) ) )

  # check equality of linearized variables
  expect_equal( attr( a1 , "linearized" ) , attr( b1 , "linearized" ) )
  expect_equal( attr( a1 , "index" ) , attr( b1 , "index" ) )

  # check equality vcov diagonals
  expect_equal( diag( vcov( a2 ) ) , suppressWarnings( diag( vcov( a2.nocov ) ) ) )
  expect_equal( diag( vcov( b2 ) ) , suppressWarnings( diag( vcov( b2.nocov ) ) ) )

  # compare with svyjdiv
  expect_equal( coef( a1 )[[1]] , coef( d1 )[[1]] )
  expect_equal( SE( a1 )[[1]] , SE( d1 )[[1]] )
  expect_equal( as.numeric( attr( a1 , "linearized" )[,1,drop=FALSE] ) , as.numeric( attr( d1 , "linearized" ) ) )
  expect_equal( as.numeric( SE( a2 )[,1] ) , as.numeric( SE( d2 ) ) )

  # compare with svygei
  expect_equal( coef( a1 ) , rowSums( sapply( list( e1 , f1 ) , coef ) ) )
  expect_equal( attr( a1 , "linearized" ) , attr( e1 , "linearized" ) + attr( f1 , "linearized" ) )

} )

### test 2: income data from eusilc --- database-backed design object

# perform tests
test_that("database svyjdivdec",{

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

  # filter postive incomes
  dbd_eusilc <- subset( dbd_eusilc , eqincome >0 )

  # calculate estimates
  c1 <- svyjdivdec( ~eqincome , dbd_eusilc , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
  c2 <- svyby( ~eqincome , ~rb090, dbd_eusilc , svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )

  # remove table and close connection to database
  dbRemoveTable( conn , 'eusilc' )
  dbDisconnect( conn )

  # peform tests
  expect_equal( coef( a1 ) , coef( c1 ) )
  expect_equal( coef( a2 ) , coef( c2 )[ match( names( coef( c2 ) ) , names( coef( a2 ) ) ) ] )
  expect_equal( SE( a1 ) , SE( c1 ) )
  expect_equal( SE( a2 ) , SE( c2 )[ 2:1 , ] )
  expect_equal( deff( a1 ) , deff( c1 ) )
  expect_equal( deff( a2 ) , deff( c2 )[ 2:1 , ] )

  # check equality of linearized variables
  expect_equal( attr( c1 , "linearized" ) , attr( a1 , "linearized" ) )
  expect_equal( attr( c1 , "index" ) , attr( a1 , "index" ) )

} )

### test 3: compare subsetted objects to svyby objects

# calculate estimates
sub_des <- svyjdivdec( ~eqincome , design = subset( des_eusilc , rb090 == "male" ) , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
sby_des <- svyby( ~eqincome, by = ~rb090, design = des_eusilc, FUN = svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )
sub_rep <- svyjdivdec( ~eqincome , design = subset( des_eusilc_rep , rb090 == "male" ) , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
sby_rep <- svyby( ~eqincome, by = ~rb090, design = des_eusilc_rep, FUN = svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )

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
  expect_lte( cv_diff , .10 )

  # check equality of linearized variables
  expect_equal( attr( sub_des , "linearized" ) , attr( sub_rep , "linearized" ) )

  # check equality of variances
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

  # filter postive incomes
  dbd_eusilc <- subset( dbd_eusilc , eqincome >0 )
  dbd_eusilc_rep <- subset( dbd_eusilc_rep , eqincome >0 )

  # calculate estimates
  sub_dbd <- svyjdivdec( ~eqincome , design = subset( dbd_eusilc , rb090 == "male" ) , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
  sby_dbd <- svyby( ~eqincome, by = ~rb090, design = dbd_eusilc, FUN = svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )
  sub_dbr <- svyjdivdec( ~eqincome , design = subset( dbd_eusilc_rep , rb090 == "male" ) , subgroup = ~db040 , deff = TRUE , linearized = TRUE )
  sby_dbr <- svyby( ~eqincome, by = ~rb090, design = dbd_eusilc_rep, FUN = svyjdivdec , subgroup = ~db040 , deff = TRUE , covmat = TRUE )

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
  expect_equal( vcov( sub_dbd ) , vcov( sub_des ) )
  expect_equal( vcov( sub_dbr ) , vcov( sub_rep ) )

  # compare equality of linearized variables
  expect_equal( attr( sub_dbd , "linearized" ) , attr( sub_dbr , "linearized" ) )
  expect_equal( attr( sub_dbd , "linearized" ) , attr( sub_des , "linearized" ) )
  expect_equal( attr( sub_dbr , "linearized" ) , attr( sub_rep , "linearized" ) )

  # compare equality of indices
  expect_equal( attr( sub_dbd , "index" ) , attr( sub_dbr , "index" ) )
  expect_equal( attr( sub_dbd , "index" ) , attr( sub_des , "index" ) )
  expect_equal( attr( sub_dbr , "index" ) , attr( sub_rep , "index" ) )


} )
