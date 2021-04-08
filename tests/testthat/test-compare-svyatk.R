# compares with vardpoor output

# load libraries
library( laeken )
library( survey )
library( convey )
library( IC2 )
library( testthat )

# test context
context( "svyatk comparison with IC2")

# collect and format data
data( eusilc )
names( eusilc ) <- tolower( names( eusilc ) )

### convey calculations

# build survey design objects
des_eusilc <- svydesign( ids = ~rb030 , strata =~db040 , weights = ~rb050 , data = eusilc )
des_eusilc_rep <- as.svrepdesign( des_eusilc , type= "bootstrap" )

# prepare survey design objects for convey
des_eusilc <- convey_prep( des_eusilc )
des_eusilc_rep <- convey_prep( des_eusilc_rep )

# filter positive incomes
des_eusilc <- subset( des_eusilc , eqincome > 0 )
des_eusilc_rep <- subset( des_eusilc_rep , eqincome > 0 )

# calculate estimates using convey
fun_atkw <- svyatk( ~eqincome , des_eusilc )
fun_atkw_rep <- svyatk( ~eqincome , des_eusilc_rep )

# collect point estimates from convey object
convest <- coef( fun_atkw )
attributes( convest ) <- NULL

# collect SE estimates from convey object
convse <- SE( fun_atkw )
attributes( convse ) <- NULL

### IC2 calculations

# compute point estimates
IC2est  <- calcAtkinson( x = eusilc$eqincome[ eusilc$eqincome > 0 ], w = eusilc$rb050[ eusilc$eqincome > 0 ] )$ineq$index[[1]]

# perform tests
test_that("compare results convey vs vardpoor",{

  # compare point estimates
  expect_equal( IC2est[[1]] , convest )

  # # compare point estimates on domains
  # expect_equal( vardestd , convestd*100 )
  #
  # # compare SE estimates
  # expect_equal( varse , convse*100 )
  # expect_equal( varsed , convsed*100 )

} )
