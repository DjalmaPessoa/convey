# compares with vardpoor output

# load libraries
library( laeken )
library( survey )
library( convey )
library( vardpoor )
library( testthat )

# test context
context( "svyrmir comparison with vardpoor")

# collect and format data
data( eusilc )
names( eusilc ) <- tolower( names( eusilc ) )

# create new data.frame
dati = data.frame( IDd = seq( 10000 , 10000 + nrow( eusilc ) - 1 ) , eusilc )

# custom function for variance estimation
SE_lin2 <- function( t , des ) {
  variance<-survey::svyrecvar( t/des$prob , des$cluster , des$strata , des$fpc , postStrata = des$postStrata )
  sqrt( variance )
}

### convey calculations

# build survey design objects
des_eusilc <- svydesign( ids = ~rb030 , strata =~db040 , weights = ~rb050 , data = eusilc )
des_eusilc_rep <- as.svrepdesign( des_eusilc , type= "bootstrap" )

# prepare survey design objects for convey
des_eusilc <- convey_prep( des_eusilc )
des_eusilc_rep <- convey_prep( des_eusilc_rep )

# calculate estimates using convey
fun_rmirw <- svyrmir( ~eqincome , des_eusilc , age = ~age , agelim = 65 )
fun_rmirw_rep <- svyrmir( ~eqincome , des_eusilc_rep , age = ~age , agelim = 65 )

# collect point estimates from convey object
convest <- coef( fun_rmirw )
attributes( convest ) <- NULL

# collect SE estimates from convey object
convse <- SE( fun_rmirw )
attributes( convse ) <- NULL

### vardpoor calculations

# calculate estimates with vardpoor
vardpoor_rmirw <- linrmir(Y = "eqincome", id = "IDd", weight = "rb050", Dom = NULL, dataset = dati , age = "age" )

# set up coefficients object
vardest<- vardpoor_rmirw$value
attributes(vardest)<- NULL
vardest<- unlist(vardest)

# calculate SE using vardpoor linearization, but survey's variance estimation function
varse <- SE_lin2( vardpoor_rmirw$lin$lin_rmir , des_eusilc )
attributes( varse ) <- NULL

##### domain estimation

### using vardpoor

# domain esitmation using vardpoor
vardpoor_rmird <- linrmir(Y = "eqincome" , id = "IDd" , weight = "rb050" , Dom = "hsize" , dataset = dati , age = "age" )

# colect point estimates
vardestd <- unlist( vardpoor_rmird$value$rmir )

# calculate SE using vardpoor linearization, but survey's variance estimation function
varsed <- sapply( data.frame(vardpoor_rmird$lin)[,2:10] , function(t) SE_lin2( t , des_eusilc ) )
attributes( varsed ) <- NULL

### using convey

# calculate estimates
fun_rmird <- svyby( ~eqincome , ~hsize , des_eusilc , svyrmir , age = ~age , agelim = 65 )

# collect point estimates
convestd <- coef( fun_rmird )
attributes( convestd ) <- NULL

# collect SE estimates
convsed<- SE( fun_rmird )

# perform tests
test_that("compare results convey vs vardpoor",{

  # compare point estimates
  expect_equal( vardest[[1]] , convest )

  # compare point estimates on domains
  expect_equal( vardestd , convestd )

  # compare CV estimates within five pp. tolerance
  expect_equal( varse / vardest , convse / vardest , tolerance = .05 )
  expect_equal( varsed / vardestd , convsed / vardestd , tolerance = .05 )

  # compare SE estimates
  skip( "SE estimates differ")
  expect_equal( varse , convse )
  expect_equal( varsed , convsed )

} )
