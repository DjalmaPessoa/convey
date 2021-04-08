# compares with vardpoor output

# load libraries
library( laeken )
library( survey )
library( convey )
library( vardpoor )
library( testthat )

# test context
context( "svygini comparison with vardpoor")

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
fun_giniw <- svygini( ~eqincome , des_eusilc )
fun_giniw_rep <- svygini( ~eqincome , des_eusilc_rep )

# collect point estimates from convey object
convest <- coef( fun_giniw )
attributes( convest ) <- NULL

# collect SE estimates from convey object
convse <- SE( fun_giniw )
attributes( convse ) <- NULL

### vardpoor calculations

# calculate estimates with vardpoor
vardpoor_giniw <- lingini(Y = "eqincome", id = "IDd", weight = "rb050", Dom = NULL, dataset = dati )

# set up coefficients object
vardest<- vardpoor_giniw$value
attributes(vardest)<- NULL
vardest<- unlist(vardest)

# calculate SE using vardpoor linearization, but survey's variance estimation function
varse <- SE_lin2( vardpoor_giniw$lin$lin_gini , des_eusilc )
attributes( varse ) <- NULL

##### domain estimation

### using vardpoor

# domain esitmation using vardpoor
vardpoor_ginid <- lingini(Y = "eqincome" , id = "IDd" , weight = "rb050" , Dom = "hsize" , dataset = dati )

# colect point estimates
vardestd <- unlist( vardpoor_ginid$value$Gini )

# calculate SE using vardpoor linearization, but survey's variance estimation function
varsed <- sapply( data.frame(vardpoor_ginid$lin)[,2:10] , function(t) SE_lin2( t , des_eusilc ) )
attributes( varsed ) <- NULL

### using convey

# calculate estimates
fun_ginid <- svyby( ~eqincome , ~hsize , des_eusilc , svygini )

# collect point estimates
convestd <- coef( fun_ginid )
attributes( convestd ) <- NULL

# collect SE estimates
convsed <- SE( fun_ginid )

# perform tests
test_that("compare results convey vs vardpoor",{

  # compare point estimates
  expect_equal( vardest[[1]] , convest*100 )

  # compare point estimates on domains
  expect_equal( vardestd , convestd*100 )

  # compare SE estimates
  expect_equal( varse , convse*100 )
  expect_equal( varsed , convsed*100 )

} )
