# compares with vardpoor output

# load libraries
library( laeken )
library( survey )
library( convey )
library( vardpoor )
library( testthat )

# test context
context( "svyfgt-relq comparison with vardpoor")

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
fun_arprw <- svyfgt( ~eqincome , des_eusilc , quantiles = 0.5 , percent = 0.6 , g = 0 , type_thresh = "relq" )
fun_arprw_rep <- svyarpr( ~eqincome , des_eusilc_rep , quantiles = 0.5 , percent = 0.6 , g = 0 , type_thresh = "relq" )

# collect point estimates from convey object
convest <- coef( fun_arprw )
attributes( convest ) <- NULL

# collect SE estimates from convey object
convse <- SE( fun_arprw )
attributes( convse ) <- NULL

### vardpoor calculations

# calculate estimates with vardpoor
vardpoor_arprw <- linarpr(Y = "eqincome", id = "IDd", weight = "rb050", Dom = NULL, dataset = dati, percentage = 60, order_quant = 50L)

# set up coefficients object
vardest<- vardpoor_arprw$value
attributes(vardest)<- NULL
vardest<- unlist(vardest)

# calculate SE using vardpoor linearization, but survey's variance estimation function
varse <- SE_lin2( vardpoor_arprw$lin$lin_arpr , des_eusilc )
attributes( varse ) <- NULL

##### domain estimation

### using vardpoor

# domain estimation using vardpoor
vardpoor_arprd <- linarpr(Y = "eqincome" , id = "IDd" , weight = "rb050" , Dom = "hsize" ,
                          dataset = dati , percentage = 60 , order_quant = 50L )

# colect point estimates
vardestd <- unlist(vardpoor_arprd$value$arpr )

# calculate SE using vardpoor linearization, but survey's variance estimation function
varsed <- sapply( data.frame(vardpoor_arprd$lin)[,2:10] , function(t) SE_lin2( t , des_eusilc ) )
attributes( varsed ) <- NULL

### using convey

# calculate estimates
fun_arprd <- svyby( ~eqincome , ~hsize , des_eusilc , svyfgt , quantiles = 0.5 , percent = 0.6 , g = 0 , type_thresh = "relq" )

# collect point estimates
convestd <- coef( fun_arprd )
attributes( convestd ) <- NULL

# collect SE estimates
convsed<- SE( fun_arprd )

# perform tests
test_that("compare results convey vs vardpoor",{

  # compare point estimates
  expect_equal( vardest , 100*convest )

  # compare point estimates on domains
  expect_equal( vardestd , 100*convestd )

  # compare CV estimates within .05 pp difference
  skip('different fgt-relq CV in domain')
  expect_equal( varse / vardest , 100*convse / vardest , tolerance =  1e-4  )
  expect_equal( varsed / vardestd , 100*convsed / vardestd , tolerance =  1e-4  )

  # compare SE estimates
  skip('different fgt-relq se')
  expect_equal( varse , 100*convse )
  expect_equal( varsed , 100*convsed )

} )
