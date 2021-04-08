context("test treatment for variables with missings")

# load libraries
library( survey )
library( convey )
library( laeken )
library( testthat )

# collect and format data
data( eusilc )
names( eusilc ) <- tolower( names( eusilc ) )

# set up survey design objects
des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 , weights = ~rb050 , data = eusilc )
des_eusilc_rep <-as.svrepdesign( des_eusilc , type= "bootstrap" , replicates = 50 )

# prepare for convey
des_eusilc <- convey_prep( des_eusilc )
des_eusilc_rep <- convey_prep( des_eusilc_rep )

# # only striclty positive incomes
# test_that( "error on income <= 0 " , expect_error( svygei( ~eqincome , des_eusilc , epsilon = .5 )  ) )

# filter positive
des_eusilc <- subset( des_eusilc , py010n > 0 | is.na( py010n ) )
des_eusilc_rep <- subset( des_eusilc_rep , py010n > 0 | is.na( py010n ) )

# test for missing

# make initial object
out <- NULL

# cycle through single-valued functions
for( this_fun in c( "svygei" , "svygpg" , "svyatk" , "svyqsr" , "svypoormed" , "svyrmpg" , "svyrmir" , "svyisq" , "svyiqalpha" , "svyarpr" , "svyarpt" , "svyfgt" , "svygini" ) ){

  # test across functions
  test_that( paste( "coef and SE matrix must return missing:" , this_fun ) , {


    # get function object
    final_fun <- FUN <- get( this_fun )

    # set up arguments
    if( identical( FUN , svyrmpg ) ) final_fun <- function( ... ) FUN( ... , thresh = TRUE )
    if( identical( FUN , svyrmir ) ) final_fun <- function( ... ) FUN( ... , age = ~ age , med_old = TRUE )
    if( identical( FUN , svyisq ) ) final_fun <- function( ... ) FUN( ... , alpha = 0.2 )
    if( identical( FUN , svyiqalpha ) ) final_fun <- function( ... ) FUN( ... , alpha = 0.5 )
    if( identical( FUN , svyfgt ) ) final_fun <- function( ... ) FUN( ... , g = 0 , abs_thresh = 10000 )
    if( identical( FUN , svygpg ) ) final_fun <- function( ... ) FUN( ... , sex = ~ rb090 )
    if( identical( FUN , svygei ) ) final_fun <- function( ... ) FUN( ... , epsilon = 0.5 )

    # evaluate function with missing
    # (these must return a missing value)
    est_lin <- final_fun( ~ py010n , des_eusilc )
    est_rep <- final_fun( ~ py010n , des_eusilc_rep )

    # evaluate function with missing
    # (these must return a missing value)
    est_lin_narm <- final_fun( ~ py010n , des_eusilc , na.rm = TRUE )
    est_rep_narm <- final_fun( ~ py010n , des_eusilc_rep , na.rm = TRUE )

    # test result in variable with missing
    expect_true( is.na( coef( est_lin ) ) )
    expect_true( is.na( SE( est_lin ) ) )
    expect_true( is.na( coef( est_rep ) ) )
    expect_true( is.na( SE( est_rep ) ) )

    # test result in variable with na.rm
    expect_true( !is.na( coef( est_lin_narm ) ) )
    expect_true( !is.na( SE( est_lin_narm ) ) )
    expect_true( !is.na( coef( est_rep_narm ) ) )
    expect_true( !is.na( SE( est_rep_narm ) ) )

  } )

}

