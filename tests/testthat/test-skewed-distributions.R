# load libraries
library(survey)
library(convey)
library(testthat)

# options( warn = 2 )
# options( error = recover )

# run tests
context( "functions work on weird distributions" )

# auxilliary function
no.na <- function( z , value = FALSE ){ z[ is.na( z ) ] <- value ; z }

# list of functions
all_funs <-
  c( "svyfgt" , "svywatts" , "svychu" , "svyarpt" , "svyarpt" ,
     "svyrmir" , "svyqsr" , "svygpg" , "svypoormed" , "svyrmpg" ,
     "svyatk" , "svygei" , "svyjdiv" ,"svygini" ,
     "svyiqalpha" , "svyisq" , "svylorenz" )

# test environment
test_that( "testing functions work on weird distributions" , {

  # skip on cran
  skip_on_cran()

  # loop through functions
  for( n in c( 50 , 1000 ) ){

    # set random numbers seed
    set.seed( n )

    # simulate variables
    dist_frame <-
      data.frame(
        wt_binom_two = rbinom( n , 1 , .5 ) ,
        wt_binom = rbinom( n , 1 , .5 ) ,
        wt_unif = runif( n ) ,
        beta_050_050 = rbeta( n , 0.5 , 0.5 ) ,
        beta_500_100 = rbeta( n , 5 , 1 ) ,
        beta_100_300 = rbeta( n , 1 , 3 ) ,
        beta_200_200 = rbeta( n , 2 , 2 ) ,
        beta_200_500 = rbeta( n , 2 , 5 ) ,
        pois5 = rpois( n , 5 ) ,
        pois500 = rpois( n , 100 ) ,
        sex = ifelse( rbinom( n , 1 , 0.5 ) , "male" , "female" ) ,
        age = sample( 0:99 , n , replace = TRUE )
      )

    # create sex variable
    dist_frame[ -which( names( dist_frame ) == 'sex' ) ] <- sapply( dist_frame[ -which( names( dist_frame ) == 'sex' ) ] , as.numeric )

    # loop through functions
    for( this_function in all_funs ){

      # get function
      FUN <- get( this_function )

      # prepare design object
      suppressWarnings( {
        unwtd_des <- convey_prep( svydesign( ~ 1 , data = dist_frame ) )
        binom_des <- convey_prep( svydesign( ~ 1 , data = dist_frame , weight = ~ wt_binom ) )
        unif_des <- convey_prep( svydesign( ~ 1 , data = dist_frame , weight = ~ wt_unif ) )
      } )

      # create replicated designs
      unwtd_rep <- convey_prep( as.svrepdesign( unwtd_des , type = 'bootstrap' ) )
      binom_rep <- convey_prep( as.svrepdesign( binom_des , type = 'bootstrap' ) )
      unif_rep <- convey_prep( as.svrepdesign( unif_des , type = 'bootstrap' ) )

      # select variables
      for( this_prefix in c( 'unwtd' , 'binom' , 'unif' ) ){

        # set up variables
        lin_des <- get( paste0( this_prefix , "_des" ) )
        rep_des <- get( paste0( this_prefix , "_rep" ) )

        # loop through variables
        for( this_formula in list( ~ beta_050_050 , ~ beta_500_100 , ~ beta_100_300 , ~ beta_200_200 , ~ beta_200_500 , ~ pois5 , ~ pois500 ) ){

          # function-specific parameters here

          if( any( unlist( lapply( list( svygei , svyatk , svywatts , svychu , svyjdiv ) , function( z ) identical( z , FUN ) ) ) ) ){
            lin_des <- subset( lin_des , no.na( lin_des$variables[ , as.character( this_formula )[2] ] > 0 ) )
            rep_des <- subset( rep_des , no.na( rep_des$variables[ , as.character( this_formula )[2] ] > 0 ) )
          }

          lin_params_list <- list( this_formula , lin_des )
          rep_params_list <- list( this_formula , rep_des )

          if( identical( FUN , svyfgt ) ){
            lin_params_list <- c( lin_params_list , list( g=0, type_thresh= "abs", abs_thresh=10000 ) )
            rep_params_list <- c( rep_params_list , list( g=0, type_thresh= "abs", abs_thresh=10000 ) )
          }

          if( identical( FUN , svychu ) ){
            lin_params_list <- c( lin_params_list , list( g=.5, type_thresh= "abs", abs_thresh=10000 ) )
            rep_params_list <- c( rep_params_list , list( g=.5, type_thresh= "abs", abs_thresh=10000 ) )
          }

          if( identical( FUN , svywatts ) ){
            lin_params_list <- c( lin_params_list , list( type_thresh= "abs", abs_thresh=10000 ) )
            rep_params_list <- c( rep_params_list , list( type_thresh= "abs", abs_thresh=10000 ) )
          }

          if( identical( FUN , svygpg ) ){
            lin_params_list <- c( lin_params_list , list( sex = ~sex ) )
            rep_params_list <- c( rep_params_list , list( sex = ~sex ) )
          }

          if( identical( FUN , svyiqalpha ) | identical( FUN , svyisq ) ){
            lin_params_list <- c( lin_params_list , list( alpha = 0.5 ) )
            rep_params_list <- c( rep_params_list , list( alpha = 0.5 ) )
          }

          if( identical( FUN , svyrmir ) ){
            lin_params_list <- c( lin_params_list , list( age = ~age ) )
            rep_params_list <- c( rep_params_list , list( age = ~age ) )
          }

          if( this_prefix == 'unwtd' ) wt_vec <- seq( nrow( dist_frame ) )
          if( this_prefix == 'binom' ) wt_vec <- which( dist_frame$wt_binom > 0 )
          if( this_prefix == 'unif' ) wt_vec <- which( dist_frame$wt_unif > 0 )

          # evaluate functions
          lin_res <- do.call( FUN , lin_params_list )
          rep_res <- do.call( FUN , rep_params_list )

          # test results
          test_that( paste( "testing functions work on weird distributions" , this_function , this_prefix , as.character( this_formula )[2] ) , {

            # result should give some non-missing number
            expect_that( coef( lin_res ) , is.numeric )
            expect_that( coef( lin_res ) , function( w ) all( !is.na( w ) ) )
            expect_that( coef( rep_res ) , is.numeric )
            expect_that( coef( rep_res ) , function( w ) all( !is.na( w ) ) )

            # coefficients should be equal
            expect_equal( coef( lin_res ) , coef( rep_res ) )

          } )

        }
      }
    }
  }
} )

