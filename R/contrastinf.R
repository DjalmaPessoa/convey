#' Generalized linearization of a smooth function of survey statistics

#'
#' @param exprlist a call
#' @param infunlist a list of lists, each having two components: value - the estimate value and lin - the linearized variable
#'
#' @details The call must use function that \code{deriv} knows how to differentiate. It allows to compute the linearized variable of a complex indicator from the linearized variables of simpler component variables, avoiding the formal derivatives calculations.
#'
#' @return a list with two components: values - the estimate value and lin - the linearized variable
#'
#' @author Djalma Pessoa and Anthony Damico
#'
#' @seealso \code{svyqsr}
#'
#' @references Guillaume Osier (2009). Variance estimation for complex indicators
#' of poverty and inequality. \emph{Journal of the European Survey Research
#' Association}, Vol.3, No.3, pp. 167-195,
#' ISSN 1864-3361, URL \url{https://ojs.ub.uni-konstanz.de/srm/article/view/369}.
#'
#'
#' @examples
#' library(survey)
#' library(laeken)
#' data(eusilc) ; names( eusilc ) <- tolower( names( eusilc ) )
#'
#' # linearized design
#' des_eusilc <- svydesign( ids = ~rb030 , strata = ~db040 ,  weights = ~rb050 , data = eusilc )
#' des_eusilc <- convey_prep(des_eusilc)
#'
#' w <- weights(des_eusilc)
#'
#' # ratio linearization
#' T1 = list(value = sum(w*eusilc$eqincome) , lin = eusilc$eqincome )
#' T2 = list(value = sum(w), lin = rep (1, nrow(eusilc)) )
#' list_all <- list( T1 = T1, T2 = T2)
#' lin_R = contrastinf (quote(T1/T2), list_all)
#'
#' # estimate of the variable eqincome mean
#' lin_R$value
#' # se estimate of the variable eqincome mean
#' SE(svytotal(lin_R$lin, des_eusilc))
#' # to check, use
#' svymean (~eqincome, des_eusilc)
#'
#' # quintile share ratio (qsr) linearization
#' S20 <- svyisq(~ eqincome, design = des_eusilc, .20 , linearized = TRUE )
#' S20_val <- coef (S20); attributes (S20_val) <- NULL
#' S20_lin <- attr(S20 , "linearized" )
#' S80 <- svyisq(~ eqincome, design = des_eusilc, .80 , linearized = TRUE )
#' S80_val <- coef (S80); attributes (S80_val) <- NULL
#' S80_lin <- attr(S80 , "linearized" )
#' SU <- list (value = S80_val, lin = S80_lin )
#' SI <- list (value = S20_val, lin = S20_lin )
#' TOT <- list(value = sum( w * eusilc$eqincome) , lin = eusilc$eqincome )
#' list_all <- list (TOT = TOT, SI = SI, SU = SU )
#' lin_QSR <- contrastinf( quote((TOT-SU)/SI), list_all)
#'
#' # estimate of the qsr
#' lin_QSR$value
#' # se estimate of the qsr:
#' SE(svytotal(lin_QSR$lin, des_eusilc))
#' # to check, use
#' svyqsr(~eqincome, des_eusilc )
#' # proportion of income below the quantile .20
#' list_all <- list (TOT = TOT, SI = SI )
#' lin_Lor <- contrastinf( quote(SI/TOT), list_all)
#' # estimate of the proportion of income below the quantile .20
#' lin_Lor$value
#' # se estimate
#' SE(svytotal(lin_Lor$lin,des_eusilc))
#'
#' @export
contrastinf <- function(exprlist, infunlist) {
    datalist <- lapply(infunlist, function(thresh) thresh$value)
    listlin <- lapply(infunlist, function(thresh) thresh$lin)
    if (!is.list(exprlist))
        exprlist <- list(contrast = exprlist)
    dexprlist <- lapply(exprlist, function(expr) deriv(expr, names(datalist))[[1]])
    value <- eval(exprlist$contrast, datalist)
    values_deriv <- lapply(dexprlist, function(dexpr) eval(do.call(substitute, list(dexpr,
        datalist))))
    matval <- attr(values_deriv$contrast, "gradient")
    matlin <- matrix(NA, length(infunlist[[1]]$lin), ncol(matval))
    for (i in 1:length(listlin)) matlin[, i] <- listlin[[i]]
    IT_lin <- matlin %*% t(matval)
    list(value = value, lin = IT_lin)
}
