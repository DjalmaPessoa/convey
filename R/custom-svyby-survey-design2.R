#' @importFrom survey svyby
#' @export
svyby.survey.design2 <-
  function (formula, by, design, FUN, ..., deff = FALSE, keep.var = TRUE,
            keep.names = TRUE, verbose = FALSE, vartype = c("se", "ci",
                                                            "ci", "cv", "cvpct", "var"), drop.empty.groups = TRUE,
            covmat = FALSE, influence = covmat, na.rm.by = FALSE, na.rm.all = FALSE,
            multicore = getOption("survey.multicore"))
  {
    if (inherits(by, "formula"))
      byfactors <- model.frame(by, model.frame(design), na.action = na.pass)
    else byfactors <- as.data.frame(by)
    if (multicore && !requireNamespace("parallel", quietly = TRUE))
      multicore <- FALSE
    if (!inherits(formula, "formula")) {
      if (NROW(formula) != NROW(byfactors))
        stop("'formula' is the wrong length")
      if (!(is.data.frame(formula) || is.matrix(formula) ||
            is.vector(formula))) {
        stop("invalid type for 'formula'")
      }
    }
    hasdeff <- is.character(deff) || deff
    byfactor <- do.call("interaction", byfactors)
    dropped <- weights(design, "sampling") == 0
    if (na.rm.by)
      dropped <- dropped | apply(byfactors, 1, function(x) any(is.na(x)))
    if (na.rm.all) {
      if (inherits(formula, "formula"))
        allx <- model.frame(formula, model.frame(design),
                            na.action = na.pass)
      else allx <- formula
      dropped <- dropped | (!complete.cases(allx))
    }
    uniquelevels <- sort(unique(byfactor[!dropped]))
    uniques <- match(uniquelevels, byfactor)
    if (missing(vartype))
      vartype <- "se"
    vartype <- match.arg(vartype, several.ok = TRUE)
    nvartype <- base::which(eval(formals(sys.function())$vartype) %in%
                              vartype)
    if (any(is.na(nvartype)))
      stop("invalid vartype")
    if (keep.var) {
      unwrap <- function(x) {
        rval <- c(coef(x))
        nvar <- length(rval)
        rval <- c(rval, c(se = SE(x), ci_l = confint(x)[,
                                                        1], ci_u = confint(x)[, 2], cv = cv(x, warn = FALSE),
                          `cv%` = cv(x, warn = FALSE) * 100, var = SE(x)^2)[rep((nvartype -
                                                                                   1) * (nvar), each = nvar) + (1:nvar)])
        if (!is.null(attr(x, "deff")))
          rval <- c(rval, DEff = deff(x))
        rval
      }
      results <- (if (multicore)
        parallel::mclapply
        else lapply)(uniques, function(i) {
          idx <- byfactor %in% byfactor[i]
          if (verbose && !multicore)
            print(as.character(byfactor[i]))
          if (inherits(formula, "formula"))
            data <- formula
          else data <- subset(formula, byfactor %in% byfactor[i])
          if (covmat || influence) {
            r <- FUN(data, design[byfactor %in% byfactor[i],
            ], deff = deff, ..., influence = influence)
          }
          else {
            r <- FUN(data, design[byfactor %in% byfactor[i],
            ], deff = deff, ...)
          }
          if ( is.null( attr(r, "index") ) ) attr( r, "index" ) <- idx
          r
        })
      rval <- t(sapply(results, unwrap))
      if ( covmat || influence) {
        infs <- lapply(results, attr, "influence")
        idxs <- lapply(results, attr, "index")
        if (!all(sapply(idxs, is.logical))) idxs <- lapply( idxs, function( id ) rownames( design$allprob ) %in% id )
        if (all(sapply(infs, is.null))) stop("FUN does not return influence functions")
        inflmats <- vector("list", length(infs))
        for (i in seq_along(infs)) {
          inflmats[[i]] <- matrix(0, ncol = NCOL(infs[[i]]), nrow = length(idxs[[i]]))
          inflmats[[i]][idxs[[i]], ] <- infs[[i]]
        }
        inflmat <- do.call(cbind, inflmats)
        rownames( inflmat ) <- rownames( design$allprob )
      }
      if ( covmat ){
        covmat.mat <- svyrecvar(inflmat, design$cluster, design$strata, design$fpc, postStrata = design$postStrata)
      } else {
        covmats <- lapply(results, vcov)
        ncovmat <- sum(sapply(covmats, ncol))
        covmat.mat <- matrix(0, ncol = ncovmat, nrow = ncovmat)
        j <- 0
        for (i in 1:length(covmats)) {
          ni <- nrow(covmats[[i]])
          covmat.mat[j + (1:ni), j + (1:ni)] <- covmats[[i]]
          j <- j + ni
        }
      }
    } else {
      unwrap2 <- function(x) {
        if (!is.null(attr(x, "deff")))
          c(statistic = unclass(x), DEff = deff(x))
        else c(statistic = unclass(x))
      }
      rval <- sapply(uniques, function(i) {
        if (verbose)
          print(as.character(byfactor[i]))
        if (inherits(formula, "formula"))
          data <- formula
        else data <- subset(formula, byfactor %in% byfactor[i])
        unwrap2(FUN(data, design[byfactor %in% byfactor[i],
        ], deff = deff, ...))
      })
      if (is.matrix(rval))
        rval <- t(rval)
    }
    nr <- NCOL(rval)
    nstats <- nr/(1 + keep.var * (length(vartype) + ("ci" %in% vartype)) + hasdeff)
    if (nr > 1)
      rval <- cbind(byfactors[uniques, , drop = FALSE], rval)
    else rval <- cbind(byfactors[uniques, , drop = FALSE], statistic = rval)
    expand.index <- function(index, reps, x = FALSE) {
      ns <- max(index)
      if (x) {
        i <- matrix(1:(ns * reps), ncol = reps)
        rval <- t(i[index, ])
      }
      else {
        i <- matrix(1:(ns * reps), ncol = reps, nrow = ns,
                    byrow = TRUE)
        rval <- i[index, ]
      }
      as.vector(rval)
    }
    if (drop.empty.groups) {
      if (keep.names)
        rownames(rval) <- paste(byfactor[uniques])
      rval <- rval[order(byfactor[uniques]), ]
      i <- expand.index(order(byfactor[uniques]), nstats)
      if (keep.var)
        covmat.mat <- covmat.mat[i, i]
    }
    else {
      a <- do.call("expand.grid", lapply(byfactors, function(f) levels(as.factor(f))))
      a <- cbind(a, matrix(NA, ncol = nr, nrow = nrow(a)))
      names(a) <- names(rval)
      a[match(byfactor[uniques], levels(byfactor)), ] <- rval
      rval <- a
      if (keep.names)
        rownames(rval) <- levels(byfactor)
      if (keep.var) {
        tmp <- matrix(ncol = nrow(a) * nstats, nrow = nrow(a) *
                        nstats)
        i <- expand.index(match(byfactor[uniques], levels(byfactor)),
                          nstats, TRUE)
        tmp[i, i] <- covmat.mat
        covmat.mat <- tmp
      }
    }
    attr(rval, "svyby") <- list(margins = 1:NCOL(byfactors),
                                nstats = nstats, vars = if (keep.var) length(vartype) else 0,
                                deffs = deff, statistic = deparse(substitute(FUN)),
                                variables = names(rval)[-(1:NCOL(byfactors))][1:nstats],
                                vartype = vartype)
    if (!keep.names)
      rownames(rval) <- 1:NROW(rval)
    if (covmat) attr(rval, "var") <- covmat.mat
    if (influence) attr(rval, "influence") <- inflmat
    attr(rval, "call") <- sys.call()
    class(rval) <- c("svyby", "data.frame")
    rval
  }
