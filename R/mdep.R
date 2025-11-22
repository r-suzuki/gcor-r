#' Estimate generalized correlation and related measures
#'
#' @description Estimate measures based on mutual dependency, which includes:
#' \itemize{
#'   \item{Generalized correlation measure (`gcor`)}
#'   \item{Dissimilarity between variables (`gdis`)}
#'   \item{Predictability score (`pscore`)}
#' }
#'
#' @param x a vector, matrix, data frame or formula. If formula, `data` should be specified.
#' `gdis` requires a matrix or data frame.
#' @param y `NULL` (default) or a vector, matrix or data frame with compatible dimensions to `x`.
#' @param k `NULL` (default) or an integer specifying the number of groups for discretization.
#' Numerical data are divided into `k` groups using `k`-quantiles.
#' If `NULL`, it is determined automatically.
#' @param data `NULL` (default) or a data frame. Required if `x` is a formula.
#' @param simplify a logical. If `TRUE`, the returned value is coerced to
#' a vector when one of its dimensions is one.
#' @param dropNA a character specifying how to handle missing values.
#' It should be one of the following:
#' \describe{
#'   \item{`"none"`}{(default) Treat missing values as observations of a single
#' categorical value, namely `NA`. Recommended for reflecting missing patterns in the analysis.}
#'   \item{`"casewise"`}{Casewise deletion; rows containing any missing value are removed.
#' Similar to `use = "complete.obs"` in \code{\link{cor}}.}
#'   \item{`"pairwise"`}{Pairwise deletion; for any pair of columns (x,y),
#' the i-th row is removed if x\[i\] or y\[i\] is missing.
#' This process is applied per pair; the same row is not removed with another pair (w,z),
#' if both w\[i\] and z\[i\] are not missing.
#' Similar to `use = "pairwise.complete.obs"` in \code{\link{cor}}.}
#' }
#' @param ... additional arguments (`diag` and `upper`) passed to `as.dist` function.
#' See \code{\link{as.dist}} for details.
#'
#' @return For `gcor` and `pscore`, a numeric matrix is returned (or a vector if `simplify = TRUE`).
#' For `gdis`, an object of class `"dist"` is returned.
#'
#' @references
#' Suzuki, R. (2025). *A generalization of correlation coefficient*. preprint.
#' \url{https://r-suzuki.github.io}
#'
#' @examples
#' # Generalized correlation measure
#' gcor(iris)
#'
#' # Predictability of Species from other variables
#' ps <- pscore(Species ~ ., data = iris)
#' dotchart(sort(ps), main = "Predictability of Species")
#'
#' # Clustering
#' gd <- gdis(iris)
#' hc <- hclust(gd, method = "ward.D2")
#' plot(hc)
#'
#' # Multidimensional scaling
#' mds <- cmdscale(gd, k = 2)
#' plot(mds, type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE,
#'      main = "cmdscale with gdis(iris)")
#' text(mds[,1], mds[,2], rownames(mds))
#' @name gcor-package
#'
#' @importFrom utils data
#' @importFrom stats as.dist complete.cases model.frame model.response setNames
NULL

# @param measure a character specifying the type of measure, one of `"cor"`, `"dist"`, `"pred"`.
# `gcor` is a wrapper for `mdep` with `measure = "cor"`.
# Similarly, `gdis` wraps `measure = "dist"`, and `pscore` wraps `measure = "pred"`.
# @param xname a character to be used as the name of `x`, when x is an atomic vector.
# @param yname a character used as the name of `y` (same as `xname` for `x`).
mdep <- function(x, y = NULL, k = NULL, data = NULL, simplify = FALSE, dropNA = "none",
                 measure,
                 xname = deparse1(substitute(x)), yname = deparse1(substitute(y)),
                 ...
                 ) {
  IS_XY_SYNMETRIC <- FALSE
  MEASURES <- c("cor", "dist", "pred")
  DROP_NA <- c("none", "casewise", "pairwise")
  xx <- yy <- kk <- ret <- NULL

  if(is.na(match(measure, MEASURES))) {
    stop(.gen_msg("measure", measure, MEASURES))
  }

  if(is.null(y) & is.atomic(x) && is.null(dim(x))) {
    stop("Supply non-NULL y if x is a non-matrix atomic vector.")
  }

  if(is.na(match(dropNA, DROP_NA))) {
    stop(.gen_msg("dropNA", dropNA, DROP_NA))
  }

  if(inherits(x, "formula")) {
    if(!is.null(y)) {
      stop("y should be NULL if x is a formula")
    }

    if(is.null(data)) {
      stop("Supply non-null data if x is a formula")
    }

    mf <- model.frame(formula = x, data = data)
    xx <- mf[, -1, drop = FALSE]
    yy <- as.data.frame(model.response(mf))

    if(ncol(yy) == 1) {
      names(yy) <- all.vars(x[[2]])
    }
  }

  if(is.null(xx)) {
    if(is.atomic(x) && is.null(dim(x))) {
      xx <- .vec2df(x, xname)
    } else {
      xx <- as.data.frame(x)
    }
  }

  if(is.null(yy)) {
    if(is.null(y)) {
      yy <- xx
      IS_XY_SYNMETRIC <- TRUE
    } else if(is.atomic(y) && is.null(dim(y))) {
      yy <- .vec2df(y, yname)
    } else {
      yy <- as.data.frame(y)
    }
  }

  stopifnot(nrow(xx) == nrow(yy))
  if(dropNA == "casewise") {
    cc <- complete.cases(xx, yy)
    xx <- subset(xx, cc)
    yy <- subset(yy, cc)
  }

  # default selection of k (10-by-2 rule)
  if(is.null(k)) k <- pmax(2, floor(nrow(xx)^log10(2)/2))
  stopifnot(length(k) == 1)

  ret <- matrix(rep(NA_real_, ncol(xx) * ncol(yy)),
                nrow = ncol(xx), ncol = ncol(yy),
                dimnames = list(names(xx), names(yy)))

  for(i in seq_len(ncol(xx))) {
    for(j in seq_len(ncol(yy))) {
      if(IS_XY_SYNMETRIC && i > j) {
        # DO NOTHING
      } else if(IS_XY_SYNMETRIC && i == j) {
        ret[i, j] <- if(measure == "dist") 0.0 else 1.0
      } else {
        m_ij <- .mdep_quantile_grid(xx[,i], yy[,j], k,
                                    useNA = (dropNA != "pairwise"))
        phi_ij <- m_ij$estimate

        # phi_ij should be greater or equal to 1, but estimated values
        # with some approximation could be less than 1. It is adjusted here.
        if(!is.na(phi_ij) && phi_ij < 1) {
          warning("Estimated mutual dependency < 1; adjusted to 1.")
          phi_ij <- 1
        }

        r2 <- 1 - 1/phi_ij

        kx <- m_ij$kx
        ky <- m_ij$ky
        kk <- sqrt(kx) * sqrt(ky)

        # If kk == 1, both x and y is constant.
        # In this case gcor(x,y) = 1 and gdis(x,y) = 0.
        if(measure == "cor") {
          ret[i, j] <- if(kk == 1) 1 else sqrt(r2 / (1 - 1/kk))
          if(IS_XY_SYNMETRIC) ret[j, i]  <- ret[i, j]
        } else if(measure == "dist") {
          ret[i, j] <- if(kk == 1) 0 else sqrt(1 - r2 / (1 - 1/kk))
          if(IS_XY_SYNMETRIC) ret[j, i] <- ret[i, j]
        } else if(measure == "pred") {
          # If ky == 1, y is constant and completely dependent on any random variable.
          # So pscore(x,y) = 1.
          ret[i, j] <- if(ky == 1) 1 else sqrt(r2 / (1 - 1/ky))
          if(IS_XY_SYNMETRIC) ret[j, i] <- if(kx == 1) 1 else sqrt(r2 / (1 - 1/kx))
        }
      }
    }
  }

  if(measure == "dist") {
    ret <- as.dist(ret, ...)
  } else if(simplify) {
    if(nrow(ret) == 1 && ncol(ret) == 1) {
      ret <- as.vector(ret)
    } else if(nrow(ret) == 1) {
      ret <- setNames(as.vector(ret), colnames(ret))
    } else if(ncol(ret) == 1) {
      ret <- setNames(as.vector(ret), rownames(ret))
    }
  }

  return(ret)
}

#' @rdname gcor-package
#' @export
gcor <- function(x, y = NULL, k = NULL, data = NULL, simplify = TRUE, dropNA = "none") {
  mdep(x = x, y = y, k = k, data = data, simplify = simplify, dropNA = dropNA, measure = "cor",
       xname = deparse1(substitute(x)), yname = deparse1(substitute(y)))
}

#' @rdname gcor-package
#' @export
gdis <- function(x, k = NULL, dropNA = "none", ...) {
  if(!is.matrix(x) && !is.data.frame(x)) {
    stop("x should be a matrix or data frame.")
  }

  mdep(x = x, y = NULL, k = k, data = data, dropNA = dropNA, measure = "dist",
       xname = deparse1(substitute(x)), yname = deparse1(substitute(y)), ...)
}

#' @rdname gcor-package
#' @export
pscore <- function(x, y = NULL, k = NULL, data = NULL, simplify = TRUE, dropNA = "none") {
  mdep(x = x, y = y, k = k, data = data, simplify = simplify, dropNA = dropNA, measure = "pred",
       xname = deparse1(substitute(x)), yname = deparse1(substitute(y)))
}
