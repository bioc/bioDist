locf2func <- function(x, ...) {
    require(locfit) || stop("can only be used if locfit library is available")
    l = locfit(~x, ...)
    den1 <- preplot(l,...)
    xvals <- den1$xev[[1]]
    yvals <- den1$trans(den1$fit)
    f <- function(w) approx(xvals, yvals, w, yleft=0, yright=0)$y
    class(f) <- "dfun"
    f
}

