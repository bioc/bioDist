"dens2func" <-
function(x, ...) {
    den1 <- density(x, ...)
    f <- function(w) approx(den1$x, den1$y, w, yleft=0, yright=0)$y
    class(f) <- "dfun"
    f
}

