"dens2func" <-
function(x, ...) {
    #rng<-range(x)
    rng<-extendrange(x[c(i,j),],r=range(x,na.rm=T),f=0.01)
    den1 <- density(x,from=rng[1],to=rng[2])
    f <- function(w) approx(den1$x, den1$y, w, yleft=0, yright=0)$y
    class(f) <- "dfun"
    f
}

