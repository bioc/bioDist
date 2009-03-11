KLD <- function(f,g, supp=c(-3,3), subdivisions=1000){
# f and g must be 1-arg funcs that evaluate
# on vector inputs
  me <- .Machine$double.eps
  wlrat <- function(x) log((f(x)+me)/(g(x)+me))*f(x)
  options(show.error.messages = FALSE)
  on.exit(options(show.error.messages = TRUE))
  xx<-try(integrate( wlrat,supp[1],supp[2], subdivisions=subdivisions,rel.tol=0.01 ))
  if(inherits(xx, "try-error") )
     return(NA)
  return(xx$value)
}
