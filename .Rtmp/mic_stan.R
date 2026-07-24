suppressMessages(devtools::load_all(".", quiet=TRUE))
make_mic_data <- function(subjects=paste0("P",1:3), base=c(.4,.6,.5), n=40L, seed=5) {
  set.seed(seed); cells<-list(HH=c(2,2),HL=c(2,1),LH=c(1,2),LL=c(1,1))
  mult<-c(HH=1.05,HL=1,LH=1,LL=.95)
  do.call(rbind, Map(function(s,b) do.call(rbind, lapply(names(cells), function(nm){
    ch<-cells[[nm]]; data.frame(Subject=s,Condition="A",RT=abs(rnorm(n,b*mult[[nm]],b*.15)),
    Correct=TRUE,Channel1=ch[1],Channel2=ch[2])})), subjects, base))
}
d <- make_mic_data()
if (requireNamespace("rstan", quietly=TRUE)) {
  a <- micGroup.bayes(d, method="HalfNormal", ndraws=100, burnin=50, chains=2, seed=1)
  cat("Stan method_code:", a$method_code, " params:", a$summary$parameter, "\n")
  cat("stan_fit non-null:", !is.null(a$stan_fit), "\n")
} else cat("rstan not available\n")
