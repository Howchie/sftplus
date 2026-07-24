suppressMessages(devtools::load_all(".", quiet=TRUE))
set.seed(1)
mk <- function(subj, base) {
  cells <- list(HH=c(2,2), HL=c(2,1), LH=c(1,2), LL=c(1,1))
  rows <- lapply(names(cells), function(nm){
    ch <- cells[[nm]]
    n <- 40
    # over-additive-ish: HH+LL slower than HL+LH, scaled by base (subject speed)
    mu <- base * switch(nm, HH=1.05, HL=1.00, LH=1.00, LL=0.95)
    data.frame(Subject=subj, Condition="A", RT=abs(rnorm(n, mu, base*0.15)),
               Correct=TRUE, Channel1=ch[1], Channel2=ch[2])
  })
  do.call(rbind, rows)
}
d <- do.call(rbind, Map(mk, paste0("P",1:4), c(0.4,0.6,0.5,0.8)))
cat("=== multi-subject analytic ===\n")
a <- micGroup.bayes(d, ndraws=500, burnin=100, chains=2, seed=42, rope=0.025)
print(a$method_code); print(a$estimand)
print(a$score)
print(a$population_summary)
print(a$posterior_probability[c("population_overadditive","population_underadditive","population_additive")])
cat("param names:", a$summary$parameter, "\n")

cat("\n=== multi-subject bootstrap ===\n")
b <- micGroup.bayes(d, ndraws=500, burnin=100, chains=2, seed=42, rope=0.025,
                    var_method="bootstrap", n_boot=200)
print(b$score[,c("subject","relative_MIC","se","se_analytic","z")])
print(b$var_method)

cat("\n=== single subject analytic ===\n")
s <- micGroup.bayes(d, Subject="P1", ndraws=500, chains=2, seed=7, rope=0.025)
print(s$method_code); print(s$summary); print(names(s$posterior_probability))

cat("\n=== ucip.bayes bootstrap var_method ===\n")
rt <- list(rexp(60,3.5), rexp(60,1.5), rexp(60,2))
u <- ucip.bayes(rt, ndraws=500, chains=2, seed=3, var_method="bootstrap", n_boot=300)
print(u$var_method); print(u$score)

cat("\n=== capacityGroup.bayes bootstrap ===\n")
rtg <- lapply(1:3, function(i) list(rexp(40,2), rexp(40,1.8), rexp(40,1.6)))
names(rtg) <- paste0("P",1:3)
g <- capacityGroup.bayes(rtg, ndraws=400, burnin=100, chains=2, seed=5, var_method="bootstrap", n_boot=200)
print(g$var_method); print(g$score[,c("subject","theta_hat","se","se_analytic")])
cat("\nDONE\n")
