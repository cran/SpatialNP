`spatial.shape` <-
function(X, score = c("sign", "symmsign", "rank", "signrank"),
location = NULL, init = NULL, steps = Inf, eps = 1e-06, maxiter = 100,
na.action = na.fail)
{
score<-match.arg(score)
switch(score,
       "sign"=tyler.shape(X,location, init, steps, eps, maxiter, na.action=na.action),
       "symmsign"=duembgen.shape(X, init, steps, eps, maxiter, na.action=na.action),
       "rank"=rank.shape(X, init, steps, eps, maxiter, na.action),
       "signrank"=signrank.shape(X, location, init, steps, eps, maxiter, na.action))
}

