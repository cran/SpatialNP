`spatial.location` <-
function(X,score=c("sign","signrank"),init=NULL,shape=TRUE,maxiter=500,eps=1e-6,na.action=na.fail)
{
 score<-match.arg(score)
 switch(score,
       "sign"=
       {
        if(is.matrix(shape))
        {
         res<-spatial.median(X%*%solve(mat.sqrt(shape)),init, maxiter,eps,na.action=na.action)
         attr(res,"shape")<-shape
         res
        }
        else if(is.logical(shape))
        {
         if(shape) 
         {
          res1<-HR.Mest(X,maxiter,eps,eps,na.action)
          res<-res1$center
          attr(res,"shape")<-res1$scatter

         }
         else 
         {
          res<-spatial.median(X,init,maxiter,eps,na.action=na.action)
          attr(res,"shape")<-diag(dim(X)[2])
         }
        }
        else stop("'shape' must be a matrix or logical")
       },
       "signrank"=
       ae.hl.estimate(X,init,shape,maxiter,eps,na.action=na.fail)
       )
}

