
# Expand matrix
#' @title Expand matrix
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export
expand_matrix<-function(Q)
{
  diag(Q)<-rep(0, nrow(Q))
  Q.vec<-c(Q)
  Q.vec.frac<-fractions(Q.vec, cycles = 10, max.denominator = 2000)
  fracs<-attr(Q.vec.frac[Q.vec.frac!=0],"fracs")
  # sometimes numbers without fraction appear like 1. add fractions
  fracs[fracs=="1"]<-"1/1"
  ####

  fracs<-strsplit(fracs, "/")
  lapply(fracs, function (x) x[2]) %>% unlist() %>% as.numeric() ->vec.denom

  lcm<-numbers::mLCM(vec.denom)
  lambda2<-1/lcm

  lapply(fracs, function (x) x[1]) %>% unlist() %>% as.numeric() ->vec.num
  theta.vec<-vec.num*(lcm/vec.denom)
  Theta<-Q
  Theta[Theta!=0]<-theta.vec

  # Expansion
  expn.by<-apply(Theta, 2, max)
  # entire row  can be 0 in theta, so change to 1
  expn.by[expn.by==0]<-1

  QD<-matrix(0, nrow = sum(expn.by), ncol = sum(expn.by), byrow = T)

  # make names
  names<-c()
  #i=2
  for (i in seq_along(expn.by))
  {
    names<-c(names, paste(i, ".", c(1:expn.by[i]), sep="") )
  }
  colnames(QD)<-rownames(QD)<-names

  # fill in lambda2
  # cumulative sum to find start of rows faster
  expn.by.cum<-cumsum(expn.by)
  expn.by.start<-c(1, expn.by.cum[-length(expn.by.cum)]+1 )

  for (i in 1:nrow(Theta))
  {
    for (j in 1:nrow(Theta))
    {
      if ( Theta[i, j]>0 )
      {
        QD[expn.by.start[i]:expn.by.cum[i],
           expn.by.start[j]:(expn.by.start[j]+(Theta[i, j]-1) )  ] <- lambda2
      }
    }
  }

  diag(QD)<-apply(QD, 1, function(x) sum(x))*-1
  return (QD)

}
######
#QD<-expand_matrix(Q)
#diag(QD)<-0


