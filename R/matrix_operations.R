### crerate initial matrix

#' @title Create intial matix for run_hifi
#' @description Creates intial matrix for evolutinary agorithm. NA's mean constant (0 or diagonal)
#' cells, 0's mean cells that can potentially change, 1's mean cells with the rate > 0.
#' @param part_scheme partitioning scheme that is a vector where element id is an observable state
#' while the element value is the number of hidden changes
#' @param set.NA a vector that sets main diagonal elements (i.k, i.l) to NA. Good to use if some observable state
#' is known to lack any hidden states. For example if observable state 1.x lacks hidden states than
#' specify set.NA=c(1)
#' @return Matrix
#' @examples
#' #part_scheme<-c(5, 2, 3)
#' #names(part_scheme)<-c("s1", "s2", "s3")
#' #make_init_matrix(part_scheme, set.NA=c(2,3) )
#' #make_init_matrix(part_scheme)
#' @export
make_init_matrix<-function(part_scheme, set.NA=NULL)
{
  st.total<-sum(part_scheme)
  M<-matrix(0, nrow = st.total, ncol = st.total)

  names<-c()
  # set names
  for (i in 1:length(part_scheme))
  {
    names<-c(names,
    paste(rep(i, part_scheme[i]), c(1:part_scheme[i]), sep="." )
    )
  }
  colnames(M)<-rownames(M)<-names

  # set 1's
  #valid.states<-paste( c(1: length(part_scheme) ), rep(1, length(part_scheme) ), sep="." )

  valid.states<-c(1, cumsum(part_scheme)+1)
  valid.states<-valid.states[-length(valid.states)]
  names(valid.states)<-names(part_scheme)

  for (i in 1:length(valid.states))
  {
    index<-valid.states[-i]
    M[valid.states[i], index]<-1
  }

  # set diag
  diag(M)<-NA

  if(!is.null(set.NA))
  {
    # add last state
    valid.states<-c(valid.states, nrow(M)+1)

    for (i in 1:length(set.NA))
    {
      M[valid.states[set.NA[i]]:(valid.states[set.NA[i]+1]-1),
        valid.states[set.NA[i]]:(valid.states[set.NA[i]+1]-1) ]<-NA
    }
  }

  M<-phyhifi(M=M)
  return(M)
}
####




##################
#
# check strong lumpability
#
###########
# part_scheme=list(c(1, 2), c(3), c(4))
# part_scheme=list(c(1), c(2), c(3,4))
# part_scheme=list(c(1,2), c(3,4))
# is_strg_lumpable(M_kr, list(c(1,2), c(3,4)))
#
# M_kr<-make_init_matrix(c(2,2))
#
# M_kr=rows2rate_matrix(M.temp, row.vec, dt_rates[,7])
###
is_strg_lumpable<-function(M_kr, part_scheme){

  # normilize diag
  diag(M_kr)<-apply(M_kr, 1, function(x) sum(x))*-1

  Nper.part=lapply(part_scheme, length)%>%unlist
  stat2M=matrix(0, nrow=length(Nper.part), ncol=ncol(M_kr))
  for (i in 1:nrow(stat2M))
    stat2M[i, part_scheme[[i]]]<-1

  M_rows=M_kr %*% t(stat2M)
  tru.vals=c()

  for (i in 1:length(part_scheme)){
    for (j in 1:ncol(M_rows)){
      tru.vals=c(tru.vals, length(unique(M_rows[part_scheme[[i]],j]))<2)
    }
  }
  # if this is false then matrix is not lumpable
  return(all(tru.vals))
}

#M3=comb2matrices(M1, M2, dependent.state=2)

#M_kr=indepen2matrices(M1, M2,  name.sep=":", diag.as=0)
#M_kr=indepen2matrices(M1, M2,  name.sep=":")
#is_strg_lumpable(M3, part_scheme)


######
#
# Summarize matrices
#
##########

# summ_matrix(best.Q)
# summ_matrix(Q.pool)

summ_matrix<-function(M)
{
  M<-lapply(M, function(x) x$M)

  M.out<-M[[1]]
  for (i in 2:length(M))
  {
    M.out<-M.out + M[[i]]
  }

  #M.out[which(M.out>1)]<-1
  #M.out<-phyhifi(M=M.out)
  return(list(M.out) )
}


#############
#
# Recode chars given hidden states and display it in a proper format for corHMM
#
###########

#sim
#part_scheme<-c(5, 4,3)

#recode.to<-list(s1=c(1:5), s2=c(6:9), s3=c(10:12))
#phytools.chars<-sim$states

#convert_chars(sim$states, list(s1=c(1:5), s2=c(1:5), s3=c(10:12)) )

convert_chars<-function(phytools.chars, recode.to)
{


  st.to<-c()
  for (i in 1:length(recode.to))
  {
    st.to<-c(st.to, paste( (recode.to[[i]] ), collapse = "&")  )

  }

  out.chars<-mapvalues(phytools.chars,
                       from =c(1:length(recode.to)), to=st.to )

  out.chars<-cbind(names(out.chars), out.chars)
  #out.chars<-cbind( out.chars)

  out.remapped<-list(from=c(1:length(recode.to)), to=st.to)

  return( list(chars=out.chars, mapping=out.remapped) )


}


