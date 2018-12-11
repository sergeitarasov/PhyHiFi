### crerate initial matrix

part_scheme=list(s1=c(1:5), s2=c(1:2))

part_scheme<-c(5, 2, 3)
names(part_scheme)<-c("s1", "s2", "s3")

make_init_matrix<-function(part_scheme, )
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

  return(M)
}



####




##################
#
# check strong lumpability
#
###########
part_scheme=list(c(1, 2), c(3), c(4))
part_scheme=list(c(1), c(2), c(3,4))
part_scheme=list(c(1,2), c(3,4))
is_strg_lumpable(M_kr, list(c(1,2), c(3,4)))

M_kr=rows2rate_matrix(M.temp, row.vec, dt_rates[,7])
###
is_strg_lumpable<-function(M_kr, part_scheme){

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

M3=comb2matrices(M1, M2, dependent.state=2)

M_kr=indepen2matrices(M1, M2,  name.sep=":", diag.as=0)
M_kr=indepen2matrices(M1, M2,  name.sep=":")
is_strg_lumpable(M3, part_scheme)
