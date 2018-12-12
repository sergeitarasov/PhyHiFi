c.dt<-convert_chars(sim$states, list(s1=c(1:4), s2=c(1:4), s3=c(5:8)) )

out<- greedy_rayDISC(phy=sim, data=c.dt$chars, rate.mat=make.cor.matrix(M),
                     root.p="maddfitz", do.ans=FALSE, do.thorough=TRUE)


out3<- greedy_rayDISC(phy=sim, data=c.dt$chars, rate.mat=make.cor.matrix(M), node.states="marginal", model="ARD",
                     root.p="maddfitz", do.ans=FALSE, do.thorough=TRUE)

out1<- rayDISC(phy=phy.data, data=char.data, rate.mat=make.cor.matrix(rhf2$Q.pool[[40]]), node.states="marginal", model="ARD",
                     root.p="maddfitz")
rhf2$Q.pool[[40]]$Cor.sol

greedy_rayDISC(phy=phy.data, data=char.data, rate.mat=make.cor.matrix(rhf2$Q.pool[[40]]), node.states="marginal", model="ARD",
               root.p="maddfitz", do.ans=FALSE, do.thorough=TRUE)

run_hifi(M.initial, phy.data, char.data, deltaAIC=10,  N.iter=5,
         do.thorough=FALSE, root.p="maddfitz",  Cor.sol=TRUE)

# argiments
# M.init
# N.iter
# deltaAIC #AIC range
# root.phylo
# root.est


# ####
# #Q.cor[which(Q.cor==1)]<-0
# M<-Q.cor
# M.r<-matrix(c(1:nrow(M)^2), nrow(M), nrow(M) )
#
#
# # set intial rates
# M[1, obs.state[[2]][1]]<-1
# M[obs.state[[2]][1], 1]<-1
# M<-list(M)
# M<-list(phyhifi(M=M[[1]]) )
# ##
# result<-matrix(0, nrow=N.iter, ncol=4)
# colnames(result)<-c("N_Q.best", "N_new_mt", "", "")
#
#
# M<-make_init_matrix(c(4,4))
# class(M)
# i=c(2,10)
# class(i)
####
##############################
#
# Serach Loop
# PhyHiFi
#
##########################
phy.data<-sim
c.dt<-convert_chars(sim$states, list(s1=c(1:4), s2=c(1:4), s3=c(5:8)) )
char.data<-c.dt$chars
M.initial<-make_init_matrix(part_scheme=c(4,4), set.NA=NULL )
result<-matrix(0, nrow=N.iter, ncol=4)
N.iter<-10

rhf<-run_hifi(M.initial, phy.data, char.data, deltaAIC=10,  N.iter=10,
                   do.thorough=TRUE, root.p="maddfitz",  Cor.sol=TRUE)

rhf2<-run_hifi(M.initial, phy.data, char.data, deltaAIC=10,  N.iter=5,
              do.thorough=TRUE, root.p="maddfitz",  Cor.sol=TRUE)

run_hifi<-function(M.initial, phy.data, char.data, deltaAIC=10, iter.remove=5, N.iter=50,
                   do.thorough=FALSE, root.p="maddfitz",  Cor.sol=TRUE)
{


  cat("************ INITIALIZATION... *******************", "\n")
  #cat("***********************************************", "\n")

  M<-list(M.initial) # all elemnts are phihifi lists
  M.r<-matrix(c(1:nrow(M[[1]]$M)^2), nrow(M[[1]]$M), nrow(M[[1]]$M) )

  # summarize results of search
  result<-matrix(NA, nrow=N.iter, ncol=8)
  colnames(result)<-c("iter", "Q_best", "Q_new", "Q_cheked",
                      "Q_HT=1", "Q_HT>1", "min(AIC)", "HT_min(AIC)")



# Initilalize
ptm <- proc.time()#
out<- greedy_rayDISC(phy=phy.data, data=char.data, rate.mat=make.cor.matrix(M[[1]]),
               root.p=root.p, do.ans=FALSE, do.thorough=do.thorough)
t.elapsed<-proc.time() - ptm #


best.Q<-phyhifi(M=M[[1]]$M, rate=out$solution, Ln=out$loglik, AIC=out$AIC, AICc=out$AICc,
                Cor.sol=Cor.sol, elapsed.time=t.elapsed[[1]])
best.Q<-list(best.Q)
Q.pool<-best.Q




#it<-1
for (it in 1:N.iter){

  cat(" \n")
  cat("***********************************************", "\n")
  #cat("-----------------------------------------------", "\n")
  cat("************ START of GENERATION", it, "************", "\n")

  #################
  # add move
  ################
  #cells.rnd<-lapply(M, function(x) add_move(x, M.r, N.trials))
  cells.rnd<-lapply(M, function(x) add_move_phf(x, M.r) )

  # construct possible matrices
  Q.add<-vector(mode = "list", length = length(cells.rnd))
  #i=1
  for (i in 1:length(cells.rnd) ){

    Q.add[[i]]<-M[[i]]
    Q.add[[i]]$M[ cells.rnd[[i]] ]<-1
    Q.add[[i]]$count<-1
  }

  #################
  # remove move
  ################
  if (it>iter.remove) # start remove move after Nth generation
  {
    cells.rnd<-lapply(M, function(x) remove_move_phf(x, M.r) )
    cells.rnd<-cells.rnd[!is.na(cells.rnd)]

    # construct possible matrices
    Q.rem<-vector(mode = "list", length = length(cells.rnd))
    #i=1
    for (i in 1:length(cells.rnd) ){
      Q.rem[[i]]<-M[[i]]
      Q.rem[[i]]$M[ cells.rnd[[i]] ]<-0
      Q.rem[[i]]$count<-1

      # Q.rem[[i]]<-M[[i]]
      # Q.rem[[i]][ cells.rnd[[i]] ]<-0
    }

  } else {
    Q.rem<-NULL
  }
  ###################################
  #################
  # two mves together
  ################

  Q.add<-c(Q.add, Q.rem)

  # add unmutated chimera
  if (length(M) > 1)
  {
    Q.add<-c(Q.add, chimera_matrix(M) )
  }

  unm<-unique_matrix_phf2(best.Q, Q.add, count=TRUE)
  best.Q<-unm[[1]]
  Q.add<-unm[[2]]

  # compare with the pool
  unm<-unique_matrix_phf2(Q.pool, Q.add, count=TRUE)
  Q.pool<-unm[[1]]
  Q.add<-unm[[2]]

  cat("Number of of newly evolved matrices: ", length(Q.add), "\n")
  cat(" \n")
  #cat("___________________________________________________________\n")
  #cat("***********************************************", "\n")


  for (i in 1:length(Q.add) ){
    #cat("___________________________________________________________\n")
    cat("Working on matrix:", i, "\n")
    cat("\n")

    ptm<-proc.time() #
    out<- greedy_rayDISC(phy=phy.data, data=char.data, rate.mat=make.cor.matrix(Q.add[[i]]),
                         root.p=root.p, do.ans=FALSE, do.thorough=do.thorough)
    t.elapsed<-proc.time() - ptm #

    sol<-c(out$solution)
    sol<-sol[which(sol>0)[1]]

    Q.add[[i]]$rate<-sol
    Q.add[[i]]$Ln<-out$loglik
    Q.add[[i]]$AIC<-out$AIC
    Q.add[[i]]$AICc<-out$AICc

    if (Cor.sol)
    Q.add[[i]]$Cor.sol<-out$solution

    Q.add[[i]]$elapsed.time<-t.elapsed[[1]] # not sure if this elsapsed time is indeed correct

    print(Q.add[[i]])
    cat("___________________________________________________________\n")


  }


  ###########

  #####
  # Select
  ####
  #AIC=2k-2Ln(L)
  # remove those AIC that do not fall within dAIC=10

  # current generation
  #Ln<-c(best.Ln, add.Ln, rem.Ln)
  #Q<-c(best.Q, Q.add, Q.rem)

  Q<-c(best.Q, Q.add)
  un.results<-unique_matrix_phf(Q, count=TRUE)
  #un.results<-unique_matrix(Q, Ln)

  # all pool of matrices
  Q.pool<- c(Q.pool, Q.add)
  Q.pool<-unique_matrix_phf(Q.pool, count=TRUE)


  aic<-lapply(un.results, function(x) x$AIC) %>% unlist()
  daic<-min(aic)+deltaAIC
  best.id<-which(aic<=daic )
  best.Q<-un.results[best.id]

  # summary
  result[it, 1]<-it
  result[it, 2]<-length(best.Q)
  result[it, 3]<-length(Q.add)
  result[it, 4]<-length(Q.pool)
  ################################################

  # slect unique matrices
  M<-best.Q

  # add chimera marix to best matrices
  if (length(M) > 1)
  {
    M<-c(M, chimera_matrix(M) )
    # leave only unqieu matrices
    M<-unique_matrix_phf(M, count=TRUE)
  }



  lapply(best.Q, function(x) x$count) %>% unlist() %>% table() -> tb.count
  q.once<-tb.count[names(tb.count)==1] %>% set_names(., NULL)
  q.Gonce<-tb.count[names(tb.count)>1] %>% set_names(., NULL) %>% sum()

  result[it, 5]<-q.once
  result[it, 6]<-q.Gonce
  result[it, 7]<-min(aic)[1]
  result[it, 8]<-length( which(aic==min(aic)) )


  # # print statistics
  #cat("************ END of GENERATION", it, "************", "\n")
  cat( "\n" )
  cat("STATISTICS OF GENERATION", it,  "\n")
  #cat("___________________________________________________________\n")

  print(result[it, 1:4])
  cat( "\n" )
  cat("Hit score, HT:", "\t", names(tb.count), "\n")
  cat("N of matrices:", "\t", set_names(tb.count, NULL), "\n" )
  cat( "\n" )
  cat("Matrices (HT=1)= ", q.once, "\t", "Matrices (HT>1)= ", q.Gonce, "\n", sep="")
  cat( "\n" )
  cat("Min(AIC)= ", min(aic)[1], " \t", "HT of Min(AIC)= ", length(min(aic)), "\n", sep="")
  cat("___________________________________________________________\n")

  #cat("*********************************************", "\n")
  #cat("*********************************************", "\n")


} # end loop

final.out<-list(best.Q=best.Q, Q.pool=Q.pool, summary=result)
return(final.out)

} #end run_hifi


