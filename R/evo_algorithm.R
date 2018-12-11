

# argiments
M.init
N.iter
deltaAIC #AIC range
root.phylo
root.est


run_hifi<-function()
{

#Q.cor[which(Q.cor==1)]<-0
M<-Q.cor
M.r<-matrix(c(1:nrow(M)^2), nrow(M), nrow(M) )


# set intial rates
M[1, obs.state[[2]][1]]<-1
M[obs.state[[2]][1], 1]<-1
M<-list(M)
M<-list(phyhifi(M=M[[1]]) )
##
result<-matrix(0, nrow=N.iter, ncol=4)
colnames(result)<-c("N_Q.best", "N_new_mt", "", "")


##############################
#
# Serach Loop
# PhyHiFi
#
##########################

# Initilalize
N.iter<-40

ptm <- proc.time()#
out<- greedy_rayDISC(tree.corr, c.dt, rate.mat=make.cor.matrix(M[[1]]), node.states="marginal", model="ARD",
                     root.p="maddfitz", do.ans=FALSE)
t.elapsed<-proc.time() - ptm #

best.Q<-phyhifi(M=M[[1]]$M, rate=out$solution, Ln=out$loglik, AIC=out$AIC, AICc=out$AICc,
                Cor.sol=T, elapsed.time=t.elapsed[[3]])
best.Q<-list(best.Q)
Q.pool<-best.Q
# best.Ln <-out$loglik
# best.Q <-M


it<-1
for (it in 1:N.iter){

  print(paste0("Generation ", it) )
  result[it, 1]<-length(best.Q)

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
  if (it>5) # start remove move after Nth generation
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

  cat("N of newly evolved matrices: ", length(Q.add), "\n")
  result[it, 2]<-length(Q.add)

  # do inference for each matrix
  #add.Ln<-c()
  #add.rate<-c()
  for (i in 1:length(Q.add) ){
    cat("Generation ", it, ";  Add move, matrix ", i, "\n", "-------------------------------------", "\n")

    t.elapsed<-proc.time() - ptm #
    out<- greedy_rayDISC(tree.corr, c.dt, rate.mat=make.cor.matrix(Q.add[[i]]), node.states="marginal", model="ARD",
                         root.p="maddfitz", do.ans=FALSE)
    t.elapsed<-proc.time() - ptm #

    sol<-c(out$solution)
    sol<-sol[which(sol>0)[1]]

    Q.add[[i]]$rate<-sol
    Q.add[[i]]$Ln<-out$loglik
    Q.add[[i]]$AIC<-out$AIC
    Q.add[[i]]$AICc<-out$AICc
    Q.add[[i]]$Cor.sol<-out$solution
    Q.add[[i]]$elapsed.time<-t.elapsed[[3]] # not sure if this elsapsed time is indeed correct

    print(Q.add[[i]])


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

  #Q<-un.results[[1]]
  #Ln<-un.results[[2]]

  #aic=2-(2*Ln)

  aic<-lapply(un.results, function(x) x$AIC) %>% unlist()
  daic<-min(aic)+deltaAIC
  best.id<-which(aic<=daic )
  best.Q<-un.results[best.id]

  #best.Ln<-Ln[ best.id]
  #best.Q<-Q[best.id]
  ################################################

  # set M to be best Q
  #M.Ln<-max(Ln)[1]

  # slect unique matrices
  M<-best.Q

  # add chimera marix to best matrices
  if (length(M) > 1)
  {
    M<-c(M, chimera_matrix(M) )
    # leave only unqieu matrices
    M<-unique_matrix_phf(M, count=TRUE)
  }



  # # print best Q
  cat("----------------------------------------------------\n")
  cat("Matrix distribution:\n")
  lapply(best.Q, function(x) x$count) %>% unlist() %>% table() %>% print()
  # cat("Generation ", it, "; Best matrix, ML=", M.Ln, " :")
  # print(Q[which(Ln==M.Ln)[1] ])
  #
  # cat("N of best matrices is ", length(M), "\n" )
  cat("----------------------------------------------------\n")
  #



} # end loop


summ_matrix(best.Q)
summ_matrix(Q.pool)

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

}#end
