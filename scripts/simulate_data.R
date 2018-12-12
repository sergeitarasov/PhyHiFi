Q<-matrix(
  c(
    0,  0.02, 0.02,
    0.02, 0, 0.02,
    0.1, 0.02,  0
  ),
  nrow = 3, ncol = 3, byrow = T
)
diag(Q)<-apply(Q, 1, function(x) sum(x))*-1
colnames(Q)<-rownames(Q)<-c("1", "2", "3")


# simulate history
require("phytools")
tree<-pbtree(n=300, scale=200)
sim<-sim.history(tree, Q, nsim=1, anc=setNames(c(1, 0, 0), colnames(Q) ) )
plot(sim)

# matrix for inference
char.state<-c(0, 1)
rate.param<-c(1, 2)
TL<-init_char_matrix(char.state, rate.param, diag.as=0)


sim.radt[[1]]$states<-mapvalues(sim.radt[[1]]$states,
                                from = c("0", "1"), to=c("0&1&2", "3&4&5") )
