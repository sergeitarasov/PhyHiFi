# package Dependencies
devtools::use_package("corHMM", type = "Depends")
#library("magrittr", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
#library("magrittr", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
#library("numbers", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
#library("MASS", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
#library("plyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

#############
#devtools::document()

# object PhyHifi

#####################
# Init obj
#####################
#' @title Create PhyHiFi object
#' @description
#' @param M matrix
#' @return vector of names.
#' @examples
#' phyhifi()
#' @export

phyhifi<-function(M=NULL, rate=NULL, Ln=NULL, AIC=NULL, AICc=NULL, elapsed.time=NULL, Cor.sol=FALSE, count=1  )
{
  # get solution from corr
  if (Cor.sol)
    Cor.sol<-rate

  # get rate
  sol<-c(rate)
  rate<-sol[which(sol>0)[1]]

  phf<- list(
    M=M,
    rate=rate,
    Ln=Ln,
    AIC=AIC,
    AICc=AICc,
    Cor.sol=Cor.sol,
    elapsed.time=elapsed.time,
    count=count
  )
  class(phf) <- append(class(phf),"phyhifi")
  return(phf)
}

#' @title Print PhyHiFi object
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export
#'
print.phyhifi<-function(x)
{
  #cat("___________________________________________________________\n")
  cat("\n")
  cat("Ln:", x$Ln, "\tAIC:", x$AIC, "\tRate:", x$rate, "\n", sep=" ")
  cat("\n")
  print(x$M)
  cat("\n")
  cat("elapsed time:", x$elapsed.time, "\tcount:", x$count, sep=" ", "\n")

  #cat("___________________________________________________________\n")
}
#print(ph)
#ph<-phyhifi(M=M[[1]], rate=out$solution, Ln=out$loglik, AIC=out$AIC, AICc=out$AICc, Cor.sol=F)

###########################
# get non absorbing states
############################
#' @title Get non-absorbing states
#' @description Get non-absorbing states
#' @param M.phf matrix
#' @return vector of names.
#' @examples
#' @export

#M.phf<-ph
get_non_absrb<-function(M.phf)
{
  M<-M.phf$M
  M[which(M==0)]<-NA
  out<-apply(M, 1, function(x) !all(is.na(x)) )
  return(out)
}
#get_non_absrb(ph)

#################
# add move
################

#' @title "Add" move
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export

add_move_phf<-function(M.phf, M.r)
{
  M<-M.phf$M
  # get cells to sample
  non.absrb<-get_non_absrb(M.phf)
  potential.cells<-c(M.r[!is.na(M) ], M.r[,non.absrb]) # non-absr and non NA,
  potential.cells<-potential.cells[!potential.cells%in%which(M==1)] #remove cells equal 1
  to.sample<-potential.cells[duplicated(potential.cells)]

  cells.rnd<-sample(to.sample, 1)
  return(cells.rnd)

}
#add_move_phf(ph, M.r)


#################
# remove move
################

#' @title "Remove" move
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export

remove_move_phf<-function(M.phf, M.r)
{
  M<-M.phf$M
  # slect those rows where a remove does not create absorbing model
  r.nonabsrb<-apply(M, 1, function(x) which(x==1) %>% length() )

  if ( length(M.r[r.nonabsrb>1,])>0 )
  {
    to.sample<-c(M.r[r.nonabsrb>1,], M.r[which(M==1)])
    to.sample<-to.sample[duplicated(to.sample) ]
    #to.sample<-M.r[which(M==1)]
    cells.rnd<-sample(to.sample, 1)
    return(cells.rnd)
  } else
  {
    return(NA)
  }

}
#remove_move_phf(ph, M.r)


#################
# select unique matrices
################

#' @title select unique matrices
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export
 # M<-Q
# M.phf.list<-M
# M.phf.list[[3]]<-Q[[2]]
# unique_matrix_phf(M.phf.list)

unique_matrix_phf<-function(M.phf.list, count=FALSE)
{
  M.trans<-lapply(M.phf.list, function(x) x$M)

  #M.trans<-M
  for (i in 1:length(M.trans)) # get rid of NAs
  {
    M.trans[[i]][is.na(M.trans[[i]])]<-0
  }

  duplic<-c()
  for (i in 1:(length(M.trans)-1) ) # get rid of NAs
  {

    for (j in (i+1):length(M.trans))
    {
      if (all(M.trans[[i]]==M.trans[[j]]))
      {
        duplic<-c(duplic, j)

        if (count)
        {
          # increment count id object is duplicatedy
          M.phf.list[[i]]$count<-M.phf.list[[i]]$count+1
        }

      }
    }

  }


  duplic<-unique(duplic)

  if (is.null(duplic))
  {
    return(M.phf.list)
  } else
  {
    return(M.phf.list[-duplic])
  }

}
#unique_matrix(M, Ln)

#################
# select unique matrices given two sets: check M.phf.list2 against M.phf.list1
################
#' @title select unique matrices given two sets: check M.phf.list2 against M.phf.list1
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export

# M.phf.list1<-M
# M.phf.list2<-c(M, M)
# unique_matrix_phf2(M.phf.list1, M.phf.list2, count=TRUE)

unique_matrix_phf2<-function(M.phf.list1, M.phf.list2, count=FALSE)
{
  list.out<-vector(mode = "list", length = 3)

  M.trans1<-lapply(M.phf.list1, function(x) x$M)
  M.trans2<-lapply(M.phf.list2, function(x) x$M)

  #M.trans<-M
  for (i in 1:length(M.trans1)) # get rid of NAs
  {
    M.trans1[[i]][is.na(M.trans1[[i]])]<-0
  }

  #M.trans<-M
  for (i in 1:length(M.trans2)) # get rid of NAs
  {
    M.trans2[[i]][is.na(M.trans2[[i]])]<-0
  }



  duplic<-c()
  for (i in 1:length(M.trans1) ) # get rid of NAs
  {

    for (j in 1:length(M.trans2))
    {
      if (all(M.trans1[[i]]==M.trans2[[j]]))
      {
        duplic<-c(duplic, j)

        if (count)
        {
          # increment count id object is duplicatedy
          M.phf.list1[[i]]$count<-M.phf.list1[[i]]$count+1
        }

      }
    }

  }


  duplic.un<-unique(duplic)

  if (is.null(duplic.un))
  {
    list.out[[1]]<-M.phf.list1
    list.out[[2]]<-M.phf.list2
    list.out[[3]]<-0

    return(list.out)
  } else
  {
    list.out[[1]]<-M.phf.list1
    list.out[[2]]<-M.phf.list2[-duplic.un]
    list.out[[3]]<-duplic

    return(list.out)
  }

}

#################
# chimera matrix
################
#' @title Make chimera matrix
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export
#'
chimera_matrix<-function(M)
{
  M<-lapply(M, function(x) x$M)

  M.out<-M[[1]]
  for (i in 2:length(M))
  {
    M.out<-M.out + M[[i]]
  }

  M.out[which(M.out>1)]<-1
  M.out<-phyhifi(M=M.out)
  return(list(M.out) )
}
# chimera_matrix(M)
# chimera_matrix(Q)

#' @title # make cor rate matrix: set 0 to NA
#' @description
#' @param
#' @return
#' @examples
#' phyhifi()
#' @export
make.cor.matrix<-function(M)
{
  M<-M$M
  M[which(M==0)]<-NA
  return(M)
}

#make.cor.matrix(M)
###
