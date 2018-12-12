# Modification of rayDIS() function from corHMM package

#' @title Mod Ray Disc
#' @description
#' @param M matrix
#' @return vector of names.
#' @examples
#' phyhifi()
#' @export
greedy_rayDISC <- function(phy, data, ntraits = 1, charnum = 1, rate.mat = NULL,
                             model = c("ARD"), node.states = c("marginal"), state.recon = c("subsequently"), lewis.asc.bias = FALSE,
                             p = NULL, root.p = NULL, ip = NULL, lb = 0, ub = 100, verbose = TRUE,
                             diagn = FALSE, do.ans=FALSE, do.thorough=FALSE){

    if (is.null(node.states)) {
      obj <- NULL
      obj$loglik <- NULL
      obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to rayDISC command for parameter 'node.states': joint, marginal, or scaled.")
      return(obj)
    }


    else {


      valid.models <- c("joint", "marginal", "scaled")
      if (!any(valid.models == node.states)) {
        obj <- NULL
        obj$loglik <- NULL
        obj$diagnostic <- paste("'", node.states, "' is not valid for ancestral state reconstruction method.  Please pass one of the following to rayDISC command for parameter 'node.states': joint, marginal, or scaled.",
                                sep = "")
        return(obj)
      }
      if (length(node.states) > 1) {
        node.states <- "marginal"
        cat("No model selected for 'node.states'. Will perform marginal ancestral state estimation.\n")
      }
    }


    if (!state.recon == "subsequently" & node.states == "marginal" |
        node.states == "scaled") {
      stop("Simultaneous estimation of rates and states using either marginal or scaled probabilities not yet implemented.",
           call. = FALSE)
    }
    if (!state.recon == "subsequently") {
      if (!is.null(phy$node.label)) {
        if (!is.na(phy$node.label[Ntip(phy) + 1])) {
          root.p <- NULL
        }
      }
    }

    ## STARTS HERE

    phy$edge.length[phy$edge.length == 0] = 1e-05
    matching <- corHMM:::match.tree.data(phy, data)
    data <- matching$data
    phy <- matching$phy

    # F
    if (nlevels(as.factor(data[, charnum + 1])) <= 1) {
      obj <- NULL
      obj$loglik <- NULL
      obj$diagnostic <- paste("Character ", charnum, " is invariant. Analysis stopped.",
                              sep = "")
      return(obj)
    } else {
      lvls <- as.factor(data[, charnum + 1])

      #F
      if (nlevels(as.factor(data[, charnum + 1])) == 2 && length(which(lvls ==
                                                                       "?"))) {
        obj <- NULL
        obj$loglik <- NULL
        obj$diagnostic <- paste("Character ", charnum, " is invariant. Analysis stopped.",
                                sep = "")
        return(obj)
      }
    }

    workingData <- data.frame(data[, charnum + 1], data[, charnum +
                                                          1], row.names = data[, 1])
    workingData <- workingData[phy$tip.label, ]
    counts <- table(workingData[, 1])
    levels <- levels(as.factor(workingData[, 1]))
    cols <- as.factor(workingData[, 1])
    if (verbose == TRUE) {
      cat("State distribution in data:\n")
      cat("States:", levels, "\n", sep = "\t")
      cat("Counts:", counts, "\n", sep = "\t")
    }

    k <- 1
    factored <- corHMM:::factorData(workingData, charnum = charnum)
    nl <- ncol(factored)
    state.names <- colnames(factored)
    bound.hit <- FALSE
    if (ub < 0) {
      ub <- 100
    }
    if (lb < 0) {
      lb <- 0
    }
    if (ub < lb) {
      ub <- 100
      lb <- 0
    }
    obj <- NULL
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    model = model
    root.p = root.p
    ip = ip
    model.set.final <- corHMM:::rate.cat.set.rayDISC(phy = phy, data = workingData,
                                                     model = model, charnum = charnum)
    if (!is.null(rate.mat)) {
      rate <- rate.mat
      model.set.final$np <- max(rate, na.rm = TRUE)
      rate[is.na(rate)] = max(rate, na.rm = TRUE) + 1
      model.set.final$rate <- rate
      model.set.final$index.matrix <- rate.mat
    }

    lower = rep(lb, model.set.final$np)
    upper = rep(ub, model.set.final$np)
    opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
                 ftol_rel = .Machine$double.eps^0.5)

    #### False, a vector of transition rates
    if (!is.null(p)) {
      if (verbose == TRUE) {
        cat("Calculating likelihood from a set of fixed parameters",
            "\n")
      }
      out <- NULL
      out$solution <- p
      if (state.recon == "subsequently") {
        out$objective <- dev.raydisc(out$solution, phy = phy,
                                     liks = model.set.final$liks, Q = model.set.final$Q,
                                     rate = model.set.final$rate, root.p = root.p,
                                     lewis.asc.bias = lewis.asc.bias)
        loglik <- -out$objective
      }
      else {
        if (lewis.asc.bias == TRUE) {
          loglik.num <- ancRECON(phy = phy, data = data,
                                 p = p, hrm = FALSE, rate.cat = NULL, rate.mat = rate.mat,
                                 ntraits = ntraits, method = node.states, model = model,
                                 charnum = charnum, root.p = root.p, get.likelihood = TRUE)
          phy.dummy <- phy
          data.dummy <- cbind(phy$tip.label, 0)
          phy.dummy$node.label <- rep(1, length(phy.dummy$node.label))
          loglik.dummy <- ancRECON(phy = phy, data = data,
                                   p = p, hrm = FALSE, rate.cat = NULL, rate.mat = rate.mat,
                                   ntraits = ntraits, method = node.states, model = model,
                                   charnum = charnum, root.p = root.p, get.likelihood = TRUE)
          loglik <- (loglik.num - log(1 - exp(loglik.dummy)))
          loglik <- out$objective
        }
        else {
          out$objective <- ancRECON(phy = phy, data = data,
                                    p = p, hrm = FALSE, rate.cat = NULL, rate.mat = rate.mat,
                                    ntraits = ntraits, method = node.states, model = model,
                                    charnum = charnum, root.p = root.p, get.likelihood = TRUE)
          loglik <- out$objective
        }
      }
      est.pars <- out$solution
    }

    #######


    else {
      if (is.null(ip)) {
        if (verbose == TRUE) {
          cat("Initializing...", "\n")
        }
        model.set.init <- corHMM:::rate.cat.set.rayDISC(phy = phy,
                                                        data = workingData, model = "ER", charnum = charnum)
        opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
                     ftol_rel = .Machine$double.eps^0.5)
        dat <- as.matrix(workingData)
        dat <- phangorn::phyDat(dat, type = "USER", levels = levels(as.factor(workingData[,
                                                                                          1])))
        #str(dat)
        par.score <- phangorn::parsimony(phy, dat, method = "fitch")
        tl <- sum(phy$edge.length)
        mean.change = par.score/tl
        if (mean.change == 0) {
          ip = 0.01 + lb
        }
        else {
          ip <- rexp(1, 1/mean.change)
        }
        if (ip < lb || ip > ub) {
          ip <- lb
        }
        lower.init = rep(lb, model.set.init$np)
        upper.init = rep(ub, model.set.init$np)

        #ptm <- proc.time() #
        init = nloptr(x0 = rep(ip, length.out = model.set.init$np), ###################### Init
                      eval_f = dev.raydisc, lb = lower.init, ub = upper.init,
                      opts = opts, phy = phy, liks = model.set.init$liks,
                      Q = model.set.init$Q, rate = model.set.init$rate,
                      root.p = root.p, lewis.asc.bias = lewis.asc.bias)

        out<-init

        if (do.thorough) # thorough search is disabled by default
        {
          if (verbose == TRUE) {
            cat("Finished. Beginning thorough search...",
                "\n")
          }

          lower = rep(lb, model.set.final$np)
          upper = rep(ub, model.set.final$np)
          if (state.recon == "subsequently") {

            out <- nloptr(x0 = rep(init$solution, length.out = model.set.final$np), ######## Inference
                          eval_f = dev.raydisc, lb = lower, ub = upper,
                          opts = opts, phy = phy, liks = model.set.final$liks,
                          Q = model.set.final$Q, rate = model.set.final$rate,
                          root.p = root.p, lewis.asc.bias = lewis.asc.bias)



          } else
          {

            out <- nloptr(x0 = rep(init$solution, length.out = model.set.final$np), ##########
                          eval_f = dev.raydisc.rates.and.states, lb = lower,
                          ub = upper, opts = opts, phy = phy, data = data,
                          hrm = FALSE, rate.cat = NULL, rate.mat = rate.mat,
                          ntraits = ntraits, method = node.states, model = model,
                          charnum = charnum, root.p = root.p, lewis.asc.bias = lewis.asc.bias,
                          get.likelihood = TRUE)

          }
        } # thorough search

        loglik <- -out$objective
        est.pars <- out$solution
      } else {############# ip specified
        if (verbose == TRUE) {
          cat("Beginning subplex optimization routine -- Starting value(s):",
              ip, "\n")
        }
        opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
                     ftol_rel = .Machine$double.eps^0.5)
        if (state.recon == "subsequently") {
          out <- nloptr(x0 = rep(init$solution, length.out = model.set.final$np),
                        eval_f = dev.raydisc, lb = lower, ub = upper,
                        opts = opts, phy = phy, liks = model.set.final$liks,
                        Q = model.set.final$Q, rate = model.set.final$rate,
                        root.p = root.p, lewis.asc.bias = lewis.asc.bias)
        } else {
          out <- nloptr(x0 = rep(init$solution, length.out = model.set.final$np),
                        eval_f = dev.raydisc.rates.and.states, lb = lower,
                        ub = upper, opts = opts, phy = phy, data = data,
                        hrm = FALSE, rate.cat = NULL, rate.mat = rate.mat,
                        ntraits = ntraits, method = node.states, model = model,
                        charnum = charnum, root.p = root.p, lewis.asc.bias = lewis.asc.bias,
                        get.likelihood = TRUE)
        }
        loglik <- -out$objective
        est.pars <- out$solution
      } ############# END ip specified

    }

    if (do.ans == TRUE) {
      ############################### Inferring ancestral states using
      if (verbose == TRUE) {
        cat("Finished. Inferring ancestral states using", node.states,
            "reconstruction.", "\n")
      }
      TIPS <- 1:nb.tip
      if (node.states == "marginal" || node.states == "scaled") {
        lik.anc <- ancRECON(phy, data, est.pars, hrm = FALSE,
                            rate.cat = NULL, rate.mat = rate.mat, ntraits = ntraits,
                            method = node.states, model = model, charnum = charnum,
                            root.p = root.p)
        pr <- apply(lik.anc$lik.anc.states, 1, which.max)
        phy$node.label <- pr
        tip.states <- lik.anc$lik.tip.states
      }
      if (!state.recon == "given") {
        if (node.states == "joint") {
          lik.anc <- ancRECON(phy, data, est.pars, hrm = FALSE,
                              rate.cat = NULL, rate.mat = rate.mat, ntraits = ntraits,
                              method = node.states, model = model, charnum = charnum,
                              root.p = root.p)
          phy$node.label <- lik.anc$lik.anc.states
          tip.states <- lik.anc$lik.tip.states
        }
      }else {
        lik.anc <- NULL
        lik.anc$lik.anc.states <- phy$node.label
        lik.anc$lik.tip.states <- workingData[, 1]
        tip.states <- lik.anc$lik.tip.states
      }

    } else
    {
      lik.anc <- NULL
      lik.anc$lik.anc.states <- NULL
      lik.anc$lik.tip.states <- NULL
      tip.states <- lik.anc$lik.tip.states
    }############################### End Inferring ancestral states using

    if (diagn == TRUE) {
      if (verbose == TRUE) {
        cat("Finished. Performing diagnostic tests.", "\n")
      }
      h <- hessian(func = dev.raydisc, x = est.pars, phy = phy,
                   liks = model.set.final$liks, Q = model.set.final$Q,
                   rate = model.set.final$rate, root.p = root.p)
      solution <- matrix(est.pars[model.set.final$index.matrix],
                         dim(model.set.final$index.matrix))
      solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix],
                            dim(model.set.final$index.matrix))
      hess.eig <- eigen(h, symmetric = TRUE)
      eigval <- signif(hess.eig$values, 2)
      eigvect <- round(hess.eig$vectors, 2)
    } else {
      solution <- matrix(est.pars[model.set.final$index.matrix],
                         dim(model.set.final$index.matrix))
      solution.se <- matrix(0, dim(solution)[1], dim(solution)[1])
      eigval <- NULL
      eigvect <- NULL
    }
    if ((any(solution == lb, na.rm = TRUE) || any(solution ==
                                                  ub, na.rm = TRUE)) && (lb != 0 || ub != 100)) {
      bound.hit <- TRUE
    }

    rownames(solution) <- rownames(solution.se) <- state.names
    colnames(solution) <- colnames(solution.se) <- state.names

    if (do.ans == TRUE) { ## Mine
      if (is.character(node.states)) {
        if (node.states == "marginal" || node.states == "scaled") {
          colnames(lik.anc$lik.anc.states) <- state.names
        }
      }
    }

    obj = list(loglik = loglik, AIC = -2 * loglik + 2 * model.set.final$np,
               AICc = -2 * loglik + (2 * model.set.final$np * (nb.tip/(nb.tip -
                                                                         model.set.final$np - 1))), ntraits = 1, solution = solution,
               solution.se = solution.se, index.mat = model.set.final$index.matrix,
               lewis.asc.bias = lewis.asc.bias, opts = opts, data = data,
               phy = phy, states = lik.anc$lik.anc.states, tip.states = tip.states,
               iterations = out$iterations, eigval = eigval, eigvect = eigvect,
               bound.hit = bound.hit)
    if (!is.null(matching$message.data)) {
      obj$message.data <- matching$message.data
      obj$data <- matching$data
    }
    if (!is.null(matching$message.tree)) {
      obj$message.tree <- matching$message.tree
      obj$data <- matching$data
    }
    class(obj) <- "raydisc"
    return(obj)
  }

