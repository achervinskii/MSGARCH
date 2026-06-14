f_CondVol <- function(object, par, data, do.its = FALSE, nahead = 1L, do.cumulative = FALSE, ctr = list(), ...) {
  object    <- f_check_spec(object)
  data_      <- f_check_y(data)
  par.check <- f_check_par(object, par)
  if (nrow(par.check) == 1) {
    ctr     <- f_process_ctr(ctr)
    nsim <- ctr$nsim
  } else {
    if(is.null(ctr$nsim)){
      nsim = 1
    } else {
      nsim = ctr$nsim
    }
  }
  ctr       <- f_process_ctr(ctr)
  variance  <- object$rcpp.func$calc_ht(par.check, data_)
  P         <- State(object = object,par = par, data = data_)
  PredProb  <- P$PredProb
  vol <- matrix(NA, nrow = dim(PredProb)[1], ncol = nrow(par.check))
  if (object$K == 1) {
    for (i in 1:nrow(par.check)) {
      vol[,i] <- variance[,i]
    }
  } else {
    for (i in 1:nrow(par.check)) {
      vol[,i] <- rowSums(PredProb[,i,] * variance[,i,])
    }
  }
  vol <- sqrt(vol)
  draw <- NULL
  vol_draw <- NULL
  if (!isTRUE(do.its)) {
    tmp    <- mean(vol[dim(PredProb)[1],])
    vol    <- vector(mode = "numeric", length = nahead)
    vol[1] <- tmp

    simul <- Sim(object = object, data = data, nahead = nahead,
                  nsim = nsim, par = par)
    draw <- simul$draw

    # save the volatility realizations from the simulations
    nsim_total <- ncol(draw)
    CondVol <- simul$CondVol
    # note that the state returned by Sim is already 1-based indexed
    state <- simul$state
    # safety checks
    stopifnot(dim(CondVol)[2] == nsim_total)
    stopifnot(dim(state)[2] == nsim_total)
    vol_draw <- matrix(data = NA, nrow = nahead, ncol =  nsim_total)
    for (h in 1:nahead) {
      for (sim in 1:nsim_total) {
        the_state <- state[h, sim]
        vol_draw[h, sim] <- CondVol[h, sim, the_state]
      }
    }

    if (nahead > 1) {
      if(isTRUE(do.cumulative)){
        draw = apply(draw, 2, cumsum)
      }
      vol[2:nahead] = apply(draw[2:nahead,, drop = FALSE], 1, sd)
    }
    names(vol) <- paste0("h=", 1:nahead)
    rownames(vol_draw) = paste0("h=",1:nahead)
    colnames(vol_draw) =  paste0("Sim #",1:nsim_total)
    if(zoo::is.zoo(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data)[length(data)]+(1:nahead))
    }
    if(is.ts(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data)[length(data)]+(1:nahead))
      vol = as.ts(vol)
    }
  } else {
    if (nrow(par.check) > 1) {
      vol <- rowMeans(vol[1:length(data_),])
      vol <- as.vector(vol)
    } else {
      vol <- as.vector(vol)
      vol <- vol[1:length(data_)]
    }
    names(vol) <- paste0("t=", 1:(length(data_)))
    if(zoo::is.zoo(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data))
    }
    if(is.ts(data)){
      vol = zoo::zooreg(vol, order.by = zoo::index(data))
      vol = as.ts(vol)
    }
  }
  out = list()
  class(vol) <- c("MSGARCH_CONDVOL",class(vol))
  out$vol = vol
  out$draw = draw
  out$vol_draw = vol_draw
  return(out)
}


