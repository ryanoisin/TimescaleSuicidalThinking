# o.ryan@uu.nl; January 2021

# -------------------------------------------------------------------------------
# Function which extracts burst measurements from a list of person-specific data
# -------------------------------------------------------------------------------

get_burst <- function(ldata, # list of length n containing observations for each ID
                      include_end = TRUE, # include the "end-points" or burst observations?
                      end_tol = 24 # if include_end = TRUE, how long must pass after the last burst to omit an "end-point"
                      ){

  bldata <- list()

  for(l in 1:length(ldata)){
    tmp <- ldata[[l]]
    ty <- tmp$Type
    bcntr <- 0
    i <- 1
    # Go through type and find pairs of "EMA" followed by "Burst"
    while(i < length(ty)){

      if(ty[i] == "EMA"){ i <- i + 1 } else {
        start <- i-1 # start Burst (coded as ema)
        end <- i + 1
        end_ind <- FALSE

        # find the "end" of the burst measurements
        while(!end_ind){
          # start counting forward from j until you find "EMA" again
          if(ty[end] == "EMA" | end == length(ty)){
            end_ind <- TRUE }else{end <- end + 1}
        }
        bcntr <- bcntr + 1
        i <- end
        if(!include_end){ # Include end point or not?
          end <- end-1
        }
        outm <- tmp[start:end,]
        nb <- nrow(outm)
        outm$Type2 <- c("Start",rep("Burst",nb-2),"End")
        outm$burst_episode <- bcntr
        outm$SeqTime[1] <- 0

        # Optional -  filter out "end" EMA segments that are above a certain threshold
        # e.g., more than the fixed interval length of regular EMA measurements
        if(include_end & !is.na(end_tol)){
         if(outm$SeqTime[nb] > end_tol){ outm <- outm[-nb,] } }

        # Collect with other bursts
        if(bcntr == 1){
          out <- outm } else { out <- rbind(out,outm) }
        i <- i + 1
      }
    }
    # write to list for each individual
    bldata[[l]] <- out
  }

  bldata
}


# -------------------------------------------------------------------------------
#------------ Function to plot Phi for a range of delta t -----------------------
# -------------------------------------------------------------------------------


get_phidt_ar <- function(drift,dts){
  p <- ncol(drift)
  out <- vapply(
    lapply(
      dts,
      FUN = function(inter)
        as.matrix(expm::expm(drift * inter))
    ),
    identity, matrix(0, p, p))
  return(out)}


# -------------------------------------------------------------------------------
# --------- Plot model-implied lagged effects (VAR)------------------------------
# -------------------------------------------------------------------------------


phiplot <- function(drift = mean_est, #  drift matrix
                    driftl = driftl, # list of drift matrices sampled from posterior
                    maxdt = 2.5, # Maximum time-interval
                    step = .01,  # decreasing this number makes the phi-dt lines smoother
                    add = FALSE, # add lines to existing plot? Default is no (add = FALSE)
                    ylim = NULL, # optional: supply a custom y-axis range
                    legtune = 0, # tuning parameter for position of the legend
                    lwd = 3, # thickness of the plotted lines
                    addlegend = TRUE, # add a legend or not?
                    main = NULL,
                    colvec = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
                    typevec = rep(1,4)
){
  p <- nrow(drift)

  # Create vector of titles
  if(is.null(main)){
      main <- c("Lagged Parameter vs Time Interval")
  }

  # vector of time-intervals
  dts <- seq(0,maxdt,step)

  # Get array of estimated dt parameters over time-intervals
  mean_ar  <- get_phidt_ar(drift,dts)

  # 4d array with dimensions dimensions [row, column,time-interval, sampled drift matrix]
  lout <- abind(lapply(driftl, get_phidt_ar,dts = dts),along = 4)

  # ok, now i want to get, for each dt and each element the 95 CIs
  lower <-  upper <-  array(dim = dim(mean_ar))
  for(i in 1:2){
    for(j in 1:2){
      for(t in 1:length(dts)){
      lower[i,j,t] <- quantile(lout[i,j,t,],probs = .025)
      upper[i,j,t] <- quantile(lout[i,j,t,],probs = .975)
    }
  }
}
  # get y-axis limits if not specified
  # if(is.null(ylim)){
  #   ylim <- pretty(range(c(blank)))
  #   ylim <- c(ylim[1],ylim[length(ylim)])
  # }

  #### Lagged parameter plot for auto-regressive effects
    plot.new()
    plot.window(xlim=c(0,maxdt), ylim=ylim)
    axis(2); title(ylab = expression(paste(Phi,"(", Delta, "t) values")), line=2.5)
    axis(1); title(xlab = expression(paste("Time-interval (", Delta, "t)", sep="")),
                   line=2.5)
    abline(h=0)
    title(main = main[1])

    lines(y= mean_ar[1,1,],x=dts,col=colvec[1],lwd=lwd, lty=typevec[1])
    lines(y= mean_ar[2,1,],x=dts,col=colvec[2],lwd=lwd, lty=typevec[2])
    lines(y=mean_ar[1,2,],x=dts,col=colvec[3],lwd=lwd, lty=typevec[3])
    lines(y=mean_ar[2,2,],x=dts,col=colvec[4],lwd=lwd, lty=typevec[4])

    lines(y=lower[1,1,],x=dts,col=colvec[1],lwd=1, lty=2)
    lines(y=lower[2,1,],x=dts,col=colvec[2],lwd=1, lty=2)
    lines(y=lower[1,2,],x=dts,col=colvec[3],lwd=1, lty=2)
    lines(y=lower[2,2,],x=dts,col=colvec[4],lwd=1, lty=2)

    lines(y=upper[1,1,],x=dts,col=colvec[1],lwd=1, lty=2)
    lines(y=upper[2,1,],x=dts,col=colvec[2],lwd=1, lty=2)
    lines(y=upper[1,2,],x=dts,col=colvec[3],lwd=1, lty=2)
    lines(y=upper[2,2,],x=dts,col=colvec[4],lwd=1, lty=2)

  # Add legend if making a new plot
    if(addlegend){
      leg <- c("Des -> Des", "Des -> Int", "Int -> Des", "Int -> Int")

      l1obj <- legend("topright",
                      leg,
                      lty = typevec,
                      col = colvec,
                      lwd = lwd)
    }

} # end of Function

# -------------------------------------------------------------------------------
## ------------------ Function for plotting Time series
# -------------------------------------------------------------------------------
plot_sts <- function(datasel, xlim= c(0,max(data$RunTime))){
  plot.new()
  plot.window(xlim = xlim, ylim = c(0,10))
  axis(1); axis(2)
  lines(datasel$RunTime, datasel$SI_DesireKill, col = "red", type = "b")
  lines(datasel$RunTime, datasel$SI_Intent, col = "blue", type = "b")
}

# -------------------------------------------------------------------------------
## -------- Function for cleaning/prepping data for ctsem -----------------------
# -------------------------------------------------------------------------------

get_ctdata <- function(datasel, standardize = TRUE){

  t1 <- as.POSIXct(datasel[,2],format="%Y-%m-%d %H:%M:%OS")
  time <- as.numeric(difftime(t1,t1[1], units="hours"))

  # Create SeqTime for later
  SeqTime <- as.numeric(difftime(t1[-1],t1[-length(t1)], units="hours"))

  # select only relevant variables
  data_esm <- datasel[,c("SI_DesireKill","SI_Intent")]

  # Scale
  if(isTRUE(standardize)){
  data_esm <- apply(data_esm, 2, scale)
  }
  data_esm <- cbind(data_esm,c(0,SeqTime))
  colnames(data_esm)[3] <- "TI"

  # create ID variable
  id <- rep(1,nrow(data_esm))

  # Create long form dataset
  data_long <- cbind(id, time, data_esm[,c("SI_DesireKill","SI_Intent")])

  out <- list(data_esm,data_long)
  names(out) <- c("data_esm","data_long")
  out
}

# -------------------------------------------------------------------------------
## --------------- FunctionS to help extract ctsem output -----------------------
# -------------------------------------------------------------------------------

get_mean <- function(post){
  matrix(c(mean(post$pop_DRIFT[,1,1]),mean(post$pop_DRIFT[,1,2]),
           mean(post$pop_DRIFT[,2,1]),mean(post$pop_DRIFT[,2,2])),
         2,2,byrow=T)
}

get_driftl <- function(post){
  driftl <- list()
  for(i in 1:1000){
    driftl[[i]] <- matrix(c(post$pop_DRIFT[i,1,1],post$pop_DRIFT[i,1,2],
                            post$pop_DRIFT[i,2,1],post$pop_DRIFT[i,2,2]),
                          2,2,byrow=T)
  }
  driftl
}


get_ctsub <- function(datasel){
  # scale before taking the subset to make predictions comparable

  data_sc <- apply(datasel[,c("SI_DesireKill","SI_Intent")],2,scale)
  data_sc <- data_sc[datasel$Type!="Burst",]

  t1sub <- as.POSIXct(datasel[datasel$Type!="Burst",2],format="%Y-%m-%d %H:%M:%OS")
  timesub <- as.numeric(difftime(t1sub ,t1sub [1], units="hours"))

  # Create SeqTime for later
  SeqTimesub <- as.numeric(difftime(t1sub[-1],t1sub[-length(t1sub)], units="hours"))

  # select only relevant variables and scale
  data_esm_sub <- cbind(data_sc,c(0,SeqTimesub))
  colnames(data_esm_sub)[3] <- "TI"

  # create ID variable
  id <- rep(1,nrow(data_esm_sub))

  # Create long form dataset
  data_long_sub <- cbind(id, timesub, data_esm_sub[,c("SI_DesireKill","SI_Intent")])
  colnames(data_long_sub)[2] <- "time"

  out <- list(data_esm_sub,data_long_sub)
  names(out) <- c("data_esm","data_long")
  out
}


# -------------------------------------------------------------------------------
## -------- Function to detect elevated leves of a variable ---------------------
# -------------------------------------------------------------------------------

detectelev <- function (x, RunTime, Type, Time) {
  library(data.table)
  x[is.na(x)] <- 0
  n <- length(x)
  ind_present <- x > 0
  counter <- 1
  ind_label <- rep(NA,n)
  ind_label[1] <- ind_present[1]
  for (p in 2:n) {
    if (!ind_present[p]) {
      ind_label[p] <- 0
      if (ind_label[p - 1] == counter)
        counter <- counter + 1
    }
    else {
      ind_label[p] <- counter
    }
  }
  if (all(!ind_present))
    counter <- 0
  elev_stats <- matrix(NA, nrow = counter, ncol = 4)
  colnames(elev_stats) <- c("id", "length","duration", "burst_tf")
  if (!all(!ind_present)) {
    elev_stats[, 1] <- 1:counter
    for (p in 1:counter) {
      ss_p <- x[ind_label == p]
      rowids <- which(ind_label==p)
      start <- RunTime[rowids[1]]
      if(max(rowids)+1 > length(x)){
        end <- RunTime[max(rowids)]
        duration <- ifelse(as.IDate(Time[max(rowids)]) == as.IDate(Time[max(rowids)-1]), end - start, NA)
      }else{
        end <- RunTime[max(rowids)+1]
        duration <- ifelse(as.IDate(Time[max(rowids+1)]) == as.IDate(Time[max(rowids)]), end - start, NA)
      }

      burst_tf <- "Burst" %in% Type[rowids]
      elev_stats[p, 2] <- length(ss_p)
      elev_stats[p, 3] <- duration
      elev_stats[p, 4] <- burst_tf
    }
  }
  elev_stats <- elev_stats[-nrow(elev_stats), ]
  outlist <- list(ind_label = ind_label, elev_stats = as.data.frame(elev_stats),
                  n_elev = nrow(elev_stats))
  return(outlist)
}

# -------------------------------------------------------------------------------
## ---------------- Helper function for plotting networks -----------------------
# -------------------------------------------------------------------------------


plotnet <- function(phimat, lower = NULL, upper = NULL, layscale =1, digits =2, maximum = 1, laymat = "2x2",  ...){
  # First determine edge labels
  if(is.null(lower) | is.null(upper)){ edge.labels = TRUE } else {
    p <- nrow(phimat)
    edge.labels <- matrix(0, p,p)
    for(r in 1:p){
      for(c in 1:p){
        edge.labels[r,c] <-  paste0(round(phimat[r,c],digits),
                                    " (", round(lower[r,c],digits), ",", round(upper[r,c],digits), ")")
      }}
  }

  if(laymat == "2x2"){
    layout <- matrix(c(-.75,.75,
                       .75,-.75),2,2, byrow = T)*layscale
  }else if (is.matrix(laymat)){
    layout = laymat*layscale
  }else{
    layout = "circle"
  }

  qgraph(phimat, fade = FALSE, edge.labels = edge.labels, maximum = maximum, layout = layout, rescale = FALSE, ...)

}

# -------------------------------------------------------------------------------
## ---------------- Helper function for computing IRFs --------------------------
# -------------------------------------------------------------------------------

# IRF function
IRFfun <- function(drift, start){
  expm::expm(drift)%*%start
}

