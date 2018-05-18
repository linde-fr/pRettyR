plot_acf_multi_withNA_gg <- function(maxlag, fitlist, datalist, Dwhole, 
                                  normalized=T, n.ahead = 1, ylimset = c(-0.06, 0.25)){
  # data in datalist should all be of form D, so has same time stamps
  # acf will be outputted with NAs placed to make it like D.state
  
  ## plots all residuals below one another, with the same y limits "reslim"
  names(fitlist) <- sapply(fitlist, function(x) x$modelname)
  n <- length(fitlist)
  
  ### set datalist0 for prediction ===============================================================
  pred <- lapply(1:n, function(x) predict(fitlist[[x]], newdata = datalist[[x]], n.ahead = n.ahead ) )
  names(pred) <- names(fitlist)
  
  ## Calculate the residuals =============================================
  yToHat <- lapply(pred, function(x) x$output$pred$yTo )
  names(yToHat) <- names(fitlist)
  residuals <- lapply(1:n, function(x) datalist[[x]]$yTo - yToHat[[x]] )
  names(residuals) <- names(fitlist)
  if(normalized){
    sds <- lapply(pred, function(x) x$output$sd$yTo )
    residuals <- lapply(1:n, function(x) residuals[[x]]/sds[[x]] )
  }

  ##----------------------------------------------------------------
  ## Plot the auto-correlation function and cumulated periodogram in a new window
  x11()
  #par(mfrow=c(n,1))
  ## The blue lines indicates the 95% confidence interval, meaning that if it is
  ##  white noise, then approximately 19 out of 20 lag correlations will be inside.
  acfs <- vector('list', length = n)
  for (j in 1:n){
    ## put residuals with NAs
    df <- data.frame(t = datalist[[1]]$t, residuals = residuals[[j]])
    DF <- merge(Dwhole, df, by = "t", all = T)
    DF <- DF[, c("t", "residuals")]
    
    # plot them
    z <- acf(DF$residuals, lag.max = maxlag, na.action = na.exclude, plot = F)
    acfs[[j]] <- z$acf[-1]
  }
  # put them in data.frame
  acf.df <- data.frame(acfs[[1]])
  for(j in 2:n){
    acf.df <- cbind(acf.df, acfs[[j]])
  }
  # add lag numbers
  acf.df <- data.frame(lag = 1:maxlag, acf.df)
  colnames(acf.df)[2:(n+1)] <- as.character(1:n)
  # compute significance levels
  siglev <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(DF$residuals)))

  ## plotmapping = aes(xend = lag, yend = 0))
  DFplot <- melt(acf.df, id.var= 'lag')
  p <- ggplot(DFplot, aes(x = lag, y = value)) + geom_segment(aes(xend = lag, yend = 0), color = 'black') + 
    geom_hline(aes(yintercept = 0)) + # scale_y_continuous(limits = c(-0.06, 0.25)) +
    facet_grid(variable ~ ., scales = "free_y") + theme(legend.position = "none") +
    ylab(NULL) + geom_hline(aes(yintercept = siglev), color = 'blue', linetype = 'dashed') +
    geom_hline(aes(yintercept = -siglev), color = 'blue', linetype = 'dashed') +
    coord_cartesian(ylim=ylimset) + xlab('lag (minutes)')
  return(p)
  
}