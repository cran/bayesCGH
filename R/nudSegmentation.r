nudSegmentation <- function(seg, factor.change = 0.75, amplified.max.width = 10, deleted.max.width = 10, amplified.magnitude = 1, deleted.magnitude = 1, cellularity)
  {
    if(missing(cellularity))
      cellularity <- rep(NA, ncol(seg$M.predicted))
    
    for(kk in 1:ncol(seg$M.predicted))
      {
        breaks <- findBreakPoints(seg,kk)
        
        breaks.start <- breaks[(1:(length(breaks) / 2) * 2 - 1)]
        breaks.end <- breaks[(1:(length(breaks) / 2) * 2)]
        
        predicted <- vector("numeric", length(breaks.start))
        chr <- vector("numeric", length(breaks.start))
        position.start <- vector("numeric", length(breaks.start))
        position.end <- vector("numeric", length(breaks.start))
        
        
        for(ii in 1:length(breaks.start))
          {
            predicted[ii] <- seg$M.predicted[breaks.start[ii],kk]
            chr[ii] <- seg$genes$Chr[breaks.start[ii]]
            position.start[ii] <- seg$genes$Position[breaks.start[ii]]
            position.end[ii] <- seg$genes$Position[breaks.end[ii]]
          }
        
        dev.from.predicted <- sum((seg$M.observed[,kk] - seg$M.predicted[,kk])^2, na.rm = TRUE) / (nrow(seg$M.observed) - 1)
        dev.from.predicted <- seg$M.observed[,kk] - seg$M.predicted[,kk]
        dev.from.predicted <- (summary(dev.from.predicted)[5] - summary(dev.from.predicted)[2])
        if(!is.na(cellularity[kk]))
          dev.from.predicted <- log2(cellularity[kk] * 2^dev.from.predicted + 1 - cellularity[kk])
        
        seg$regions[[kk]] <- data.frame(region.start = I(breaks.start), region.end = I(breaks.end), chr = I(chr), position.start = position.start, position.end = position.end, predicted = I(predicted), class = NA, color = NA, spot = NA)
        
        seg$regions[[kk]]$class[which(abs(seg$regions[[kk]]$predicted) <= factor.change * dev.from.predicted)] <- "Normal"
        seg$regions[[kk]]$class[which(seg$regions[[kk]]$predicted < -factor.change * dev.from.predicted)] <- "Down"
        seg$regions[[kk]]$class[which(seg$regions[[kk]]$predicted > factor.change * dev.from.predicted)] <- "Up"
        seg$regions[[kk]]$class[which(seg$regions[[kk]]$predicted > amplified.magnitude & 
                                      seg$regions[[kk]]$position.end - seg$regions[[kk]]$position.start < amplified.max.width)] <- "Amplified"
        seg$regions[[kk]]$class[which(seg$regions[[kk]]$predicted < -deleted.magnitude & 
                                      seg$regions[[kk]]$position.end - seg$regions[[kk]]$position.start < deleted.max.width)] <- "Deleted"
        
        seg$regions[[kk]]$color[which(seg$regions[[kk]]$class == "Normal")] <- "orange"
        seg$regions[[kk]]$color[which(seg$regions[[kk]]$class == "Down")] <- "green"
        seg$regions[[kk]]$color[which(seg$regions[[kk]]$class == "Up")] <- "red"
        seg$regions[[kk]]$color[which(seg$regions[[kk]]$class == "Amplified")] <- "darkred"
        seg$regions[[kk]]$color[which(seg$regions[[kk]]$class == "Deleted")] <- "darkgreen"
        seg$regions[[kk]]$spot[which(seg$regions[[kk]]$class == "Amplified")] <- 1
        seg$regions[[kk]]$spot[which(seg$regions[[kk]]$class == "Deleted")] <- 1
      }
    seg <- new("ClassifySegList", seg)
    seg <- classifyProbes(seg)
  }
