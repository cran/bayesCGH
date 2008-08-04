plotProportions <- function(seg, arrays, chrom.info = chrominfo.Mb, ylimit=c(-1,1),
                             main="Frequency of genetic changes across tumours", num.chr = 24,
                             amps.dels = TRUE, errors = TRUE, select = FALSE, colors = c("green", "red", "dark green", "dark red"))
  {
    overlay = FALSE
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")

    probe.class <- seg$probes.class
    gene.list <- seg$genes
    
    props <- proportions(seg, arrays)
    
    dat <- props$proportions
    no.samples <- props$no.samples
    
    positions <- gene.list$Position
    for(ii in 1:num.chr)
      positions[which(gene.list$Chr == ii)] <- positions[which(gene.list$Chr == ii)] + sum(chrom.info$length[which(chrom.info$chrom < ii)])
    
    genes.to.plot <- 1:nrow(dat)
    if(is.na(select[1]))
      {
        genes.to.plot <- NULL
      } else if(select[1]) {
        genes.to.plot <- select
      }
    
    if(errors)
      {
        dirichlet.prior <- MLDirichletPriors(seg, arrays, list(c("Normal"), c("Up", "Amplified"), c("Down", "Deleted")))
        up.errors <- matrix(ncol = 2, nrow = nrow(dat))
        down.errors <- matrix(ncol = 2, nrow = nrow(dat))
        up.extremes <- matrix(ncol = 2, nrow = nrow(dat))
        down.extremes <- matrix(ncol = 2, nrow = nrow(dat))

        cases <- unique(dat[genes.to.plot,])
        
        for(ii in 1:nrow(cases))
          {
            x <- rdirichlet(10^5, dirichlet.prior+no.samples * cases[ii,])
            datrows <- which(apply(t(dat) == cases[ii,], 2, sum) == ncol(dat))
            up.errors[datrows,] <- matrix(c(quantile(x[,2], 0.25, na.rm = TRUE), quantile(x[,2], 0.75, na.rm = TRUE)), ncol = 2, nrow = length(datrows), byrow = TRUE)
            down.errors[datrows,] <- matrix(c(quantile(x[,3], 0.25, na.rm = TRUE), quantile(x[,3], 0.75, na.rm = TRUE)), ncol = 2, nrow = length(datrows), byrow = TRUE)
            up.extremes[datrows,] <- matrix(c(quantile(x[,2], 0.05, na.rm = TRUE), quantile(x[,2], 0.95, na.rm = TRUE)), ncol = 2, nrow = length(datrows), byrow = TRUE)
            down.extremes[datrows,] <- matrix(c(quantile(x[,3], 0.05, na.rm = TRUE), quantile(x[,3], 0.95, na.rm = TRUE)), ncol = 2, nrow = length(datrows), byrow = TRUE)
          }
          plot(x = NA, y = NA, ylim = ylimit, xlim = c(0,max(positions)),
               xlab="Chromosome",ylab="Proportion of Genomic changes",
               main=paste(main, "  {number of samples: ", no.samples, "}", sep = ""), axes = FALSE, type = "n")
        
        for(ii in 1:nrow(cases))
          {
            datrows <- which(apply(t(dat) == cases[ii,], 2, sum) == ncol(dat))
            splitdatrows <- c(0, which(datrows[-length(datrows)] - datrows[-1] < -1), length(datrows))
            rect(xleft = positions[datrows[splitdatrows[-length(splitdatrows)] + 1]], xright = positions[datrows[splitdatrows]],
                 ybottom = up.extremes[datrows[1],1], ytop = up.extremes[datrows[1],2], col = colors[3], density = NA)
            rect(xleft = positions[datrows[splitdatrows[-length(splitdatrows)] + 1]], xright = positions[datrows[splitdatrows]],
                 ybottom = up.errors[datrows[1],1], ytop = up.errors[datrows[1],2], col = colors[1], density = NA)
            rect(xleft = positions[datrows[splitdatrows[-length(splitdatrows)] + 1]], xright = positions[datrows[splitdatrows]],
                 ybottom = -down.extremes[datrows[1],1], ytop = -down.extremes[datrows[1],2], col = colors[4], density = NA)
            rect(xleft = positions[datrows[splitdatrows[-length(splitdatrows)] + 1]], xright = positions[datrows[splitdatrows]],
                 ybottom = -down.errors[datrows[1],1], ytop = -down.errors[datrows[1],2], col = colors[2], density = NA)
          }
      }  else {
        suppressWarnings(par(new = overlay))
          plot(y = c(dat[genes.to.plot,2], -dat[genes.to.plot,3]), x = rep(positions[genes.to.plot], 2), type = "h",
               col = c(rep(colors[1], length(genes.to.plot)), rep(colors[2], length(genes.to.plot))),
               xlab="Chromosome",ylab="Proportion of Genomic changes", ylim = ylimit, xlim = c(0,max(positions)),
               main=paste(main, "  {number of samples: ", no.samples, "}", sep = ""), axes = FALSE)
        par(new = FALSE)
        
      if(amps.dels)
        {
          prop.amp.dels <- proportions(seg, arrays, list(gains = "Amplified", losses = "Deleted"))
          dat.amp.dels <- prop.amp.dels$proportions
          par(new = TRUE)
          plot(y = c(dat.amp.dels[genes.to.plot,1], -dat.amp.dels[genes.to.plot,2]), x = rep(positions[genes.to.plot], 2), type = "h",
               col = c(rep(colors[3], length(genes.to.plot)), rep(colors[4], length(genes.to.plot))),
               xlab = "",ylab = "", main = "",
               ylim = ylimit, xlim = c(0,max(positions)), axes = FALSE)
          par(new = FALSE)
        }

    }

  axis(side = 2)

  abline(h = 0)
  
  kk <- 0
  abline(v = kk, col = "blue")
  for(ii  in 1:num.chr)
    {
      abline(v = kk + chrom.info$centromere[ii], col = "orange", lty = 2)
      kk <- kk + chrom.info$length[ii]
      abline(v = kk, col = "blue")
      text(kk - (chrom.info$length[ii] / 2), 0.7, as.character(ii), col = "blue")
    }
}
