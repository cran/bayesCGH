classifyProbes <- function(seg)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")

    probe.class <- matrix(NA, nrow = nrow(seg$genes), ncol = length(seg$regions))
    for(ii in 1:length(seg$regions))
      for(kk in 1:nrow(seg$regions[[ii]]))
        probe.class[seg$regions[[ii]]$region.start[kk]:seg$regions[[ii]]$region.end[kk], ii] <- seg$regions[[ii]]$class[kk]

    seg$probes.class <- probe.class
    
    seg
  }


MLDirichletPriors <- function(seg, arrays, probe.classes, X = FALSE, Y = FALSE)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")
    chrs <- 1:22
    if(X) chrs <- c(chrs, 23)
    if(Y) chrs <- c(chrs, 24)

    eq.arrays <- arrays
    boots <- NULL
    for(bb in 1:10^3)
      {
        sel.gen <- NULL
        for(cc in chrs)
          sel.gen <- c(sel.gen, sample(which(seg$genes$Chr == cc), 1))
        x <- NULL
        for(ii in 1:length(probe.classes))
          x <- cbind(x, apply(seg$probes.class[sel.gen, eq.arrays], 1, rowCount, probe.classes[[ii]]))
        boots <- rbind(boots, optim(rep(1,ncol(x)), MLalphas, us = x, control = list(fnscale = -1, trace = FALSE))$par)
      }
    dirichlet.prior <- apply(boots, 2, median)
}        


rdirichlet <- function(n, alphas)
  {
    y <- matrix(NA, nrow = n, ncol = length(alphas))
    for(ii in 1:length(alphas))
      y[,ii] <- rgamma(n, alphas[ii], 1)
    x <- y / apply(y, 1, sum)
    return(x)
  }


proportions <- function(seg, sel.arrays, classes = list(normals = c("Normal"), gains = c("Up", "Amplified"), losses = c("Down", "Deleted")))
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")

    probe.matrix <- as.matrix(seg$probes.class[,sel.arrays])
    proportions.probe.matrix <- matrix(0, nrow = nrow(probe.matrix), ncol = length(classes), dimnames = list(c(), names(classes)))
    for(jj in 1:length(classes))
      for(ii in 1:length(classes[[jj]]))
        proportions.probe.matrix[,jj] <- proportions.probe.matrix[,jj] + apply(probe.matrix == classes[[jj]][ii], 1, sum, na.rm = TRUE)
    proportions.probe.matrix <- proportions.probe.matrix / ncol(probe.matrix)
    list(proportions = proportions.probe.matrix, no.samples = ncol(probe.matrix))
  }

find.MCR <- function(summary.set, plotID.col)
  {
    if(nrow(summary.set) == 0) return(character(0)) else {
      MCR <- vector("character", nrow(summary.set))
      for(ii in 1:nrow(summary.set))
        {
          section.arrays <- as.numeric(strsplit(summary.set[ii, plotID.col], ", ")[[1]])
          if(length(section.arrays) == 0)
            {
              MCR[ii] <- "*"
            } else {
              section.arrays.left <- section.arrays.right <- NULL
              if(ii > 1)
                if(summary.set$begin.clone.no[ii] - 1 == summary.set$end.clone.no[ii - 1])
                  section.arrays.left <- as.numeric(strsplit(summary.set[ii - 1, plotID.col], ", ")[[1]])
              if(ii < nrow(summary.set))
                if(summary.set$end.clone.no[ii] + 1 == summary.set$begin.clone.no[ii + 1])
                  section.arrays.right <- as.numeric(strsplit(summary.set[ii + 1, plotID.col], ", ")[[1]])
              
              if(sum(!section.arrays %in% section.arrays.left) > 0 & sum(!section.arrays %in% section.arrays.right) > 0)
                MCR[ii] <- "*"
            }
        }
      return(MCR)
    }
  }

beta <- function(alpha)
  {
    prod(gamma(alpha)) / gamma(sum(alpha))
  }


findChangePoints <- function(sel.arrays, seg)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")

    probe.class <- seg$probes.class
    probe.class[which(probe.class == "Deleted")] <- "Down"
    probe.class[which(probe.class == "Amplified")] <- "Up"
    
    changepoints <- 1
    for(ii in 2:nrow(probe.class))
      if(sum(probe.class[ii, sel.arrays] != probe.class[ii - 1, sel.arrays]) > 0 | seg$genes$Chr[ii] != seg$genes$Chr[ii - 1]) changepoints <- c(changepoints, ii)

    changepoints <- c(changepoints, nrow(probe.class) + 1)
    return(changepoints)
  }

equivalenceClasses <- function(x)
  {
    eq.classes <- list()
    if(length(x) == 1) {
      eq.classes <- list(list(x))
    } else {
      sub.eq.classes <- equivalenceClasses(x[-1])
      for(ii in 1:length(sub.eq.classes))
        {
          eq.classes <- c(eq.classes, list(c(sub.eq.classes[[ii]], x[1])))
          for(jj in 1:length(sub.eq.classes[[ii]]))
            {
              new.class <- sub.eq.classes[[ii]]
              new.class[[jj]] <- c(new.class[[jj]], x[1])
              eq.classes <- c(eq.classes, list(new.class))
            }
        }
    }
    eq.classes
  }

propsList <- function(seg, sel.groups)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")
    
    props <- list()
    for(ii in 1:length(sel.groups))
      {
        props[[ii]] <- proportions(seg, sel.groups[[ii]])
        props[[ii]]$proportions[,2] <- props[[ii]]$proportions[,2]
      }
    props
  }

rowCount <-function(class.row, probe.classes)
  {
    prop <- sum(class.row %in% probe.classes)
    return(prop)
  }

MLalphas <- function(alphas, us)
  {
    if(sum(alphas <= 0) > 0)
      return(-Inf)
    sum(apply(us, 1, logBetaOverBeta, alphas))
  }

logsubgamma <- function(x)
  sum(log(0:(x[1]-1) + x[2]))

logBetaOverBeta <- function(us, alpha)
  {
    if(sum(us) > 0)
      {
        ws <- which(us > 0)
        sum(apply(cbind(us[ws], alpha[ws]), 1, logsubgamma)) - logsubgamma(c(sum(us), sum(alpha)))
      } else return(0)
  }

cytobandIds <- function(locations, cytobands)
  {
    if(nrow(locations) > 0)
      {
        bands <- vector("character", nrow(locations))
        for(ii in 1:nrow(locations))
          {
            chr.bands <- which(cytobands$chr == locations$chr[ii])
            starts <- union(chr.bands[cytobands$start[chr.bands] <= locations$position.start[ii] * 10^6], chr.bands[which.min(cytobands$end[chr.bands])])
            starts <- starts[which.max(cytobands$start[starts])]
            ends <- union(chr.bands[cytobands$end[chr.bands] >= locations$position.end[ii] * 10^6], chr.bands[which.max(cytobands$end[chr.bands])])
            ends <- ends[which.min(cytobands$end[ends])]
            if(starts == ends) bands[ii] <- paste(locations$chr[ii], cytobands$cytoband[starts], sep = "") else
            bands[ii] <- paste(locations$chr[ii], cytobands$cytoband[starts], " - ", locations$chr[ii], cytobands$cytoband[ends], sep = "")
          }
      } else bands <- NULL
    bands
  }


summariseRegions <- function(regions)
  {
    summary <- data.frame(statistic = NA, begin.clone.no = NA, end.clone.no = NA, chr = NA, position.start = NA, position.end = NA, length = NA, genes = NA)
    sel.reg <- regions$genes
    
    group.summary <- list()
    for(jj in 1:length(regions$group.genes))
      {
        group.summary[[jj]] <- as.data.frame(matrix(NA, ncol = ncol(regions$group.genes[[jj]])))
        colnames(group.summary[[jj]]) <- colnames(regions$group.gene[[jj]])
      }
    
    if(nrow(sel.reg) > 0)
      for(ii in 1:nrow(sel.reg))
        {
          contig.reg <- which(summary$end.clone.no + 1 == sel.reg$Clone.No[ii] &
                              summary$chr == sel.reg$Chr[ii] & summary$statistic == sel.reg$statistic[ii])[1]
          if(!is.na(contig.reg))
            for(jj in 1:length(regions$group.genes))
              if(sum(group.summary[[jj]][contig.reg, grep("Hybs", colnames(regions$group.gene[[jj]]))] != regions$group.genes[[jj]][ii, grep("Hybs", colnames(regions$group.gene[[jj]]))]) != 0)
                {
                  contig.reg <- NA
                  break()
                }
          
          if(is.na(contig.reg))
            {
              temp.new <- data.frame(sel.reg$statistic[ii], sel.reg$Clone.No[ii], sel.reg$Clone.No[ii],
                                     sel.reg$Chr[ii], NA, NA, NA,
                                     sel.reg$GeneName[ii])                  
              colnames(temp.new) <- colnames(summary)
              
              new.positions <- strsplit(sel.reg$SystematicName[ii], ":")[[1]][2]
              temp.new$position.start <- as.numeric(strsplit(new.positions, "-")[[1]][1])
              temp.new$position.end <- as.numeric(strsplit(new.positions, "-")[[1]][2])
              temp.new$length <- temp.new$position.end - temp.new$position.start
              
              if(length(grep("chr", as.character(temp.new$genes))) > 0)
                temp.new$genes = ""
              summary <- rbind(summary, temp.new)
              
              for(jj in 1:length(regions$group.genes))
                {
                  temp.group <- regions$group.genes[[jj]][ii,]
                  group.summary[[jj]] <- rbind(group.summary[[jj]], temp.group)
                }
            } else {
              summary$end.clone.no[contig.reg] <- sel.reg$Clone.No[ii]
              summary$position.end[contig.reg] <- as.numeric(strsplit(sel.reg$SystematicName[ii], "-")[[1]][2])
              summary$length[contig.reg] <- summary$position.end[contig.reg] - summary$position.start[contig.reg]
              
              for(jj in 1:length(regions$group.genes))
                {
                  group.summary[[jj]]$min[contig.reg] <- min(group.summary[[jj]]$min[contig.reg], regions$group.genes[[jj]]$min[ii], na.rm = TRUE)
                  group.summary[[jj]]$max[contig.reg] <- max(group.summary[[jj]]$max[contig.reg], regions$group.genes[[jj]]$max[ii], na.rm = TRUE)
                  group.summary[[jj]]$seg.min[contig.reg] <- min(group.summary[[jj]]$seg.min[contig.reg], regions$group.genes[[jj]]$seg.min[ii], na.rm = TRUE)
                  group.summary[[jj]]$seg.max[contig.reg] <- max(group.summary[[jj]]$seg.max[contig.reg], regions$group.genes[[jj]]$seg.max[ii], na.rm = TRUE)
                }
              if(length(grep(sel.reg$GeneName[ii], summary$genes[contig.reg])) == 0 & length(grep("chr", sel.reg$GeneName[ii])) == 0)
                if(summary$genes[contig.reg] == "") summary$genes[contig.reg] <- sel.reg$GeneName[ii] else
              summary$genes[contig.reg] <- paste(summary$genes[contig.reg], sel.reg$GeneName[ii], sep = ", ")              
            }
        }

    summary <- summary[-1,]
    if(nrow(summary) > 0)
      rownames(summary) <- 1:nrow(summary)

    for(jj in 1:length(group.summary))
      {
        group.summary[[jj]] <- group.summary[[jj]][-1,]
        if(nrow(group.summary[[jj]]) > 0)
          rownames(group.summary[[jj]]) <- 1:nrow(group.summary[[jj]])
      }
    
    summary$position.start <- round(summary$position.start / 10^6, 2)
    summary$position.end <- round(summary$position.end / 10^6, 2)
    summary$length <- round(summary$length / 10^6, 2)

    summary <- list(regions = summary, group.summary = group.summary)

    invisible(summary)
  }
