summarySelection <- function(seg, regions.clone.no, stats, group.arrays, group.names,
                             probe.classes = list(c("Up", "Amplified"), c("Down", "Deleted"), c("Normal")),
                             probe.class.name = c("Normal", "Gain", "Loss"), local.maximum = FALSE, cytobands = cytobands)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")

    props <- propsList(seg, group.arrays)
    
    probe.class <- seg$probes.class
    
    if(length(regions.clone.no) > 0)
      {
        regions <- list(group.genes = list())
        for(ii in 1:length(group.arrays))
          {
            regions$group.genes[[ii]] <- data.frame(min = rep(NA, length(regions.clone.no)), max = NA, seg.min = NA, seg.max = NA)
            for(kk in length(probe.classes):1)
              {
                temp.df <- data.frame(Hybs = I(""), ID = I(""), plot.ID = I(""), proportion = rep(0, length(regions.clone.no)))
                colnames(temp.df) <- paste(probe.class.name[kk], ":", colnames(temp.df), sep = "")
                regions$group.genes[[ii]] <- cbind(temp.df, regions$group.genes[[ii]])
              }
            for(jj in 1:length(regions.clone.no))
              {
                for(kk in 1:length(probe.classes))
                  {
                    sel.samples <- intersect(which(probe.class[regions.clone.no[jj],] %in% probe.classes[[kk]]), group.arrays[[ii]])
                                        #if(length(sel.samples) != 0)
                   
#                    regions$group.genes[[ii]][jj, 1 + (kk - 1) * 4]<- paste(acgh.table$Hyb[sel.samples], sep = "", collapse = ", ")
#                    regions$group.genes[[ii]][jj, 2 + (kk - 1) * 4] <- paste(acgh.table$ID[sel.samples], sep = "", collapse = ", ")
                    regions$group.genes[[ii]][jj, 3 + (kk - 1) * 4] <- paste(sel.samples, sep = "", collapse = ", ")
                    regions$group.genes[[ii]][jj, 4 + (kk - 1) * 4] = signif(length(sel.samples) / props[[ii]]$no.samples, 2)
                  }
                regions$group.genes[[ii]]$min[jj] <- signif(min(seg$M.observed[regions.clone.no[jj], group.arrays[[ii]]], na.rm = TRUE), 2)
                regions$group.genes[[ii]]$max[jj] <- signif(max(seg$M.observed[regions.clone.no[jj], group.arrays[[ii]]], na.rm = TRUE), 2)
                regions$group.genes[[ii]]$seg.min[jj] <- signif(min(seg$M.predicted[regions.clone.no[jj], group.arrays[[ii]]], na.rm = TRUE), 2)
                regions$group.genes[[ii]]$seg.max[jj] <- signif(max(seg$M.predicted[regions.clone.no[jj], group.arrays[[ii]]], na.rm = TRUE), 2)
              }
          }
      
        regions$genes <- cbind(Clone.No = regions.clone.no, seg$genes[regions.clone.no,], statistic = signif(stats[regions.clone.no], 2))

        summary.regions <- summariseRegions(regions)        
        summary.regions$regions <- cbind(cytoband = cytobandIds(summary.regions$regions, cytobands), summary.regions$regions)
        
        if(length(probe.classes) == 1)
          {
            MCR.groups <- matrix(NA, nrow = nrow(summary.regions$regions), ncol = length(summary.regions$group.summary))
            for(ii in 1:length(summary.regions$group.summary))
              {
                summary.set <- cbind(summary.regions$regions, summary.regions$group.summary[[ii]])
                MCR.groups[,ii] <- find.MCR(summary.set, grep("plot.ID", colnames(summary.set)))
              }
            
            MCR <- rep("", nrow(MCR.groups))
            MCR[which(apply(MCR.groups != "*", 1, sum) == 0)] <- "*"
            
            summary.regions$regions <- cbind(MCR = MCR, summary.regions$regions)
          }
        summary.regions.table <- summary.regions$regions
        for(ii in 1:length(summary.regions$group.summary))
          {
            if(!is.na(group.names[1]))
              colnames(summary.regions$group.summary[[ii]]) <- paste(colnames(summary.regions$group.summary[[ii]]), " (", group.names[ii], ")", sep = "")
            summary.regions.table <- cbind(summary.regions.table, summary.regions$group.summary[[ii]])
          }

        if(local.maximum)
          {
            contig.group <- c(1)
            local.max <- c()
            for(ss in 1:(nrow(summary.regions.table)-1))
              {
                if(summary.regions.table$chr[ss+1] == summary.regions.table$chr[ss] & summary.regions.table$begin.clone.no[ss+1] == summary.regions.table$end.clone.no[ss] + 1)
                  {
                    contig.group <- c(contig.group, ss + 1)
                  } else {
                    print(contig.group)
                    local.max <- c(local.max, contig.group[which.max(summary.regions.table$statistic[contig.group])])
                    contig.group <- c(ss+1)
                  }
              }
            local.max <- c(local.max, contig.group[which.max(summary.regions.table$statistic[contig.group])])
     
            summary.regions.table = summary.regions.table[local.max,]
            regions.clone.no <- c()
            for(ss in 1:nrow(summary.regions.table))
              regions.clone.no <- c(regions.clone.no, summary.regions.table$begin.clone.no[ss]:summary.regions.table$end.clone.no[ss])
            }

        summary <- list(summary = summary.regions.table, genes = seg$genes[regions.clone.no,])
      } else summary <- list(summary = NULL, genes = NULL)
    
    return(summary)
  }
