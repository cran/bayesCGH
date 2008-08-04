selectAmpsDels <- function(seg, arrays, group.name = "")
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")

    probe.class <- seg$probes.class

    props <- proportions(seg, arrays, classes = list("Amplified", "Deleted"))
    prop.class <- data.frame(unique(props$proportions))
    
    data <- list()
    data[[1]] <- list()
    data[[1]]$group.name <- group.name
    data[[1]]$plot.ids <- arrays
    data[[1]]$number.samples <- props$no.samples
    
    selection.table <- list(data = data, amplifications = props$proportions[,1], deletions = props$proportions[,2]) 

    selection.table$summary <- list(amplifications = summarySelection(seg, which(selection.table$amplifications > 0), selection.table$amplifications, list(arrays), NA,
                                      list(c("Amplified"), c("Up", "Amplified")), c("amplifications", "gains"), list(props)),
                                    deletions = summarySelection(seg, which(selection.table$deletions > 0),
                                      selection.table$deletions, list(arrays), NA, list(c("Deleted"), c("Down", "Deleted")), c("deletions", "losses"), list(props)))

    invisible(selection.table)
  }
