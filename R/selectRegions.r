selectRegions <- function(seg, arrays, plot = FALSE, summarise = TRUE, prob.prop = 0.9, prop.greater.than = 0.2, group.name = "", ...)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")
    
    probe.class <- seg$probes.class
    
    props <- proportions(seg, arrays)
    prop.class <- data.frame(unique(props$proportions))
    probabilities <- data.frame(clone.no = 1:nrow(props$proportions), gains = NA, losses = NA)

    dat <- props$proportions
    no.samples <- props$no.samples

    dirichlet.prior <- MLDirichletPriors(seg, arrays, list(c("Normal"), c("Up", "Amplified"), c("Down", "Deleted")))
    cases <- unique(dat)

    for(ii in 1:nrow(cases))
      {
        x <- rdirichlet(10^5, dirichlet.prior+no.samples * cases[ii,])
        datrows <- which(apply(t(dat) == cases[ii,], 2, sum) == ncol(dat))
        probabilities$gains[datrows] <- sum(x[,2] > prop.greater.than) / nrow(x)
        probabilities$losses[datrows] <- sum(x[,3] > prop.greater.than) / nrow(x)
      }

    data <- list()
    data[[1]] <- list()
    data[[1]]$group.name <- group.name
    data[[1]]$plot.ids <- arrays
    data[[1]]$number.samples <- props$no.samples

    selection.table <- list(data = data, gains = probabilities$gains, losses = probabilities$losses)
    selected.regions <- union(which(probabilities$gains > prob.prop), which(probabilities$losses > prob.prop))
    
    if(summarise)
      selection.table$summary <- list(gains = summarySelection(seg, which(probabilities$gains > prob.prop), probabilities$gains, list(arrays), NA, list(c("Up", "Amplified")), "gains", list(props)),
                                      losses = summarySelection(seg, which(probabilities$losses > prob.prop), probabilities$gains, list(arrays), NA, list(c("Down", "Deleted")), "losses", list(props)))
    
    if(plot)
      plotProportions(seg = seg, arrays = arrays, select = selected.regions, ...)
    
    invisible(selection.table)
  }
