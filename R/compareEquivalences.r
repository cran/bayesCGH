compareEquivalences <- function(seg, group.arrays, group.names, probe.classes = list(c("Up", "Amplified"), c("Down", "Deleted"), c("Normal")), classes = NA, class.names = NULL, class.priors = NA, X.distinct = TRUE, Y.distinct = TRUE, same.priors = TRUE)
  {
    if(class(seg) != "ClassifySegList")
      stop("'seg' argument must be of class 'ClassifySegList'")
    
    probe.class <- seg$probes.class

    num.groups <- length(group.arrays)
    
    if(class(classes) != "list")
      {
        eq.classes <- equivalenceClasses(1:num.groups)
      } else eq.classes <- classes
    
    props <- propsList(seg, group.arrays)
    
    data <- list()
    for(dd in 1:num.groups)
      {
        data[[dd]] <- list()
        data[[dd]]$group.name <- group.names[dd]
        data[[dd]]$plot.ids <- group.arrays[[dd]]
        data[[dd]]$number.samples <- props[[dd]]$no.samples
      }

    compare <- list(data = data, equivalence.classes = eq.classes, probabilities = matrix(NA, nrow = length(eq.classes), ncol = nrow(seg$probes.class)))
    
    all.arrays <- NULL
    for(ii in 1:length(group.arrays))
      all.arrays <- union(all.arrays, group.arrays[[ii]])
    
    changepoints <- findChangePoints(all.arrays, seg)


    dirichlet.prior <- list()
    for(ee in 1:length(eq.classes))
      {
        dirichlet.prior[[ee]] <- list()
        for(qq in 1:length(eq.classes[[ee]]))
            dirichlet.prior[[ee]][[qq]] <- list()
      }

    if(!same.priors)
      {
        for(ee in 1:length(eq.classes))
          for(qq in 1:length(eq.classes[[ee]]))
              {
                found.existing <- FALSE
                if(ee > 1)
                  for(eeee in 1:(ee - 1))
                    for(qqqq in 1:length(eq.classes[[eeee]]))
                      if(sum(suppressWarnings(eq.classes[[ee]][[qq]] != eq.classes[[eeee]][[qqqq]])) == 0)
                        {
                          dirichlet.prior[[ee]][[qq]] <- dirichlet.prior[[eeee]][[qqqq]]
                          found.existing <- TRUE
                        }
                if(!found.existing)
                  dirichlet.prior[[ee]][[qq]] <- MLDirichletPriors(seg, unlist(group.arrays[eq.classes[[ee]][[qq]]]), probe.classes, !X.distinct, !Y.distinct)
              }
      }  else {
        standard.dirichlet.prior <- MLDirichletPriors(seg, unlist(group.arrays), probe.classes, !X.distinct, !Y.distinct)
        for(ee in 1:length(eq.classes))
          for(qq in 1:length(eq.classes[[ee]]))
              dirichlet.prior[[ee]][[qq]] <- standard.dirichlet.prior
            }

    print(dirichlet.prior)
    if(is.na(class.priors))
       {
         prior <- rep(1, length(eq.classes))
         prior <- prior / sum(prior)
       } else prior <- class.priors

    
    for(kk in 1:(length(changepoints) - 1))
      {
        sets.obs <- matrix(NA, ncol = length(probe.classes), nrow = num.groups)
        for(gg in 1:num.groups)
          for(pp in 1:length(probe.classes))
            sets.obs[gg,pp] <- sum(seg$probes.class[changepoints[kk], group.arrays[[gg]]] %in% probe.classes[[pp]])

        Pr.eq.class <- c()
           
        for(ee in which(prior != 0))
          {
            eq.class.prob <- 1
            for(qq in 1:length(eq.classes[[ee]]))
              {
                D <- matrix(sets.obs[eq.classes[[ee]][[qq]],], ncol = length(probe.classes))
                eq.class.prob <- eq.class.prob * (prod(factorial((apply(D, 1, sum)))) / prod(factorial(D))) * beta(dirichlet.prior[[ee]][[qq]] + apply(D, 2, sum)) / beta(dirichlet.prior[[ee]][[qq]])
              }
            Pr.eq.class[ee] <- eq.class.prob * prior[ee]
          }
        Pr.eq.class <- Pr.eq.class / (sum(Pr.eq.class, na.rm = TRUE))
        for(ee in 1:length(eq.classes))
            compare$probabilities[ee, (changepoints[kk]):(changepoints[kk+1] - 1)] <-  Pr.eq.class[ee]        
      }

    rownames(compare$probabilities) <- class.names
    return(compare)
  }
