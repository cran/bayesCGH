writeSummaries <- function(data, dir = ".", table.name)
  {
    for(nn in 1:length(data$summary))
        {
          class.name <- names(data$summary)[nn]
          filename = paste(dir, "/", table.name, ".", class.name, ".txt", sep = "")
          for(ii in 1:length(data$data))
            {
              write(paste(data$data[[ii]]$group.name, " - Number of samples: ", data$data[[ii]]$number.samples), file = filename, append = (ii > 1))
              write(paste(data$data[[ii]]$group.name, " - Hybs: ", paste(data$data[[ii]]$Hybs, collapse = ", ", sep = "")), file = filename, append = TRUE)
              write(paste(data$data[[ii]]$group.name, " - PlotIDs: ", paste(data$data[[ii]]$plot.ids, collapse = ", ", sep = "")), file = filename, append = TRUE)
              write("\n", file = filename, append = TRUE)
            }
          write("\n\n", file = filename, append = TRUE)
          suppressWarnings(write.table(data$summary[[nn]]$summary, file = filename, sep = "\t", row.names = FALSE, append = TRUE))
          write("\n\n\n\n\nGenes\n\n", file = filename, append = TRUE)
          suppressWarnings(write.table(data$summary[[nn]]$genes, file = filename, sep = "\t", row.names = FALSE, append = TRUE))
        }
  }
