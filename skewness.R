library(moments)
library(gridExtra)

remove_outliers = function(target, level = 0.999) {
  target <- target[which(target < quantile(target, level))]
  return(target)
}

  remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
  low_threshold <- quantile(target, low_percentile / 100)
  high_threshold <- quantile(target, high_percentile / 100)
  target <- target[target >= low_threshold & target <= high_threshold]
  return(target)
}

 # Define skew_gate function
skew_gate <- function(x, alpha) {
  sk <- moments::skewness(x)
  
  n <- length(x)
  a <- min(x) + (min(x) + median(x)) / 10
  b <- max(x)
  
  if (sk < 0) {
    message("The skewness is negative!")
    
    return(list(
    skewness = sk,
    cutoff = b,
    N = n,
    N_removed = 0,
    percentage_removed = 0,
    returnvalplot = x[which(x> a+(b-a)/2)]

  ))
  }

  iteration <- 0
  
  while (abs(sk) > alpha & iteration <= 100) {
    # print(iteration)
    if (sk >= 0) {
      b <- a + (b - a) / 2
    } else {
      a <- a + (b - a) / 2
    }
    a <- min(a, b)
    b <- max(a, b)
    
    sk <- skewness(x[which(x < b)])
    if (is.nan(sk)) {
      message("Warning: skewness is NaN")
      break
    }
    iteration <- iteration + 1
  }
  
  n_removed <- sum(x > b)
  perc_removed <- round(n_removed / n, 3)
  list(
    skewness = sk,
    cutoff = b,
    N = n,
    N_removed = n_removed,
    percentage_removed = perc_removed,
    returnvalplot = x[which(x> a+(b-a)/2)]

  )
  }

# plot_gating_individual = function(target, output, marker, do_ggplot = TRUE) {
#   if (do_ggplot == TRUE) {
#     plot <- ggplot(data.frame(x = log(target)), aes(x = x)) +
#       geom_histogram(bins = 40,
#                      fill = "lightblue",
#                      color = "black") +
#       geom_vline(xintercept = min(log(output$cutoff)), color = "red") +
#       labs(
#         x = "",
#         y = "",
#         title = paste(
#           "N+",
#           output$N_removed,
#           "+R",
#           output$percentage_removed,
#           sep = "="
#         )
#       ) +
#       # print(marker)
#       xlab(marker)
#     return(plot)
#   } else{
#     hist(
#       log(target),
#       breaks = 200,
#       main = "",
#       xlab = "",
#       ylab = "",
#       col = "lightblue"
#     )
#     abline(v = min(log(output$cutoff)), col = "red")
#     title(paste(
#       paste("N+", output$N_removed, sep = "="),
#       paste("+R", output$percentage_removed, sep = "="),
#       sep = "; "
#     ),
#     xlab = marker,
#     ylab = "") +
#       theme(plot.title = element_text(hjust = 0.5))
#   }
# }

  # summary_text <- paste(
  #   "Crevasse Cutoff:", round(cutoff_value,2),
  #   "\nSkewness:", round(skewness_value,3),
  #   "\nN:", round(N_val,2),
  #   "\nN Removed (N+):", length(truncated),
  #   "\nPercentage Removed (+R):", round(length(truncated)/length(targetforplot), 3),
  #   "\n\nGMM Cutoff:", round(gmm_gate_value,2)

  # )


#  truncated <- result_to_plot$returnvalplot
#   # targetforplot <- gate_vliner[which(gate_vliner<quantile(gate_vliner, 0.999))]
#   targetforplot <- remove_outliers2(gate_vliner)

plot_gating_individual = function(target, output, marker, do_ggplot = TRUE) {
  truncated <- output$returnvalplot
  targetforplot <- remove_outliers2(target)
  
  if (do_ggplot == TRUE) {
    plot <- ggplot(data.frame(x = target), aes(x = x)) +
      geom_histogram(bins = 40,
                     fill = "lightblue",
                     color = "black") +
      geom_vline(xintercept = output$cutoff, color = "red") +
      labs(
        x = "",
        y = "",
        title = paste(
          "N+",
          length(truncated),
          "+R",
          round(length(truncated)/length(targetforplot), 3),
          sep = "="
        )
      ) +
      # print(marker)
      xlab(marker)
    return(plot)
  } else{
    hist(
      target,
      breaks = 200,
      main = "",
      xlab = "",
      ylab = "",
      col = "lightblue"
    )
    abline(v = output$cutoff, col = "red")
    title(paste(
      paste("N+", length(truncated), sep = "="),
      paste("+R", round(length(truncated)/length(targetforplot), 3), sep = "="),
      sep = "; "
    ),
    xlab = marker,
    ylab = "") +
      theme(plot.title = element_text(hjust = 0.5))
  }
}
# plot_gating_grid = function(df, unique_id, column_indices) {
#   # names_columns = names(df[, column_index])
#   plot_list = list()

#   names_columns = names(df)
  
#   for (j in column_indices) {
#     target <- df[df$imageid == unique_id, j]
#     target <- remove_outliers2(target)
#     target <- target[!is.infinite(target)]
#     # print(target)
#     marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
#     # print(marker)
#     output <- skew_gate(target, 0.5)
#     # output["marker"] = marker
#     p_temp = plot_gating_individual(target, output, marker)
#     plot_list = c(plot_list, list(p_temp))
#   }
#   p = gridExtra::grid.arrange(grobs = plot_list, ncol = 5, do_ggplot=FALSE)
#   return(p)
# }
# par(mfrow=c(length(unique_ids),length()), mar = c(5.1, 3, 4.1, 2))
    # par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))  # Adjust these parameters as needed

plot_gating_grid = function(df, unique_ids, column_indices) {
  plot_list = list()

  names_columns = names(df)
  # Calculate the number of rows and columns for mfrow
  n_plots_per_row = length(column_indices)  # Adjust this as needed
  n_plots = length(unique_ids) * length(column_indices)
  n_rows = ceiling(n_plots / n_plots_per_row)

  print(n_rows)
  print(n_plots)
  print(n_plots_per_row)

  for (unique_id in unique_ids) {

  # for (a in c(1:2)) {
    # par(mfrow = c(n_rows, n_plots_per_row), mar = c(5.1, 4.1, 4.1, 2.1))
    # par(mfrow = c(1, 1))

    for (j in column_indices) {
    target <- df[df$imageid == unique_id, j]


    # Calculate the maximum value for each column (excluding imageid)
    max_values <- sapply(target, max)
    max_values <- max_values[-1]  # Remove the imageid column


    target <- remove_outliers2(target)


    if (any(max_values > 20)) {
        cat("Working with raw data\n")
        # Apply log transformation and other processing steps here
        # sub_data <- log(sub_data + 1e-10)  # Adding 1 to avoid log(0)
        target <- log(target)  # Adding 1 to avoid log(0)
        # sub_data[is.infinite(sub_data)] <- 1e-10  # Replace -Inf values
        } else {
          cat("Working with logged values\n")
        }

    target <- target[!is.infinite(target)]

    # print(target)

    marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
    output <- skew_gate(target, 0.5)
    p_temp = plot_gating_individual(target, output, marker)
    plot_list = c(plot_list, list(p_temp))
  }
  }

  # Limit the number of plots to a maximum of 180
  plot_list <- plot_list[1:min(length(plot_list), 180)]
  
  # p = gridExtra::grid.arrange(grobs = plot_list, ncol = n_plots_per_row, nrow= n_rows)
  p = gridExtra::grid.arrange(grobs = plot_list)

  return(p)
}

plot_gating_grid2 = function(df, unique_ids, column_indices_list) {

  plot_list = list()


  names_columns = names(df)

  for (imageid in unique_ids){
  # par(mfrow=c(5,4), mar = c(5.1, 3, 4.1, 2))
  # par(mfrow=c(5,4), mar = c(5.1, 3, 4.1, 2))
  par(mfrow = c(5, 4), mar = c(5.1, 3, 6, 2.1))


  for (j in column_indices_list){
    target <- df[df$imageid == unique_id, j]
    target <- remove_outliers2(target)
    target <- target[!is.infinite(target)]
    marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
    print(marker)
    truncated<-Sgate(target, 0.8)
    N[j,a]<-length(target)
    n[j,a]<-length(truncated)
    hist(log(target), breaks=200, main = "", xlab = "",  ylab = "", col = "lightblue")
    abline(v=min(log(truncated)), col="red")
    title(paste(paste("N+",length(truncated), sep="="), paste("+R", round(length(truncated)/length(target), 3), sep = "="), sep = "; "), xlab = paste(paste(strsplit(files[j], "-")[[1]][1], sep=""), marker, sep="; "), ylab = "")
  }
}
}

# files<- list.files()
# pdf("Gating_all.pdf")
# for (a in c(1:9)){
#   par(mfrow=c(5,4), mar = c(5.1, 3, 4.1, 2))
#   for (j in c(1:4, 6:8, 10:12, 14:21, 23:24)){
#     target<- as.numeric(unlist(GlSnorm[[j]][a]))
#     target<- target[which(target<quantile(target, 0.999))]
#     #target<- target[which(target>quantile(target, 0.0001))]
#     marker<- names(GlSnorm[[j]][a])
#     truncated<-Sgate(target, 0.8)
#     N[j,a]<-length(target)
#     n[j,a]<-length(truncated)
#     hist(log(target), breaks=200, main = "", xlab = "",  ylab = "", col = "lightblue")
#     abline(v=min(log(truncated)), col="red")
#     title(paste(paste("N+",length(truncated), sep="="), paste("+R", round(length(truncated)/length(target), 3), sep = "="), sep = "; "), xlab = paste(paste(strsplit(files[j], "-")[[1]][1], sep=""), marker, sep="; "), ylab = "")
#   }
# }

# plot_gating_grid = function(df, unique_ids, column_indices) {
#   # names_columns = names(df[, column_index])
#   plot_list = list()


#   names_columns = names(df)
#   print(names_columns)

#   # print()

#   # print(column_index)
#   print(column_indices)
#   print("unq ids inside plot gate grid")
#   print(unique_ids)
  
#   for (unique_id in unique_ids) {
#     print("inside first ")
#     print(unique_id)
#     for (j in column_indices) {
#     print(names_columns[j])
#     target <- df[df$imageid == unique_id, j]
#     print(colnames(target))
#     target <- remove_outliers2(target)
#     target <- target[!is.infinite(target)]
#     # print(target)
#     marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
#     # print(marker)
#     output <- skew_gate(target, 0.5)
#     # output["marker"] = marker
#     p_temp = plot_gating_individual(target, output, marker)
#     plot_list = c(plot_list, list(p_temp))
#   }
#   p = gridExtra::grid.arrange(grobs = plot_list, ncol = 5)
#   return(p)
# }
#   }

# plot_gating_grid <- function(df, unique_ids, column_indices) {
#   plot_lists <- list()  # Create a list to store plots for each unique ID

#   names_columns <- names(df)

#   for (unique_id in unique_ids) {
#     unique_plots <- list()  # Create a list to store plots for the current unique ID
#     for (j in column_indices) {
#       print(unique_id)
#       print(names_columns[j])
#       target <- df[df$imageid == unique_id, j]
#       target <- remove_outliers2(target)
#       target <- target[!is.infinite(target)]
#       marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
#       output <- skew_gate(target, 0.5)
#       p_temp <- plot_gating_individual(target, output, marker)
#       unique_plots <- c(unique_plots, list(p_temp))
#     }
#     plot_lists <- c(plot_lists, list(unique_plots))
#   }

#   # Now 'plot_lists' is a list of lists, where each sublist contains plots for one unique ID
#   # You can arrange these plots as needed

#   # For example, if you want to arrange plots for each unique ID in separate rows:
#   p_list <- lapply(plot_lists, function(unique_plots) {
#     gridExtra::grid.arrange(grobs = unique_plots, ncol = 5)
#   })

#   return(p_list)
# }

  

# generatePlot <- function() {
#     req(uploaded_df())
#     df <- uploaded_df()

#     selected_columns <- input$selected_columns

#     # Get the column indices for the selected columns
#     column_indices <- sapply(selected_columns, function(col_name) {
#       which(names(df) == col_name)
#     })

#     unique_ids <- unique(df$imageid)  # You can select any specific unique_id here

#     # Get the names of the selected columns based on their indices
#     names_columns <- names(df)[column_indices]

#     # Plot
#     plot_gating_grid(df, unique_ids, column_indices, names_columns)
# }

# plot_gating_individual = function(target, output, marker, do_ggplot = TRUE) {
#   if (do_ggplot == TRUE) {
#     plot <- ggplot(data.frame(x = target), aes(x = x)) +
#       geom_histogram(bins = 40,
#                      fill = "lightblue",
#                      color = "black") +
#       geom_vline(xintercept = output$cutoff, color = "red") +
#       labs(
#         x = marker,  # Use marker for x-axis label
#         y = "",
#         title = paste(
#           "N+",
#           output$N_removed,
#           "+R",
#           output$percentage_removed,
#           sep = "="
#         )
#       )
#     return(plot)
#   } else{
#     hist(
#       log(target),
#       breaks = 200,
#       main = "",
#       xlab = "",
#       ylab = "",
#       col = "lightblue"
#     )
#     abline(v = min(log(output$cutoff)), col = "red")
#     title(paste(
#       paste("N+", output$N_removed, sep = "="),
#       paste("+R", output$percentage_removed, sep = "="),
#       sep = "; "
#     ),
#     xlab = marker,  # Use marker for x-axis label
#     ylab = "") +
#       theme(plot.title = element_text(hjust = 0.5))
#   }
# }

# plot_gating_grid = function(df, unique_ids, column_indices, names_columns) {
#   plot_list = list()

#   # print(unique_ids)

#   # print(unique_id)
  
#   for (unique_id in unique_ids) {
#     for (j in column_indices) {
#       # print(unique_id)
#       # print(j)
#       target <- df[df$imageid == unique_id, j]
#       # print(target)
#       target <- remove_outliers2(target)
#       marker <- unique_id
#       target <- target[!is.infinite(target)]
#       output <- skew_gate(target, 0.5)
#       p_temp = plot_gating_individual(target, output, marker)
#       plot_list = c(plot_list, list(p_temp))
#     }
#   }
#   p = gridExtra::grid.arrange(grobs = plot_list, ncol = 5)
#   return(p)
# }

