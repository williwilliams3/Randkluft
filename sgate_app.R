# Load required library
library(moments)

get_gates_csv <- function(dataframe, csv_save_file_path) {

  # #  Define remove_outliers function
  # remove_outliers <- function(target, level = 0.999) {
  #   target <- target[which(target < quantile(target, level))]
  #   return(target)
  # }

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
      return(b
      )
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
    return(b)
  }

  # Record start time for CSV reading
  csv_read_start_time <- Sys.time()

  data <- dataframe # Replace with your actual CSV file path

  # print(data)

  # Record end time for CSV reading
  csv_read_end_time <- Sys.time()

  # Calculate elapsed time for CSV reading
  csv_read_elapsed_time <- csv_read_end_time - csv_read_start_time

  # Get unique values from the imageid column
  unique_imageids <- unique(data$imageid)

  # Define skew_gate function above

  alpha <- 0.5  # Set your desired alpha value

  # Create an empty data frame to store the results
  results_df <- data.frame(Patient = character(),
                          Marker = character(),
                          Gate = numeric(),
                          stringsAsFactors = FALSE)


  # Record start time
  start_time <- Sys.time()
  # Loop through each unique imageid
  for (imageid in unique_imageids) {
    # sub_data <- data[data$imageid == imageid, 1:46]  # Select columns 1 to 46 2 10 
    # sub_data <- data[data$imageid == imageid, c(4:6, 8:10, 12:14, 16:18, 20:22, 24:26, 28:30, 32:34, 36:38, 40:42)]  # Select columns 1 to 46 2 10  3:42 for mouseagg
    # sub_data <- data[data$imageid == imageid, 3:42]  
    sub_data <- data[data$imageid == imageid, -1]  


    # Calculate the maximum value for each column (excluding imageid)
    max_values <- sapply(sub_data, max)
    max_values <- max_values[-1]  # Remove the imageid column

    
    # if (any(max_values > 20)) {
    #   cat("Working with raw data\n")
    #   # Apply log transformation and other processing steps here
    #   # sub_data <- log(sub_data + 1e-10)  # Adding 1 to avoid log(0)
    #   sub_data <- log(sub_data)  # Adding 1 to avoid log(0)
    #   # sub_data[is.infinite(sub_data)] <- 1e-10  # Replace -Inf values
    # } else {
    #   cat("Working with logged values\n")
    # }


    # # Remove rows with -Inf values
    # sub_data <- sub_data[!apply(sub_data, 1, function(row) any(row == -Inf)), ]
    
    # Loop through each column in sub_data
    for (col_idx in 1:ncol(sub_data)) {
      col_name <- colnames(sub_data)[col_idx]
      values <- sub_data[, col_idx]

      values <- remove_outliers2(values)

      if (any(max_values > 20)) {
        cat("Working with raw data\n")
        # Apply log transformation and other processing steps here
        # sub_data <- log(sub_data + 1e-10)  # Adding 1 to avoid log(0)
        values <- log(values)  # Adding 1 to avoid log(0)
        # sub_data[is.infinite(sub_data)] <- 1e-10  # Replace -Inf values
        } else {
          cat("Working with logged values\n")
        }

      # values <- log(values)

      # # Replace -Inf values with a small number (e.g., 1e-10)
      # values[values == -Inf] <- 1e-10

      # if (quality_control) {
      #   values <- remove_outliers(values)
      # }
      
      # values <- remove_outliers(values)

      # values[values == -Inf] <- 1e-10


      # values <- log(values)

      # # Replace -Inf values with a small number (e.g., 1e-10)
      # values[values == -Inf] <- 1e-10

      values <- values[!is.infinite(values)]


      result <- skew_gate(values, alpha)
      
      # Append row to results_df
      results_df <- rbind(results_df, data.frame(Patient = imageid, Marker = col_name, Gate = result))
      
      cat("Image ID:", imageid, "Marker:", col_name, "Result:", result, "\n")
    }
  }
  # Record end time
  end_time <- Sys.time()

  # Calculate elapsed time for the loop
  loop_elapsed_time <- end_time - start_time

  # Record start time for CSV creation
  csv_start_time <- Sys.time()
  # write.csv(results_df, "DEV_crevasse_gating_results.csv", row.names = FALSE)  # Adjust the file name as needed
  # write.csv(results_df, paste0("Crev_", save_file_name, "_gating_results.csv"), row.names = FALSE)
  return(results_df)

  # Record end time for CSV creation
  csv_end_time <- Sys.time()

  # Calculate elapsed time for CSV creation
  csv_elapsed_time <- csv_end_time - csv_start_time

  # Print elapsed times
  cat("CSV Reading Elapsed Time:", csv_read_elapsed_time, "\n")
  cat("Loop Elapsed Time:", loop_elapsed_time, "\n")
  cat("CSV Creation Elapsed Time:", csv_elapsed_time, "\n")


  }
  


