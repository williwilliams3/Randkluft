# Load required library
library(moments)


remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
low_threshold <- quantile(target, low_percentile / 100)
high_threshold <- quantile(target, high_percentile / 100)
target <- target[target >= low_threshold & target <= high_threshold]
return(target)
}

skew_gate_plot <- function(dataframe) {

  # #  Define remove_outliers function
  # remove_outliers <- function(target, level = 0.999) {
  #   target <- target[which(target < quantile(target, level))]
  #   return(target)
  # }

  # Define skew_gate function
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
      percentage_removed = 0
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
      percentage_removed = perc_removed
    )
    }

  # Record start time for CSV reading
  csv_read_start_time <- Sys.time()

  data <- dataframe # Replace with your actual CSV file path

  # Record end time for CSV reading
  csv_read_end_time <- Sys.time()

  # Calculate elapsed time for CSV reading
  csv_read_elapsed_time <- csv_read_end_time - csv_read_start_time

  # # Get unique values from the imageid column
  # unique_imageids <- unique(data$imageid)

  # Define skew_gate function above

  alpha <- 0.5  # Set your desired alpha value

  # Record start time
  start_time <- Sys.time()
# Loop through each unique imageid
  # Calculate the maximum value for each column (excluding imageid)
  max_values <- sapply(data, max)

  max_values <- max_values[-1]  # Remove the imageid column

  
  # Loop through each column in sub_data
  
  values <- data

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

  values <- values[!is.infinite(values)]

  # values <- log(values)

  # # Replace -Inf values with a small number (e.g., 1e-10)
  # values[values == -Inf] <- 1e-10

  result <- skew_gate(values, alpha)
  
  # Append row to results_df      
  # cat("Image ID:", imageid, "Marker:", col_name, "Result:", result, "\n")

  return(result)

}