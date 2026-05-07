# Load required library
library(moments)

get_gates_csv_single <- function(dataframe, csv_save_file_path) {

  # #  Define remove_outliers function
  # remove_outliers <- function(target, level = 0.999) {
  #   target <- target[which(target < quantile(target, level))]
  #   return(target)
  # }

  remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
    print("wowzer")
  print(target)
  low_threshold <- quantile(target, low_percentile / 100)
  high_threshold <- quantile(target, high_percentile / 100)
  print("barr")
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


  data <- dataframe # Replace with your actual CSV file path

  # # Calculate the maximum value for each column (excluding imageid)
  max_values <- sapply(data[-1], max)


  # List of markers
  markers <- colnames(data)[-1]

  # Create a list to store ggplot objects
  histogram_list <- list()

   # Define skew_gate function above

  alpha <- 0.5  # Set your desired alpha value

  # Create an empty data frame to store the results
  results_df <- data.frame(Patient = character(),
                          Marker = character(),
                          Gate = numeric(),
                          stringsAsFactors = FALSE)

  # Loop over unique imageids
  unique_imageids <- unique(data$imageid)
  for (imageid in unique_imageids) {
    for (marker in markers) {
      # Subset data for the current imageid and marker
      subset_data <- data[data$imageid == imageid, c("imageid", marker)]

    
      if (any(max_values > 20)) {
          cat("Working with raw data\n")
          if (input$GMM_gate_on_off) {
            cat("gmm on")
            subset_data_hist <- log(subset_data[-1])
          } else {
            cat("gmm off")
            subset_data_hist <- log(subset_data[-1])
          }
        } else {
          cat("Working with logged values\n")
          subset_data_hist <- subset_data
        }


  
      # Calculate gate value using skew_gate function
      result_to_plot <- skew_gate_plot(subset_data[[marker]])  # Replace with your desired alpha
      skewness_value <- result_to_plot$skewness
      cutoff_value <- result_to_plot$cutoff

       # Append row to results_df
      results_df <- rbind(results_df, data.frame(Patient = imageid, Marker = marker, Gate = cutoff_value))
      
      cat("Image ID:", imageid, "Marker:", marker, "Result:", cutoff_value, "\n")
    }}




  return(results_df)


  }
  


