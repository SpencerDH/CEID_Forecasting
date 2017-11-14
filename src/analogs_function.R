## The method of analogs in functional form.

scale_01 <- function(x) {
  scale <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  return(scale)
}

# Function to convert the 26 week datasets into 25 week datasets.
# convert_df <- function(df) {
#   if (nrow(df) == 26) {
#     new_pi <- (df$PI[13] + df$PI[14]) / 2
#     df$PI[13] <- new_pi
#     df <- df[-14, ]
#     df$PI <- (20 - 0)*(df$PI - min(df$PI, na.rm = TRUE)) / (max(df$PI, na.rm = TRUE) - min(df$PI, na.rm = TRUE))
#   }
#   return(df)
# }

convert_df <- function(df) {
  if (nrow(df) == 26) {
    new_pi <- (df$PI[13] + df$PI[14]) / 2
    df$PI[13] <- new_pi
    df <- df[-14, ]
  }
  return(df)
}

# Function to find NAs and replace them with some arbitrary high value.
repeal_and_replace <- function(value) {
  if (value == 0) {
    value <- 1000
  }
  return(value)
}

# Function to get the difference between the max and min of a vector of numbers.
range_diff <- function(x) {
  return(max(x$PI, na.rm = TRUE) - min(x$PI, na.rm = TRUE))
}

##################################
#        get_analogs()           #
##################################

get_analogs <- function(train, start_year, data) {
  # Fits a time series model based on Mandel's
  # method of analogs. The model splits all past
  # data into candidates for analogs, where every
  # city / flu season combination is stored as a
  # small (25 rows) dataframe in a large list
  # that contains all such candidates.
  # The function then evaluates the candidates
  # to see which ones most closely match the
  # training set given different metrics.
  # Currently the function supports two such
  # metrics - the SSE of the time series values
  # themselves, and the absolute difference of
  # the standard deviations of the training set
  # and candidate analog. The next iteration of the
  # function will allow the user to specify additional
  # covariates to be used in evaluating the potential
  # analogs.
  # Args:
  #   train: The training set, inputted as a vector.
  #          Generally this will be the first 4 - 12
  #          weeks of a given flu season, following
  #          the standard convention that seasons begin
  #          in week 40 of a given year.
  #   start_year: The year in which the flu season starts.
  #   data: The entire flu season dataframe.
  # Returns:
  #   The average of the best analogs as a vector.
  #   This vector's length will be equal to the length
  #   of the entire season (considered here to be 25
  #   weeks) minus the length of the training vector.

  # Function to get the classes.
  class_fun <- function(n) {
    vec <- unname(unlist(nm_years[n, ]))
    new_vector <- c(as.numeric(vec[1]), vec[2])
    return(new_vector)
  }
  
  # Function to get SSE.
  get_sse <- function(train, city) {
    diff <- train - city$PI[1:length(train)]
    return(sum(sapply(diff, function(x) x^2), na.rm = TRUE))
  }
  
  # Scale the training set.
 # train <- (20 - 0)*(train - min(train, na.rm = TRUE)) / (max(train, na.rm = TRUE) - min(train, na.rm = TRUE))
  .GlobalEnv$train_list[[length(.GlobalEnv$train_list) + 1]] <- train
  
  sub_data <- filter(data, (Year <= (start_year - 1)) | (Year == start_year & Week.Number <= 12))
  city_names <- as.character(rep(unique(sub_data$City.Name), each = (start_year - 2000)))
  years <- rep(2000:(start_year - 1), length(unique(sub_data$City.Name)))
  
  nm_years <- as.data.frame(cbind(years, city_names), stringsAsFactors = FALSE)
  damian <- lapply(1:nrow(nm_years), class_fun)
  raw_split_cities <- lapply(damian, function(x) filter(sub_data, 
                                                        ((Year == as.numeric(x[1]) & Week.Number >= 40) 
                                                         | (Year == (as.numeric(x[1]) + 1) & Week.Number <= 12)), 
                                                        City.Name == x[2]))
  
  # Now we scale all of the city PI data to the 0 to 20 range.
  split_cities <- lapply(raw_split_cities, convert_df)
  .GlobalEnv$total_pi[[length(.GlobalEnv$total_pi) + 1]] <- lapply(split_cities, function(x) return(x$PI)) %>% 
    unlist(.)
  
  
  # Get the SSE for each potential analog.
  .GlobalEnv$train_head[[length(.GlobalEnv$train_head) + 1]] <- head(train)
  sse_list <- sapply(split_cities, get_sse, train = train) %>% sapply(., repeal_and_replace)
  .GlobalEnv$sse_archive[[length(.GlobalEnv$sse_archive) + 1]] <- sse_list
  
  # Get the 95th percentile for the SSEs.
  percentile <- sort(sse_list, decreasing = TRUE)[round(0.95*length(sse_list))]
  .GlobalEnv$percentile_vector <- c(.GlobalEnv$percentile_vector, percentile)
  
  # Get the standard deviation for each potential analog, as well as for the training set.
  analog_sd <- sapply(split_cities, function(x) sd(x$PI))
  train_sd <- sd(train, na.rm = TRUE)
  
  criterion <- integer(0)
  n <- 1
  while (length(criterion) == 0) {
    criterion <- which(sse_list <= percentile & 
                         abs(analog_sd - train_sd) <= n)
    matches <- split_cities[criterion] %>% 
      lapply(., function(x) x$PI)
    divide <- function(x) return(x / sum(x, na.rm = TRUE))
    weights <- sse_list[criterion] %>% divide(.)
    avg_matches <- mapply(function(x, y) return(unlist(x)*y), matches, weights)
    n <- n + 1
    if (n > 1000) {
      break
    }
  }
  
  # analogs <- data.frame(matrix(unlist(matches), ncol = length(matches))) %>%
  #   rowSums(., na.rm = TRUE) / length(matches)
  
  analogs <- matches[[sample(length(matches), size = 1)]]

  avg_range <- lapply(matches, function(x) return(sd(x))) %>%
    unlist(.) %>%
    mean(.)
  
  # Here we're randomly getting the standard deviation from a normal distribution
  # with the mean avg_range and standard deviation 1.
  
  analog_sd <- rnorm(n = 1, mean = avg_range, sd = 1)
  
  analogs_upper <- analogs + 1.96*analog_sd
  analogs_lower <- analogs - 1.96*analog_sd
  
  result <- list(analogs_lower, analogs, analogs_upper)
  names(result) <- c("lower", "mean", "upper")

  return(result)
}

# plot_analogs()

plot_analogs <- function(real, analogs) {
  real_data <- data.frame(weeks = 1:25, PI = real)
  pred_data <- data.frame(weeks = (25 - length(analogs) + 1):25, PI = analogs)
  ggplot(data = real_data, aes(x = weeks, y = PI, group = 1)) + geom_point(colour = "green") + 
    geom_line(colour = "green") +
    geom_point(data = pred_data, aes(x = weeks, y = PI, group = 1), colour = "blue") +
    geom_line(data = pred_data, aes(x = weeks, y = PI, group = 1), colour = "blue") +
    xlab("Weeks") + ylab("PI")
}






