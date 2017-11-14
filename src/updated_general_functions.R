## Updating all the general functions for the new general_fit() function.

# general_fit()
# general_fit <- function(train, model_type, n, ...) {
  # Fits a time-series model using the model type specified
  # and the training set, and then predicts n steps ahead.
  # This function is agnostic to the specific type of model
  # it is called upon to fit, and can accommodate both
  # traditional time-series models like AR, MA, or ARIMA,
  # as well as unconventional models like method of analogs.
  # Args:
  #   train: The training set, inputted as a single vector.
  #   model_type: The type of model to fit. This should be
  #               inputted as the name of the function as
  #               it is defined in R, but without the 
  #               parentheses; e.g., "lm" or "auto.arima"
  #   n: The number of steps ahead to predict.
  #   ... : These are the miscellaneous arguments passed to
  #         the function specified in model_type. These arguments
  #         are of two kinds: strings which will be converted into
  #         arguments for the model function, and values for any
  #         variables that may be contained in these arguments.
  #         For instance, to fit the model auto.arima(train, max.q = q),
  #         where q = 5, the arguments to ... would be:
  #         "max.q = q", q
  #         q would need to have already been assigned the value of 5;
  #         else, the function will return an error message showing that
  #         it does not know the value of q.
  # Returns:
  #   The n predicted values from the model function.
#   args <- unlist(list(...), recursive = FALSE)
#   which_strings <- unlist(lapply(args, function(x) ifelse(class(x) == "character", TRUE, FALSE)))
#   fit_model_args <- unlist(args[which(which_strings)])
#   fit_model <- paste0(model_type, "(train, ", paste(fit_model_args, collapse = ", "), ")")
#   g <- eval(parse(text = fit_model))
#   forecast(g, h = n)
# }

unearth_levels <- function(l) {
  next_level <- unlist(l, recursive = FALSE)
  if (length(next_level) == 1 & class(next_level) == "list") {
    unearth_levels(next_level)
  }
  else {
    return(next_level)
  }
}

general_fit <- function(train, model_type, n, ...) {
  args <- list(...)
  
  # We have to see whether args represents a list of arguments, or whether
  # it is a list containing a list of arguments, as it will be if args is
  # passed from another function rather than as a series of arguments.
  
  args_list <- unearth_levels(args)
  
  # Now we actually fit the model and make predictions.
  
  fit_model_args <- lapply(args_list, function(x) ifelse(class(x) == "character", x, NA)) %>%
    unlist(.) %>%
    sapply(., function(x) if (is.na(x) == FALSE) {return(x)}) %>% unlist(.)
  fit_model <- paste0(model_type, "(train, ", paste(fit_model_args, collapse = ", "), ")")
  g <- eval(parse(text = fit_model))
  if (is.numeric(g$mean)) {
    result <- list(lower = g$lower,
                   mean = g$mean,
                   upper = g$upper)
    return(result)
  }
  else {
    forecast(g, h = n)
  }
}

# predict_and_plot()
predict_and_plot <- function(pi_vec, model_type, k, with_bar = FALSE, ...) {
  # Computes the absolute error of a model's prediction of PI deaths k weeks ahead of time,
  # where k is an initially specified value that increments by 1 for a set number of times,
  # determined by the length of pi_vec.
  # The error for each number of weeks out from the target week is plotted in a single graph.
  # Args:
  #   pi_vec: The PI death data being inputted that will be split into training and
  #           testing sets.
  #   model_type: The name of the time-series model inputted as a string, e.g., "auto.arima"
  #   k: The length of the training set.
  #   with_bar: A boolean that allows the user to choose whether the plot will contain 95%
  #             confidence interval bars.
  # Returns:
  #   The plot of the absolute errors in predicting the value at the last week - that is,
  #   the last value of pi_vec - from k weeks out from the beginning of the training set,
  #   k + 1 weeks out, k + 2 weeks out, and so on to 1 week out from the last week.
  model_args <- list(...)
  train <- pi_vec[1:k]
  unseen <- pi_vec[(k + 1):length(pi_vec)]
  
  mean_error <- vector()
  upr_error <- vector()
  lwr_error <- vector()
  
  while (length(unseen) > 1) {
    test_predict <- unseen[1] %>% ifelse(is.na(.), mean(train, na.rm=TRUE), .)
    pred <- general_fit(train = train, model_type = model_type, n = 1, model_args)
    
    if (class(pred) == "forecast") {
      mean_error <- c(mean_error, abs(pred$mean - test_predict))
      upr_error <- c(upr_error, abs(pred$upper - test_predict))
      lwr_error <- c(lwr_error, abs(pred$lower - test_predict))
    }
    else {
      idx <- length(train)
      mean_error <- c(mean_error, abs(pred$mean[idx] - test_predict))
      upr_error <- c(upr_error, abs(pred$upper[idx] - test_predict))
      lwr_error <- c(lwr_error, abs(pred$lower[idx] - test_predict))
    }

    train <- c(train, unseen[1])
    unseen <- unseen[-1]
  }
  
  # Create the data frame
  d.frame <- as.data.frame(cbind(1:length(mean_error), mean_error, upr_error, lwr_error))
  names(d.frame) <- c("Weeks.Ahead", "MAE", "Upper", "Lower")
  
  # Plot the declining prediction error, with or without confidence interval bars,
  # as specified in the with_bar argument
  
  if (with_bar == TRUE) {
    k_weeks_plot <- ggplot(data=d.frame, aes(x=Weeks.Ahead, y=MAE, group=1)) + geom_line() + geom_point() +
      geom_errorbar(data=d.frame, aes(ymax=Upper, ymin=Lower)) +
      xlab("Weeks Ahead") + ylab("Absolute Error") +
      scale_x_reverse() + ggtitle(paste("Predicting", k, "Weeks Ahead"))
  }
  else {
    k_weeks_plot <- ggplot(data=d.frame, aes(x=Weeks.Ahead, y=MAE, group=1)) + geom_line() + geom_point() +
      xlab("Weeks Ahead") + ylab("Absolute Error") +
      scale_x_reverse() + ggtitle(paste("Predicting", k, "Weeks Ahead"))
  }
  
  return(k_weeks_plot)
}


# Updated predict_with_noise()
predict_with_noise <- function(pi_vec, model_type, steps_ahead, sigma, seed = 1995, ...) {
  # Computes noisy predictions for specified number of steps ahead based on the
  # original training set and model, with new "observations" added to the training
  # set from the noisy predictions at each step. The results are then plotted with
  # a 95% confidence band.
  # Args:
  #   pi_vec: The PI death data being inputted that will be split into training and
  #           testing sets.
  #   model_type: The name of the time-series model inputted as a string, e.g., "auto.arima"
  #   steps_ahead: Number of noisy predictions that will be made.
  #   sigma: A scale term to increase or decrease the noise from the model residuals.
  #   train: The initial time-series training set.
  #   seed: The random number seed to control the random sampling from the model residuals.
  #         Can be set to NULL if reproducibility is not desired.
  # Returns:
  #   A ggplot object with the plot of the model's predictions for PI deaths for the number
  #   of weeks specified in steps.ahead, with the 95% confidence band surrounding the
  #   predicted points.
  
  model_args <- list(...)
  cutoff <- length(pi_vec) - steps_ahead
  train <- pi_vec[1:cutoff]
  test <- pi_vec[(cutoff + 1):length(pi_vec)]

  set.seed(seed)
  pred_interval <- data.frame(pred_mean = c(rep(NA, length(train)), rep(0, steps_ahead)),
                              test = c(rep(NA, length(train)), test),
                              train_vec = c(train, rep(NA, length(test))),
                              weeks = 1:length(pi_vec))
  band_interval <- data.frame(upper = rep(0, steps_ahead), lower = rep(0, steps_ahead),
                     weeks = (length(train) + 1):(length(train) + steps_ahead))

  train_length <- length(train)
  for (i in 1:steps_ahead) {
    pred <- general_fit(train = train, model_type = model_type, n = 1, model_args)
    
    if (class(pred) == "forecast") {
      res_part <- sigma*sample(pred$residuals, 1)
      fixed_part <- pred_interval$pred_mean[i + train_length] <- ifelse(pred$mean > 0, pred$mean, 0)
      band_interval$upper[i] <- ifelse(pred$upper[2] > 0, pred$upper[2], 0)
      band_interval$lower[i] <- ifelse(pred$lower[2] > 0, pred$lower[2], 0)
    }
    else if (class(pred) != "forecast") {
      res_part <- sigma*rnorm(n = 1, mean = 0, sd = 4)
      fixed_part <- pred_interval$pred_mean[i + train_length] <- ifelse(pred$mean[i + train_length] > 0,
                                                                        pred$mean[i + train_length],
                                                                        0)
      band_interval$upper[i] <- ifelse(pred$upper[i + train_length] > 0, 
                                       pred$upper[i + train_length],
                                       0)
      band_interval$lower[i] <- ifelse(pred$lower[i + train_length] > 0,
                                       pred$upper[i + train_length],
                                       0)
    }
    
    add_noise <- res_part + fixed_part
    train <- c(train, add_noise)
  }
  
  pred_interval <- melt(pred_interval, id = "weeks")
  print("band_interval$lower")
  print(band_interval$lower)

  rib_plot <- ggplot(data = pred_interval, aes(x = weeks, y = value, colour = variable)) +
    geom_point() + geom_line() +
    scale_colour_manual(values = c("skyblue", "magenta", "orange"),
                        name = "Series", breaks = c("pred_mean", "test", "train_vec"),
                        labels = c("Predicted", "Real", "Training Set")) +
     geom_ribbon(data = band_interval, aes(x = weeks, ymin = lower, ymax = upper), alpha = 0.3, inherit.aes = FALSE) +
     ggtitle("Prediction Plot - 95% Confidence Band") + scale_x_continuous(breaks=seq(1, length(train) + steps_ahead, 1)) +
     xlab("Weeks") + ylab("Predicted PI Deaths")
  return(rib_plot)
}

# ts_cv()

ts_cv <- function(pi_vec, train_length, model_type, k, type, change_by, type_transform, ...) {
  # Runs Hyndman's rolling cross-validation on the time-series data and returns a plot of the 
  # residuals.
  # Args:
  #   pi_vec: The vector of PI deaths for the entire flu season in question.
  #   train_length: The length of the initial training set, inputted as an integer.
  #   model_type: The model type name used inputted as a string, e.g., "HoltWinters" or "auto.arima"
  #   k: The number of weeks by which the size of the training set increases on each iteration,
  #      if an option for cross-validation is chosen in which the training set's size increases.
  #   type: The method to be used for cross-validation. Available options are moving_train,
  #         where the training set remains the same size but is moved forward by one week
  #         on each iteration, exp_fixed_k, where the training set always begins at the original start
  #         but increases by k weeks in size (where k is the argument passed to ts_cv()),
  #         and exp_moving_k, which works the same way as exp_fixed_k, except that k
  #         decreases by the reciprocal of change_by (which is passed as an argument
  #         to ts_cv()), until k decays to less than 1. At that point, the training set
  #         simply increases by 1 on each iteration until the completion of the validation.
  #   change_by: The argument mentioned in type; it is the parameter by which k decays on
  #              each iteration of the cross-validation.
  #   type_transform: The type of transformation to be applied to the time-series. Options
  #                   are "log", "box.cox", or "none."
  # Returns:
  #   Global variables "pred" and "res," which contain the predicted values and residuals
  #   from the cross-validation, and the plot of residuals by week.
  
  
  # We first test to make sure that there is a gap of at least one week between the end of the
  # training set and the season's end; if there is, we proceed with cross-validation. If not,
  # the cross-validation is finished, and we skip to the end of the function.
  
  model_args <- list(...)
  
  if ((length(pi_vec) - train_length) > 0) {
    # This block checks to see what transformation, if any, should be used, and completes it
    # if one is specified. It then specifies the training set and testing set (a single
    # week which follows the training set).
    if (type_transform == "log") {
      to.each <- function(x) ifelse(is.na(x) == FALSE, x + 0.00001, mean(pi_vec, na.rm = TRUE) + 0.00001)
      log_pi <- pi_vec %>% sapply(., to.each) %>% log(.)
      train <- log_pi[1:train_length]
      test <- pi_vec[train_length + 1]
    }
    else if (type_transform == "box.cox") {
      this_lambda <<- BoxCox.lambda(pi_vec)
      bc_transform <- BoxCox(pi_vec, lambda = this_lambda)
      train <- bc_transform[1:train_length]
      test <- bc_transform[train_length + 1]
    }
    else if (type_transform == "none") {
      train <- pi_vec[1:train_length]
      test <- pi_vec[train_length + 1]
    }
    
    fit <- general_fit(train = train, model_type = model_type, n = 1, model_args)
    pred <- fit$mean
    
    pred_values <- c(pred_values, pred)
    .GlobalEnv$res <- c(.GlobalEnv$res, (test - pred))
    # At this point the residual for this prediction has been determined and assigned to res,
    # and likewise the predicted value itself has been assigned to pred.values. Both are
    # global variables so that they can be incremented within this function and called directly
    # outside of the function.
    
    
    # The remaining conditionals determine the type of cross-validation preferred, and recursively
    # call ts.cv with parameters that are updated based on the type of cross-validation.
    if (type == "move_train") {
      train_length <- train_length + 1
      ts_cv(pi_vec = pi_vec, train_length = train_length, model_type = model_type, k = NULL,
            type = "move_train", change_by = NULL, type_transform = type_transform, model_args)
    }
    else if (type == "exp_fixed_k") {
      train_length <- train_length + k
      ts_cv(pi_vec = pi_vec, train_length = train_length, model_type = model_type, k = k, 
            type = "exp_fixed_k", change_by = NULL, type_transform = type_transform, model_args)
    }
    else if (type == "exp_moving_k") {
      train_length <- ifelse(floor(train_length + k / change_by) == length(pi_vec), 
                    train_length + 1, floor(train_length + k / change_by))
      ts_cv(pi_vec = pi_vec, train_length = train_length, model_type = model_type, k = (k / change_by),
            type = "exp_moving_k", change_by = change_by, type_transform = type_transform,
            model_args)
    }
  }
  else if (length(pi_vec) - train_length == 0) {
    # Finally, we complete the function by plotting the residuals by week. No return() statement is needed,
    # since res and pred.values are both global variables that can be accessed directly without having
    # to be called through this function's output.
    if (type_transform == "log") {
      for (i in 1:length(res)) {
        res[i] <- log.reverse(res[i])
      }
    }
    if (type_transform == "box.cox") {
      for (i in 1:length(res)) {
        res[i] <- box.cox.reverse(res[i], this_lambda)
      }
    }
    res.dat <- data.frame(residuals=res)
    ggplot(res.dat, aes(x=1:length(residuals), y=residuals)) + geom_point(shape=1) + xlab("Intervals") + 
      ylab("Residuals") + ggtitle("Residuals")
  }
}

get_city <- function(city.name, dataset = data) {
  # Extracts the PI deaths data from the larger CDC dataset by city.
  # Args:
  #   city.name: The name of the city inputted as a string with ordinary 
  #                   capitalization and spacing, e.g., "New Providence"
  #   dataset: The data from which the city data is extracted.
  #			 By default this is the CDC dataset.
  # Returns:
  #   A data frame with four columns, whose rows each represent one week: 
  #     1. Year: The year as an int
  #     2. Week.Number: The number of the week in the year, e.g., 
  #				the first week of the year is 1, second is 2, etc., until
  #                     week 52, when a new year begins and Week.Number 
  #			resets to 1.
  #     3. week: The week value for the week as a time series object which
  #			 begins at 2000.019 and increments by 0.019
  #              for each week of the year. Unlike Week.Number, 
  #		these values do not reset at the beginning of each year,
  #              but continue incrementing by the same amount each time.
  #     4. PI: The number of pneumonia and influenza (PI) deaths 
  #		for the week represented by the row.
  new.data <- filter(dataset, City.Name == city.name) %>% 
    select(., Year, Week.Number, P, I)
  time <- round(as.numeric(as.character(new.data$Year + 
                                          (new.data$Week.Number)/53)), digits=3)
  new.data$week <- time
  new.data$PI <- new.data$P + new.data$I
  new.data <- select(new.data, Year, Week.Number, week, PI)
  return(new.data)
}





