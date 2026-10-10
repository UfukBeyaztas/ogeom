og_center_data <- function(data, mu_fun) {
  list(t = data$t, w = Map(function(t, w) w - mu_fun(t), data$t, data$w))
}