#' Standard Error bars
#'
#' @description Plots error bars, taken from https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
#'
#' @param x your x variable
#' @param y is your y variable 
#' @param sdev your standard error (will do y-sdev)
#' @param length the length of the bar (default = 0.05)
#' @param angle the angle (default = 90)
#' @param code the end of the bars (default = 3)
#'
#' @return conditional density object
#' 
error.bar <- function(x,y,sdev,length = 0.05, angle = 90, code = 3){
  arrows(x, y-sdev, x, y+sdev, length= length, angle=angle, code=code)
}
