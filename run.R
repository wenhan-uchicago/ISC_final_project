library(mvtnorm)
library(ggplot2)
#library(snow)

source("./functions_codes.R")
##location <- "./seq_data_1.txt"
##location <- c("./seq_data_1.txt", "./seq_data_2.txt", "./seq_data_3.txt")
location <- c("./seq_data_1.txt", "./seq_data_2.txt", "./seq_data_3.txt", "./seq_data_4.txt", "./seq_data_5.txt")

################################################################


if (length(location) == 1) {
    ## seq_data is a matrix now
    seq_data <- read_data(location)
    seq_data <- process_seq_data(seq_data)
} else if (length(location) > 1) {
    ## seq_data is a list of matrix now
    seq_data <- lapply(X = location, FUN = read_data)
    seq_data <- lapply(X = seq_data, FUN = process_seq_data)

} else {
    stop("the location of seq data is no valid.")
}


result <- calculate_log_likelihood(0.1)
result$likelihood <- exp(result$log_likelihood)
result$likelihood <- result$likelihood / sum(result$likelihood)

ggplot(data = result, aes(x = s, y = likelihood)) + geom_point() + geom_line() + labs(title = "Likelihood of Different s Values", x = "selection coefficients", y = "Likelihood") + theme(plot.title = element_text(lineheight=.8, face="bold"))
