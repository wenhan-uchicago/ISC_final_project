library(mvtnorm)
library(ggplot2)
#library(snow)

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


################################################################
## defining funcitons below; final function is read_data() & reads_p2_given_reads_p1()


## iterating possible s values, and save log-likelihood in data.frame log_results
calculate_log_likelihood <- function(step_size = 0.1) {
    log_results <- data.frame(s = numeric(0), log_likelihood = numeric(0))
    for (s in seq(-0.9, 1.9, by = step_size)) {
        log <- log_likelihood_per_s_h(c(s, 0), seq_data)
        log_results <- rbind(log_results, data.frame(s = s, log_likelihood = log))
    }

    return(log_results)
}

## read data as data.frame
read_data <- function(location) {
    seq_data <- read.table(location, col.names = c("depth_1", "reads_1", "generations_1"))
    seq_data$depth_2 <- c(seq_data$depth[2:length(seq_data$depth)], -Inf)
    seq_data$reads_2 <- c(seq_data$reads_1[2:length(seq_data$reads_1)], -Inf)
    seq_data$generations_2 <- c(seq_data$generations_1[2:length(seq_data$generations_1)], -Inf)
    seq_data$generations_in_between <- seq_data$generations_2 - seq_data$generations_1

    ##    return(seq_data)
    return(seq_data[, c("depth_2", "reads_2", "depth_1", "reads_1", "generations_in_between")])
}

process_seq_data <- function(seq_data_single_df) {
    seq_data_matrix <- as.matrix(seq_data_single_df)
    ## since last line has no depth_2, reads_2; truncate it
    seq_data_matrix <- seq_data_matrix[1:(dim(seq_data_matrix)[1] - 1), ]
    return(seq_data_matrix)
    ##apply(seq_data_matrix, MARGIN)
}

################################################################
################################################################


## P(p2 | p1)
p2_given_p1 <- function(absolute_p2, absolute_p1, s, h, N) {
    if (length(absolute_p2) != 1) {
        stop("length of absolute_p2 should be 1.")
    }
    if (length(absolute_p1) != 2*N + 1) {
        stop("Have to iterate all possible p1, from 0 to 2*N")
    }

    i <- absolute_p1; j <- absolute_p2
    ## calculate eta
    eta = ((1 + s) * (i^2) + (1 + s*h) * i * (2 * N - i)) / ((1 + s) * (i^2) + 2 * (1 + s*h) * i * (2 * N - i) + (2*N-i)^2)



    prob <- dbinom(x = absolute_p2, size = 2*N, prob = eta, log = FALSE)

    return(prob)

}

## P(p1 | reads_1, depth_1)
p1_given_depth_1_reads_1 <- function(absolute_p1, depth_p1, reads_p1, N) {
    ## suppose uniform prior P(p1)
    if (length(absolute_p1) != 2*N + 1) {
        stop("Have to iterate all possible p1, from 0 to 2*N")
    }
    if (length(depth_p1) != 1 | length(reads_p1) != 1) {
        stop("depth_p1 and reads_p1 have to be size 1.")
    }

    p1 <- absolute_p1 / (2*N)
    prob <- dbinom(x = reads_p1, size = depth_p1, prob = p1, log = FALSE)

    return(prob / sum(prob))
}

## for calculating single absolute_p2 value

p2_given_depth_p1_reads_p1_singular <- function(absolute_p2, depth_p1, reads_p1, s, h, N) {
    absolute_p1 <- 0:(2*N)

    ## probability of single p2 value
    if (length(absolute_p2) == 1) {
        temp_prob <- sum(p2_given_p1(absolute_p2, absolute_p1, s, h, N) * p1_given_depth_1_reads_1(absolute_p1, depth_p1, reads_p1, N))
        return(temp_prob)
    } else {
        stop("this step is only used to calculate the else if clause; as this is not probability; they won't sum to 1. Thus, this clause is only recruited to calcaulate all of them and divide by sum(prob).")
    }
}




## P(p2 | reads_1, depth_1)
## P(p2 | Model, depth_p1, reads_p1) = sum_all_p1{P(p2 | Model, p1) * P(p1 | Model, depth_p1, reads_p1)}
p2_given_depth_p1_reads_p1 <- function(absolute_p2, depth_p1, reads_p1, s, h, N) {

    ## probability of single p2 value
    if (length(absolute_p2) == 1) {
        stop("should use the singular version.")
    } else if (length(absolute_p2) == 2*N + 1) {
        ## could use for (j in 0:(2*N)), but I think maybe use apply is better for parallel computing
        ## or maybe not use parallel here??
        ##prob <- parApply(cl, X = matrix(absolute_p2, nrow = length(absolute_p2)), MARGIN = 1, FUN = p2_given_depth_p1_reads_p1_singular, depth_p1 = depth_p1, reads_p1 = reads_p1, s = s, h = h, N = N)
        prob <- apply(X = matrix(absolute_p2, nrow = length(absolute_p2)), MARGIN = 1, FUN = p2_given_depth_p1_reads_p1_singular, depth_p1 = depth_p1, reads_p1 = reads_p1, s = s, h = h, N = N)

        return(prob)

        ## no standarizing constant is necessary!! Since this sum_all_B P(A, B) = P(A) is true probability
##        return(prob / sum(prob))
    } else {
        stop("the return value for length(absolute_p2) != 2*N + 1 will not be probability. since they need to be divided by sum(prob).")
    }

}


## P(reads_p2 | p2, depth_p2)
reads_p2_given_p2 <- function(reads_p2, depth_p2, absolute_p2, N) {
##    if (length(reads_p2) != depth_p2 + 1) {
##        stop()
##    }

    p2 <- absolute_p2 / (2*N)
    ## no standarizing constant is necessary
    prob <- dbinom(x = reads_p2, size = depth_p2, prob = p2)
    return(prob)
}

## singular version
reads_p2_given_reads_p1_singular <- function(reads_p2, depth_p2, reads_p1, depth_p1, s, h, N) {
    if (length(reads_p2) != 1) {
        stop("this is only the singular version, to calculate the values for sum(prob), the standarizing constant; the output is not probability.")
    }
    absolute_p2 <- 0:(2*N)

    temp_prob <- sum(reads_p2_given_p2(reads_p2, depth_p2, absolute_p2, N) * p2_given_depth_p1_reads_p1(absolute_p2, depth_p1, reads_p1, s, h, N))

    return(temp_prob)
}

## P(reads_p2 | reads_p1, depth_p1, depth_p2, s, h)
## reads_p2_given_reads_p1 <- function(reads_p2, depth_p2, reads_p1, depth_p1, s, h, N) {

## reads_p2_given_reads_p1 <- function(depth_p2, reads_p2, depth_p1, reads_p1, generations_in_between, s, h, N) {
reads_p2_given_reads_p1 <- function(data, s, h, N) {
    if (length(data) != 5) {
        stop("not enough data in reads_p2_given_reads_p1().")
    }
    depth_p2 <- data[1]; reads_p2 <- data[2]; depth_p1 <- data[3]; reads_p1 <- data[4]
    generations_in_between <- data[5]

    ## just for calculating the standarizing constant
##    temp_reads_p2 <- 0:(depth_p2)

    ##    if(length(reads_p2) == 1) {
    if(max(length(reads_p2), length(depth_p2), length(reads_p1), length(depth_p1)) == 1) {
        ## seems no standarizing constant is needed for this
##        total_prob <- apply(X = matrix(temp_reads_p2, nrow = length(temp_reads_p2)), MARGIN = 1, FUN = reads_p2_given_reads_p1_singular, depth_p2 = depth_p2, reads_p1 = reads_p1, depth_p1 = depth_p1, s = s, h = h, N = N)

        prob <- reads_p2_given_reads_p1_singular(reads_p2, depth_p2, reads_p1, depth_p1, s, h, N)

        ## similarly, no standrizing constant is necessary
        ## return the log value
        return(log(prob))
        ##return(prob / sum(total_prob))
    } else {
        stop("need to make sure length(all parameters) == 1, since reads_p2 is always a KNOWN value.")
    }

}


log_likelihood_per_s <- function(s) {

    stop("probably should use log_likelihood_per_s_h() now.")

    ##    sum(apply(X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = s, h = 1, N = 100))
    sum(parApply(cl, X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = s, h = 1, N = 100))
}

################################################################
## 2d version

log_likelihood_per_s_h <- function(parameters, seq_data) {
    if (is.matrix(seq_data) == TRUE) {
        ## no replications
        ## seq_data = seq_data_matrix, is a matrix
        s <- parameters[1]; h <- parameters[2]

        sum(apply(X = seq_data, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = s, h = h, N = 100))
        ##        sum(parApply(cl, X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = s, h = h, N = 100))



#        sum(parApply(cl, X = seq_data, MARGIN = 1, FUN = function(data) reads_p2_given_reads_p1(data, s = s, h = h, N = 100)))

    } else if (is.list(seq_data) == TRUE) {
        ## with replicate sequencing data; now, restricted to be sequenced at the same generations
        ## seq_data = list(seq_data_matrix), is a list of seq_data_matrix
        sum(sapply(X = seq_data, FUN = function(seq_data_matrix) log_likelihood_per_s_h(parameters, seq_data_matrix)))
    } else {
        stop("Something wrong here. seq_data should either be a matrix or a list.")
    }
}



################################################################

## step <- 101

## new method!!

## s_matrix <- matrix(seq(-0.95, 0.95, length = step), ncol = 1)
## out <- parApply(cl, X = s_matrix, MARGIN = 1, FUN = log_likelihood_per_s)

################################################################
################################################################

## using MH-algorithm

## use mu = log(s+1), so that mu is in range R; s = exp(mu) - 1;
## fmu(mu) = fs(exp(mu) - 1) * abs(exp(mu)) = fs(exp(mu) - 1) * exp(mu)


met_norm_proposal_1d <- function(iters) {

    logxvec <- numeric(iters)

    likelihood_fun <- function(logx) {
        ## sum of log probabilities
        ##        p <- sum(parApply(cl, X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = x, h = 1, N = 100))
        ##        p <- exp(log_likelihood_per_s(exp(logx) - 1)) * exp(logx)
        p <- log_likelihood_per_s(exp(logx) - 1) + logx
        ## print(p)
        return(p)
    }


    accept <- 0

    logx <- -1

    for (i in 1:iters) {
        logxs <- logx + 0.5 * rnorm(1, mean = 0, sd = 1)
        ## likelihood
        ##        A <- likelihood_fun(xs) / likelihood_fun(x)
        print(paste("#iteration = ", i))
##p        print(exp(logxs) - 1)
        if (exp(logxs) - 1 > -1) {
            ##            A <- exp(likelihood_fun(logxs)) / exp(likelihood_fun(logx))
            A <- exp(likelihood_fun(logxs) - likelihood_fun(logx))

            if (runif(1) < A) {
                logx <- logxs
                accept = accept + 1
            }
        }

        logxvec[i] <- logx
    }
    print(accept / iters)
    return(logxvec)
}



################################################################
################################################################
## 2d version of mh method

met_norm_proposal_2d <- function(iters) {

    logx_matrix <- matrix(0, nrow = iters, ncol = 2)

    likelihood_fun <- function(logx) {
        ## sum of log probabilities
        ##        p <- sum(parApply(cl, X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = x, h = 1, N = 100))
        ##        p <- exp(log_likelihood_per_s(exp(logx) - 1)) * exp(logx)
        s <- exp(logx[1]) - 1; h <- logx[2]
        p <- log_likelihood_per_s_h(c(s, h), seq_data) + logx[1]
        ## print(p)
        return(p)
    }


    accept <- 0

    logx <- c(-1, 1)

    for (i in 1:iters) {
        ##        logxs <- logx + 0.5 * rnorm(1, mean = 0, sd = 1)
##        print(logx)
        logxs <- logx + 0.5*rmvnorm(1, mean = c(0, 0), sigma = matrix(c(1, 0, 0, 1), byrow = T, nrow = 2))
        ## likelihood
        ##        A <- likelihood_fun(xs) / likelihood_fun(x)
        print(paste("#iteration = ", i))
##p        print(exp(logxs) - 1)
        ##        if (exp(logxs[1]) - 1 > -1) {
        temp_s <- exp(logxs[1]) - 1; temp_h <- logxs[2]

        ## set h to be (-1, 2); s to be (-1, 2); because if s = 100, h = 1, then the ratio is 1+100 : 1+100*1 : 1....
        ##        if ((temp_s > -1) & (temp_s*temp_h > -1) & (temp_h > -1) & (temp_h < 2)) {
        if ((temp_s > -1) & (temp_s*temp_h > -1) & (temp_h > -1) & (temp_h < 2) & (temp_s < 2)) {
            ##            A <- exp(likelihood_fun(logxs)) / exp(likelihood_fun(logx))

            A <- exp(likelihood_fun(logxs) - likelihood_fun(logx))
##            print(A)
            if (runif(1) < A) {
                logx <- logxs
                accept = accept + 1
            }
        }

        logx_matrix[i,] <- logx
    }
    print(accept / iters)
    return(logx_matrix)
}




