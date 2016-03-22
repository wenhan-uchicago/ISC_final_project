library(mvtnorm)
library(ggplot2)
library(snow)

## P(p2 | p1)
p2_given_p1_single_gen <- function(absolute_p2, absolute_p1, s, h, N) {
    stop("use the new one!")

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

    ## this does NOT need a standarizing constant!
    ## also, I was wrong before, by making sum_p1(P(p2|p1)) == 1, this should actually be sum_p2(P(p2|p1)) == 1!!!!
                                        #    return(prob/sum(prob))
}

## P(p_generation_j | p_generation_1); allow between generations other than 1
## p2_given_p1 <- function(absolute_p2, absolute_p1, s, h, N, generations_in_between) {
p2_given_p1_only_once_no_use <- function(absolute_p2, absolute_p1, s, h, N, generations_in_between) {
    ## this gives the probablity matrix of all p2, given all possible p1
    ## each column is all p2 | given single p1
    ## each row is single p2 | given all p1
    if (length(absolute_p2) != 2*N + 1) {
        ##        stop("length of absolute_p2 should be 1.")
        stop("length of absolute_p2 should be 2*N+1.")
    }
    if (length(absolute_p1) != 2*N + 1) {
        stop("Have to iterate all possible p1, from 0 to 2*N")
    }

    stop("I am called. Which should not happen. Should only be called in calculate_slim.R")
##    return(matrix(1, nrow = 201, ncol = 201))

    temp_all_p2_given_single_p1 <- function(absolute_p1, s, h, N) {
        if (length(absolute_p1) != 1) {
            stop("absolute_p1 here should be of length 1.")
        }

        absolute_p2 <- 0:(2*N)
        i <- absolute_p1; j <- absolute_p2

        ## calculate eta
        eta = ((1 + s) * (i^2) + (1 + s*h) * i * (2 * N - i)) / ((1 + s) * (i^2) + 2 * (1 + s*h) * i * (2 * N - i) + (2*N-i)^2)
        if (length(eta) != 1) {
            stop("Wrong")
        }
        ## probabilities of p2 = 0-2N, given a single p1 value
        prob <- dbinom(x = absolute_p2, size = 2*N, prob = eta, log = FALSE)

        return(prob)
    }

    ## sapply(X = 0:(2*N), FUN = function(p1) temp_all_p2_given_single_p1(p1, s, h, N))
    ##    prob_at_p1 <- rep(x = 1/(2*N+1), times = 2*N + 1)
    ##    single_p2_given_single_p1 <- function(absolute_p2, absolute_p1, s, h, N, generations_in_between) {
    single_p2_given_all_p1 <- function(absolute_p2, s, h, N, generations_in_between) {
        if (length(absolute_p2) != 2*N + 1) {
            stop("only when 2*N + 1 absolute_p2")
        }

        absolute_p1 <- 0:(2*N)

        ##        prob_at_p2 <- temp_all_p2_given_single_p1(absolute_p1, s, h, N)
        prob_at_p2_matrix <- sapply(X = 0:(2*N), FUN = function(absolute_p1) temp_all_p2_given_single_p1(absolute_p1, s, h, N))

        prob_at_p1_matrix <- prob_at_p2_matrix
        generations_in_between <- generations_in_between - 1

        while (generations_in_between > 0) {
            transition_matrix <- sapply(X = 0:(2*N), FUN = function(p1) temp_all_p2_given_single_p1(p1, s, h, N))
            prob_at_p2_matrix <- transition_matrix %*% prob_at_p1_matrix

            prob_at_p1_matrix <- prob_at_p2_matrix
            generations_in_between <- generations_in_between - 1
        }
        ##return(prob_at_p2[absolute_p2+1])


        ## the row that is composed of prob. p2, given all 0:2N p1 values
        ##        return(prob_at_p2_matrix[absolute_p2+1,])
        return(prob_at_p2_matrix)
    }

    ##    prob <- sapply(X = absolute_p1, FUN = function(p1) single_p2_given_single_p1(absolute_p2, p1, s, h, N, generations_in_between))
    prob <- single_p2_given_all_p1(absolute_p2 = absolute_p2, s = s, h = h, N = N, generations_in_between = generations_in_between)
    return(prob)
}

p2_given_p1 <- function(absolute_p2, absolute_p1, s, h, N, generations_in_between, list_of_transition_matrix, para) {
    if (! (length(list_of_transition_matrix) == nrow(para))) {
        stop("Incorrect! These 2 lengths have to be equal!")
    }
    if (length(absolute_p2) != 2*N + 1) {
        ##        stop("length of absolute_p2 should be 1.")
        stop("length of absolute_p2 should be 2*N+1.")
    }
    if (length(absolute_p1) != 2*N + 1) {
        stop("Have to iterate all possible p1, from 0 to 2*N")
    }
    if (N != 500) {
        stop("N should be 500 now!")
    }
    row.is.a.match <- apply(para, 1, identical, c(s, h))
    if (sum(row.is.a.match) != 1) {
        print(s)
        print(h)
        print(sum(row.is.a.match))
        print(head(para))
        stop("each s&h combination should only appear in para once and only once!")
    }
    row_index <- which(row.is.a.match)
    prob <- list_of_transition_matrix[[row_index]][[generations_in_between]]
    if (is.null(prob)) {
        stop("mistake! the list_of_transition_matrix does not have this generations_in_between, thus an error.")
    }
    return(prob)
}

## P(p1 | reads_1, depth_1)
##p1_given_depth_1_reads_1 <- function(absolute_p1, depth_p1, reads_p1, N) {
p1_given_depth_1_reads_1 <- function(absolute_p1, depth_p1, reads_p1, N, sample_method) {
    ## suppose uniform prior P(p1)
    if (length(absolute_p1) != 2*N + 1) {
        stop("Have to iterate all possible p1, from 0 to 2*N")
    }
    if (length(depth_p1) != 1 | length(reads_p1) != 1) {
        stop("depth_p1 and reads_p1 have to be size 1.")
    }

    if (sample_method == "binomial") {
        p1 <- absolute_p1 / (2*N)
        prob <- dbinom(x = reads_p1, size = depth_p1, prob = p1, log = FALSE)
    } else if (sample_method == "hypergeometric") {
        p1 <- absolute_p1
        prob <- dhyper(x = reads_p1, m = p1, n = 2*N - p1, k = depth_p1)
    } else {
        stop("you need to specify the correct sampling method")
    }

    ## should I use prob/sum(prob) or not???
    ## YES! because I am pretending I have used Bayes' Theorem w/ uniform priors!!!
    ## so, the return should be P(p1|depth_p1, reads_p1), which is a probability distribution of p1; sum_p1(P(p1|reads_1..)) == 1
    return(prob / sum(prob))
}

## for calculating single absolute_p2 value

p2_given_depth_p1_reads_p1_singular <- function(absolute_p2, depth_p1, reads_p1, s, h, N, generations_in_between) {
    stop("no use")
    absolute_p1 <- 0:(2*N)

    ## probability of single p2 value
    if (length(absolute_p2) == 1) {
        temp_prob <- sum(p2_given_p1(absolute_p2, absolute_p1, s, h, N, generations_in_between) * p1_given_depth_1_reads_1(absolute_p1, depth_p1, reads_p1, N, sample_method = sample_method))
        return(temp_prob)
    } else {
        stop("this step is only used to calculate the else if clause; as this is not probability; they won't sum to 1. Thus, this clause is only recruited to calcaulate all of them and divide by sum(prob).")
    }
}




## P(p2 | reads_1, depth_1)
## P(p2 | Model, depth_p1, reads_p1) = sum_all_p1{P(p2 | Model, p1) * P(p1 | Model, depth_p1, reads_p1)}
p2_given_depth_p1_reads_p1 <- function(absolute_p2, depth_p1, reads_p1, s, h, N, generations_in_between) {
    stop("no use")

    ## probability of single p2 value
    if (length(absolute_p2) == 1) {
        stop("should use the singular version.")
    } else if (length(absolute_p2) == 2*N + 1) {
        ## could use for (j in 0:(2*N)), but I think maybe use apply is better for parallel computing
        ## or maybe not use parallel here??
        ##prob <- parApply(cl, X = matrix(absolute_p2, nrow = length(absolute_p2)), MARGIN = 1, FUN = p2_given_depth_p1_reads_p1_singular, depth_p1 = depth_p1, reads_p1 = reads_p1, s = s, h = h, N = N)
        prob <- apply(X = matrix(absolute_p2, nrow = length(absolute_p2)), MARGIN = 1, FUN = p2_given_depth_p1_reads_p1_singular, depth_p1 = depth_p1, reads_p1 = reads_p1, s = s, h = h, N = N, generations_in_between = generations_in_between)

        return(prob)

        ## no standarizing constant is necessary!! Since this sum_all_B P(A, B) = P(A) is true probability
        ##        return(prob / sum(prob))
    } else {
        stop("the return value for length(absolute_p2) != 2*N + 1 will not be probability. since they need to be divided by sum(prob).")
    }

}

p2_given_depth_p1_reads_p1_matrix <- function(absolute_p2, depth_p1, reads_p1, s, h, N, generations_in_between, sample_method, list_of_transition_matrix, para) {

    ## probability of single p2 value
    if (length(absolute_p2) == 1) {
        stop("should use the singular version.")
    } else if (length(absolute_p2) == 2*N + 1) {

        absolute_p1 <- 0:(2*N)
        prob_p1 <- p1_given_depth_1_reads_1(absolute_p1, depth_p1, reads_p1, N, sample_method = sample_method)

        ## prob of p1 is in each column
        ##    prob_p1_given_reads_matrix <- matrix(rep(prob_p1, time = 2*N+1), nrow = 2*N+1, byrow = F)

        ## this is already a matrix, with each row representing P(single p2 | given 0:2N p1);
        ## thus, sum(prob_p2_given_p1_matrix[1,] * prob_p1_given_reads_matrix[,1]) is P(single p2 | given reads_p1, depth_p1)
        prob_p2_given_p1_matrix <- p2_given_p1(0:(2*N), 0:(2*N), s = s, h = h, N = N, generations_in_between = generations_in_between, list_of_transition_matrix, para)

        ## P(p2 | reads_p1);
        prob_matrix <- prob_p2_given_p1_matrix %*% prob_p1

        ##        prob <- apply(X = matrix(absolute_p2, nrow = length(absolute_p2)), MARGIN = 1, FUN = p2_given_depth_p1_reads_p1_singular, depth_p1 = depth_p1, reads_p1 = reads_p1, s = s, h = h, N = N, generations_in_between = generations_in_between)

        return(as.numeric(prob_matrix))

        ## no standarizing constant is necessary!! Since this sum_all_B P(A, B) = P(A) is true probability
        ##        return(prob / sum(prob))
    } else {
        stop("the return value for length(absolute_p2) != 2*N + 1 will not be probability. since they need to be divided by sum(prob).")
    }

}


## P(reads_p2 | p2, depth_p2)
##reads_p2_given_p2 <- function(reads_p2, depth_p2, absolute_p2, N) {
reads_p2_given_p2 <- function(reads_p2, depth_p2, absolute_p2, N, sample_method) {
    if (length(reads_p2) != 1) {
        stop("the length of reads_p2 should be 1")
    }
    if (length(absolute_p2) != 2*N + 1) {
        stop("WRONG here")
    }

    if (sample_method == "binomial") {
        p2 <- absolute_p2 / (2*N)
        ## no standarizing constant is necessary
        prob <- dbinom(x = reads_p2, size = depth_p2, prob = p2)
    } else if (sample_method == "hypergeometric") {
        p2 <- absolute_p2
        prob <- dhyper(x = reads_p2, m = p2, n = 2*N - p2, k = depth_p2)
    } else {
        stop("you need to specify the correct sampling method")
    }

    ## I think should not use prob/sum(prob) here, as the returned is P(reads_p2|p2), not P(p2|reads_p2)
    ## this is different from p1_given_depth_1_reads_1(), since that returns a "posterior" after Bayes' Theorm w/ uniform priors
    ## but this is not a "posterior", no uniform prior is used here
    return(prob)
}

## singular version
reads_p2_given_reads_p1_singular <- function(reads_p2, depth_p2, reads_p1, depth_p1, s, h, N, generations_in_between, sample_method, list_of_transition_matrix, para) {
    if (length(reads_p2) != 1) {
        stop("this is only the singular version, to calculate the values for sum(prob), the standarizing constant; the output is not probability.")
    }
    absolute_p2 <- 0:(2*N)

    if (((reads_p2 == depth_p2) & (reads_p1 == depth_p1)) | ((reads_p2 == 0) & (reads_p1 == 0))) {
        ## should be very IMPORTANT
        ## 2 consecutive reads = 0 or reads = depth; return 1 because these actually gives NO effective information
        return(1)
    } else {
        ##    temp_prob <- sum(reads_p2_given_p2(reads_p2, depth_p2, absolute_p2, N) * p2_given_depth_p1_reads_p1(absolute_p2, depth_p1, reads_p1, s, h, N, generations_in_between))
        temp_prob <- sum(reads_p2_given_p2(reads_p2, depth_p2, absolute_p2, N, sample_method = sample_method) * p2_given_depth_p1_reads_p1_matrix(absolute_p2, depth_p1, reads_p1, s, h, N, generations_in_between, sample_method = sample_method, list_of_transition_matrix, para = para))

        return(temp_prob)
    }
}

## P(reads_p2 | reads_p1, depth_p1, depth_p2, s, h)
## reads_p2_given_reads_p1 <- function(reads_p2, depth_p2, reads_p1, depth_p1, s, h, N) {

## reads_p2_given_reads_p1 <- function(depth_p2, reads_p2, depth_p1, reads_p1, generations_in_between, s, h, N) {
reads_p2_given_reads_p1 <- function(data, s, h, N, sample_method, list_of_transition_matrix, para) {
    if (length(data) != 5) {
        stop("not enough data in reads_p2_given_reads_p1().")
    }
    depth_p2 <- data[1]; reads_p2 <- data[2]; depth_p1 <- data[3]; reads_p1 <- data[4]
    generations_in_between <- data[5]
##    print(N)
    ## just for calculating the standarizing constant
    ##    temp_reads_p2 <- 0:(depth_p2)

    ##    if(length(reads_p2) == 1) {
    if(max(length(reads_p2), length(depth_p2), length(reads_p1), length(depth_p1)) == 1) {
        ## seems no standarizing constant is needed for this
        ##        total_prob <- apply(X = matrix(temp_reads_p2, nrow = length(temp_reads_p2)), MARGIN = 1, FUN = reads_p2_given_reads_p1_singular, depth_p2 = depth_p2, reads_p1 = reads_p1, depth_p1 = depth_p1, s = s, h = h, N = N)

        prob <- reads_p2_given_reads_p1_singular(reads_p2, depth_p2, reads_p1, depth_p1, s, h, N, generations_in_between, sample_method = sample_method, list_of_transition_matrix, para = para)

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

call_log_2d_singular <- function(seq_data_vector, N, sample_method) {
    ## get seq_data_vector from calculate_rep.R, as a vector of (position, generation, reads, depth, ...)
    ## only for a single experiment; multiple experiments are combined in calculate_rep.R
    if (is.numeric(seq_data_vector) == TRUE) {
        reads_1_index <- seq.int(from = 3, to = length(seq_data_vector), by = 3)
        depth_1_index <- seq.int(from = 4, to = length(seq_data_vector), by = 3)
        generation_1_index <- seq.int(from = 2, to = length(seq_data_vector), by = 3)

        seq_data <- data.frame(depth_1 = seq_data_vector[depth_1_index], reads_1 = seq_data_vector[reads_1_index], generations_1 = seq_data_vector[generation_1_index])

        seq_data$depth_2 <- c(seq_data$depth[2:length(seq_data$depth)], -Inf)
        seq_data$reads_2 <- c(seq_data$reads_1[2:length(seq_data$reads_1)], -Inf)
        seq_data$generations_2 <- c(seq_data$generations_1[2:length(seq_data$generations_1)], -Inf)
        seq_data$generations_in_between <- seq_data$generations_2 - seq_data$generations_1

        seq_data <- seq_data[, c("depth_2", "reads_2", "depth_1", "reads_1", "generations_in_between")]


        seq_data <- as.matrix(seq_data)
        ## since last line has no depth_2, reads_2; truncate it
        seq_data <- seq_data[1:(dim(seq_data)[1] - 1), ]

        ## needs to apply par[] matrix
        ##        log_likelihood_per_s_h(par, seq_data)
        results_in_vec <- numeric(0)
        for (t in 1:length(para_list)) {
            para <- para_list[[t]]
            load(file_name_vec[t])
            ##            print(length(list_of_transition_matrix))
##            print(list_of_transition_matrix[[1]][[10]][4,])
            ##            results_in_vec <- apply(para, 1, function(parameters) log_likelihood_per_s_h(parameters, seq_data, N, sample_method))
            results_in_vec <- c(results_in_vec, apply(para, 1, function(parameters) log_likelihood_per_s_h(parameters, seq_data, N, sample_method, list_of_transition_matrix, para = para)))
        }


        ## maybe I should return 95% C.I., mode etc instead of raw results, just for the saving memory

        ## this print seems useless in parallel computing
##        print("Finished");
##        return(summarize_results(cbind(para, results_in_vec)));
        return(results_in_vec)
    } else {
        print(class(seq_data_vector))
        stop("only used for single-seq data");
    }

}

################################################################

call_log_2d <- function(seq_data_frame, N, sample_method) {
    ## get a data.frame of values, call call_log_2d_singular(), and sum results
    ##    if (is.list(seq_data_frame) == FALSE) {
    if (is.data.frame(seq_data_frame) == FALSE) {
        stop("Only accept a data.frame for input.\n");
    } else if (is.data.frame(seq_data_frame) == TRUE) {
        if (dim(seq_data_frame)[1] >= 2) {
            result <- apply(seq_data_frame, 1, FUN = function(seq_data_vector) call_log_2d_singular(seq_data_vector, N, sample_method));
            sum_result <- rowSums(result);
            ##        return(sum_result);

            ## this will return the prob matrix
##                    return(cbind(para, sum_result));

            return(summarize_results(cbind(big_para, sum_result)));
        } else if (dim(seq_data_frame)[1] == 1) {
            result <- call_log_2d_singular(as.numeric(seq_data_frame), N, sample_method);
            ## return(result)
            return(summarize_results(cbind(big_para, result)));
        } else {
            stop("The dim() is not correct.\n");
        }
    } else {
        stop("Wrong data input time")
    }
}


################################################################

log_likelihood_per_s_h <- function(parameters, seq_data, N, sample_method, list_of_transition_matrix, para) {
    if (is.matrix(seq_data) == TRUE) {
        ## no replications
        ## seq_data = seq_data_matrix, is a matrix
        s <- parameters[1]; h <- parameters[2]

        ##    sum(apply(X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = s, h = 1, N = 100))
        ##        sum(parApply(cl, X = seq_data_matrix, MARGIN = 1, FUN = reads_p2_given_reads_p1, s = s, h = h, N = 100))
        ## N is hard coded here; and seems only hard coded here
        ##        sum(parApply(cl, X = seq_data, MARGIN = 1, FUN = function(data) reads_p2_given_reads_p1(data, s = s, h = h, N = 100)))
        ## just a test

        ##        sum(apply(X = seq_data, MARGIN = 1, FUN = function(data) reads_p2_given_reads_p1(data, s = s, h = h, N = 100)))
                sum(apply(X = seq_data, MARGIN = 1, FUN = function(data) reads_p2_given_reads_p1(data, s = s, h = h, N = N, sample_method = sample_method, list_of_transition_matrix, para = para)))

    } else if (is.list(seq_data) == TRUE) {
        stop("should not be here as this time seq_data is just a matrix")
        ## with replicate sequencing data; now, restricted to be sequenced at the same generations
        ## seq_data = list(seq_data_matrix), is a list of seq_data_matrix
        sum(sapply(X = seq_data, FUN = function(seq_data_matrix) log_likelihood_per_s_h(parameters, seq_data_matrix)))
    } else {
        stop("Something wrong here. seq_data should either be a matrix or a list.")
    }
}



################################################################
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




summarize_results <- function(results) {
    ## define the C.I. here
    ##    conf_interval <- c(0.025, 0.975)

    ## if h is no longer a single value

##    results_s_only <- matrix(, nrow = length(unique(para[,1])), ncol = 3)
##    if (length(unique(results[,2])) != 1) {
        ## then, need to transform log_prob to prob = exp(log_prob)
        ## and p(s) = sum_all_h(p(s, h))

        ## transform to data.frame for dplyr
        results_df <- as.data.frame(results)
        colnames(results_df) <- c("s", "h", "log_prob")
        results_df$prob <- exp(results_df$log_prob)

        ## IMPORTANT; maybe not necessary?
        results_df$prob <- results_df$prob / sum(results_df$prob)

        ## p(s) = sum_h(p(s,h))
        results_df_s <- dplyr:::group_by(results_df, s)
        results_df_s <- dplyr:::summarise(results_df_s, prob = sum(prob))
        ## now, the log_prob is in 3rd column again
        results_df_s$log_prob <- log(results_df_s$prob)
        ## now, the results is a matrix, with [,1] = s, [,2] = prob, [,3] = log_prob
##        results_s_only <- as.matrix(results_df_s)
##        colnames(results_s_only) <- NULL

        ## p(h) = sum_s(p(s,h))
        results_df_h <- dplyr:::group_by(results_df, h)
        results_df_h <- dplyr:::summarise(results_df_h, prob = sum(prob))
        ## now, the log_prob is in 3rd column again
        results_df_h$log_prob <- log(results_df_h$prob)
        ## now, the results is a matrix, with [,1] = s, [,2] = prob, [,3] = log_prob
##        results_h_only <- as.matrix(results_df_h)
##        colnames(results_h_only) <- NULL

##    }


    conf_interval <- c(0.005, 0.995)

    ##    total_results <- data.frame(lower_bound = rep(0, 1), upper_bound = rep(0, 1), mean = rep(0, 1), mode = rep(0, 1), CI = rep(0, 1), BF = rep(0, 1))
    total_results <- data.frame(lower_bound_of_s = rep(0, 1), upper_bound_of_s = rep(0, 1), mean_of_s = rep(0, 1), mode_of_s = rep(0, 1), mode_of_h = rep(0, 1), grand_mode_s = rep(0, 1), grand_mode_h = rep(0, 1), CI_of_s = rep(0, 1), BF = rep(0, 1))

    trans <- data.frame(prob = results[,3], s = results[,1], h = results[,2])
    ##    p <- exp(trans$prob) / sum(exp(trans$prob))
    ## p <- data.frame(s = trans$s, prob = p)
    trans$prob <- exp(trans$prob) / sum(exp(trans$prob))


    ## IMPORTANT!!!! to prevent multiple max values; pick the smalled ones
    ##    mode <- p$s[p$prob == max(p$prob)][1]
    ## the mode of s, and mode of h, seprated; p(s), p(h)
    total_results[1,]$grand_mode_s <- trans[trans$prob == max(trans$prob),]$s[1]
    total_results[1,]$grand_mode_h <- trans[trans$prob == max(trans$prob),]$h[1]


    ## now, turn to s_only and h_only

    ## calculte mean here
##    total_results[1,]$mean <- sum(p$prob * p$s)
    total_results[1,]$mode_of_h <- results_df_h[results_df_h$log_prob == max(results_df_h$log_prob),]$h[1]
    total_results[1,]$mode_of_s <- results_df_s[results_df_s$log_prob == max(results_df_s$log_prob),]$s[1]


    ## NEW: calculate Bayes' Factor; just a test
    ## I think BF = p() / p(s = 0), instead of BF = p() / p(s = 0, h = 0), since we are rejecting s = 0, not (s = 0 & h = 0)

    ## since p == results_df_s, this 2 equations equal
    ##    total_results[1,]$BF <- mean(p$prob) / p[p$s == 0,]$prob
    total_results[1,]$BF <- mean(results_df_s$prob) / results_df_s[results_df_s$s == 0,]$prob


    p <- results_df_s
    p$prob <- p$prob / sum(p$prob)
    total_results[1,]$mean_of_s <- sum(p$prob * p$s)

    boundary_by_cdf <- c(0, 0, 0)

################################################################
                                        #cdf_method()
    ## what's inside cdf_method()

    ## !!!!!!!!!!!
    ## COULD be calculated simply using quantile() function!!! or maybe not


    discrete_cdf <- numeric(0)
    discrete_cdf[1] <- p$prob[1]
    for (k in (2:length(p$prob))) {
                                        #        print(t)
        discrete_cdf[k] = discrete_cdf[k-1] + p$prob[k]
    }

    ## generate the cdf
    cdf <- approxfun(x = p$s, y = discrete_cdf)
    ## find cdf(left) = 0.025 & cdf(right) = 0.975
    found_left = FALSE
    found_right = FALSE
    step = 0.005
    for (s in seq.int(min(p$s), max(p$s), by = step)) {
                                        #    print(cdf(s))

        ## just to make sure the lower_bound exist; since cdf(min(p$s)) could be > conf_interval[1]
        ## but cdf(max(p$s)) = 1 always true
        if (cdf(min(p$s)) >= conf_interval[1] & found_left == FALSE) {
            ## even cdf(min(p$s)) is larger than conf_interval[1]
            ## then, set the lower_bound to be min(p$s)
            boundary_by_cdf[1] <- min(p$s)
            found_left = TRUE
        }

        if (cdf(s) >= conf_interval[1] & cdf(s-step) <= conf_interval[1] & found_left == FALSE) {
            ## left = F(2.5%)
            ## because cdf(s) already counds pmf(s), thus, need to add step
            boundary_by_cdf[1] <- s + step
            found_left = TRUE
        }
                                        #    if (all.equal(cdf(s), 0.975, tol = 0.025) == TRUE & found_right == FALSE) {
        if (cdf(s) >= conf_interval[2] & cdf(s-step) <= conf_interval[2] & found_right == FALSE) {
            boundary_by_cdf[2] <- s + step
            found_right = TRUE

            ## calculate the true percent
            ## to make sure it is within cdf function's range, otherwise NA
            if (boundary_by_cdf[2] >= max(p$s)) {
                boundary_by_cdf[3] <- 1 - cdf(boundary_by_cdf[1])
            } else {
                boundary_by_cdf[3] <- cdf(boundary_by_cdf[2]) - cdf(boundary_by_cdf[1])
            }

            break
        }
    }

    total_results[1,]$lower_bound_of_s <- boundary_by_cdf[1]
    total_results[1,]$upper_bound_of_s <- boundary_by_cdf[2]
    total_results[1,]$CI_of_s <- boundary_by_cdf[3]

    ##    return(total_results)
    return(as.matrix(total_results))
}
