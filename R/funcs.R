

# this is the main function, executing both new_pi and log_likelihood
gl_to_gp <- function(matrix_A, pivec = c(1 / 3, 1 / 3, 1 / 3), tol = 10^-5, intermax = 100) {
  if (is.null(matrix_A)) {
    return(1)
  }
  err <- 100000
  iter <- 0
  curr_pivec <- pivec
  err_buff <- abs(log_like(matrix_A))
  while (err > tol && iter < intermax) {
    curr_pivec <- new_pi(matrix_A, curr_pivec)
    err <- err_buff - abs(log_like(matrix_A, curr_pivec))
    err_buff <- abs(log_like(matrix_A, curr_pivec))
    iter <- iter + 1
  }
  curr_pivec
}


# this function calculates the pi (new) values,
# the function starts with a for loop for each k (1-3)
# it then initializes a variable called local_pi which is just going to hold
# each new pik we find through the other iterations. this local_pi is then added
# at the end of each iteration to a vector called pik_new witch at the end of the program
# is supposed to contain new pik values for each k choice, then the pik_new vector is
# returned.
new_pi <- function(A, piks = c(1 / 3, 1 / 3, 1 / 3)) {
  n <- nrow(A)
  pik_new <- c()
  for (k in 1:3) { # this is really 0 to 2
    local_pi <- 0
    for (i in 1:n) {
      bottom <- 0
      for (j in 1:3) {
        bottom <- bottom + (piks[j] * A[i, j])
      }
      local_pi <- local_pi + ((piks[k] * A[i, k]) / bottom)
    }
    pik_new <- append(pik_new, (local_pi * 1 / n))
  }
  pik_new
}

# this function calculates the log likelihood. This function iterates over all n
# and creates a variable called the inner_sum witch is the sum of pik * aik for
# each individual (1 - n). The function then adds the log of the inner sum to the
# total likelihood variable and outputs a double.
log_like <- function(A, piks = c(1 / 3, 1 / 3, 1 / 3)) {
  n <- nrow(A)
  total_like <- 0
  for (i in 1:n) {
    inner_sum <- 0
    for (k in 1:3) {
      inner_sum <- inner_sum + (piks[k] * A[i, k])
    }
    total_like <- total_like + log(inner_sum)
  }
  total_like
}
