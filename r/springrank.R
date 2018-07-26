library(igraph)
library(Matrix)
library(Rlinsolve)

spring_rank <- function(A, alpha = 0, l0 = 1.0, l1 = 1.0, shift = TRUE) {
  #' Core function for calculating SpringRank.
  #' Default parameters follow stanadard model.
  #'
  #' @param A The adjacency matrix of the graph. This should be a sparse dgCMatrix;
  #' if it isn't it will be coerced to a dgCMatrix.
  #' @param alpha Controls the impact regularization term.
  #' @param l0 Regularization spring's rest length.
  #' @param l1 Interaction springs' rest length.
  #' @param shift (Optional, default TRUE) normalize such that the lowest-ranked
  #'  node has a SpringRank value of zero.
  #'
  #' @return A vector of SpringRank scores for each node. Sort or order the
  #' vector for ordinal rankings of each node.

  if (class(A) == "matrix") {
    # coerce dense matrix to sparse matrice so the user doesn't have to.
    A <- as(Matrix(A, sparse = TRUE, doDiag = FALSE), "dgCMatrix")
  } else {
    # confirm it's the right kind of sparse matrix. might throw an error
    # if it's one of the less common sparse matrix types.
    A <- as(A, "dgCMatrix")
  }

  N <- dim(A)[1]
  k_in <- colSums(A)
  k_out <- rowSums(A)
  One <- as(Matrix(rep(1, N)), "dgCMatrix")
  C <- A + t(A)

  D1 <- as(Matrix(0, ncol = N, nrow = N), "dgCMatrix")
  D2 <- as(Matrix(0, ncol = N, nrow = N), "dgCMatrix")

  for (i in 1:N) {
    D1[i, i] <- k_out[i] + k_in[i]
    D2[i, i] <- l1 * (k_out[i] - k_in[i])
  }

  if (alpha != 0) {
    print("assuming invertible matrix")

    B = One * l0 + D2 %*% One
    A_ = alpha * diag(nrow = N, ncol = N) + D1 - C
  } else {
    print("fixing a rank degree of freedom")

    C <- C +
      matrix(rep(A[N,] , times = N),
             ncol = N,
             nrow = N,
             byrow = T) +
      matrix(rep(A[, N], times = N),
             ncol = N,
             nrow = N,
             byrow = T)

    D3 <- as(Matrix(0, ncol = N, nrow = N), "dgCMatrix")
    for (i in 1:N) {
      D3[i, i] <- l1 * (k_out[N] - k_in[N])
    }

    B <- D2 %*% One + D3 %*% One
    A_ <- D1 - C
  }

  rank <- lsolve.bicgstab(A_, B, verbose = F)
  rank <- rank$x

  if (shift) {
    rank <- rank - min(rank)
  }

  # coerce matrix to vector, so we can use names()
  rank <- rank[,1]
  names(rank) <- colnames(A)
  return(rank)
}


spring_rank_network <- function(N, beta, alpha, K, l0 = 0.5, l1 = 1.0) {
  #' Generative model for SpringRank. Builds a weighted graph with self-loops.
  #' The model first generates node scores using a normal distribution, then
  #' generates edges (and edge weights) from a Poisson distribution with mean
  #' equal to the energy of the system.
  #'
  #' @param N The number of nodes.
  #' @param beta Inverse temperature, a noise parameter.
  #' @param alpha Variance of the Normal prior.
  #' @param K The average degree of the network.
  #' @param l0 Prior spring's rest length.
  #' @param l1 Interaction spring's rest length.
  #'
  #' @return A directed graph (potentially weighted, potentially containing
  #' self-loops)

  scores <- rnorm(N, l0, 1/sqrt(alpha*beta))
  Z <- 0
  for (i in 1:N) {
    for (j in 1:N) {
      Z <- Z + exp(-0.5 * beta * (scores[i] - scores[j] - l1)^2)
    }
  }
  C <- (K*N)/Z

  # for loops are slow in R so make a matrix of element-wise subtractions,
  # each element i, j being scores[i]-scores[j]
  # basically, trading off increased memory usage (dense matrix) for speed.
  scores_mat <- matrix(1, length(scores), 1) %*% t(scores)
  scores_mat <- scores_mat - scores - l1
  H <- .5 * scores_mat^2
  lambda <- C * exp(-1*beta*H)
  A <- rpois(length(scores)^2, lambda) %>% matrix(nrow = dim(lambda)[1])

  return(graph_from_adjacency_matrix(A, mode = "directed", weighted = "weight"))
}
