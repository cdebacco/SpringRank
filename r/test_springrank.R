# This just runs the test script, same as the Python version, but a
# bit more barebones.

library(readr)
library(igraph)
library(Matrix)

source("springrank.R")

df <- read_delim("../data/US_CS_adjacency.dat",
                 delim = " ",
                 col_names = c("i", "j", "weight"))

G <- graph.data.frame(df)
A <- as.matrix(as_adjacency_matrix(G, attr = "weight"))
A <- Matrix(A) # sparsify A

alpha <- 0.0
l0 <- 1.0
l1 <- 1.0

r <- spring_rank(A, alpha, l0, l1)
print(round(rev(r[order(r)]), 3))