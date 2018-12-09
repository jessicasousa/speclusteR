apply_gaussian_similarity <- function(xi, xj, alpha = 1){
  exp(- alpha * norm(as.matrix(xi - xj), type = "F"))
}

create_similarity_matrix <- function(m) {
  n <- nrow(m)
  s <- matrix(NA, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      s[i, j] <- apply_gaussian_similarity(m[i, ], m[j, ])
    }
  }
  s
}
