construct_similarity_graph <- function(A, neighboors = 3){
  n <- nrow(A)
  if (neighboors >= n) {
    S <- A
  } else {
    S <- matrix(0, ncol = n, nrow = n)
    for(i in seq_len(n)) {
      best_similarities <- sort(A[i, ], decreasing = TRUE)[1:neighboors]
      indices <- which(A[i, ] %in% best_similarities)
      S[i, indices] <- A[i, indices]
      S[indices, i] <- A[i, indices]
    }
  }
  S
}

#Unnormalized spectral clustering
apply_spectral_clustering <- function(A, k){
  n <- ncol(A)
  #Construct a similarity graph
  A <- construct_similarity_graph(A)
  #Compute the unnormalized Laplacian L.
  degrees <- rowSums(A)
  L <- diag(degrees) - A
  #Compute the first k eigenvectors u1 , ... , uk of L.
  eigenvl <- eigen(L, symmetric = TRUE)
  eigenvectors <- eigenvl$vectors
  eigenvectors <- eigenvectors[ , (n - k + 1):n]
  #Cluster the points
  kmeans(eigenvectors, centers = k)
}
