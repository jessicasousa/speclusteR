#' speclusteR: Spectral Clustering Algorithms in R
#'
#' Trabalho para disciplina Álgebra Linear para Data Science Machine Learning - 2018.2,
#' implementação dos três principais Spectral Clustering Algorithms.
#'
#' @docType package
#' @name speclusteR
#' @author Jessica Cardoso
#' @references
#' \itemize{
#' \item https://bit.ly/2B6qFnn
#' \item https://bit.ly/2zRckeK
#' \item https://bit.ly/2G5UmL0
#' }
NULL


#' apply_squared_exponential
#'
#' Métrica de similaridade utilizada para construir o grafo de similaridade.
#' A métrica utilizada é dada pela seguinte equação:
#' \deqn{S(i, j) = ||x_i - x_j||^2 / 2 x sigma1 x sigma 2}
#'
#' @return A função retorna um numérico correspondendo a distância entre duas linhas da matriz.
#' @param i Corresponde ao número que representa a i-ésima linha matriz
#' @param j Corresponde ao número que representa a j-ésima linha matriz
#' @param A matriz numérica que terá suas linhas comparadas
#' @param sig corresponde a um real
#' @param sig2 corresponde a um número real
#' @export
#' @examples
#'mat <- matrix(sample(10,9), nrow = 3)
#'print(mat)
#'res <- apply_squared_exponential(1, 3, mat)
#'print(res)
apply_squared_exponential <- function(i, j, A, sig = 0.8, sig2 = 1){
  xi <- A[i, ]
  xj <- A[j, ]
  norm_squared <- function(x, y) norm(as.matrix(x - y))^2
  exp( - norm_squared(xi, xj) / (2 * sig * sig2) )
}

#' build_similarity_graph
#'
#' Essa função realiza a transformação da matriz em um grafo de similaridade.
#' Usando a função apply_squared_exponential para obter a similaridade.
#'
#' @return A função retorna uma matriz correspondendo a grafo de similaridade.
#' @param A matriz numérica quadrada
#' @param sig1 corresponde a um real utilizado na função de similaridae
#' @param sig2 corresponde a um real utilizado na função de similaridae
#' @export
#' @examples
#' A <- matrix(sample(100), nrow = 10)
#' print(A)
#' res <- build_similarity_graph(A)
#' print(res)
#'
build_similarity_graph <- function(A, sig1 = 0.8, sig2 = 1){
  n <- nrow(A)
  stopifnot(n == ncol(A))
  indices <- seq_len(n)
  vectFunc <- Vectorize(apply_squared_exponential, vectorize.args = c('i','j'))
  outer(indices, indices, FUN = vectFunc, A, sig1, sig2)
}


#' create_graph_laplacian
#'
#' Essa função realiza o cálculo do grafo Laplaciano, o qual pode ser não normalizado, ou
#' normalizado segundo Shi e Malik (2000) ou Ng, Jordan, e Weiss (2002).
#'
#' @return A função retorna uma matriz correspondendo a grafo de laplaciano.
#' @param similarity_matrix corresponde a matriz de similaridade calculada.
#' @param type Corresponde a um inteiro indicando qual tipo de normalização deseja
#' executar.
#' \itemize{
#' \item 1: Não normalizado.
#' \item 2: Normalizado segundo Shi e Malik (2000).
#' \item 3: Normalizado segundo Ng, Jordan, e Weiss (2002).
#' }
#' @export
#' @examples
#' set.seed(2018)
#' x <- sample(100, 10)
#' A <- as.matrix(dist(x))
#' S <- build_similarity_graph(A)
#' L <- create_graph_laplacian(S)
#'
create_graph_laplacian <- function(similarity_matrix, type = 1){
  #Calcula a matriz de graus
  degrees <- rowSums(similarity_matrix) # graus de vertice
  degrees[degrees == 0] <- .Machine$double.eps
  D <- diag(degrees)

  #computa o Laplaciano não-normalizado
  L <- D - similarity_matrix
  #No tipo 1, retornar versao não-normalizada
  if(type == 1){
    return(L)
  }
  #Variações de Laplaciano normalizado
  if(type == 2){
    #Relacionado ao Random Walk
    Di <- diag(1 / degrees)
    NL <- Di %*% similarity_matrix
  }
  else if(type == 3){
    #Matriz simétrica
    Di <- diag(1 / sqrt(degrees))
    NL <- Di %*% similarity_matrix %*% Di
  }
  return(NL)
}

#' apply_spectral_clustering
#'
#' Essa função realiza o cálculo do spectral clustering de acordo com o tipo de normalização.
#'
#' @return Retorna uma lista contendo o resultado do kmeans sobre o
#' @param A corresponde a matriz de similaridade calculada.
#' @param k corresponde ao número de grupos a ser considerado
#' @param sig1 corresponde a um real
#' @param sig2 corresponde a um real
#' @param type Corresponde a um inteiro indicando qual tipo de normalização deseja
#' executar.
#' \itemize{
#' \item 1: Não normalizado.
#' \item 2: Normalizado segundo Shi e Malik (2000).
#' \item 3: Normalizado segundo Ng, Jordan, e Weiss (2002).
#' }
#' @export
#' @examples
#' \dontrun{
#' set.seed(2018)
#' n <- 150
#' r <- rnorm(n, 5, .25)
#' theta <- runif(n, 0, 2 * pi)
#' c1 <- data.frame(x = rnorm(n), y = rnorm(n))
#' c2 <- data.frame(x = r * cos(theta), y = r * sin(theta))
#' points1 <- rbind(c1, c2)
#' A <- as.matrix(dist(points1))
#' S <- build_similarity_graph(A)
#' L <- create_graph_laplacian(S)
#'}
apply_spectral_clustering <- function(A, k, sig1 = 0.8, sig2 = 1, type = 1){

  if(!(type %in% 1:3)){
    warning("Tipo deve ser um dos seguintes: 1, 2 or 3. Tipo continua com valor padrão 1")
    type = 1
  }
  similarity_graph <- build_similarity_graph(A, sig1, sig2)
  L <- create_graph_laplacian(similarity_graph, type)
  ei <- eigen(L, symmetric = TRUE)
  U <- ei$vectors[ , 1:k]

  if(type == 3){
    U <- t(apply(U, 1, function(x) x / sqrt(sum(x^2))))
  }

  km <- stats::kmeans(U, centers = k)
  return(km)
}
