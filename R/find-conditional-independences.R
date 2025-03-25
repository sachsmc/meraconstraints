#' D separation
#'
#' Is node X d-separated from Y by Z in the graph?
#'
#' @param graph An igraph object representing a finite directed graph with named vertices
#' @param X Name of a node in graph
#' @param Y Name of a node in graph
#' @param Z Name of a node in graph
#'
#' @returns Logical, TRUE if is node X d-separated from Y by Z in the graph
#' @export
#' @references Dethlefsen C, Højsgaard S (2005). “A Common Platform for Graphical Models in R: The gRbase Package.” Journal of Statistical Software, 14(17), 1–12. https://www.jstatsoft.org/v14/i17/.
#' @examples
#'
#' dag0 <- graph_from_literal(a -+ b:c:e, b -+ c, c -+ d, e -+ d)
#' d_separates("c", "e", "a", dag0)
#'
d_separates <- function(X, Y, Z, graph) {

  ag <- gRbase::ancestralGraph(c(X, Y, Z), graph)
  gRbase::separates(X, Y, Z, gRbase::moralize(ag))

}


#' D separation
#'
#' Is set of nodes X d-separated from Y by Z in the graph?
#'
#' @param graph An igraph object representing a finite directed graph with named vertices
#' @param X Vector of names of a node in graph
#' @param Y Vector of names of a node in graph
#' @param Z Vector of names of a node in graph
#'
#' @returns Logical, TRUE if is node X d-separated from Y by Z in the graph
#' @export
#' @references Dethlefsen C, Højsgaard S (2005). “A Common Platform for Graphical Models in R: The gRbase Package.” Journal of Statistical Software, 14(17), 1–12. https://www.jstatsoft.org/v14/i17/.
#' @examples
#'
#' dag0 <- graph_from_literal(a -+ b:c:e, b -+ c, c -+ d, e -+ d)
#' d_separate_sets(c("d", "c"),  c("a", "b"), "e", dag0)
#'
d_separate_sets <- function(X, Y, Z, graph) {

  totest <- cbind(rep(X, each = length(Y)),
        rep(Y, times = length(X)))
  for(i in 1:nrow(totest)) {

    chk <- d_separates(totest[i,1], totest[i, 2], Z, graph)
    if(!chk) return(FALSE)

  }

  TRUE

}


#' List all conditional independence relations
#'
#' Get all independence relationships in a DAG, among the observed variables
#'
#' @param dag A dag with optional vertex attributes "latent"
#' @returns A list, each element representing an independence, given 0 or more other variables
#' @export
#' @examples
#' dag <- initialize_graph(graph_from_literal(V1 -+ V2, V2 -+ V3, V3 -+ V4, U -+ V2:V4))
#' dag_list_independences(dag)
#'
dag_list_independences <- function(dag) {

  if("latent" %in% vertex_attr_names(dag)) {
    obs <- V(dag)[!latent]$name
  } else obs <- V(dag)$name

  pairs <- combn(obs, 2)

  res <- vector(mode = "list", length = ncol(pairs))
  for(j in 1:ncol(pairs)) {

    # check for d-separation given nothing
    nulchk <- d_separate_sets(pairs[1,j], pairs[2,j], NULL, dag)
    if(nulchk) {
      res[[j]] <- list(indep = pairs[,j], given = NULL)
      next
    }

    Zset <- setdiff(obs, pairs[, j])
        # check for d-separation given everything
    #maxset_chk <- d_separate_sets(pairs[1,j], pairs[2,j], Zset, dag)
    # find the minimal set
      breakout <- FALSE
      for(zi in 1:length(Zset)) {
        czin <- combn(Zset, zi, simplify = FALSE)
        for(k in 1:length(czin)) {

          inchk <- d_separate_sets(pairs[1,j], pairs[2,j], czin[[k]], dag)
          if(inchk) {
            res[[j]] <- list(indep = pairs[,j], given = czin[[k]])
            breakout <- TRUE
            break
          }

        }
        if(breakout) break
      }



  }

  res[!sapply(res, is.null)]

}


