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
  gRbase::separates(X, Y, Z, moralize(ag))

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

  obs <- V(dag)[!latent]$name
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


#' Partition into districts
#'
#' Districts are defined for Acyclic directed mixed graphs, but they are basically the same as c-components
#'
#' @param dag A directed acyclic graph. It should have vertex attributed called "latent", that is logical and indicates if a vertex is latent/unobserved/hidden.
#' @returns A list of the districts
#' @export
#' @examples
#' dag <- initialize_graph(graph_from_literal(V1 -+ V2, V2 -+ V3, V3 -+ V4, U -+ V2:V4))
#' dag_get_districts(dag)
#'

dag_get_districts <- function(dag) {

  if(is.null(V(dag)$latent)){
    warning("Graph vertices do not have attribute latent, all are in their own district:")
    return(V(dag))
  }

  dagu <- delete_edges(dag, E(dag)[V(dag)[latent == 0] %--% V(dag)[latent == 0]])
  dists <- components(dagu)$membership

  lapply(sort(unique(dists)), \(i) {

    V(dag)[names(dists[dists == i])]

  })

}

#' Get the nodes that the given node depends on
#'
#' @param dag A dag with vertex attribute "latent"
#' @param node The name of a node to query
#' @returns A character vector of the names of observed nodes that node depends on
#' @export
#' @examples
#' dag <- initialize_graph(graph_from_literal(V1 -+ V2, V2 -+ V3, V3 -+ V4, U -+ V2:V4))
#' dag_get_dependencies(dag, "V4")
#' dag_get_dependencies(dag, "V3")

dag_get_dependencies <- function(dag, node) {

  districts <- dag_get_districts(dag)
  dist_of_node <- districts[[which(sapply(districts, \(d) node %in% names(d)))]][!latent]
  nd_node_in_dist <- difference(dist_of_node, V(dag)[descendants(node, dag)])
  pa_dist <- adjacent_vertices(dag, dist_of_node, mode = "in")

  little_set <- V(dag)[latent == TRUE | name == node]
  difference(do.call("c", c(list(nd_node_in_dist), pa_dist)),
                   little_set)$name

}

#' Find all descendants of a set in graph
#'
#' @param set Vector of node names
#' @param graph A dag
#' @returns a vector of the descendents of set
#' @export
#' @examples
#' g1 <- graph_from_literal(X1 -+ X2, X2 -+ X3, X3 -+ X4,
#'      U1 -+ X1:X3, U2 -+ X2:X4)
#' descendants(c("U1", "U2"), g1)

descendants <- function(set, graph) {

  amat <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
  An <- setorig <- set
  #amat <- amat[match(set, rownames(amat)), , drop = FALSE]
  repeat {
    set2 <- colSums(amat[set, , drop = FALSE])
    set <- names(which(set2 > 0))
    if (!length(set))
      (break)()
    An <- c(An, set)
  }
  setdiff(An, setorig)

}

