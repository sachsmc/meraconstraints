#' Get the moebius parameters
#'
#' @param dag An igraph object representing a dag
#' @export
#' @returns A list
#' @examples
#' dag <- initialize_graph(graph_from_literal(V1 -+ V2, V2 -+ V3, U -+ V1:V3))
#' dag_get_moebius_params(dag)

dag_get_moebius_params <- function(dag) {

  am <- as_adjacency_matrix(dag, sparse = FALSE)  # get adjacency matrix
  if (igraph::is_named(dag)) {
    vnames <- igraph::V(dag)$name
  } else vnames <- NULL
  class(am) <- c(class(am), "adjMatrix")

  mdag <- MixedGraphs::mixedgraph(n=length(igraph::V(dag)), vnames=vnames,
                          edges=list(directed=am))

  admg <- ADMGs2::latentProject(mdag, latent = match(names(V(dag)[latent == TRUE]), mdag$vnames))

  ADMGs2:::maps(admg, dims = V(dag)[latent == FALSE]$nvals)
  mobin <- ADMGs2::moebius(admg)

  ## get out and cond for each district
  pforms <- NULL

  for(i in seq_along(mobin$q)) {

    for(j in seq_along(mobin$q[[i]])) {

      pforms[[length(pforms) + 1]] <- list(out = vnames[mobin$heads[[i]][[j]]],
                                           cond = vnames[mobin$tails[[i]][[j]]])

    }

  }

  pforms

}

#' Get linear probabilities according to district factorization
#'
#' @param dag An igraph object representing a dag
#' @returns A list of probabilities in character form
#' @export
#' @examples
#' dag <- initialize_graph(graph_from_literal(V1 -+ V2, V2 -+ V3, U -+ V1:V3))
#' dag_get_district_params(dag)

dag_get_district_params <- function(dag) {

  districts <- dag_get_districts(dag)

  res <- lapply(districts, \(vs) {

    num_latent <- sum(vs$latent == 1)
    if(num_latent <= 1) {

      candconds <- V(dag)[setdiff(V(dag)[latent == 0], vs[vs$latent == 0])]
      confconds <- sapply(adjacent_vertices(dag, candconds, mode = "out"), \(ch) any(ch %in% vs))

      list(out = names(vs[latent == 0]), cond = names(candconds)[confconds])

    } else NULL

  })

  res[!sapply(res, is.null)]

}


