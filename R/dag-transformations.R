#' Transform a regular dag into canonical partitions
#'
#' Replaces latent common causes with correlated response function variables
#'
#' @param dag An igraph object representing a dag
#' @returns An igraph object with response function variables for each observed variable
#' @export
#' @examples
#'
#' # classic IV model
#'
#' ivg <- initialize_graph(graph_from_literal(Z -+ X, X -+ Y, U -+ X:Y))
#' dag_transform_canonical_partition(ivg)
#'


dag_transform_canonical_partition <- function(dag) {

  vs <- names(V(dag)[latent == 0])
  newvs <- paste0("R_", vs)
  newnvals <- sapply(vs, \(v) {

    pa_nvals <- neighbors(dag, v, mode = "in")[latent == 0]$nvals

    # if(length(pa_nvals) == 0) {
    #    V(dag)[v]$nvals
    #  } else {
    V(dag)[v]$nvals^prod(pa_nvals)
    # }
  })

  dag2 <- add_vertices(dag, nv = length(newvs),
                       latent = rep(1, length(newvs)),
                       name = newvs,
                       nvals = newnvals,
                       respfvar = rep(TRUE, length(newvs))
  ) |>
    add_edges(unlist(lapply(vs, \(v) c(paste0("R_", v), v))))

  us <- V(dag)[latent == 1]
  for(u in us) {

    tocon <- paste0("R_", names(neighbors(dag, u, mode = "out")))
    tocon2 <- as.matrix(subset(expand.grid(tocon, tocon), Var1 != Var2))
    dag2 <- add_edges(dag2, c(t(tocon2)))

  }

  V(dag2)$respfvar[is.na(V(dag2)$respfvar)] <- FALSE

  delete_vertices(dag2, V(dag2)[latent == 1 & !respfvar])

}


#' Add response functions to a canonical partitioned dag
#'
#' @param dag A dag igraph object as returned by \link{dag_transform_canonical_partition}
#' @export
#' @returns An igraph object with response functions as vertex attributes
#' @examples
#' # classic IV example
#' ivg <- initialize_graph(graph_from_literal(Z -+ X, X -+ Y, U -+ X:Y))
#' dag_transform_response_functions(dag_transform_canonical_partition(ivg))
#'
dag_transform_response_functions <- function(dag) {

  if(is.null(V(dag)$respfvar)) {
    stop("Dag must be in canonical partition form, see dag_transform_canonical_partition")
  }

  vs <- V(dag)[latent == 0]

  rtabs <- vector(mode = "list", length = length(vs))
  for(i in 1:length(vs)) {

    v <- vs[i]
    pa_v <- neighbors(dag, v, mode = "in")
    resptable <- do.call(expand.grid, lapply(pa_v$nvals, \(x) 0:(x-1)))
    colnames(resptable) <- names(pa_v)

    ntmp <- prod(pa_v[!respfvar]$nvals)
    resptable[[names(v)]] <- c(t(as.matrix(expand.grid(lapply(1:ntmp,
                                                              \(x) 0:(v$nvals - 1))))))

    rtabs[[i]] <- resptable
  }

  V(dag)[vs]$response_function <- rtabs

  dag

}


#' Merge latent variables within districts
#'
#' Results in a linear dag that has only has 1 latent variable per district
#'
#' @param dag A igraph object representing a dag
#' @returns An igraph object representing the modified dag
#' @export
#' @examples
#' duargr <- initialize_graph(graph_from_literal(A -+ B:E, B -+ C, C -+ D, D -+ F, E -+ C,
#' U1 -+ B:D, U2 -+ A:C:E, U3 -+ D:F))
#' dag_transform_merge_latents(dag)
#'
dag_transform_merge_latents <- function(dag) {

  dists <- dag_get_districts(dag)

  latent_count <- sapply(dists, function(vs) {
    sum(vs$latent)
  })

  tomerge <- which(latent_count > 1)

  for(i in seq_along(tomerge)) {

    thisto <- dists[[tomerge[i]]][!latent]
    thisfrom <- dists[[tomerge[i]]][latent == TRUE][1]
    thisremove <- dists[[tomerge[i]]][latent == TRUE][-1]
    thisvs <- paste0(thisfrom$name, "|", thisto$name)
    newedges <- setdiff(thisvs, as_ids(E(dag)))


    dag <- add_edges(dag, unlist(strsplit(newedges, "|", fixed = TRUE)))
    dag <- delete_vertices(dag, thisremove$name)

  }

  dag
}

#' Add attributes to an igraph
#'
#' @param graph an igraph object with named vertices
#' @param latent.pattern The regular expression pattern used to identify latent variables
#' @param nvals A named vector that maps node names to the number of categories for that node. If NULL defaults to 2.
#' @returns An igraph with additional attributes latent (logical) and nvals (integer)
#' @export
#' @examples
#' initialize_graph(graph_from_literal(Z -+ X, X -+ Y, U -+ X:Y), nvals = c(Z = 3, X = 3))
#'

initialize_graph <- function(graph, latent.pattern = "^U", nvals = NULL) {

  nvars <- length(V(graph))

  if(is.null(V(graph)$latent)) {
    Unodes <- grep(latent.pattern, names(V(graph)))
    V(graph)$latent <- rep(FALSE, nvars)
    V(graph)$latent[Unodes] <- TRUE
  }


  if(is.null(V(graph)$nvals)) {
    V(graph)$nvals <- rep(2, nvars)
    if(!is.null(nvals)) {
      thisnvals <- names(nvals)
      V(graph)[thisnvals]$nvals <- nvals
    }

  }


  graph

}


#' Latent projection for use with causaleffect
#'
#' @param graph an igraph object with names vertices and latent variables
#' @returns An igraph with no latent variables but bidirected edges for use with causaleffect
#' @export
#' @examples
#' graph_from_literal(V1 -+ V2, V2 -+ V3, V3 -+ V4,
#'   U1 -+ V1:V3, U2 -+ V2:V4) |> initialize_graph() |>
#'   dag_transform_latent_project()

dag_transform_latent_project <- function(graph) {


  for(vs in V(graph)[latent == TRUE]) {

    chofu <- names(neighbors(graph, vs, mode = "out"))

    graph <- add_edges(graph, c(combn(chofu, 2)), attr = list(description = "U"))

  }

  newgraph <- delete_vertices(graph, V(graph)[latent == TRUE])

  newgraph


}
