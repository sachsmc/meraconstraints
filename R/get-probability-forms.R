#' Get the moebius parameters
#'
#' @param dag An igraph object representing a dag
#' @export
#' @returns A list
#' @examples
#' dag <- initialize_graph(graph_from_literal(V1 -+ V2, V2 -+ V3, U -+ V1:V3))
#' dag_get_moebius_params(dag)

dag_get_moebius_params <- function(dag) {

  am <- as_adjacency_matrix(dag, sparse = FALSE)
  if (igraph::is_named(dag)) {
    vnames <- igraph::V(dag)$name
  } else vnames <- NULL
  class(am) <- c(class(am), "adjMatrix")

  mdag <- MixedGraphs::mixedgraph(n=length(igraph::V(dag)), vnames=vnames,
                          edges=list(directed=am))

  admg <- ADMGs2::latentProject(mdag, latent = match(names(V(dag)[latent == TRUE]), mdag$vnames))

  amap <- ADMGs2:::maps(admg, dims = V(dag)[latent == FALSE]$nvals)
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
  dvars <- lapply(districts, \(vs) names(vs[latent == FALSE]))
  admg <- dag_transform_latent_project(dag)
  dord <- dag_order(dag)[[1]]

  res <- character(length(districts))
  remove <- NULL
  vterms <- poterms <- distrs <- resterms <- vector(mode = "list", length = length(districts))
  for(i in 1:length(districts)) {

    paofdist <- adjacent_vertices(dag, districts[[i]][latent == FALSE], mode = "in")
    intervenset <- setdiff(unique(unlist(lapply(paofdist, \(vs) names(vs[latent == FALSE])))),
                           names(districts[[i]][latent == FALSE]))

    if(length(intervenset) == 0 & length(dvars[[i]] == 1)){
      remove <- c(remove, i)
      next
    }

    ## product of each variable given it's parents

      padistudist <- union(dvars[[i]], unlist(lapply(paofdist, \(pa) names(pa[latent == FALSE]))))
      prodterm <- lapply(dvars[[i]], \(di) {
        setdiff(intersect(padistudist, dord[1:which(di == dord)]), di)
      })

    names(prodterm) <- dvars[[i]]
    resterms[[i]] <- prodterm

    res[i] <- Reduce(\(x, y) paste0(x, y), lapply(names(prodterm), \(y) {
      sprintf("P(%s = {%s}|%s)", y, y, paste(paste0(prodterm[[y]], " = {", prodterm[[y]], "}"), collapse = ","))
    }))
    #res[i] <- causal.effect(y = dvars[[i]],
    #              x = intervenset, G = admg, prune = TRUE)

    poterms[[i]] <- list(out = dvars[[i]], do = intervenset)
    vterms[[i]] <- c(dvars[[i]], intervenset)
    distrs[[i]] <- paste0("R_", dvars[[i]])
  }

  if(!is.null(remove)) {
    res <- res[-remove]
    distrs <- distrs[-remove]
    vterms <- vterms[-remove]
    resterms <- resterms[-remove]
    poterms <- poterms[-remove]
  }


  list(vterms = vterms,
       pforms = res, distrs = distrs,
       plist = resterms, poterms = poterms)

}


