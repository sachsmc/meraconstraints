#' Get all topological orders of observed vertices
#'
#' @param graph An igraph object
#' @export
#' @returns A list of vectors of vertex names

obs_tops <- function (graph) {

  obsVs <- V(graph)[latent == FALSE]
  PAs <- lapply(obsVs, function(x) neighbors(graph,x, mode = "in")[latent == FALSE])
  CHs <- lapply(obsVs, function(x) neighbors(graph,x, mode = "out")[latent == FALSE])

  inDeg <- lengths(PAs)
  sq <- seq_along(obsVs)

  out <- list()
  vst <- rep(FALSE, length(CHs))

  ## get vertices with no parents
  zeroIn <- which(inDeg == 0)
  if (length(zeroIn) == 0) stop("Graph is cyclic")

  for (k in seq_along(zeroIn)) {
    inDeg[sq[CHs[[sq[k]]]]] <- inDeg[sq[CHs[[sq[k]]]]] - 1
    vst[sq[k]] <- TRUE

    tmp <- ADMGs2:::to_recurse(sapply(CHs, names), inDeg, vst)
    tmp <- lapply(tmp, function(x) c(k, x))

    inDeg[sq[CHs[[sq[k]]]]] <- inDeg[sq[CHs[[sq[k]]]]] + 1
    vst[sq[k]] <- FALSE

    out <- c(out, tmp)
  }

  lapply(out, \(i) names(obsVs[i]))

}
