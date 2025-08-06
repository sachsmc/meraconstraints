#' Get all topological orders of observed vertices
#'
#' @param graph An igraph object
#' @export
#' @returns A list of vectors of vertex names

dag_order <- function (graph) {

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


#' Sample a distribution from a causal model
#'
#' @param dag A causal model possibly with latent variables, currently only supports binary observed variables
#' @export
#' @returns A list of probability tables
#' @examples
#' # example code
#' ivdag <- graph_from_literal(Z -+ X, X -+ Y, U -+ X:Y) |> initialize_graph()
#' dag_sample_distribution(ivdag)
#'

dag_sample_distribution <- function(dag, parmsamp = \(k) rnorm(k)) {


  latents <- names(V(dag)[latent == TRUE])
  V(dag)$latent <- rep(FALSE, length(V(dag)))
  revtop <- rev(dag_order(dag)[[1]])
  nvals <- V(dag)[revtop]$nvals

  pmats <- vector(mode = "list", length = length(revtop) - length(latents))
  names(pmats) <- setdiff(revtop, latents)
  for(var in names(pmats)) {
    pa <- names(neighbors(dag, var, mode = "in"))

    if(length(pa) == 0) {
      prob <- runif(1)
      outmat <- data.frame(0:1, probs = c(prob, 1 - prob))
      names(outmat)[1] <- var

      pmats[[var]] <- outmat
      next
    }
    pr.var <- function(dmat, u) {

      iform <- as.formula(paste("~", paste(pa, collapse = "*")))
      inmat <- cbind(dmat, u)
      unms <- intersect(pa, latents)
      colnames(inmat)[ncol(inmat)] <- if(length(unms) == 0) "U" else unms
      dmat <- model.matrix(iform, inmat)
      parm <- parmsamp(ncol(dmat))
      plogis(dmat %*% parm)[, 1]

    }

    outmat <- do.call(expand.grid, sapply(pa, \(v) {
      if(v %in% latents) -1 else 0:(nvals[which(revtop == v)]-1)
    }, simplify = FALSE))[, setdiff(pa, latents), drop = FALSE]

    probs <- suppressWarnings(sapply(1:nrow(outmat), \(i){
      integrate(\(u) u * pr.var(outmat[i,, drop = FALSE], u), 0, 1)$value
      }))
    outmat1 <- outmat0 <- outmat
    outmat1[[var]] <- 1
    outmat1$probs <- probs
    outmat0[[var]] <- 0
    outmat0$probs <- 1-probs

    pmats[[var]] <- rbind(outmat1, outmat0)

  }

  pmats

}

