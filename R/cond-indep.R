#' List all conditional independence relations
#'
#' Get the conditional independence relations among the observed variables that are
#' implied by a DAG, as the ordered local Markov statements for every ancestral set
#' (Richardson, 2003). Every implied d-separation over disjoint sets follows from these
#' by the semi-graphoid rules.
#'
#' @param dag A object of class `dhvm` i.e., a caugi dag with node attributes "latent"
#' @returns A list, each element representing an independence of two variables, given 0 or
#'   more other variables. Together, the elements are equivalent to the statements for sets
#'   of variables, which are in the attribute `"sets"`.
#' @export
#' @examples
#' dag <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, U %-->% V2 + V4))
#' dhvm_list_independences(dag)
#'
dhvm_list_independences <- function(dag) {
  sets <- ordered_local_markov(dag)
  out <- list()
  for (s in sets) {
    v <- s$indep[[1]]; w <- s$indep[[2]]
    for (k in seq_along(w)) {
      z <- c(s$given, w[seq_len(k - 1L)])
      out[[length(out) + 1L]] <- list(indep = c(v, w[k]), given = if (length(z)) z else NULL)
    }
  }
  key <- vapply(out, \(x) paste(sort(x$indep), collapse = "~") |> paste(paste(sort(x$given), collapse = ","), sep = "|"), "")
  out <- out[!duplicated(key)]
  attr(out, "sets") <- sets
  out
}

## The ordered local Markov statements (Richardson, 2003) in set form: for every observed V and
## every ancestral set A of observed variables with V in A and A contained in pre(V) + V,
##   V _||_ A \ (V + mb(V, A)) | mb(V, A),
## where mb(V, A) is the district of V in the subgraph on A together with the parents of that
## district, without V. Statements implied by another one by decomposition are dropped.
ordered_local_markov <- function(dag) {
  obs   <- subset(dag$nodes_data, latent == FALSE)$name
  unobs <- subset(dag$nodes_data, latent == TRUE)$name
  ord   <- intersect(topological_sort(dag$caugi_graph), obs)
  adm   <- if (length(unobs)) latent_project(dag$caugi_graph, latents = unobs) else dag$caugi_graph
  e     <- as.data.frame(edges(adm))
  dir   <- e[e$edge == "-->", , drop = FALSE]
  bi    <- e[e$edge == "<->", , drop = FALSE]
  anc   <- lapply(setNames(obs, obs), \(v) {
    a <- v
    repeat { n <- union(a, dir$from[dir$to %in% a]); if (length(n) == length(a)) return(a); a <- n }
  })
  res <- list()
  for (i in seq_along(ord)) {
    v     <- ord[i]
    pre   <- ord[seq_len(i - 1L)]
    must  <- setdiff(anc[[v]], v)
    free  <- setdiff(pre, must)
    extra <- c(list(character()),
               unlist(lapply(seq_along(free), \(k) combn(free, k, simplify = FALSE)), recursive = FALSE))
    for (x in extra) {
      A <- c(must, x, v)
      if (!all(unlist(anc[A]) %in% A)) next                       # not ancestral
      comp <- v
      repeat {
        b <- bi[bi$from %in% A & bi$to %in% A, , drop = FALSE]
        n <- union(comp, c(b$to[b$from %in% comp], b$from[b$to %in% comp]))
        if (length(n) == length(comp)) break
        comp <- n
      }
      mb   <- intersect(ord, setdiff(union(comp, dir$from[dir$to %in% comp]), v))
      rest <- intersect(ord, setdiff(A, c(v, mb)))
      if (length(rest)) res[[length(res) + 1L]] <- list(indep = list(v, rest), given = mb)
    }
  }
  ## drop duplicates and statements implied by decomposition
  keep <- vapply(seq_along(res), \(i) !any(vapply(seq_along(res), \(j) {
    j != i && res[[j]]$indep[[1]] == res[[i]]$indep[[1]] && setequal(res[[j]]$given, res[[i]]$given) &&
      all(res[[i]]$indep[[2]] %in% res[[j]]$indep[[2]]) &&
      (length(res[[j]]$indep[[2]]) > length(res[[i]]$indep[[2]]) || j < i)
  }, logical(1))), logical(1))
  res[keep]
}
