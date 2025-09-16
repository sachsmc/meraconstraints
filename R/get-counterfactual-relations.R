#' Get response patterns by recursive substitution
#'
#' @param tab A data frame of values for the observed variables
#' @param dag An igraph object representing a dag as returned by \link{dag_transform_response_functions}
#' @param target A character string for the variable of interest
#' @returns A data frame
#' @export

upmerge <- function(tab, dag, target) {

  pa_d <- neighbors(dag, target, mode = "in")[!respfvar]
  if(length(pa_d) == 0) {

    merge(V(dag)[target]$response_function[[1]], tab, all = FALSE)

  } else {

    inres <- vector(mode = "list", length = length(pa_d))
    for(i in 1:length(pa_d)) {
      v <- pa_d[i]
      inres[[i]] <- upmerge(v$response_function[[1]], dag,
                            target = names(v))

    }

    Reduce(\(x, y) merge(x, y, all = FALSE), init = tab,
           inres, simplify = FALSE)


  }

}



#' Get table of response patterns matching a set of observables
#'
#' @param dag A dag as returned by \link{dag_transform_response_functions}
#' @param event.table data.frame with the observable events of interest
#' @returns A data.frame including the response variable events
#' @export
#'
dag_get_cf_tables <- function(dag, event.table){

  if(missing(event.table)) {
    args <- lapply(names(V(dag)[latent == 0]),
                   \(v) 0:(V(dag)[v]$nvals - 1))
    names(args) <- names(V(dag)[latent == 0])
    event.table <- do.call(expand.grid, args)
  }

  if(!all(sapply(colnames(event.table), match, names(V(dag)), nomatch = 0) > 0)) {
    stop("Variables in event.table that are not in dag")
  }

  thesevars <- paste0(rep(c("", "R_"), each = ncol(event.table)), colnames(event.table))
  subdag <- delete_vertices(dag, setdiff(names(V(dag)), thesevars))
  nochildren <- adjacent_vertices(subdag, V(subdag), mode = "out")
  nochildren <- intersect(names(nochildren[sapply(nochildren, length) == 0]),
                          colnames(event.table))

  toporder <- dag_order(subdag)[[1]]

  thistab2 <- event.table
  for(vs in toporder) {
    thistab2 <- upmerge(thistab2, dag, vs)
  }
  ## merge back in any variables without children
  if(!is.null(nochildren)) {
    for(v in nochildren){
      thistab2 <- merge(V(dag)[v]$response_function[[1]], thistab2, all = FALSE)
    }
  }

  unique(thistab2[, thesevars])


}

#' For testing

dag_get_cf_tables2 <- function(dag) {

  obsVars <- V(dag)[latent == FALSE]
  qlist <- lapply(obsVars$response_function, \(vf) unique(vf[, intersect(names(vf), paste0("R_", names(obsVars)))]))
  names(qlist) <- paste0("R_", names(obsVars))
  q.vals <- do.call(expand.grid, qlist)
  respvars <- obsVars$response_function
  names(respvars) <- names(obsVars)

  PAs <- lapply(obsVars, function(x) names(neighbors(dag, x, mode = "in")[latent == FALSE]))

  gee_r <- function(r, i) {

    parents <- PAs[[i]]
    if(length(parents) == 0) {

      respvars[[i]][respvars[[i]][[paste0("R_", i)]] == r[, paste0("R_", i)], i]

    } else {

      lookin <- lapply(parents, function(gu) {
        as.numeric(gee_r(r, gu))
      })
      names(lookin) <- parents
      lookind <- as.data.frame(lookin)
      lookind[[paste0("R_", i)]] <- r[, paste0("R_", i)]

      inres <- merge(respvars[[i]], lookind)[, i]
      inres

    }

  }

  res.mat <- matrix(NA, nrow = nrow(q.vals), ncol = length(obsVars))
  colnames(res.mat) <- names(obsVars)
  for (k in 1:nrow(q.vals)) {
    for (j in names(obsVars)) {
      res.mat[k, j] <- gee_r(r = q.vals[k, , drop = FALSE], i = j)
    }
  }

  cbind(res.mat, q.vals)

}

#' Get relations between probabilities and response variables
#'
#' @param dag A causal dag
#' @returns A list of character vectors
#' @export

dag_get_linear_relations <- function(dag) {

  d1r <- dag_transform_canonical_partition(dag) |>
    dag_transform_response_functions()
  vcombs <- lapply(V(dag)[latent == FALSE]$nvals, \(v) seq_len(v) - 1)
  names(vcombs) <- names(V(dag)[latent == FALSE])

  ptabs <- dag_get_cf_tables(d1r, do.call(expand.grid, vcombs))

  dparms <- dag_get_district_params(dag)

  relation <- prob.fancy <- problevels <- vector(mode = "list", length = length(dparms$pforms))
  for(j in 1:length(dparms$pforms)) {

    splitptab <- split(ptabs, ptabs[, dparms$vterms[[j]]])
    relation[[j]] <- lapply(splitptab, \(df) {
      unique(make_prob_names(list(out = dparms$distrs[[j]], cond = NULL),
                             df[, dparms$distrs[[j]], drop = FALSE]))
    })
    prob.fancy[[j]] <- lapply(splitptab, \(df) {
      unique(glue::glue_data(df, dparms$pforms[[j]]))
    })
    problevels[[j]] <- unique(do.call(rbind, splitptab)[, dparms$vterms[[j]]])

  }

  list(relation = relation,
       probform = prob.fancy,
       problevels = problevels,
       problist = dparms$plist)

}




