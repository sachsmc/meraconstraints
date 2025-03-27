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

  thistab2 <- event.table
  for(vs in names(event.table)) {
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

