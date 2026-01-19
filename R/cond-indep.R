#' List all conditional independence relations
#'
#' Get all independence relationships in a DAG, among the observed variables
#'
#' @param dhvm A object of class `dhvm` i.e., a caugi dag with node attributes "latent"
#' @returns A list, each element representing an independence, given 0 or more other variables
#' @export
#' @examples
#' dag <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, U %-->% V2 + V4))
#' dag_list_independences(dag)
#'
dhvm_list_independences <- function(dag) {

  obs <- subset(dag$nodes_data, latent == FALSE)$name
  pairs <- combn(obs, 2)

  res <- vector(mode = "list", length = ncol(pairs))
  for(j in 1:ncol(pairs)) {

    # check for d-separation given nothing
    nulchk <- d_separated(dag$caugi_graph, X = pairs[1,j], Y = pairs[2,j], Z = NULL)
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

          inchk <- d_separated(dag$caugi_graph, X = pairs[1,j], Y = pairs[2,j], Z = czin[[k]])
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


