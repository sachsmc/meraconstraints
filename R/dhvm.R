#' @import caugi
NULL

#' Create a Discrete Hidden Variables Model (dhvm)
#'
#' A DAG over some unobserved variables and some observed variables, where
#' all the observed variables are categorical, and the unobserved variables
#' have no restriction.
#'
#' @param caugi_graph A caugi DAG
#' @param latent.pattern A regular expression that is used to identify latent variables
#' @param latent.variable A character vector indicating which variables are latent. If not NULL, then `latent.pattern` is ignored
#' @param nvals A named numeric vector specifying the number of categories for each of the observed variables. If not specified, then 2 categories are assumed
#' @export
#' @returns An object of class `dhvm`, which is a list with the caugi object, and a data.frame with the node information
#' @examples
#' ivd1 <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y))
#' ivd2 <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y), latent.variables = "U")
#' ivd3 <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y), nvals = c(Z = 3))
#'

dhvm <- function(caugi_graph, latent.pattern = "^U",
                                             latent.variables = NULL, nvals = NULL) {

                        nodes_latent <- setNames(rep(FALSE, nrow(nodes(caugi_graph))),
                                                 nodes(caugi_graph)$name)

                        if(!is.null(latent.variables)) {
                          mnm <- match(latent.variables, names(nodes_latent))
                          if(any(is.na(mnm))) stop("Cannot find variable(s): ", latent.variables[is.na(mnm)] )
                          nodes_latent[mnm] <- TRUE
                        } else {
                          Unodes <- grep(latent.pattern, nodes(caugi_graph)$name)
                          nodes_latent[Unodes] <- TRUE
                        }

                        nodes_nvals <- setNames(rep(2L, length(nodes(caugi_graph)$name)),
                                                nodes(caugi_graph)$name)
                        if(!is.null(nvals)) {

                          mnn <- match(names(nvals), nodes(caugi_graph)$name)
                          if(any(is.na(mnn))) stop("Cannot find variable(s): ", names(nvals)[is.na(mnn)] )
                          nodes_nvals[names(nvals)] <- nvals

                        }
                        nodes_data <- nodes(caugi_graph)
                        nodes_data$latent <- nodes_latent
                        nodes_data$nvals <- nodes_nvals

                       structure(list(caugi_graph = caugi_graph, nodes_data = nodes_data), class = "dhvm")
                      }


#' @export
print.dhvm <- function(x,
                                    max_nodes = getOption("caugi.max_nodes"),
                                    max_edges = getOption("caugi.max_edges"),
                                    ...) {

  cat("A Discrete Hidden Variables Model (dhvm)\n")
  cat("* graph data: \n")
  print(x$caugi_graph, max_nodes, max_edges)
  cat("* node data: ", "\n")
  print(head(x$nodes_data, max_nodes %||% Inf))
  invisible(x)

}

#' @export
plot.dhvm <- function(x,
                      layout = "auto",
                      latent_fill = "grey75",
                      fill = "skyblue",
                      edge_style = list(),
                      label_style = list(),
                      main = NULL,
                      title_style = list(),
                      outer_margin = grid::unit(2, "mm"),
                      title_gap = grid::unit(1, "lines"),
                      ...) {
  filllist <- lapply(1:nrow(x$nodes_data), \(i) {
    list(fill = if(x$nodes_data$latent[i]) latent_fill else fill)
  })
  names(filllist) <- x$nodes_data$name
  plot(x$caugi_graph,
       node_style = list(by_node = filllist),
       edge_style = edge_style,
       label_style = label_style,
       main = main,
       title_style = title_style,
       outer_margin = outer_margin,
       title_gap = title_gap,
       layout = layout,
       ...)
}


obs_parents <- function(dhvm, node) {

  patmp <- parents(dhvm$caugi_graph, nodes = node)
  intersect(patmp, subset(dhvm$nodes_data, latent == FALSE)$name)

}

nvals <- function(dhvm, node) {

  subset(dhvm$nodes_data, name == node)$nvals

}


#' Get the districts of a dhvm
#'
#' The districts of a DHVM are the sets of nodes that remain connected after removing all edges that have tails at observed nodes. The set of districts forms a partition of the variables in the model.
#'
#' @param dhvm A object of class dhvm
#' @returns A list of character vectors
#' @export
#' @examples
#' dhvm_districts(dhvm(caugi(U1 %-->% X + Y, U2 %-->% Z + X)))
#' dhvm_districts(dhvm(caugi(Z %-->% X %-->% Y, U1 %-->% X + Y)))
#'
dhvm_districts <- function(dhvm) {

  unobs <- subset(dhvm$nodes_data, latent == TRUE)$name
  lpdists <- districts(latent_project(dhvm$caugi_graph, latents = unobs))
  lapply(lpdists, \(i) union(i, intersect(unique(unlist(parents(dhvm$caugi_graph, i))), unobs)))

}


#' Get the c-degree of a dhvm
#'
#' The c-degree of a district is the number of latent variables in it. A district
#' without latent variables (a single observed variable) has c-degree 0.
#'
#' @param dhvm A object of class dhvm
#' @returns An integer vector
#' @export
#' @examples
#' dhvm_c_degree(dhvm(caugi(U1 %-->% X + Y, U2 %-->% Z + X)))
#' dhvm_c_degree(dhvm(caugi(Z %-->% X %-->% Y, U1 %-->% X + Y)))
#'

dhvm_c_degree <- function(dhvm) {

  unobs <- subset(dhvm$nodes_data, latent == TRUE)$name
  ccomps <- dhvm_districts(dhvm)

  cdegs <- rep(NA, length(ccomps))
  unobs_dist <- vector(mode = "list", length = length(ccomps))
  nmstmp <- rep(NA_character_, length(ccomps))
  for(i in seq_along(ccomps)) {

    unobs_dist <- sort(intersect(ccomps[[i]], unobs))
    obs_dist <- sort(setdiff(ccomps[[i]], unobs))
    cdegs[[i]] <- length(unobs_dist)
    nmstmp[i] <- paste0("{", paste(c(obs_dist,
                        if(length(unobs_dist) == 0) NULL else paste0(unobs_dist, "[u]")), collapse = ", "), "}")

  }

  setNames(cdegs, nmstmp)

}


