#' @import caugi
NULL

#' @export
dhvm <- function(caugi_graph, latent.pattern = "^U",
                                             latent.variables = NULL, nvals = NULL) {

                        nodes_latent <- rep(FALSE, nrow(nodes(caugi_graph)))
                        names(nodes_latent) <- nodes(caugi_graph)$name

                        if(!is.null(latent.variables)) {
                          mnm <- match(latent.variables, names(nodes_latent))
                          if(any(is.na(nmn))) stop("Cannot find variable(s): ", latent.variables[is.na(mnm)] )
                          nodes_latent[mnm] <- TRUE
                        } else {
                          Unodes <- grep(latent.pattern, nodes(caugi_graph)$name)
                          nodes_latent[Unodes] <- TRUE
                        }

                        nodes_nvals <- rep(2L, length(nodes(caugi_graph)))
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
       node_style = list(by_node = filllist))
}


obs_parents <- function(dhvm, node) {

  patmp <- parents(dhvm$caugi_graph, nodes = node)
  intersect(patmp, subset(dhvm$nodes_data, latent == FALSE)$name)

}

nvals <- function(dhvm, node) {

  subset(dhvm$nodes_data, name == node)$nvals

}

#' Get the c-degree of a dhvm
#'
#' The c-degree is the maximum of 1 and the numbers of latent variables in each district
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
  ccomps <- districts(latent_project(dhvm$caugi_graph, latents = unobs))

  cdegs <- rep(NA, length(ccomps))
  unobs_dist <- vector(mode = "list", length = length(ccomps))
  nmstmp <- rep(NA_character_, length(ccomps))
  for(i in seq_along(ccomps)) {

    unobs_dist <- sort(intersect(unique(unlist(parents(dhvm$caugi_graph, ccomps[[i]]))), unobs))
    cdegs[[i]] <- max(1, length(unobs_dist))
    nmstmp[i] <- paste0("{", paste(c(ccomps[[i]],
                        if(length(unobs_dist) == 0) NULL else paste0(unobs_dist, "[u]")), collapse = ", "), "}")

  }

  setNames(cdegs, nmstmp)

}

#' Transform a regular dag into canonical partitions and add response function tables
#'
#' Replaces latent common causes with correlated response function variables
#'
#' @param dag An object of class "dhvm" representing a dag
#' @returns A modified object with response function variables for each observed variable
#' @export
#' @examples
#'
#' # classic IV model
#'
#' ivd <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y))
#' ivdcp <- dhvm_canonical_partition(ivd)
#'

dhvm_canonical_partition <- function(dag) {

  stopifnot(inherits(dag, "dhvm"))
  ## check for multiple latents

  thiscdeg <- dhvm_c_degree(dag)
  if(!all(thiscdeg == 1)) {
    stop("c-degree is greater than 1, cannot do the canonical partitioning. Check district(s): ",
         names(thiscdeg[thiscdeg > 1])
         )
  }

  vs <- subset(dag$nodes_data, latent == FALSE)$name
  newvs <- paste0("R_", vs)
  newnvals <- sapply(vs, \(v) {

    obspa <- obs_parents(dag, node = v)
    pa_nvals <- if(length(obspa) == 0) 1 else sapply(obspa, \(v2) nvals(dag, v2))
    subset(dag$nodes_data, name == v)$nvals^prod(pa_nvals)

  })

  newcg <- caugi::remove_nodes(dag$caugi_graph, name = subset(dag$nodes_data, latent == TRUE)$name) |>
    caugi::add_edges(from = newvs, edge = rep("-->", length(vs)), to = vs) |>
    mutate_caugi("PDAG")

  for(u in subset(dag$nodes_data, latent == TRUE)$name) {
    chu <- children(dag$caugi_graph, nodes = u)
    addud <- combn(paste0("R_", chu), 2)
    newcg <- caugi::add_edges(newcg, from = addud[1, ],
                              edge = rep("---", ncol(addud)), to = addud[2, ], inplace = TRUE)
  }

  newcg <- build(newcg)

  newnd <- data.frame(name = c(vs, newvs), latent = rep(c(FALSE, TRUE),
                                                        each = length(vs)),
                      nvals = c(sapply(vs, \(v) nvals(dag, v)), newnvals))

  newdhvm <- structure(list(caugi_graph= newcg, nodes_data = newnd), class = "dhvm")


  rtabs <- vector(mode = "list", length = length(vs))
  for(i in 1:length(vs)) {

    v <- vs[i]
    pa_v <- parents(newcg, v)
    resptable <- do.call(expand.grid, lapply(pa_v, \(x) 0:(nvals(newdhvm, x)-1)))
    colnames(resptable) <- pa_v

    opa_v <- obs_parents(newdhvm, v)

    ntmp <- if(length(opa_v) == 0) 1 else prod(sapply(opa_v, \(op) nvals(newdhvm, op)))
    resptable[[v]] <- c(t(as.matrix(expand.grid(lapply(1:ntmp,
                           \(x) 0:(nvals(newdhvm, v) - 1))))))

    rtabs[[i]] <- resptable
  }
  names(rtabs) <- vs

  newdhvm$response_function <- rtabs
  class(newdhvm) <- c("canpart", class(newdhvm))
  newdhvm

}



dhvm_upmerge <- function(tab, dhvm, target) {

  pa_d <- obs_parents(dhvm, target)
  if(length(pa_d) == 0) {

    merge(dhvm$response_function[[target]], tab, all = FALSE)

  } else {

    inres <- vector(mode = "list", length = length(pa_d))
    for(i in 1:length(pa_d)) {
      v <- pa_d[i]
      inres[[i]] <- dhvm_upmerge(dhvm$response_function[[v]], dhvm,
                            target = v)

    }

    Reduce(\(x, y) merge(x, y, all = FALSE), init = tab,
           inres, simplify = FALSE)


  }

}



#' Get linear probabilities according to district factorization
#'
#' @param dhvm An object of class "dhvm" representing a dag
#' @returns A list of probabilities in character form
#' @export
#' @examples
#' ivd <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y))
#' ivdcp <- dhvm_district_params(ivd)
#'

dhvm_district_params <- function(dhvm) {

  dvars <- districts(latent_project(dhvm$caugi_graph, subset(dhvm$nodes_data, latent == TRUE)$name))
  dord <- topological_sort(dhvm$caugi_graph)

  res <- character(length(dvars))
  remove <- NULL
  vterms <- poterms <- distrs <- resterms <- vector(mode = "list", length = length(dvars))
  for(i in 1:length(dvars)) {

    intervenset <- setdiff(unique(unlist(lapply(dvars[[i]], \(dv) obs_parents(dhvm, dv)))),
                           dvars[[i]])

    if(length(intervenset) == 0 & length(dvars[[i]]) == 1){
      remove <- c(remove, i)
      next
    }

    ## product of each variable given it's parents

    padistudist <- union(dvars[[i]], intervenset)
    prodterm <- lapply(dvars[[i]], \(di) {
      setdiff(intersect(padistudist, dord[1:which(di == dord)]), di)
    })

    names(prodterm) <- dvars[[i]]
    resterms[[i]] <- prodterm

    res[i] <- Reduce(\(x, y) paste0(x, y), lapply(names(prodterm), \(y) {
      if(length(prodterm[[y]] > 0)) {
        sprintf("P(%s = {%s}|%s)", y, y,
                paste(paste0(prodterm[[y]], " = {", prodterm[[y]], "}"), collapse = ","))
      } else {
        sprintf("P(%s = {%s})", y, y)
      }

    }))

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



#' Get table of response patterns matching a set of observables
#'
#' @param dhvm An object of class "dhvm" representing a dag with response functions
#' @param event.table data.frame with the observable events of interest
#' @returns A data.frame including the response variable events
#' @export
#'
dhvm_cf_tables <- function(dhvm, event.table){

  if(missing(event.table)) {
    args <- lapply(subset(dhvm$nodes_data, latent == FALSE)$name,
                   \(v) 0:(nvals(dhvm, v) - 1))
    names(args) <- subset(dhvm$nodes_data, latent == FALSE)$name
    event.table <- do.call(expand.grid, args)
  }

  if(!all(sapply(colnames(event.table), match, dhvm$nodes_data$name, nomatch = 0) > 0)) {
    stop("Variables in event.table that are not in dag")
  }

  thesevars <- paste0(rep(c("", "R_"), each = ncol(event.table)), colnames(event.table))
  subdag <- if(length(setdiff(dhvm$nodes_data$name, thesevars)) == 0) {
    dhvm$caugi_graph
    } else caugi::subgraph(dhvm$caugi_graph, nodes = setdiff(dhvm$nodes_data$name, thesevars))

  nochildren <- unlist(lapply(nodes(subdag)$name, \(an) {
    if(length(children(subdag, nodes = an)) == 0) an else NULL
  }))
  nochildren <- intersect(nochildren,
                          colnames(event.table))

  toporder <- topological_sort(mutate_caugi(caugi::subgraph(subdag, nodes = subset(dhvm$nodes_data, latent == FALSE)$name),
                                            "DAG"))

  thistab2 <- event.table
  for(vs in toporder) {
    thistab2 <- dhvm_upmerge(thistab2, dhvm, vs)
  }
  ## merge back in any variables without children
  if(!is.null(nochildren)) {
    for(v in nochildren){
      thistab2 <- merge(dhvm$response_function[[v]], thistab2, all = FALSE)
    }
  }

  unique(thistab2[, thesevars])


}


#' Get relations between probabilities and response variables
#'
#' @param dhvm A causal dag
#' @returns A list of character vectors
#' @export

dhvm_linear_relations <- function(dhvm) {

  d1r <- dhvm_canonical_partition(dhvm)

  vcombs <- lapply(subset(dhvm$nodes_data, latent == FALSE)$name, \(v) seq_len(nvals(dhvm, v)) - 1)
  names(vcombs) <- subset(dhvm$nodes_data, latent == FALSE)$name

  ptabs <- dhvm_cf_tables(d1r, do.call(expand.grid, vcombs))

  dparms <- dhvm_district_params(dhvm)

  relation <- prob.fancy <- prob.star <- problevels <- vector(mode = "list", length = length(dparms$pforms))
  for(j in 1:length(dparms$pforms)) {

    splitptab <- split(ptabs, ptabs[, dparms$vterms[[j]]])
    relation[[j]] <- lapply(splitptab, \(df) {
      unique(make_prob_names(list(out = dparms$distrs[[j]], cond = NULL),
                             df[, dparms$distrs[[j]], drop = FALSE]))
    })
    prob.fancy[[j]] <- lapply(splitptab, \(df) {
      unique(glue::glue_data(df, dparms$pforms[[j]]))
    })

    startemplate <- if(length(dparms$poterms[[j]]$do > 0)) {
      sprintf("P*(%s|%s)", paste(paste0(dparms$poterms[[j]]$out, " = {", dparms$poterms[[j]]$out, "}"), collapse = ","),
              paste(paste0(dparms$poterms[[j]]$do, " = {", dparms$poterms[[j]]$do, "}"), collapse = ",")
              )
    } else {
      sprintf("P(%s)", paste(paste0(dparms$poterms[[j]]$out, " = {", dparms$poterms[[j]]$out, "}"), collapse = ","))
    }
    prob.star[[j]] <- lapply(splitptab, \(df) {
      unique(glue::glue_data(df,  startemplate))
    })
    problevels[[j]] <- unique(do.call(rbind, splitptab)[, dparms$vterms[[j]]])

  }

  structure(list(relation = relation,
       probform = prob.fancy,
       probstarform = prob.star,
       problevels = problevels,
       problist = dparms$plist,
       poterms = dparms$poterms),
       class = "fnrels")

}

