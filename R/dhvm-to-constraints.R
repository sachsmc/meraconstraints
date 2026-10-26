
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
  if(any(thiscdeg > 1)) {
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

  newcg <- caugi::remove_nodes(dag$caugi_graph,
                               name = subset(dag$nodes_data, latent == TRUE)$name) |>
    caugi::add_edges(from = newvs, edge = rep("-->", length(vs)), to = vs) |>
    mutate_caugi("PDAG")

  for(u in subset(dag$nodes_data, latent == TRUE)$name) {
    chu <- children(dag$caugi_graph, nodes = u)
    addud <- combn(paste0("R_", chu), 2)
    newcg <- caugi::add_edges(newcg, from = addud[1, ],
                              edge = rep("---", ncol(addud)),
                              to = addud[2, ])
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

    ## product over district members of P(V_1i | C_i), where C_i is the minimal
    ## conditioning set: the earlier members of the district together with the
    ## external observed parents of V_1i and of every earlier member
    ## (Proposition 1 of the paper; Tian and Pearl, 2002).

    padistudist <- union(dvars[[i]], intervenset)
    prodterm <- lapply(dvars[[i]], \(di) {
      upto    <- dvars[[i]][match(dvars[[i]], dord) <= match(di, dord)]
      earlier <- setdiff(upto, di)
      extpa   <- intersect(intervenset,
                           unique(unlist(lapply(upto, \(dv) obs_parents(dhvm, dv)))))
      intersect(padistudist, c(earlier, extpa))
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
                 district_parms = dparms,
                 dhvm_partitioned = d1r
  ),
  class = "fnrels")

}

#' Get constraints implied by a discrete hidden variables model
#'
#' @param dhvm An object of class `dhvm`, see [dhvm()]
#' @returns An object of class `meraconstraints`, with some methods.
#' @export
#' @examples
#' ivd <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y))
#' meraconstraints(ivd)
#'

meraconstraints <- function(dhvm) {

  c_degs <- dhvm_c_degree(dhvm)
  dhrels <- dhvm_linear_relations(dhvm)
  indeps <- dhvm_list_independences(dhvm)
  dhtope <- dhvm_relations_to_polytope(dhrels)
  trivial_if_false <- prob_simplex_screen(dhrels, dhtope$Hrep)

  cpart <- dhrels$dhvm_partitioned
  dhrels$dhvm_partitioned <- NULL

  names(dhtope$Vrep) <- names(dhtope$Hrep) <- names(trivial_if_false) <-
    sapply(dhrels$district_parms$poterms, \(x) paste(x$out, collapse = "_"))

  structure(list(
    model_name = deparse1(substitute(dhvm)),
    district_info = c_degs,
    dhvm_partitioned = cpart,
    functional_relations = dhrels,
    independences = indeps,
    polytope = dhtope,
    trivial_if_false = trivial_if_false
  ), class = "meraconstraints")

}

trunc_text <- function(x, width) {

  x2 <- substr(x, 1, width - 20)
  x2[nchar(x) > width - 20] <- paste(x2[nchar(x) > width - 20], "...[truncated]")
  x2

}

#' Print a meraconstraints object
#'
#' @param x As returned by [meraconstraints()]
#' @param ... Used later
#' @returns x invisibly
#' @export
print.meraconstraints <- function(x, max_constraints = 4, width = getOption("width"), ...) {

  cat("MERA constraints on the dhvm:", x$model_name, "\n")
  cat("The model has districts (unobserved variables marked with [u]:\n")
  cat("    ", paste(paste0(names(x$district_info), " with c-degree ", x$district_info), collapse = ", "), "\n")

  cat("The canonical partitioned model is available in x$dhvm_partitioned\n")
  #print(x$dhvm_partitioned)

  cat("\nConstraints:\n",
      "------------\n", sep = "")
  if(length(x$independence) == 0) {
    cat("There are no conditional independence relations.\n")
  } else {
    cat("* Conditional independence relations:\n")
    indeptext <- paste0(sapply(x$independences, \(ii) {
      paste0(paste0(ii$indep, collapse = " \u27C2 "), " | ", ii$given)
    }), collapse = ", ")

    cat(indeptext, "\n")
  }

  cat("* Equality and inequality constraints by district:\n")
  for(i in 1:length(x$polytope$Hrep)) {
    dnm <- names(x$polytope$Hrep)[i]
    thisHrep <- x$polytope$Hrep[[i]][x$trivial_if_false[[i]], ]
    n_ineq <- sum(thisHrep[, 1] == 0)
    n_eq <- sum(thisHrep[, 1] == 1)

    cat("There are ", nrow(thisHrep), " possibly nontrivial constraints involving ", dnm, ", ", n_ineq,
        " inequalities and ", n_eq, " equalities.\n", sep = "")
    if(nrow(thisHrep) > 0) {
      cat("The first", max_constraints, "expressions are\n")
      cat(trunc_text(hrep_to_text(thisHrep, x$functional_relations$probform[[i]])[1:min(nrow(thisHrep), max_constraints)],
                   width = width), sep = "\n")
    }
  }

  invisible(x)

}



