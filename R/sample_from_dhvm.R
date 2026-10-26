#' @import data.table
NULL

#' Sample a distribution from a causal model
#'
#' @param dhvm A causal model possibly with latent variables, currently only supports binary observed variables
#' @param parmsamp A function that takes a integer as the single argument, and returns the parameters used in the linear predictors
#' @param nsamp_latent The number of random draws to take from each unobserved variable
#' @export
#' @returns A data frame where each row represents a distinct observation vector plus the joint probability
#' @examples
#' # example code
#' ivdag <- caugi(Z %-->% X %-->% Y, U %-->% X + Y) |> dhvm()
#' dhvm_sample_distribution(ivdag)
#'

dhvm_sample_distribution <- function(dhvm, parmsamp = \(k) rnorm(k), nsamp_latent = 1e3) {


  latents <- subset(dhvm$nodes_data, latent == TRUE)$name
  revtop <- rev(topological_sort(dhvm$caugi_graph))
  nvals_lookup <- setNames(dhvm$nodes_data$nvals, dhvm$nodes_data$name)

  pmats <- vector(mode = "list", length = length(revtop) - length(latents))
  names(pmats) <- setdiff(revtop, latents)

  if(length(latents) == 0) {
    for(var in names(pmats)) {
      pa <- parents(dhvm$caugi_graph, nodes = var)

      if(length(pa) == 0) {
        prob <- diff(sort(c(0, 1, runif(nvals_lookup[var] - 1))))
        outmat <- data.table(as.factor(0:(nvals_lookup[var] - 1)), probs = prob)
        names(outmat)[1] <- var

        pmats[[var]] <- outmat
        colnames(pmats[[var]])[which(colnames(pmats[[var]]) == "probs")] <- paste0("prob.", var)
        next
      }


        iform <- as.formula(paste("~", paste(pa, collapse = "*")))

        outmat <- data.table(do.call(expand.grid, sapply(pa, \(var) {
          as.factor(0:(nvals_lookup[var]-1))
        }, simplify = FALSE)))

        dmat <- model.matrix(iform, outmat)
        ## need nval - 1 paramters per column
        probs <- lapply(1:(nvals_lookup[var] - 1), \(i){
          exp(dmat %*% parmsamp(ncol(dmat)))
        })

        denom <- cbind(1 + rowSums(do.call(cbind, probs)))
        oprobs <- c(list(1 / denom), lapply(probs, \(p) p / denom))
        pmats[[var]] <- do.call(rbind, lapply(1:length(oprobs), \(i){

          df <- outmat
          df[[var]] <- factor(i - 1, levels = 0:nvals_lookup[var])
          df[[paste0("prob.", var, "_", paste(setdiff(pa, latents), collapse = "."))]] <- oprobs[[i]]
          df

        }))

    }

    jdf <- Reduce(\(x, y) merge.data.table(x, y, sort = FALSE, allow.cartesian = TRUE), pmats)
    jdf[, prob.joint := apply(as.matrix(.SD), 1, prod), .SDcols = grep("prob.", colnames(jdf), fixed = TRUE)]

    return(as.data.frame(jdf))

  } else {

    ocmat <- data.table(do.call(expand.grid, sapply(dhvm$nodes_data$name, \(var) {
      if(var %in% latents) {
        runif(nsamp_latent)
        } else {
          as.factor(0:(nvals_lookup[var]-1))
        }
    }, simplify = FALSE)))

    for(var in names(pmats)) {
      pa <- parents(dhvm$caugi_graph, nodes = var)

      if(length(pa) == 0) {
        prob <- diff(sort(c(0, 1, runif(nvals_lookup[var] - 1))))
        outmat <- data.table(as.factor(0:(nvals_lookup[var] - 1)), probs = prob)
        names(outmat)[1] <- var

        pmats[[var]] <- outmat
        colnames(pmats[[var]])[which(colnames(pmats[[var]]) == "probs")] <- paste0("prob.", var)
        next
      }

        iform <- as.formula(paste("~", paste(pa, collapse = "*")))
        thisocmat <- unique(ocmat[, ..pa])
        dmat <- model.matrix(iform, thisocmat)
        ## need nval - 1 paramters per column
        probs <- lapply(1:(nvals_lookup[var] - 1), \(i){
          exp(dmat %*% parmsamp(ncol(dmat)))
        })

        denom <- cbind(1 + rowSums(do.call(cbind, probs)))
        oprobs <- c(list(1 / denom), lapply(probs, \(p) p / denom))
        pmats[[var]] <- do.call(rbind, lapply(1:length(oprobs), \(i){

          df <- thisocmat
          df[[var]] <- factor(i - 1, levels = 0:nvals_lookup[var])
          df[[paste0("prob.", var, "_", paste(setdiff(pa, latents), collapse = "."))]] <- oprobs[[i]]
          unique(df)

        }))


    }

    jdf <- Reduce(\(x, y) merge.data.table(x, y, sort = FALSE, allow.cartesian = TRUE), pmats)
    jdf[, prob.joint := apply(as.matrix(.SD), 1, prod), .SDcols = grep("prob.", colnames(jdf), fixed = TRUE)]

    ## average over latents

    jdf[, lapply(.SD, mean), by = eval(subset(dhvm$nodes_data, latent == FALSE)$name),
        .SDcols = grep("prob.", colnames(jdf), fixed = TRUE)] |>
      as.data.frame()




  }


}

