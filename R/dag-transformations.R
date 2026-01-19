#' Merge latent variables within districts
#'
#' Results in a linear dag that has only has 1 latent variable per district
#'
#' @param dhvm A object representing a dag of as returned by [dhvm()]
#' @returns A `dhvm` object representing the modified dag
#' @export
#' @examples
#' duargr <- dhvm(caugi(A %-->% B+E, B %-->% C, C %-->% D, D %-->% F, E %-->% C,
#' U1 %-->% B+D, U2 %-->% A+C+E, U3 %-->% D+F))
#' dhvm_merge_latents(duargr)
#'
dhvm_merge_latents <- function(dhvm) {

  latent_count <- dhvm_c_degree(dhvm)
  dists <- dhvm_districts(dhvm)
  unobs <- subset(dhvm$nodes_data, latent == TRUE)$name

  tomerge <- which(latent_count > 1)

  newcg <- dhvm$caugi_graph

  for(i in seq_along(tomerge)) {

    thisto <- setdiff(dists[[tomerge[i]]], unobs)
    thisfrom <- intersect(dists[[tomerge[i]]], unobs)[1]
    thisremove <- intersect(dists[[tomerge[i]]], unobs)[-1]

    newcg <- add_edges(newcg, from = rep(thisfrom, length(thisto)),
              edge = rep("-->", length(thisto)), to = thisto)

    newcg <- remove_nodes(newcg, name = thisremove)

  }

  dhvm(build(newcg))

}
