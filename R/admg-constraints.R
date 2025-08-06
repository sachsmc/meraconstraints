admg_find_constraints <- function(admg) {

  L <- NULL
  admgord <- topologicalOrder(admg)
  for(v in admgord) {

    prebarv <- admgord[seq(1, match(v, admgord))]
    mbv <- union(markov_blanket(admg, v), v)
    indepcheck <- setdiff(prebarv, mbv)

    if(length(indepcheck) > 0) {
      L <- c(L, list(indep = list(admg$vnames[v], admg$vnames[indepcheck]),
           given = admg$vnames[setdiff(mb(admg, v, A = prebarv), v)]))
    }

  }
  L

}

markov_blanket <- function(graph, v) {

  for(ds in districts(graph)) {
    if(v %in% ds) {
      S <- ds
      break
    }
  }
  union(pa(graph, S), setdiff(S, v))


}






