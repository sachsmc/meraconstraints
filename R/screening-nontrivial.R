#' Identify nontrivial constraints
#'
#' Exact screen of the constraints of each district (Proposition 6 of the
#' paper). An equality constraint is trivial if it holds for every strictly
#' positive distribution that satisfies the conditional independencies implied
#' by the graph. These are the linear combinations of the normalization
#' constraints and of the constraints saying that, for every set S of members
#' that contains all members that are ancestors of its members, the margin of S
#' under the intervention does not depend on the external parents of the other
#' members. An inequality constraint is trivial if it is implied by the
#' probability simplex together with the equality constraints of its district.
#' This is the case if and only if some probability is zero at every vertex of
#' the district polytope on its facet, and then the constraint is a
#' nonnegativity constraint written in a different form. A nontrivial inequality
#' constraint is violated by some strictly positive distribution that satisfies
#' the conditional independencies and every other constraint.
#'
#' Unlike [prob_simplex_screen()], the result does not depend on how the rows
#' of the H-representation are written, and it flags no trivial constraints.
#'
#' @param dhvm An object of class `dhvm`, see [dhvm()], with c-degree at most 1
#' @param drels A list of relations as returned by [dhvm_linear_relations()]
#'   applied to `dhvm`
#' @param dpoly The polytopes as returned by [dhvm_relations_to_polytope()]
#'   applied to `drels`. Which vertices lie on which facet is computed in exact
#'   rational arithmetic from its `HrepQ` element, or, if that is missing, from
#'   `Hrep` converted with [rcdd::d2q()].
#' @param tol Numeric tolerance for the ranks in the screen of the equality
#'   constraints
#' @returns A named list of logical vectors, one for each district, with one
#'   element per row of the H-representation. For an inequality constraint it
#'   is `TRUE` if the constraint is nontrivial. For the equality constraints,
#'   the rows marked `TRUE` form a basis of the nontrivial ones: together with
#'   the conditional independencies they imply all equality constraints of the
#'   district.
#' @export
#' @examples
#' xab <- dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B))
#' rels <- dhvm_linear_relations(xab)
#' poly <- dhvm_relations_to_polytope(rels)
#' nt <- nontrivial_screen(xab, rels, poly)
#' ## the 8 CHSH inequalities
#' hrep_to_text(poly$Hrep[[1]], rels$probstarform[[1]])[nt[[1]]]
nontrivial_screen <- function(dhvm, drels, dpoly, tol = 1e-9) {

  obs <- subset(dhvm$nodes_data, latent == FALSE)$name
  anc <- lapply(setNames(obs, obs), \(v) intersect(ancestors(dhvm$caugi_graph, nodes = v), obs))

  checks <- lapply(seq_along(dpoly$Hrep), \(j) {
    H  <- dpoly$Hrep[[j]]
    Hq <- if (is.null(dpoly$HrepQ)) rcdd::d2q(H) else dpoly$HrepQ[[j]]
    eq <- H[, 1] == 1
    out <- logical(nrow(H))
    out[!eq] <- screen_inequalities(Hq[!eq, , drop = FALSE], dpoly$Vrep[[j]])
    out[eq]  <- screen_equalities(H[eq, , drop = FALSE], dhvm, drels, j, anc, tol)
    out
  })
  names(checks) <- names(dpoly$Hrep)
  checks

}

## An inequality row is trivial iff some entry of p is zero at every vertex on its facet.
## Hq holds rational rows (0, b, -A) of b - A p >= 0, so the slack of vertex v is
## b + (-A) v, computed exactly.
screen_inequalities <- function(Hq, Vrep) {

  if (nrow(Hq) == 0) return(logical(0))
  Vq <- Vrep[, -(1:2), drop = FALSE]
  if (!is.character(Vq)) Vq <- rcdd::d2q(Vq)
  Av <- rcdd::qmatmult(Hq[, -(1:2), drop = FALSE], t(Vq))
  slack <- rcdd::qpq(rep(Hq[, 2], times = ncol(Av)), as.vector(Av))
  on <- matrix(slack == "0", nrow = nrow(Hq))   # facet-vertex incidence
  V <- rcdd::q2d(Vq)                              # 0/1, so the product below is exact
  rowSums((on %*% V) == 0) == 0                   # every entry positive at some vertex on the facet

}

## Equality rows that, added one at a time to (N) and (O), increase the rank.
screen_equalities <- function(He, dhvm, drels, j, anc, tol) {

  if (nrow(He) == 0) return(logical(0))
  lev <- district_levels(drels, j)
  mem <- drels$district_parms$poterms[[j]]$out
  ext <- drels$district_parms$poterms[[j]]$do
  ext_pa <- lapply(setNames(mem, mem), \(m) intersect(obs_parents(dhvm, m), ext))

  rk <- \(m) if (is.null(m) || nrow(m) == 0) 0L else qr(m, tol = tol)$rank
  E  <- -He[, -(1:2), drop = FALSE]
  N  <- indicator_rows(lev, ext)
  O  <- do.call(rbind, lapply(closed_member_sets(mem, anc), margin_rows,
                             lev = lev, ext = ext, ext_pa = ext_pa))
  M  <- rbind(N, O)
  if (rk(rbind(M, E)) != rk(E))
    stop("normalization and margin constraints do not hold on the polytope of district ",
         names(drels$probform)[j])

  new <- logical(nrow(E))
  for (r in seq_len(nrow(E))) {
    if (rk(rbind(M, E[r, ])) > rk(M)) {
      M <- rbind(M, E[r, ])
      new[r] <- TRUE
    }
  }
  new

}

## Levels of the variables of district j, one row per entry of p, in the order of the
## columns of the H-representation.
district_levels <- function(drels, j) {

  pn  <- names(unlist(drels$probform[[j]]))
  lev <- drels$problevels[[j]]
  lev[match(pn, do.call(paste, c(as.list(lev), sep = "."))), , drop = FALSE]

}

## Normalization: one indicator row for each value of cols.
indicator_rows <- function(lev, cols) {

  if (length(cols) == 0) return(matrix(1, nrow = 1, ncol = nrow(lev)))
  key <- do.call(paste, c(as.list(lev[, cols, drop = FALSE]), sep = "."))
  t(vapply(unique(key), \(k) as.numeric(key == k), numeric(nrow(lev))))

}

## Nonempty sets of members that contain every member that is an ancestor of one of
## their members.
closed_member_sets <- function(mem, anc) {

  sets <- unlist(lapply(seq_along(mem), \(k) combn(mem, k, simplify = FALSE)), recursive = FALSE)
  Filter(\(S) all(intersect(unlist(anc[S]), mem) %in% S), sets)

}

## The margin of S depends on the external parents only through those of S: one row for
## each value of (S, external parents) in which the other external parents are not all at
## their lowest level, minus the matching row with them at their lowest level.
margin_rows <- function(lev, S, ext, ext_pa) {

  early <- intersect(ext, unlist(ext_pa[S]))
  late  <- setdiff(ext, early)
  if (length(late) == 0) return(NULL)
  key  <- \(cols) do.call(paste, c(as.list(lev[, cols, drop = FALSE]), sep = "."))
  cell <- key(c(S, ext))
  same <- key(c(S, early))
  base <- Reduce(`&`, lapply(late, \(w) lev[[w]] == min(lev[[w]])))
  t(vapply(unique(cell[!base]), \(cl) {
    i0 <- which(cell == cl)[1]
    as.numeric(cell == cl) - as.numeric(same == same[i0] & base)
  }, numeric(nrow(lev))))

}
