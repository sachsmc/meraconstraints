# meraconstraints 0.3.0

* All per-district lists are now named by their district members (e.g. `"X_Y"`): the elements of `dhvm_district_params()`, the `relation`, `probform`, `probstarform` and `problevels` elements of `dhvm_linear_relations()`, the `Vrep` and `Hrep` elements of `dhvm_relations_to_polytope()`, and the results of `prob_simplex_screen()` and `prob_simplex_extremes()`. Previously only the components of `meraconstraints()` output were named, so `functional_relations$problevels[["X_Y"]]` returned `NULL`.
* `combinations()` is no longer exported.
* screening of potentially nontrivial constraints improved (thanks to anonymous reviewer)
* New `nontrivial_screen()` identifies exactly the nontrivial constraints of each district (Proposition 4 of the paper). An inequality is trivial if and only if it is a nonnegativity constraint, possibly rewritten using the equality constraints; an equality is trivial if and only if it follows from the conditional independencies. `meraconstraints()` now uses it, stores the result in `nontrivial` (with `trivial_if_false` kept as an alias), and no longer flags trivial constraints such as those of the front door model. `prob_simplex_screen()` and `prob_simplex_extremes()` are kept but superseded.
* `dhvm_relations_to_polytope()` also returns `HrepQ`, the H-representations in exact rational arithmetic, which `nontrivial_screen()` uses to decide exactly which vertices lie on each facet.
* dedupe columns of V-rep

# meraconstraints 0.2.1

* Bugfix in check independence when there are only 2 variables. 

# meraconstraints 0.2.0

* Initial GitHub release.
