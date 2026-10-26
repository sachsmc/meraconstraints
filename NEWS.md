# meraconstraints 0.3.0

* All per-district lists are now named by their district members (e.g. `"X_Y"`): the elements of `dhvm_district_params()`, the `relation`, `probform`, `probstarform` and `problevels` elements of `dhvm_linear_relations()`, the `Vrep` and `Hrep` elements of `dhvm_relations_to_polytope()`, and the results of `prob_simplex_screen()` and `prob_simplex_extremes()`. Previously only the components of `meraconstraints()` output were named, so `functional_relations$problevels[["X_Y"]]` returned `NULL`.
* `combinations()` is no longer exported.

# meraconstraints 0.2.1

* Bugfix in check independence when there are only 2 variables. 

# meraconstraints 0.2.0

* Initial GitHub release.
