## Deriving Complete Constraints in Hidden Variable Models

Hidden variable graphical models can sometimes imply constraints on the observable distribution that are more complex than simple conditional independence relations. These observable constraints can allow for falsification of assumptions of the model that would otherwise be untestable due to the unobserved variables and can be used to constrain estimation procedures to improve efficiency. Knowing the complete set of observable constraints is thus ideal, but this can be difficult to determine in many settings.  In models with categorical observed variables and a joint distribution that is completely characterized by linear relations to the unobservable response function variables, we develop a systematic method for deriving the complete set of observable constraints. This package implements that method and provides some examples in the vignette. 

This is the companion R package to a paper soon to appear on Arxiv. 

This package is in early development. The documentation is incomplete, and the API is subject to change. 

### Installation

For now it is only available on github.

```r
remotes::install_github("sachsmc/meraconstraints")
```

### Basic usage

See vignette. 

