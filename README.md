# lineage

`lineage` is an R package to help infer lineage graph from single cell transcriptomic 
data where a pseudotime and a celltype score matrix is available.


# Installation

``` r
devtools::install_github("BioinfoSupport/lineage/lineage")
```


# Usage

``` r
library(lineage)

# Generate a random pseudotime and identity score matrix
set.seed(123)
pseudotime <- runif(1000)
identity_scores <- matrix(runif(10000),1000,dimnames=list(NULL,LETTERS[1:10]))

# Find lineage
g <- lineage_graph_build(pseudotime,identity_scores) |>
  lineage_graph_prune()
plot_lineage_graph(g)
```


