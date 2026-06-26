
library(tidyverse)
library(tidygraph)
library(ggraph)


#' Internal helper function to build a cells tibble from pseudotime and identity_score matrix.
#'
#' It is useful to check parameter values and produce uniformized output.
#'
#' @param pseudotime a numeric vector of pseudotime value for each cell
#' @param identity_scores a numeric matrix of cell identity scores with number of
#'		row matching length(pseudotime). Each row give identity scores of the cell
#'		in each cluster. The higher the score the more likely the cell belong to the cluster.
#' @return a tibble with combined informations
get_cells <- function(pseudotime,identity_scores) {
	stopifnot("pseudotime and identity_scores should both be specified or not" = is.null(pseudotime) == is.null(identity_scores))
	if (is.null(pseudotime)) return(NULL)
	pseudotime <- as.numeric(pseudotime)
	identity_scores <- as.matrix(identity_scores)
	stopifnot("Dimensions of pseudotime and identity scores don't match" = nrow(identity_scores)==length(pseudotime))
	cells <- tibble(pseudotime,identity_scores) |>
		mutate(
			identity = max.col(identity_scores),
			identity_label = colnames(identity_scores)[identity] %||% as.character(identity),
			identity_score = identity_scores[cbind(seq_along(identity),identity)]
		)
	cells
}


#' Build complete lineage graph with associated metrics to infer lineage
#'
#' Fit linear models between identity_scores and pseudotime, and compute
#' metrics to help identify edges with target identity increasing with time
#' and source indentity decreasing with time.
#'
#' @param pseudotime a numeric vector of pseudotime value for each cell
#' @param identity_scores a numeric matrix of cell identity scores with number of
#'		row matching length(pseudotime). Each row give identity scores of the cell
#'		in each cluster.
#' @return a fully connected tbl_graph
#' @export
lineage_graph_build <- function(pseudotime,identity_scores) {
	#pseudotime <- runif(1000);identity_scores <- matrix(runif(10000),1000)
	cells <- get_cells(pseudotime,identity_scores)
	nodes <- cells |>
		group_by(identity,identity_label) |>
		summarize(
			med_pseudotime = median(pseudotime),
			lm_coefs = list(lsfit(pseudotime,identity_scores) |> coef())
		)

	edges <- expand_grid(from=nodes$identity,to=nodes$identity) |>
		filter(from!=to)
	g <- tbl_graph(nodes,edges) %>%
		activate(edges) %>%
		mutate(
			target_cells.source_id.slope     = map2_dbl(.N()$lm_coefs[to],from,~.x[2L,.y]),
			target_cells.source_id.intercept = map2_dbl(.N()$lm_coefs[to],from,~.x[1L,.y]),
			target_cells.target_id.slope     = map2_dbl(.N()$lm_coefs[to],to,~.x[2L,.y]),
			target_cells.target_id.intercept = map2_dbl(.N()$lm_coefs[to],to,~.x[1L,.y]),
			target_cells.branching.pseudotime = (target_cells.target_id.intercept - target_cells.source_id.intercept) / (target_cells.source_id.slope-target_cells.target_id.slope)
		) %>%
		activate(nodes)
	g
}

#' Prune a dense lineage graph to produce a lineage tree
#'
#' Filter edges of the provided graph. First remove edges where median source node
#' pseudotime is later than median target node pseudotime. Then rank edges according
#' to lineage metrics and keep for each target node the best ranking parent.
#'
#' @param g a lineage graph typically obtained with `lineage_graph_build`
#' @return the filter input graph
#' @export
lineage_graph_prune <- function(g) {
	g <- g %>%
		activate(edges) %>%
		filter(target_cells.source_id.slope < target_cells.target_id.slope) %>%
		arrange(desc(target_cells.branching.pseudotime)) %>%
		group_by(to) %>%
		slice_head(n=1) %>%
		activate(nodes)
}


#' Show pseudotime/identity relationship
#'
#' @param g a lineage graph typically obtained with `lineage_graph_build`
#' @param pseudotime a numeric vector of pseudotime value for each cell
#' @param identity_scores a numeric matrix of cell identity scores with number of
#'		row matching length(pseudotime). Each row give identity scores of the cell
#'		in each cluster. The higher the score the more likely the cell belong to the cluster.
#' @return a ggplot graph
#' @export
plot_lineage_incidence_matrix <- function(g,pseudotime=NULL,identity_scores=NULL) {
	#pseudotime <- runif(1000);identity_scores <- matrix(runif(10000),1000,dimnames=list(NULL,LETTERS[1:10]));g <- lineage_graph_build(pseudotime,identity_scores)
	g <- as_tbl_graph(g)
	cells <- get_cells(pseudotime,identity_scores)

	E <- g %>%
		activate(edges) %>%
		mutate(
			from_identity_label =.N()$identity_label[from],
			to_identity_label   =.N()$identity_label[to]
		) %>%
		as_tibble("edges")

	if (!is.null(cells)) {
		cells <- E %>%
			select(from_identity_label,to_identity_label) %>%
			left_join(cells,by=c("to_identity_label"="identity_label"),relationship = "many-to-many") %>%
			mutate(from_identity_score = identity_scores[cbind(seq_along(from_identity_label),match(from_identity_label,colnames(identity_scores)))]) %>%
			select(-identity_scores,-identity)
	}

	p <- ggplot(E) +
		facet_grid(to_identity_label~from_identity_label)
		if (!is.null(cells)) {
			p <- p +
				geom_point(aes(x=pseudotime,y=from_identity_score,color="SOURCE identity (should decrease with time)"),data=cells,size=0.3) +
				geom_point(aes(x=pseudotime,y=identity_score,color="TARGET identity (should increase with time)"),data=cells,size=0.3)
		} else {
			p <- p + scale_x_continuous(limits = g %>% activate(nodes) %>% pull(med_pseudotime) %>% range())
		}
		p <- p +
			geom_abline(aes(slope=target_cells.source_id.slope,intercept = target_cells.source_id.intercept),linewidth=1,color="red") +
			geom_abline(aes(slope=target_cells.target_id.slope,intercept = target_cells.target_id.intercept),linewidth=1,color="blue") +
			xlab("pseudotime / SOURCE-identity") +
			ylab("identity score / TARGET-identity") +
			labs(colour="") +
			ggtitle("Identity evolution of TARGET cells over time") +
			theme_bw() +
			theme(legend.position="top")

	p
}


#' Display lineage tree
#'
#' @param g a lineage graph
#' @return a gggraph plot
#' @export
plot_lineage_graph <- function(g) {
	g |>
		ggraph() +
		geom_edge_link(arrow = grid::arrow(angle=10,type="closed")) +
		geom_node_label(aes(label=identity_label))
}


