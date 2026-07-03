

#' @importFrom tibble tibble
#' @importFrom dplyr mutate filter pull select left_join ungroup slice_max group_by
#' @importFrom tidyr expand_grid
#' @importFrom tidygraph tbl_graph activate .N select left_join ungroup slice_max group_by local_members
#' @importFrom magrittr %>%
#' @importFrom purrr map2_dbl
#' @importFrom stringr str_c
NULL



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
	if (is.null(colnames(identity_scores))) {
		colnames(identity_scores) <- as.character(seq(ncol(identity_scores)))
	}
	identity <- max.col(identity_scores)
	cells <- tibble(pseudotime,identity_scores) |>
		mutate(
			cell_id = seq_along(pseudotime),
			identity_label = factor(colnames(identity_scores),colnames(identity_scores))[identity],
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
	nodes <- cells %>%
		dplyr::group_by(identity_label) %>%
		dplyr::summarise(
			min_pseudotime = min(pseudotime),
			med_pseudotime = median(pseudotime),
			max_pseudotime = max(pseudotime),
			lm_coefs = list(lsfit(pseudotime,identity_scores) |> coef())
		)

	edges <- tidyr::expand_grid(from=nodes$identity_label,to=nodes$identity_label) %>%
		dplyr::filter(from!=to)
	g <- tidygraph::tbl_graph(nodes,edges) %>%
		tidygraph::activate(edges) %>%
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

#' Prune a lineage graph to keep one parent per node
#'
#' Filter edges of the provided graph. First remove edges where median source node
#' pseudotime is later than median target node pseudotime. Then rank edges according
#' to lineage metrics and keep for each target node the best ranking parent.
#'
#' @param g a lineage graph typically obtained with `lineage_graph_build`
#' @param n number of parent to keep for each node
#' @return the filter input graph
#' @export
lineage_graph_prune <- function(g,n=1L) {
	g <- g %>%
		activate(edges) %>%
		filter(target_cells.source_id.slope < target_cells.target_id.slope) %>%
		group_by(to) %>%
		slice_max(order_by=target_cells.branching.pseudotime,n=n,with_ties = FALSE) %>%
		ungroup() %>%
		activate(nodes)
	g
}

#' Compute ancestor table
#'
#' @param g a graph
#' @return a sparse logical matrix where each row list all ancestors of a node
#' @return a 2-column tibble of ancestors
#' @export
lineage_ancestor_tbl <- function(g) {
	g %>%
		mutate(ancestor_label = local_members(mode="in",mindist = 0,order = +Inf) %>% map(~.N()$identity_label[.])) %>%
		as_tibble("nodes") %>%
		select(node_label=identity_label,ancestor_label) %>%
		unnest_longer(ancestor_label)
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
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_point facet_grid geom_abline xlab ylab labs ggtitle theme theme_bw scale_x_continuous
#' @importFrom tidygraph as_tbl_graph as_tibble
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
			mutate(
				facet_x = str_c(to_identity_label,"\nTARGET cells"),
				facet_y = str_c("vs ",from_identity_label,"\nSOURCE id")
			) %>%
			select(-identity_scores)
	}

	p <- E %>%
		mutate(
			facet_x = str_c(to_identity_label,"\nTARGET cells"),
			facet_y = str_c("vs ",from_identity_label,"\nSOURCE id")
		) %>%
		ggplot() +
			facet_grid(facet_y ~ facet_x)
		if (!is.null(cells)) {
			p <- p +
				geom_point(aes(x=pseudotime,y=from_identity_score,color="SOURCE identity (should decrease with time)"),data=cells,size=0.3) +
				geom_point(aes(x=pseudotime,y=identity_score,color="TARGET identity (should increase with time)"),data=cells,size=0.3)
		} else {
			xlim <- range(
				g %>% activate(nodes) %>% pull(min_pseudotime),
				g %>% activate(nodes) %>% pull(max_pseudotime)
			)
			p <- p + scale_x_continuous(limits = xlim)
		}
		p <- p +
			geom_abline(aes(slope=target_cells.source_id.slope,intercept = target_cells.source_id.intercept),linewidth=1,color="red") +
			geom_abline(aes(slope=target_cells.target_id.slope,intercept = target_cells.target_id.intercept),linewidth=1,color="blue") +
			xlab("pseudotime") +
			ylab("identity score") +
			labs(colour="") +
			ggtitle("Evolution of identity scores over pseudotime on TARGET cells") +
			theme_bw() +
			theme(legend.position="top")
	p
}

#' Display lineage tree
#'
#' @param g a lineage graph
#' @return a gggraph plot
#' @export
#' @importFrom ggraph ggraph geom_edge_link geom_node_label
plot_lineage_graph <- function(g) {
	g |>
		ggraph() +
			geom_edge_link(arrow = grid::arrow(angle=10,type="closed")) +
			geom_node_label(aes(label=str_c(seq_along(identity_label),"-",identity_label)))
}

#' Compute cells coordinates in a lineage graph from their identity_matrix
#'
#' @param g a graph
#' @param identity_scores a numeric matrix of cell identity scores with number of
#'		row matching length(pseudotime). Each row give identity scores of the cell
#'		in each cluster.
#' @return a tibble
#' @export
lineage_coords <- function(g,pseudotime,identity_scores) {
	#set.seed(123);pseudotime <- runif(1000);identity_scores <- matrix(runif(10000),1000,dimnames=list(NULL,LETTERS[1:10]));g <- lineage_graph_build(pseudotime,identity_scores) %>% lineage_graph_prune()
	ancestors <- lineage_ancestor_tbl(g)
	cells <- get_cells(pseudotime,identity_scores)
	lineages <- cells %>%
		inner_join(ancestors,by = join_by(identity_label==ancestor_label),relationship = "many-to-many") %>%
		dplyr::rename(lineage_label=node_label)
	lineages
}




