
# Code parking


library(tidyverse)
library(tidygraph)
library(ggraph)


#' Return number of leaf seen so far in a DFS search, to build tree layout
node_rank_dfs_leaf <- function(g) {
	#pseudotime <- runif(1000);identity_scores <- matrix(runif(10000),1000)
	g <- lineage_graph_build(pseudotime,identity_scores) %>%
		lineage_graph_prune()
	plot_lineage(g)

	g <- g %>%
		activate(nodes) %>%
		select(!any_of("dfs_visited_leaf")) %>%
		mutate(
			is_sink = node_is_sink(),
			dfs_rk = dfs_rank(mode="out",unreachable = TRUE,root=which.min(centrality_degree(mode="in")))
		)

	g <- g %>%
		left_join(by="identity",
							select(g,identity,dfs_rk,is_sink) %>%
								as_tibble("nodes") %>%
								arrange(dfs_rk) %>%
								mutate(dfs_visited_leaf = cumsum(c(FALSE,is_sink)[-length(identity)])) %>%
								select(identity,dfs_visited_leaf)
		)

	ggraph(g,layout="manual",x=med_pseudotime,y=dfs_visited_leaf) +
		geom_edge_link(arrow = grid::arrow(angle=10,type="closed")) +
		geom_node_label(aes(label=identity_label))


	ggraph(g) +
		geom_edge_link(arrow = grid::arrow(angle=10,type="closed")) +
		geom_node_label(aes(label=identity_label))
}



predict_pseudotime <- function(counts,model) {
	w <- readRDS(model)
	w <- w[names(w) %in% rownames(counts)]
	counts <- counts[names(w),]
	logcounts <- scuttle::normalizeCounts(counts,center.size.factors=FALSE,size.factors=colSums(counts)/100)
	as.vector(w %*% logcounts)
}

