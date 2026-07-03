#include <Rcpp.h>
#include <queue>
#include <vector>
#include <utility>
#include <functional>
using namespace Rcpp;


// [[Rcpp::export]]
List edge_mst_hclust_impl(IntegerVector from, IntegerVector to, NumericVector weight, int n_nodes) {
	int m = from.size();
	if (to.size() != m || weight.size() != m) stop("from, to and weight must have the same length");

	// Define and fill priority queue
	std::priority_queue<
		std::pair<double, int>,
		std::vector< std::pair<double, int> >,
		std::greater< std::pair<double, int> >
	> pq;
	for (int e = 0; e < m; ++e) pq.push(std::make_pair(weight[e], e));

	// union-find storage
	IntegerVector node_parent(n_nodes);
	for(size_t i=0;i<node_parent.size();++i) node_parent[i]=i;
	IntegerVector node_rank(n_nodes, NA_INTEGER);
	IntegerVector node_edge(n_nodes, NA_INTEGER);

	int r = 0;
	while (!pq.empty()) {
		int e = pq.top().second;
		pq.pop();
		int a = from[e] - 1;
		int b = to[e] - 1;
		while (node_parent[a] != a) { a = node_parent[a]; };
		while (node_parent[b] != b) { b = node_parent[b]; };
		if (a == b) continue;
		node_parent[a] = b;
		node_rank[a] = ++r;
		node_edge[a] = e;
	}

	return List::create(node_rank,node_parent+1,node_edge+1);
}


/*** R

node_mst_hierarchy <- function(weights = NULL) {
	if (!.graph_context$free() && .graph_context$active() != "nodes") {
		cli::cli_abort("This call requires nodes to be active", call. = FALSE)
	}
	weights <- (enquo(weights) |> rlang::eval_tidy(.E())) %||% rep(1,graph_size())
	edge_mst_hclust_impl(.E()$from, .E()$to, as.numeric(weights), graph_order()) |>
		as_tibble(.name_repair=~c("rank","parent","edge_idx"))
}

# Create a complete graph with random weights
set.seed(1234)
g <- create_complete(30) |>
	activate(edges) |>
	mutate(weight=runif(graph_size())) |>
	activate(nodes) |>
	mutate(node_hmst(weights = weight),)


# Display original graph
ggraph(g) + geom_node_point() + geom_edge_link()

# Convert to hclust
H <- g |>
	mutate(height = .E()$weight[edge_idx]) |>
	as_tibble("nodes") |>
	rowid_to_column() |>
	select(from=parent,to=rowid,rank,height,edge_idx) |>
	as_tbl_graph() |>
	activate(edges) |>
	filter(!edge_is_loop()) |>
	activate(nodes)

H |>
	ggraph("dendrogram") +
	geom_node_point() +
	geom_edge_link(aes(label=rank),arrow = grid::arrow(angle=10,type="closed")) +
	geom_node_point(aes(filter=node_is_root()),color="red",size=3)
*/
