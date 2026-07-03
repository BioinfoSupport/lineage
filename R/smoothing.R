
#' Compute weight matrices to smooth gene expression signal along a pseudotime
#'
#' Weights are generated using a gaussian kernel.
#'
#' @param t cells pseudotime
#' @param at timepoints at which signal is evaluated
#' @param sd standard deviation of the gaussian kernel
#' @return a list of 2 matrices `fwd` and `rev`. That can be used to smooth gene expression signals.
#'         If m is a (gene x cells) matrix, the smoothed signal estimate at given timepoints is obtained
#'         with: m |> tcrossprod(w$fwd). The smoothed (gene x cells) expression matrix can be obtained
#'         with: m |> tcrossprod(w$fwd) |> tcrossprod(w$rev)
#' @export
#' @examples
#' t <- runif(100)
#' m <- matrix(runif(1000*100),ncol=100)
#' w <- waves_weights_1d(t)
#' m |> tcrossprod(w$fwd) |> tcrossprod(w$rev)
waves_weights_1d <- function(t,at=seq(min(t),max(t),length.out=30),sd=1) {
	t <- as.vector(t)
	at <- as.vector(at)
	w <- outer(at,t,FUN="-") |> abs() |> dnorm(sd=sd)
	list(fwd = w/rowSums(w),rev = t(w)/colSums(w))
}



#' Generate a uniform 2D grid around space defined by x,y points
#'
#' @param x x coordinate of the points in 2D space
#' @param y y coordinate of the points in 2D space
#' @param nx number of grid point on x-axis
#' @param ny number of grid point on y-axis
#' @param filter a boolean, when TRUE remove grid point that are not the direct
#'               nearest neighbor of a given point
#' @export
#' @import FNN
grid_2d <- function(x,y,nx=30,ny=30,filter=TRUE) {
	grid <- expand_grid(
		x = seq(min(x),max(x),length.out=nx),
		y = seq(min(y),max(y),length.out=ny)
	)
	if (filter) {
		knn <- FNN::get.knnx(as.matrix(grid),cbind(x,y),k=1)
		n <- tabulate(knn$nn.index,nbins = nrow(grid))
		grid <- grid[n>0,]
	}
	grid
}

#' Compute weight matrices to smooth gene expression signal along a 2D landscape
#'
#' Weights are generated using a gaussian kernel on euclidian distances to
#' nearest neighbors.
#'
#' @param x x-axis position of the cells in the 2D space
#' @param y y-axis position of the cells in the 2D space
#' @param grid a 2-column tibble with x,y columns defining the grid where signal is estimated.
#' @param k number of nearest neighboring grid points to consider for each cell
#' @param sd standard deviation of the gaussian kernel
#' @return a list of 2 matrices `fwd` and `rev`. That can be used to smooth gene expression signals.
#'         If m is a (gene x cells) matrix, the smoothed signal estimate at given timepoints is obtained
#'         with: m |> tcrossprod(w$fwd). The smoothed (gene x cells) expression matrix can be obtained
#'         with: m |> tcrossprod(w$fwd) |> tcrossprod(w$rev)
#' @export
#' @import Matrix
#' @import FNN
#' @examples
#'   x <- runif(1000); y <- runif(1000)
#'   m <- matrix(runif(1000*100),ncol=1000)
#'   w <- landscape_weights_2d(x,y,k=50,sd=0.05)
#'   M <- m |> tcrossprod(w$fwd) |> tcrossprod(w$rev)
#'   tibble(x,y,z=M[1,]) |> ggplot(aes(x=x,y=y,color=z)) + geom_point(size=3) + coord_equal()
landscape_weights_2d <- function(x,y,grid=grid_2d(x,y),k=5,sd=1) {
	# Look for k closest grid points of each cell
	knn <- FNN::get.knnx(as.matrix(grid),cbind(x,y),k=k)

	# Compute smooth weights
	knn$w <- dnorm(knn$nn.dist,mean = 0,sd = sd)

	W <- Matrix::sparseMatrix(
		i = knn$nn.index,
		j = row(knn$nn.index),
		x = as.vector(knn$w),
		dims = c(nrow(grid),length(x))
	)

	list(
		fwd = W / pmax(Matrix::rowSums(W),1e-6),   # Weight matrix to compute avg signal at each grid point
		rev = Matrix::t(W) / pmax(Matrix::colSums(W),1e-6) # Weight matrix to compute avg signal from grid values
	)
}



