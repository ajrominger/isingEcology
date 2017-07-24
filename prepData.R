library(raster)
library(socorro)
library(igraph)

setwd('~/Dropbox/Research/isingEcology')


## function to generate data for testing/fitting/understanding ising model 
## as applied to spatial occupancy
#' @param path the path to the raw data file to be loaded and processed
#' @return a list with components:
#' \describe{
#'   \item{\code{Lambda}}{a data.frame describing the Ising latice by columns `spp`, `cell`, `spin`}
#'   \item{\code{edgeList}}{a matrix describing the edges connecting neighbors in the latice}
#'   \item{\code{corFun}}{the empirical correlation function defined as c(d) = <x_i * x_j> - <x>^2}
#' }

prepDataIsing <- function(path) {
    browser()
    ## take only most recent census
    x <- read.csv(path, as.is = TRUE)
    x <- x[x$year == max(x$year), ]
    
    ## make raster for plot, for PASO and UCSC, min distance was approx 2m, so
    ## i'm just going to use a resolution of 1m
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## create a df with columns of spp, cellID, spin
    cellz <- cellFromXY(r, xy = x[, c('x', 'y')])
    cellz <- factor(cellz, levels = 1:ncell(r))
    LambdaMat <- tidy2mat(cellz, x$spp, x$count)
    LambdaMat[LambdaMat > 1] <- 1 # just incase some cells have more than 1 individ...this should be rare
    LambdaMat[LambdaMat == 0] <- -1
    Lambda <- data.frame(spp = rep(colnames(LambdaMat), each = nrow(LambdaMat)), 
                         cell = rep(rownames(LambdaMat), ncol(LambdaMat)), 
                         spin = as.vector(LambdaMat))
    
    ## create edge list that indicates IDs of adjacent cells; use adjacent
    edges <- adjacent(r, 1:ncell(r), directions = 4, sorted = TRUE)
    
    ## calcualte correlation function for each spp: c(d) = < x_i * x_j > - < x >< x >
    bla <- traverceEL(edges, 2)
    head(bla[order(bla[, 1]), ])
    
}

## function to traverse an edge list, returning all nodes within `d` links
## as applied to spatial occupancy
##
## DOESN'T WORK, SEE BELOW USING `EGO`
##
#' @param el the edge list to be traversed
#' @param d the desired number of edges between nodes
#' @return a new edge list, now with all nodes that are `d` away from the nodes in column `from`
traverceEL <- function(el, d) {
    if(d == 1) {
        return(el)
    } else {
        newEL <- lapply(1:nrow(el), function(i) cbind(el[i, 1], el[el[i, 2] == el[, 1], 2]))
        newEL <- do.call(rbind(newEL))
        newEL <- newEL[newEL[, 1] != newEL[, 2]]
        newEL <- newEL[!duplicated(newEL), ]
        newEL <- rbind(el, newEL)
        return(traverceEL(newEL, d - 1))
    }
}

bla <- traverceEL(edges, 2)


matrix(1:9, nrow = 3, byrow = TRUE)
el <- matrix(c(1, 2,
               1, 4, 
               2, 1, 
               2, 3, 
               2, 5, 
               3, 2, 
               3, 6, 
               4, 1,
               4, 5, 
               4, 7, 
               5, 2, 
               5, 4, 
               5, 6,
               5, 8,
               6, 3, 
               6, 5, 
               6, 9,
               7, 4, 
               7, 8, 
               8, 7, 
               8, 5, 
               8, 9, 
               9, 6, 
               9, 8), ncol = 2, byrow = TRUE)

foo <- graph_from_edgelist(el)
plot(foo)
as_edgelist(graph_from_adj_list(ego(foo, 2)))



prepDataIsing('../data/stri/UCSC.csv')