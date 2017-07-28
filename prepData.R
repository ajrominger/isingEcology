library(raster)
library(socorro)
library(Matrix)
library(parallel)

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
    LambdaMat[LambdaMat > 1] <- 1 # just incase some cells have more than 1 individ
    LambdaMat[LambdaMat == 0] <- -1
    Lambda <- data.frame(spp = rep(colnames(LambdaMat), each = nrow(LambdaMat)), 
                         cell = rep(rownames(LambdaMat), ncol(LambdaMat)), 
                         spin = as.vector(LambdaMat), 
                         stringsAsFactors = FALSE)
    
    ## create adjacency matrix that indicates adjacent cells
    adj <- bandSparse(ncell(r), k = c(-ncol(r), -1, 0, 1, ncol(r)))
    adj[(1:nrow(r)) * ncol(r), (1:(nrow(r))-1) * ncol(r) + 1] <- FALSE
    adj[(1:(nrow(r) - 1)) * ncol(r) + 1, (1:nrow(r)) * ncol(r)] <- FALSE
    
    ## `ego` contains the ego graph for successive distances from the focal node
    ego <- adj
    
    ## calculate empirical correlation function defined as c(d) = <x_i * x_j> - <x>^2
    maxD <- round(sqrt(ncell(r)) / 2)
    cfun <- lapply(1:maxD, function(i) {
        ## use the current ego matrix to find neighbors
        groups <- which(ego >= 1, arr.ind = TRUE)
        print(head(groups))
        
        ## calculate cor fun for each spp
        out <- mclapply(unique(Lambda$spp), mc.cores = 6, FUN = function(sp) {
            dat <- Lambda[Lambda$spp == sp, ]
            spinMean <- mean(dat$spin)
            
            cr <- tapply(dat$spin[match(groups[, 1], dat$cell)], groups[, 2], 
                         function(spin) {
                             temp <- outer(spin, spin, '*')
                             mean(temp[lower.tri(temp)])
                         })
            cr <- cr - spinMean^2
            
            return(c(m = mean(cr), 
                     lo = quantile(cr, 0.025, names = FALSE), 
                     hi = quantile(cr, 0.975, names = FALSE)))
        })
        out <- do.call(rbind, out)
        
        
        ## update the ego matrix in the parent frame for the next iteration
        ego <<- ego %*% adj + ego
        
        ## return output
        return(out)
    })
    
    cfun <- data.frame(scale = rep(1:maxD, each = length(unique(Lambda$spp))), 
                       spp = unique(Lambda$spp),
                       do.call(rbind, cfun))
    
    
    return(list(Lambda = Lambda, cfun = cfun))
}

isingUCSC <- prepDataIsing('../data/stri/UCSC.csv')
isingPASO <- prepDataIsing('../data/stri/PASO.csv')
