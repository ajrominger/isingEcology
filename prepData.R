library(raster)
library(socorro)
library(Matrix)
library(parallel)

setwd('~/Dropbox/Research/isingEcology')


## function to generate data for testing/fitting/understanding ising model 
## as applied to spatial occupancy
#' @param path the path to the raw data file to be loaded and processed
#' @param writePath the path where processed data re to be written
#' @return writes two datafiles:
#' \describe{
#'   \item{\code{*Lambda}}{a data.frame describing the Ising latice by columns `spp`, `cell`, `spin`}
#'   \item{\code{*corFun}}{the empirical correlation function defined as c(d) = <x_i * x_j> - <x>^2}
#' }

prepDataIsing <- function(path, writePath) {
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
    
    browser()
    LambdaMat <- tidy2mat(cellz, x$spp, x$count)
    LambdaMat[LambdaMat > 1] <- 1 # just incase some cells have more than 1 individ
    LambdaMat[LambdaMat == 0] <- -1
    Lambda <- data.frame(spp = rep(colnames(LambdaMat), each = nrow(LambdaMat)), 
                         cell = as.numeric(rep(rownames(LambdaMat), ncol(LambdaMat))), 
                         spin = as.vector(LambdaMat), 
                         stringsAsFactors = FALSE)
    
    
    ## add xy coords for each cell
    Lambda <- cbind(Lambda, xyFromCell(r, Lambda$cell))
    
    ## create adjacency matrix that indicates adjacent cells
    adj <- bandSparse(ncell(r), k = c(-ncol(r), -1, 0, 1, ncol(r)))
    adj[(1:nrow(r)) * ncol(r), (1:(nrow(r))-1) * ncol(r) + 1] <- FALSE
    adj[(1:(nrow(r) - 1)) * ncol(r) + 1, (1:nrow(r)) * ncol(r)] <- FALSE
    
    ## `ego` contains the ego graph for successive distances from the focal node
    ego <- adj
    
    ## `theseCells` contains a subset of cells from which to grow out the neighborhoods
    ## selected to be toward the center of the plot to avoid boundary effects
    theseCells <- Lambda$cell[Lambda$x <= 0.75 * max(Lambda$x) &
                                  Lambda$x >= 0.25 * max(Lambda$x) &
                                  Lambda$y <= 0.75 * max(Lambda$y) & 
                                  Lambda$x >= 0.25 * max(Lambda$y)]
    theseCells <- unique(theseCells)
    if(length(theseCells) > 200) theseCells <- sample(theseCells, 200)
    
    ## calculate empirical correlation function defined as c(d) = <x_i * x_j> - <x>^2
    maxD <- round(sqrt(ncell(r)) / 2)
    cfun <- lapply(1:maxD, function(i) {
        ## use the current ego matrix to find neighbors
        groups <- which(ego >= 1, arr.ind = TRUE)
        groups <- groups[groups[, 2] %in% theseCells, ]

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
    
    ## write output
    outPath <- file.path(writePath, gsub('.*/|.csv', '', path))
    
    write.csv(Lambda, file = paste0(outPath, '_Lambda.csv'), row.names = FALSE)
    write.csv(cfun, file = paste0(outPath, '_corFun.csv'), row.names = FALSE)
}

## a faster (parallel) version of tapply
tapply2 <- function(X, INDEX, FUN = NULL, ..., default = NA, simplify = TRUE) {
    nI <- length(INDEX)
    namelist <- lapply(INDEX, levels)
    extent <- lengths(namelist, use.names = FALSE)
    cumextent <- cumprod(extent)
    storage.mode(cumextent) <- "integer"
    ngroup <- cumextent[nI]
    group <- as.integer(INDEX[[1L]])
    for (i in 2L:nI) {
        group <- group + cumextent[i - 1L] * (as.integer(INDEX[[i]]) - 1L)
    }
    levels(group) <- as.character(seq_len(ngroup))
    class(group) <- "factor"
    ans <- split(X, group)
    names(ans) <- NULL
    index <- as.logical(lengths(ans))
    
    ans <- mclapply(X = ans[index], mc.cores = 6, FUN = FUN, ...)
    ans <- unlist(ans, recursive = FALSE, use.names = FALSE)
    
    ansmat <- array(vector(typeof(ans)), dim = extent, dimnames = namelist)
    
    ansmat[index] <- ans
    
    return(ansmat)
}

## set tapply to parallel version
foo <- base::tapply
tapply <- tapply2

## run and then re-set tapply
prepDataIsing('../data/stri/UCSC.csv', 'data')
prepDataIsing('../data/stri/PASO.csv', 'data')

tapply <- foo
