## install packages
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}
packages <- c("RandomFields","raster","devtools","landscapetools","RColorBrewer",)
ipak(packages)
## install blockCV
library(devtools)
devtools::install_github("rvalavi/blockCV")

## load 
library(RandomFields)
library(raster)
library(blockCV)
library(landscapetools)
library(RColorBrewer)

## create dir to save csv files
# 2 sampling strategies
dir.create("coords")
dir.create("coords/coords_even")
dir.create("coords/coords_clmp")
# 4 combinations: data type of response variable and distribution of sample points
dir.create("samples")
dir.create("samples/samples_dsc_even")
dir.create("samples/samples_dsc_clmp")
dir.create("samples/samples_ctn_even")
dir.create("samples/samples_ctn_clmp")
# 3 block sizes: splitting based on coords, range, not considering the response variable's values
dir.create("blockcv")
dir.create("blockcv/blockcv_l_even")
dir.create("blockcv/blockcv_m_even")
dir.create("blockcv/blockcv_s_even")
dir.create("blockcv/blockcv_l_clmp")
dir.create("blockcv/blockcv_m_clmp")
dir.create("blockcv/blockcv_s_clmp")
##
dir.create("plots")

## init
n.sims <- 100

## set landscape extent
gridsize = c(50L, 50L)
Xvec <- seq(0, 50, len = gridsize[1])
Yvec <- seq(0, 50, len = gridsize[2])
grd <- expand.grid(Y = Yvec, X = Xvec)
lon <- grd$X
lat <- grd$Y

for(i in 1:n.sims){

  ############## STEP 1 - Simulate the Landscapes #################
  RFoptions(seed=8*i-7)
  expCov <- RMexp(var = 5, scale = 5)
  x.1 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.1 <- scale(x.1)
  
  RFoptions(seed=8*i-6)
  expCov <- RMexp(var = 15, scale = 5)
  x.2 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.2 <- scale(x.2)
  
  RFoptions(seed=8*i-5)
  expCov <- RMgauss(var = 5, scale = 15) 
  x.3 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.3 <- scale(x.3)
  
  p.d <- 0.05
  t.ps <- -floor(min(x.2, x.3))
  x.4 <- rep(1, prod(gridsize))
  x.4[((x.2 + t.ps)/(x.3 + t.ps)) < quantile(((x.2 + t.ps)/(x.3 + t.ps)), (1 - p.d))] <- 0
  
  x.5 <- (x.1 + x.2 + x.3 + x.2*x.3)
  x.5 <- scale(x.5)

  RFoptions(seed=8*i-4)
  expCov <- RMexp(var = 5, scale = 5)
  x.6 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.6 <- scale(x.6)

  RFoptions(seed=8*i-3)
  expCov <- RMexp(var = 5, scale = 5)
  x.7 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.7 <- scale(x.7)

  RFoptions(seed=8*i-2)
  expCov <- RMexp(var = 5, scale = 5)
  x.8 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.8 <- scale(x.8)

  RFoptions(seed=8*i-1)
  expCov <- RMgauss(var = 5, scale = 15)
  x.9 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.9 <- scale(x.9)

  RFoptions(seed=8*i)
  expCov <- RMgauss(var = 5, scale = 15)
  x.10 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.10 <- scale(x.10)

  x.11 <- x.2/(x.3 + 273)
  x.11 <- scale(x.11)

  x.12 <- exp(-(x.3^2/2))/(sqrt(2*pi))
  x.12 <- scale(x.12)
  
  x.13 <- exp(-(x.2^2/2))/(sqrt(2*pi))
  x.13 <- scale(x.13)

  y.0 <- x.1 + x.5 + x.12 + x.13 + x.6
  y.0 <- scale(y.0)
  y.0[y.0>x.11] <- x.11[y.0>x.11]
  y.0[x.4 == 1] <- min(y.0)
  
  ## rasterization
  y.0 = raster(matrix(y.0,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.2 = raster(matrix(x.2,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.3 = raster(matrix(x.3,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.6 = raster(matrix(x.6,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.7 = raster(matrix(x.7,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.8 = raster(matrix(x.8,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.9 = raster(matrix(x.9,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  x.10 = raster(matrix(x.10,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
  
  ############## STEP 2 - Extract the Sample Points #################
  ## generate evenly distributed coordinates
  
  ## CASE 1: continuous response variable, sample points evenly distributed
  set.seed(i)
	vals <- sample.int(50 ^ 2, 200)
	coord_x = vals %/% 50 
	coord_y = vals %% 51 
	coords_even = data.frame(long = coord_x, lat = coord_y) 
	myfile1.1 <- file.path(getwd(), paste0("coords/coords_even/coord", "_", i, ".csv"))
	write.csv(coords_even, file = myfile1.1)

  y <- extract(y.0, coords_even) 
	y <- data.frame(y) 
	x2 <- extract(x.2, coords_even) 
	x2 <- data.frame(x2)
	x3 <- extract(x.3, coords_even)
	x3 <- data.frame(x3)
	x6 <- extract(x.6, coords_even)
	x6 <- data.frame(x6)
	x7 <- extract(x.7, coords_even)
	x7 <- data.frame(x7)
  x8 <- extract(x.8, coords_even)
	x8 <- data.frame(x8)
	x9 <- extract(x.9, coords_even)
	x9 <- data.frame(x9)
	x10 <- extract(x.10, coords_even)
	x10 <- data.frame(x10)

  simdata <- do.call("cbind", list(y,x2,x3,x6,x7,x8,x9,x10))
	myfile2.1 <- file.path(getwd(), paste0("samples/samples_ctn_even/sample", "_", i, ".csv"))
	write.csv(simdata, file = myfile2.1)

  ############## STEP 3 - BlockCV Splitting #################
  species_even.1 <- SpatialPointsDataFrame(coords = coords_even, data = y)

	bcv_l <- spatialBlock(speciesData = species_even.1, rasterLayer = y.0, theRange = 2783125, k = 4, selection = "systematic", iteration = 100, xOffset = 0, yOffset = 0, showBlocks = FALSE) #2*2=4 blocks, 50*111,325/2=2,783,125(round) #2783124->9blks
  myfile3.1 <- file.path(getwd(), paste0("blockcv/blockcv_l_even/bcv_l", "_", i, ".csv"))
  write.csv(bcv_l[["foldID"]], file = myfile3.1)

  bcv_m <- spatialBlock(speciesData = species_even.1, rasterLayer = y.0, theRange = 1113251, k = 4, selection = "systematic", iteration = 100, xOffset = 0, yOffset = 0, showBlocks = FALSE) #5*5=25 blocks,50*111,325/5=1,113,250(round)
  myfile4.1 <- file.path(getwd(), paste0("blockcv/blockcv_m_even/bcv_m", "_", i, ".csv"))
  write.csv(bcv_m[["foldID"]], file = myfile4.1)

  bcv_s <- spatialBlock(speciesData = species_even.1, rasterLayer = y.0, theRange = 556624, k = 4, selection = "systematic", iteration = 100, xOffset = 0, yOffset = 0, showBlocks = FALSE) #10*10=100 blocks,50*111,325/10=556,625(round)
  myfile5.1 <- file.path(getwd(), paste0("blockcv/blockcv_s_even/bcv_s", "_", i, ".csv"))
  write.csv(bcv_s[["foldID"]], file = myfile5.1)

  ## CASE 2: discrete response variable, sample points evenly distributed
  y[y>=1.5*mean(y$y)] <- 1
  y[y<1.5*mean(y$y)] <- 0
  ## add noise
  set.seed(i)
	rows <- sample.int(200,80) # not with 0
  y[rows,1] <- 0

  simdata <- do.call("cbind", list(y,x2,x3,x6,x7,x8,x9,x10))
	myfile2.2 <- file.path(getwd(), paste0("samples/samples_dsc_even/sample", "_", i, ".csv"))
	write.csv(simdata, file = myfile2.2)
  species_even.2 <- SpatialPointsDataFrame(coords = coords_even, data = y)

  ## CASE 3: continuous response variable, sample points in clumps
  set.seed(i)
  vals1 <- sample.int(20 ^ 2, 150)
	coord_x1 = vals1 %/% 20 
	coord_y1 = vals1 %% 21 
	coords1 = data.frame(long = coord_x1, lat = coord_y1) 
  set.seed(i)
	vals2 <- sample.int(15 ^ 2, 50)
	coord_x2 = (vals2 %/% 15) +35 
	coord_y2 = (vals2 %% 16) +35
	coords2 = data.frame(long = coord_x2, lat = coord_y2) 
	coords_clmp <- rbind(coords1, coords2)
	myfile1.2 <- file.path(getwd(), paste0("coords/coords_clmp/coord", "_", i, ".csv"))
	write.csv(coords_clmp, file = myfile1.2)

  y <- extract(y.0, coords_clmp) 
	y <- data.frame(y) 
	x2 <- extract(x.2, coords_clmp) 
	x2 <- data.frame(x2)
	x3 <- extract(x.3, coords_clmp)
	x3 <- data.frame(x3)
	x6 <- extract(x.6, coords_clmp)
	x6 <- data.frame(x6)
	x7 <- extract(x.7, coords_clmp)
	x7 <- data.frame(x7)
  x8 <- extract(x.8, coords_clmp)
	x8 <- data.frame(x8)
	x9 <- extract(x.9, coords_clmp)
	x9 <- data.frame(x9)
	x10 <- extract(x.10, coords_clmp)
	x10 <- data.frame(x10)

  simdata <- do.call("cbind", list(y,x2,x3,x6,x7,x8,x9,x10))
	myfile2.3 <- file.path(getwd(), paste0("samples/samples_ctn_clmp/sample", "_", i, ".csv"))
	write.csv(simdata, file = myfile2.3)

  species_clmp.1 <- SpatialPointsDataFrame(coords = coords_clmp, data = y)

	bcv_l <- spatialBlock(speciesData = species_clmp.1, rasterLayer = y.0, theRange = 2783124, k = 4, selection = "systematic", iteration = 100, xOffset = 0, yOffset = 0, showBlocks = FALSE) #2*2=4 blocks, 50*111,325/2=2,783,125(round)
  myfile3.2 <- file.path(getwd(), paste0("blockcv/blockcv_l_clmp/bcv_l", "_", i, ".csv"))
  write.csv(bcv_l[["foldID"]], file = myfile3.2)

  bcv_m <- spatialBlock(speciesData = species_clmp.1, rasterLayer = y.0, theRange = 1113251, k = 4, selection = "systematic", iteration = 100, xOffset = 0, yOffset = 0, showBlocks = FALSE) #5*5=25 blocks,50*111,325/5=1,113,250(round)
  myfile4.2 <- file.path(getwd(), paste0("blockcv/blockcv_m_clmp/bcv_m", "_", i, ".csv"))
  write.csv(bcv_m[["foldID"]], file = myfile4.2)

  bcv_s <- spatialBlock(speciesData = species_clmp.1, rasterLayer = y.0, theRange = 556624, k = 4, selection = "systematic", iteration = 100, xOffset = 0, yOffset = 0, showBlocks = FALSE) #10*10=100 blocks,50*111,325/10=556,625(round)
  myfile5.2 <- file.path(getwd(), paste0("blockcv/blockcv_s_clmp/bcv_s", "_", i, ".csv"))
  write.csv(bcv_s[["foldID"]], file = myfile5.2)

  ## CASE 4: discrete response variable, sample points in clumps
  y[y>=1.5*mean(y$y)] <- 1
  y[y<1.5*mean(y$y)] <- 0
  set.seed(i)
	rows <- sample.int(200,80) # not with 0
  y[rows,1] <- 0

  simdata <- do.call("cbind", list(y,x2,x3,x6,x7,x8,x9,x10))
	myfile2.4 <- file.path(getwd(), paste0("samples/samples_dsc_clmp/sample", "_", i, ".csv"))
	write.csv(simdata, file = myfile2.4)
  species_clmp.2 <- SpatialPointsDataFrame(coords = coords_clmp, data = y)
} 


############## STEP 4 - Visualization (only for the last sim) #################
## plot all variables
x.1 <- raster(matrix(x.1,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
x.4 <- raster(matrix(x.4,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
x.5 <- raster(matrix(x.5,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
x.11 <- raster(matrix(x.11,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
x.12 <- raster(matrix(x.12,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
x.13 <- raster(matrix(x.13,nrow=gridsize[1],ncol=gridsize[2],byrow=TRUE),xmn=0, xmx=50, ymn=0, ymx=50)
varstack <- stack(x.1,x.2,x.3,x.4,x.5,x.6,x.7,x.8,x.9,x.10,x.11,x.12,x.13,y.0)
names(varstack) <- c("x.01","x.02","x.03","x.04","x.05","x.06","x.07","x.08","x.09","x.10","x.11","x.12","x.13","y.0")
png("plots/varstack.png", width=800)
show_landscape(varstack, n_col=7, n_row=2)
dev.off()

## plot sample species data points - continuous
point.size <- 0.7
my.palette <- brewer.pal(n = 5, name = "YlGnBu")
edges = species_even.1@bbox
edges[1, 1] = -2
edges[1, 2] = 52
edges[2, 2] = 52
edges[2, 1] = -2

png("plots/species_ctn_even.png", width=640)
spplot(species_even.1, col.regions = my.palette, cex = point.size, xlim = edges[1, ], ylim = edges[2, ],scales =list(draw = TRUE))
dev.off()

png("plots/species_ctn_clmp.png", width=640)
spplot(species_clmp.1, col.regions = my.palette, cex = point.size, xlim = edges[1, ], ylim = edges[2, ],scales =list(draw = TRUE))
dev.off()

## plot sample species data points - discrete
my.palette <- c("blue","red")

species_even.2@data[["y"]] = as.factor(species_even.2@data[["y"]])
png("plots/species_dsc_even.png", width=410)
spplot(species_even.2, col.regions = my.palette, cex = point.size, xlim = edges[1, ], ylim = edges[2, ], scales =list(draw = TRUE))
dev.off()

species_clmp.2@data[["y"]] = as.factor(species_clmp.2@data[["y"]])
png("plots/species_dsc_clmp.png", width=410)
spplot(species_clmp.2, col.regions = my.palette, cex = point.size, xlim = edges[1, ], ylim = edges[2, ], scales =list(draw = TRUE))
dev.off()
