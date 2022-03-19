library(NLMR)
library(raster)
library(sdmpredictors)
library(leaflet)
library(raster)
library(gdistance)
library(vegan)
library(ResistanceGA)
library(spdep)
library(rmaxent)
library(rJava)
library(dismo)
library("coenocliner")
library(resGF)

setwd("~/Desktop/SAMARCH/Paper1/submission")
out <- '~/Desktop/SAMARCH/Paper1/submission/test/multi_surface'
dir.create(out, showWarnings = F)
setwd(out)
out_scenario0 <- paste0(out, "/test0/")
dir.create(out_scenario0, showWarnings = F)
out_time <- paste0(out_scenario0, "/out_time/")
dir.create(out_time, showWarnings = F)
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


iters <- 50
# dataframe to save
final_cor <- NULL;
GA_GF_cor <- NULL;
mydf_top_GF  <- NULL;
mydf_top_GA  <- NULL;
cor_ld <- NULL;
final_dfIBE <- NULL;
GF1GF2_df <- NULL;
imp_gf <- NULL;
mydf_top_GF0 <- NULL;
eval_test_df <- NULL;

# Parameters
nb_ind <- 10000 # inital number of ind
per_cent <- 0.01 # % of ind kept
M <- 200 # inital snps

n1 <- 0.5 # % of SNP in raster1
n2 <- 0.3 # % of SNP in raster2
n3 <- 0.2
final_snp <- 50 # final number of snps

for (i in 1:iters) {
  library(adegenet)
  library(pegas)
  library(poppr)
  out_scenario <- paste0(out_scenario0, "/iter_", i, '/')
  dir.create(out_scenario, showWarnings = F)
  setwd(out_scenario)
  
  print("Building rasters")
  crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  
  e <- c(-5.9, 2.8, 47, 53.5)
  # use the following raster to have be in the same spatial environment
  bathy <- raster("../../../BO_bathymean_lonlat.tif")
  # can be loaded as followed:
  # bathy <- load_layers( layercodes = 'BO_bathymean', equalarea=FALSE, rasterstack=TRUE)
  bathy.crop <- crop(bathy, e, snap="out") 
  
  raster1 <- nlm_gaussianfield(ncol = 100, nrow = 100, autocorr_range = 70, mag_var = 25, nug = 0, resolution = res(bathy.crop)[1])
  crs(raster1) <- crs.wgs
  extent(raster1) <- extent(bathy.crop)
  raster2 <-  nlm_fbm(ncol = 100, nrow = 100, fract_dim = 0.5, resolution= res(bathy.crop)[1])
  crs(raster2) <- crs.wgs
  extent(raster2) <- extent(bathy.crop)
  raster4 <- nlm_randomcluster(ncol = 100, nrow = 100,
                               p = 0.4,
                               ai = c(0.25, 0.25, 0.5),
                               resolution = res(bathy.crop)[1])
  extent(raster4) <- extent(bathy.crop)
  crs(raster4) <- crs.wgs
  raster3 <- nlm_distancegradient(ncol = 100, nrow = 100, origin = c(20, 30, 10, 15), resolution = res(bathy.crop)[1]) # distance_gradient
  crs(raster3) <- crs.wgs
  extent(raster3) <- extent(bathy.crop)
  raster5 <-  nlm_mosaicfield(ncol = 100, nrow = 100,n = NA,infinit = TRUE,
                              collect = FALSE,  resolution = res(bathy.crop)[1])
  extent(raster5) <- extent(bathy.crop)
  crs(raster5) <- crs.wgs
  clipped_stack <- stack(raster1, raster2, raster3, raster4, raster5)
  r1resampled <- projectRaster(clipped_stack, bathy.crop, method = 'bilinear', crs = crs.wgs)
  names(r1resampled) <- cbind("raster1", "raster2", "raster3", "raster4", "raster5")

  r1resampled <- stack(r1resampled)
  mylist <- rep(list("A"),length(r1resampled@layers))
  
  raster1 <- r1resampled$raster1
  names(raster1) <- 'raster1'
  raster1 <- rescale0to1(raster1)
  raster2 <- r1resampled$raster2
  names(raster2) <- 'raster2'
  raster2 <- rescale0to1(raster2)
  raster3 <- r1resampled$raster3
  names(raster3) <- 'raster3'
  raster3 <- rescale0to1(raster3)
  raster4 <- r1resampled$raster4
  names(raster4) <- 'raster4'
  raster4 <- rescale0to1(raster4)
  raster5 <- r1resampled$raster5
  names(raster5) <- 'raster5'
  raster5 <- rescale0to1(raster5)
  
  clipped_stack <- stack(raster1, raster2, raster3, raster4, raster5)
  
  print("Bernoulli distribution (occurrence)")
  
  
  ## Bernoulli distribution (occurrence)
  ## ===================================
  ##############################################
  # raster 1
  x <- seq(from = 0, to = 1, length = nb_ind)
  
  ming <- 0                                # gradient minimum...
  maxg <- 1 
  opt  <- runif(M, min = ming, max = maxg)   # species optima
  tol  <- rep(0.25, M)  
  h <- c(1,3,5,7,9) / 10
  y <- coenocline(x, responseModel = "gaussian",
                  params = cbind(opt = opt, tol = tol, h = h),
                  countModel = "bernoulli")
  # making it into populations
  locs <- as.numeric(x)
  locs0 <- as.data.frame(locs)

  library(pegas)
  yy <- as.data.frame(y)
  t1 <- as.loci(yy,ploidy=1, col.loci = c(1:(dim(yy)[2])))
  x <- loci2genind(t1, ploidy = 1) # no point in the column name
  
  ## get 50% of individual
  smp_size <- floor(per_cent * nrow(locs0))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs0)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  spam.train <- locs0[train_ind, ]
  locs1 <- data.matrix(spam.train)
  
  # snp
  myalleles <- x$tab
  t.sub <- t1[train_ind, ]
  snp.sub <- myalleles[train_ind, ]
  snp1 <- snp.sub[order(row.names(snp.sub)),]
  dim(snp1)
  
  x1 <- loci2genind(t.sub, ploidy = 1) # no point in the column name
  dps1 <- propShared(x1) # using adegenet
  # mydps <- dps
  mydps1 <- 1- dps1
  dps.out1 <- mydps1
  
  
  
  # using NLM
  library(NLMR)
  library(sf)
  library(exactextractr)
  # these are mines
  
  #nlm_distancegradient
  r.pts <- rasterToPoints(raster1, spatial=TRUE)
  names(raster1) <- "raster1"
  # get xy value from raster
  library(rgdal)
  spts <- rasterToPoints(raster1, spatial = TRUE)
  llpts <- spTransform(spts, CRS(crs.wgs))
  x_raster <- as.data.frame(llpts)
  x_raster <- x_raster[order(x_raster$raster1), ]
  locs1 <- as.data.frame(locs1)
  
  # individuals coordinates
  ind_coord_df <- NULL;
  for(m in as.matrix(locs1)){
    temp_df <- NULL
    new_value <- x_raster[which.min(abs(m-x_raster$raster1)),]
    temp_df <- rbind(temp_df, c(m, new_value))
    ind_coord_df <- rbind(ind_coord_df, temp_df)
  }
  ind_coord_df <- as.data.frame(ind_coord_df)
  colnames(ind_coord_df) <- c("average_value", "nearest_value", "x", "y")
  xcoords <- as.numeric(ind_coord_df$x)
  ycoords <- as.numeric(ind_coord_df$y)
  ind_coords1 <- cbind(xcoords, ycoords)
  sample.ind1 <- SpatialPoints(ind_coords1)
  plot(raster1)
  plot(sample.ind1, add = T)
  x1@other$xy <- ind_coords1
  coords1 <- ind_coords1
  pca1 <- dudi.pca(x1, scale = FALSE, scannf = FALSE, nf = 3)
  sample.pts1 <- coords1
  sample.k1 <- SpatialPoints(coords1)
  plot(raster1)
  colorplot(ind_coords1, pca1$li, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  
  print("distribution raster1 - done")
  
  
  
  ##############################################
  # raster 2
  x <- seq(from = 0, to = 1, length = nb_ind)
  
  ming <- 0                                # gradient minimum...
  maxg <- 1 
  opt  <- runif(M, min = ming, max = maxg)   # species optima
  tol  <- rep(0.25, M)  
  
  h <- c(1,3,5,7,9) / 10
  y <- coenocline(x, responseModel = "gaussian",
                  params = cbind(opt = opt, tol = tol, h = h),
                  countModel = "bernoulli")
  # making it into populations
  locs <- as.numeric(x)
  locs0 <- as.data.frame(locs)

  library(pegas)
  yy <- as.data.frame(y)
  t2 <- as.loci(yy,ploidy=1, col.loci = c(1:(dim(yy)[2])))
  x2 <- loci2genind(t2, ploidy = 1) # no point in the column name
  
  ## get 1% of individual
  smp_size <- floor(per_cent * nrow(locs0))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs0)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  spam.train <- locs0[train_ind, ]
  locs2 <- data.matrix(spam.train)
  
  # snp
  myalleles <- x2$tab
  colSums(myalleles)
  t.sub <- t2[train_ind, ]
  snp.sub <- myalleles[train_ind, ]
  snp2 <- snp.sub[order(row.names(snp.sub)),]
  dim(snp2)
  
  x2 <- loci2genind(t.sub, ploidy = 1) # no point in the column name
  
  r.pts <- rasterToPoints(raster2, spatial=TRUE)
  
  # get xy value from raster
  names(raster2) <- "raster2"
  names(raster3) <- "raster3"
  library(rgdal)
  spts <- rasterToPoints(raster2, spatial = TRUE)
  llpts <- spTransform(spts, CRS(crs.wgs))
  x_raster <- as.data.frame(llpts)
  x_raster <- x_raster[order(x_raster$raster2), ]
  locs2 <- as.data.frame(locs2)
  
  # individuals coordinates
  ind_coord_df <- NULL;
  for(m in as.matrix(locs2)){
    temp_df <- NULL
    new_value <- x_raster[which.min(abs(m-x_raster$raster2)),]
    temp_df <- rbind(temp_df, c(m, new_value))
    ind_coord_df <- rbind(ind_coord_df, temp_df)
  }
  ind_coord_df <- as.data.frame(ind_coord_df)
  colnames(ind_coord_df) <- c("average_value", "nearest_value", "x", "y")
  xcoords <- as.numeric(ind_coord_df$x)
  ycoords <- as.numeric(ind_coord_df$y)
  ind_coords2 <- cbind(xcoords, ycoords)
  sample.ind2 <- SpatialPoints(ind_coords2)
  plot(raster2)
  plot(sample.ind2, add = T)
  x2@other$xy <- ind_coords2
  coords2 <- ind_coords2
  pca2 <- dudi.pca(x2, scale = FALSE, scannf = FALSE, nf = 3)
  # gab <- chooseCN(X$other$xy,ask=FALSE,type=2, check.duplicates =F)
  sample.pts2 <- coords2
  sample.k2 <- SpatialPoints(coords2)

  print("# distribution raster2 - done")
  
  
  ##############################################
  # raster 3
  x <- seq(from = 0, to = 1, length = nb_ind)
  
  ming <- 0                                # gradient minimum...
  maxg <- 1 
  opt  <- runif(M, min = ming, max = maxg)   # species optima
  # locs <- seq(ming, maxg, length = 100) 
  # tol <- rep(0.25, 5)
  tol  <- rep(0.25, M)  
  
  h <- c(1,3,5,7,9) / 10
  # h    <- ceiling(rlnorm(M, meanlog = 3)) /100
  y <- coenocline(x, responseModel = "gaussian",
                  params = cbind(opt = opt, tol = tol, h = h),
                  countModel = "bernoulli")
  # making it into populations
  locs <- as.numeric(x)
  locs0 <- as.data.frame(locs)
  # matplot(locs, y, lty = "solid", type = "l", xlab = "pH", ylab = "Abundance")
  
  library(pegas)
  yy <- as.data.frame(y)
  t3 <- as.loci(yy,ploidy=1, col.loci = c(1:(dim(yy)[2])))
  x3 <- loci2genind(t3, ploidy = 1) # no point in the column name
  
  ## get 1% of individual
  smp_size <- floor(per_cent * nrow(locs0))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs0)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  spam.train <- locs0[train_ind, ]
  # spam.train2 <- spam.train[order(spam.train$group),]
  locs3 <- data.matrix(spam.train)
  
  # snp
  myalleles <- x3$tab
  t.sub <- t3[train_ind, ]
  # snp.sub <- myalleles[train_ind, ]
  # snp3 <- snp.sub[order(row.names(snp.sub)),]
  # dim(snp3)
  
  x3 <- loci2genind(t.sub, ploidy = 1) # no point in the column name
  
  r.pts <- rasterToPoints(raster3, spatial=TRUE)
  
  # get xy value from raster
  spts <- rasterToPoints(raster3, spatial = TRUE)
  llpts <- spTransform(spts, CRS(crs.wgs))
  x_raster <- as.data.frame(llpts)
  x_raster <- x_raster[order(x_raster$raster3), ]
  locs3 <- as.data.frame(locs3)
  
  # individuals coordinates
  ind_coord_df <- NULL;
  for(m in as.matrix(locs3)){
    temp_df <- NULL
    new_value <- x_raster[which.min(abs(m-x_raster$raster3)),]
    temp_df <- rbind(temp_df, c(m, new_value))
    ind_coord_df <- rbind(ind_coord_df, temp_df)
  }
  ind_coord_df <- as.data.frame(ind_coord_df)
  colnames(ind_coord_df) <- c("average_value", "nearest_value", "x", "y")
  xcoords <- as.numeric(ind_coord_df$x)
  ycoords <- as.numeric(ind_coord_df$y)
  ind_coords3 <- cbind(xcoords, ycoords)
  sample.ind3 <- SpatialPoints(ind_coords3)
  plot(raster3)
  plot(sample.ind3, add = T)
  x3@other$xy <- ind_coords3
  coords3 <- ind_coords3
  pca3 <- dudi.pca(x3, scale = FALSE, scannf = FALSE, nf = 3)
  # gab <- chooseCN(X$other$xy,ask=FALSE,type=2, check.duplicates =F)
  sample.pts3 <- coords3
  sample.k3 <- SpatialPoints(coords3)

  par(mfrow=c(1,3))
  plot(raster1)
  colorplot(ind_coords1, pca1$li, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  plot(raster2)
  colorplot(ind_coords2, pca2$li, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  plot(raster3)
  colorplot(ind_coords3, pca2$li, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  par(mfrow=c(1,1))
  
  
  ##########################
  # mixing isolated
  # raster1
  # sample indi
  smp_size <- floor(n1 * nrow(locs1))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs1)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  ind_coords1b <- ind_coords1[train_ind,]
  # reduce number of allele
  myalleles1 <- x1$tab
  snp.sub <- myalleles1[train_ind, ]
  snp <- snp.sub[order(row.names(snp.sub)),]
  
  myalleles1 <- snp
  allele_sampled <- sample(seq_len(ncol(myalleles1)), size = final_snp)
  snp.sub <- myalleles1[, allele_sampled]
  snp1b <- snp.sub[order(row.names(snp.sub)),]
  dim(snp1b)
  colnames(snp1b) <- seq(1:dim(snp1b)[2])
  x1b <- df2genind(snp1b, sep = "\t")

  
  # raster2
  # sample indi
  smp_size <- floor(n2 * nrow(locs2))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs2)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  ind_coords2b <- ind_coords2[train_ind,]
  
  # reduce number of allele
  myalleles2 <- x2$tab
  snp.sub <- myalleles2[train_ind, ]
  snp <- snp.sub[order(row.names(snp.sub)),]
  
  myalleles2 <- snp
  allele_sampled <- sample(seq_len(ncol(myalleles2)), size = final_snp)
  snp.sub <- myalleles2[, allele_sampled]
  snp2b <- snp.sub[order(row.names(snp.sub)),]
  dim(snp2b)
  colnames(snp2b) <- seq(1:dim(snp2b)[2])
  x2b <- df2genind(snp2b, sep = "\t")
  
  
  # raster3
  # sample indi
  smp_size <- floor(n3 * nrow(locs3))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs3)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  ind_coords3b <- ind_coords3[train_ind,]
  
  # reduce number of allele
  myalleles3 <- x3$tab
  snp.sub <- myalleles3[train_ind, ]
  snp <- snp.sub[order(row.names(snp.sub)),]
  
  myalleles3 <- snp
  allele_sampled <- sample(seq_len(ncol(myalleles3)), size = final_snp)
  snp.sub <- myalleles3[, allele_sampled]
  snp3b <- snp.sub[order(row.names(snp.sub)),]
  dim(snp3b)
  colnames(snp3b) <- seq(1:dim(snp3b)[2])
  x3b <- df2genind(snp3b, sep = "\t")
  
  
  
  ind_coords_final <- rbind(ind_coords1b, ind_coords2b, ind_coords3b)
  newDat <- repool(x1b, x2b, x3b)
  pca_final <- dudi.pca(newDat, scale = FALSE, scannf = FALSE, nf = 3)
  
  
  plot(raster1)
  colorplot(ind_coords_final, pca_final$li, axes = 1:2, transp = F, add = TRUE, cex = 2)
  
  
  # gradient forest
  library(pegas)
  t <- as.loci(newDat)
  mysamples <- as.data.frame(ind_coords_final)
  colnames(mysamples) <- c('x','y')
  rownames(mysamples)
  lat <- cbind(ID = rownames(mysamples), mysamples) 
  p1 <- lat
  p2 <- cbind(ID = rownames(p1), p1) 
  p2[2] <- NULL
  tail(p2)
  sample.coord.sp <- SpatialPointsDataFrame(p1[,c('x','y')], proj4string=CRS(crs.wgs), data=p2)
  # clipped_stack <- stack(gaus,fbm, combine_surface)
  clim.points <- extract(clipped_stack, sample.coord.sp) 
  clim.points <- cbind(p2, clim.points) 
  # snp <- simp_fre[c(2:(dim(simp_fre)[2]))]
  # thingy to prepare
  sample.coord <- clim.points[,c('x', 'y')]
  id <- To.From.ID(nrow(sample.coord))
  # myfst <- pairwise.fstb(x)
  pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  #this generates the PCNMs, you could stop here if you want all of them
  keep <- round(length(which(pcnm$value > 0))/2)
  pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
  head(pcnm.keep)
  # colnames(clim.points) <- c("ID", "x", "y", "gaus")
  nbvar <- length(colnames(clim.points))
  env.gf <- cbind(clim.points[ , c(seq(4, nbvar, by=1)) ], pcnm.keep)
  library(gradientForest)
  maxLevel <- log2(0.368*nrow(env.gf)/2)
  env.gf <- as.data.frame(env.gf)
  # colnames(env.gf)[1] <- "gaus"
  snp <- as.data.frame(newDat$tab)
  print(paste0("## performing gradient forest"))
  start_time <- Sys.time()
  gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  
  
  pdf(paste0(out_scenario, "Split_Density.pdf"))
  plot(gf,plot.type="Split.Density", imp.vars="raster1")
  dev.off()
  pdf(paste0(out_scenario, "importance.pdf"))
  plot(gf, plot.type = "O")
  dev.off()
  by.importance <- names(extendedForest::importance(gf))
  pdf(paste0(out_scenario, "cummulative_curves.pdf"))
  plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
  dev.off()
  
  
  single_r <- resGF(gf, clipped_stack, save.image = T, results_dir = out_scenario)
  colorplot(ind_coords3, pca2$li, axes = 1:2, transp = F, add = TRUE, cex = 1.75)
  
  dps <- propShared(newDat) # using adegenet
  # mydps <- dps
  mydps <- 1- dps
  dps.out <- mydps
  
  ###
  print("evaluation")
  library(AICcmodavg)
  library(cAIC4)
  sample.k <-SpatialPoints(coords)
  # mydps <- as.matrix(dps.out)
  id <- To.From.ID(nrow(coords))
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance')
  
  print("starting Run_gdistance")
  gd.gaus <- Run_gdistance(gdist.inputs = gdist.inputs, r = gaus_gf)
  gd.fbm_gf <- Run_gdistance(gdist.inputs = gdist.inputs, r = fbm_gf)
  gd.com_surface <- Run_gdistance(gdist.inputs = gdist.inputs,  r = com_surface_gf)
  print("end Run_gdistance")
  
  
  #mod0
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.gaus), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod0 <- mlpe_rga(formula = g.dist  ~  resist.dist + (1 | pop),  data = df)
  #mod1
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.fbm_gf), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod1 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  #mod2
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.com_surface), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod2 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  
  print("model design")
  library(performance)
  perf0 <- model_performance(mod0, metric = "all")
  perf1 <- model_performance(mod1, metric = "all")
  perf2 <- model_performance(mod2, metric = "all")
  AICc_mod0 <- AICcmodavg::AICc(mod0)
  AICc_mod1 <- AICcmodavg::AICc(mod1)
  AICc_mod2 <- AICcmodavg::AICc(mod2)

  
  
  eval_df0 <- NULL;
  final_mod0 <- cbind(perf0, AICc_mod0)
  colnames(final_mod0) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df0 <- rbind(eval_df0, final_mod0)
  final_mod1 <- cbind(perf1, AICc_mod1, r2_mod1$Rsq[1])
  colnames(final_mod1) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df0 <- rbind(eval_df0, final_mod1)
  final_mod2 <- cbind(perf2, AICc_mod2, r2_mod2$Rsq[1])
  colnames(final_mod2) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df0 <- rbind(eval_df0, final_mod2)
  print(eval_df0)
  
  write.table(eval_df0, paste0(out_scenario, "/df_GF0_", i, ".txt"))
  
  rownames(eval_df0) <- c("gaus", "fbm", "combined")
  
  
  AICdf <- eval_df0[order(eval_df0["AIC"]),]
  topAIC <-rownames(AICdf)[1]
  BICdf <- eval_df0[order(eval_df0["BIC"]),]
  topBIC <-rownames(BICdf)[1]
  myR2_conditional <- eval_df0[order(-eval_df0["R2_conditional"]),]
  topR2c <-rownames(myR2_conditional)[1]
  myR2_marginal <- eval_df0[order(-eval_df0["R2_marginal"]),]
  topR2m <-rownames(myR2_marginal)[1]
  myAICc <- eval_df0[order(eval_df0["AICc"]),]
  topAICc <-rownames(myAICc)[1]
  resultsGF0 <- c(topAIC, topBIC, topR2c, topR2m, topAICc)
  mydf_top_GF0 <- rbind(mydf_top_GF0, resultsGF0)
  
  # combined surface
  names(combine_surface) <- 'combine_surface'
  # com_surface_gf0 <- resGF(gf, combine_surface, save.image = F, results_dir = out_scenario)
  clipped_stack <- stack(combine_surface)
  clim.points <- extract(clipped_stack, sample.coord.sp) 
  clim.points <- cbind(p2, clim.points) 
  # snp <- simp_fre[c(2:(dim(simp_fre)[2]))]
  # thingy to prepare
  sample.coord <- clim.points[,c('x', 'y')]
  id <- To.From.ID(nrow(coords))
  # myfst <- pairwise.fstb(x)
  pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  #this generates the PCNMs, you could stop here if you want all of them
  keep <- round(length(which(pcnm$value > 0))/2)
  pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
  head(pcnm.keep)
  colnames(clim.points) <- c("ID", "x", "y", "combine_surface")
  nbvar <- length(colnames(clim.points))
  env.gf <- cbind(clim.points[ , c(seq(4, nbvar, by=1)) ], pcnm.keep)
  library(gradientForest)
  maxLevel <- log2(0.368*nrow(env.gf)/2)
  env.gf <- as.data.frame(env.gf)
  colnames(env.gf)[1] <- "combine_surface"
  snp <- as.data.frame(snp)
  print(paste0("## performing gradient forest - combine_surface"))
  start_time <- Sys.time()
  gf2 <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  com_surface_gf <- resGF(gf2, combine_surface, save.image = F, results_dir = out_scenario)
  comb_surface_imp <- get_imp(gf2, combine_surface)
  
  # dealing with fbm
  names(fbm) <- "fbm"
  # fbm_gf0 <- resGF(gf, fbm, save.image = F, results_dir = out_scenario)
  clipped_stack <- stack(fbm)
  clim.points <- extract(clipped_stack, sample.coord.sp) 
  clim.points <- cbind(p2, clim.points) 
  # snp <- simp_fre[c(2:(dim(simp_fre)[2]))]
  # thingy to prepare
  sample.coord <- clim.points[,c('x', 'y')]
  id <- To.From.ID(nrow(coords))
  # myfst <- pairwise.fstb(x)
  pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  #this generates the PCNMs, you could stop here if you want all of them
  keep <- round(length(which(pcnm$value > 0))/2)
  pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
  head(pcnm.keep)
  colnames(clim.points) <- c("ID", "x", "y", "fbm")
  nbvar <- length(colnames(clim.points))
  env.gf <- cbind(clim.points[ , c(seq(4, nbvar, by=1)) ], pcnm.keep)
  library(gradientForest)
  maxLevel <- log2(0.368*nrow(env.gf)/2)
  env.gf <- as.data.frame(env.gf)
  colnames(env.gf)[1] <- "fbm"
  snp <- as.data.frame(snp)
  print(paste0("## performing gradient forest - fbm"))
  start_time <- Sys.time()
  gf3 <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  fbm_gf <- resGF(gf3, fbm, save.image = F, results_dir = out_scenario)
  fbm_imp <- get_imp(gf3, fbm)
  
  print("fbm_gf - done")
  
  final_gf_res <- stack(gaus_gf, fbm_gf, com_surface_gf)
  plot(final_gf_res)
  cor(values(final_gf_res), values(gaus_gf))
  
  iter_gf <- cbind(gaus_imp, fbm_imp, comb_surface_imp)
  imp_gf <- rbind(imp_gf, iter_gf)
  
  ###
  print("evaluation")
  library(AICcmodavg)
  library(cAIC4)
  sample.k <-SpatialPoints(coords)
  # mydps <- as.matrix(dps.out)
  id <- To.From.ID(nrow(coords))
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance')
  
  print("starting Run_gdistance")
  gd.gaus <- Run_gdistance(gdist.inputs = gdist.inputs, r = gaus_gf)
  gd.fbm_gf <- Run_gdistance(gdist.inputs = gdist.inputs, r = fbm_gf)
  gd.com_surface <- Run_gdistance(gdist.inputs = gdist.inputs,  r = com_surface_gf)
  print("end Run_gdistance")
  
  
  mycor_gf <- cor(lower(mydps), gd.gaus)
  
  # e.dist <- lower(as.matrix(dist(coords)))
  #mod0
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.gaus), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod0 <- mlpe_rga(formula = g.dist  ~  resist.dist + (1 | pop),  data = df)
  #mod1
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.fbm_gf), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod1 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  #mod2
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.com_surface), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod2 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  print("model design")
  
  library(performance)
  perf0 <- model_performance(mod0, metric = "all")
  perf1 <- model_performance(mod1, metric = "all")
  perf2 <- model_performance(mod2, metric = "all")
  AICc_mod0 <- AICcmodavg::AICc(mod0)
  AICc_mod1 <- AICcmodavg::AICc(mod1)
  AICc_mod2 <- AICcmodavg::AICc(mod2)

  
  
  eval_df <- NULL;
  final_mod0 <- cbind(perf0, AICc_mod0)
  colnames(final_mod0) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df <- rbind(eval_df, final_mod0)
  final_mod1 <- cbind(perf1, AICc_mod1)
  colnames(final_mod1) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df <- rbind(eval_df, final_mod1)
  final_mod2 <- cbind(perf2, AICc_mod2)
  colnames(final_mod2) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df <- rbind(eval_df, final_mod2)
  print(eval_df)
  
  write.table(eval_df, paste0(out_scenario, "/df_GF_", i, ".txt"))
  
  
  save.image(paste0(out_scenario0, "GF-vs-GA.R"))
  
  
  
  ########################
  ##################
  ## resistanceGA ##
  ##################
  out_basic_scenario <- paste0(out_scenario, "/basicRasterGA/")
  dir.create(out_basic_scenario, showWarnings = F)
  setwd(out_scenario)
  
  # resistanceGA
  print("## resistanceGA")
  
  start_time <- Sys.time()
  # mystack <- stack(clipped_stack)
  # crs(sample.k) <- crs(gaus)
  
  ## Create Rga objects
  mydps <- as.data.frame(mydps)
  rownames(mydps) <- rownames(ind_coord_df)
  colnames(mydps) <- rownames(ind_coord_df)
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance') # commuteDistance
  # writeRaster(gaus, paste0(mydir, "/asc/gaus.asc"))
  
  folder <- paste0(out_scenario, "/resistanceGA/")
  dir.create(folder, showWarnings = F)
  GA.inputs <- GA.prep(ASCII.dir = gaus,
                       Results.dir = paste0(out_scenario,"/resistanceGA/"),
                       select.trans = list("A"),
                       seed = 123,
                       parallel = 2
  )
  SS_RESULTS <- SS_optim(gdist.inputs = gdist.inputs,
                         GA.inputs = GA.inputs)
  # Make combined, optimized resistance surface.
  # optim.resist <- Combine_Surfaces(PARM = Multi.Surface_optim$GA.summary@solution,
  # gdist.inputs = gdist.inputs,
  # GA.inputs = GA.inputs,
  # rescale = TRUE)
  optim.resist <-raster(paste0(out_scenario,"/resistanceGA/Results/gaus.asc"))
  # crs(optim.resist) <- crs.wgs
  optim.resist@data@min <- min(values(optim.resist))
  optim.resist@data@max <- max(values(optim.resist))
  optim.resist <- rescale0to1(optim.resist)
  plot(optim.resist)
  writeRaster(optim.resist, paste0(out_scenario, "/basicRasterGA.tif"), overwrite = T)
  end_time <- Sys.time()
  time_GA1 <- end_time - start_time
  
  
  # fbm
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance') # commuteDistance
  # writeRaster(gaus, paste0(mydir, "/asc/gaus.asc"))
  names(fbm) <- "fbm"
  folder <- paste0(out_scenario, "/fbm/")
  dir.create(folder, showWarnings = F)
  GA.inputs <- GA.prep(ASCII.dir = fbm,
                       #Results.dir = 'all.comb',
                       Results.dir = paste0(out_scenario,"/fbm/"),
                       select.trans = list("A"),
                       seed = 123,
                       parallel = 2
  )
  SS_RESULTS <- SS_optim(gdist.inputs = gdist.inputs,
                         GA.inputs = GA.inputs)
  optim.fbm <-raster(paste0(out_scenario,"/fbm/Results/fbm.asc"))
  # crs(optim.fbm) <- crs.wgs
  optim.fbm@data@min <- min(values(optim.fbm))
  optim.fbm@data@max <- max(values(optim.fbm))
  optim.fbm <- rescale0to1(optim.fbm)
  plot(optim.fbm)
  # writeRaster(optim.fbm, paste0(out_scenario, "/optim.fbm.tif"), overwrite = T)
  
  
  # combine_surface
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance') # commuteDistance
  # writeRaster(gaus, paste0(mydir, "/asc/gaus.asc"))
  names(combine_surface) <- "combine_surface"
  folder <- paste0(out_scenario, "/combine_surface/")
  dir.create(folder, showWarnings = F)
  GA.inputs <- GA.prep(ASCII.dir = combine_surface,
                       #Results.dir = 'all.comb',
                       Results.dir = paste0(out_scenario,"/combine_surface/"),
                       select.trans = list("A"),
                       seed = 123,
                       parallel = 2
  )
  SS_RESULTS <- SS_optim(gdist.inputs = gdist.inputs,
                         GA.inputs = GA.inputs)
  optim.combine_surface <-raster(paste0(out_scenario,"/combine_surface/Results/combine_surface.asc"))
  # crs(optim.combine_surface) <- crs.wgs
  optim.combine_surface@data@min <- min(values(optim.combine_surface))
  optim.combine_surface@data@max <- max(values(optim.combine_surface))
  optim.combine_surface <- rescale0to1(optim.combine_surface)
  plot(optim.combine_surface)
  # writeRaster(optim.combine_surface, paste0(out_scenario, "/optim.combine_surface.tif"), overwrite = T)
  
  
  ###
  print("evaluation")
  library(AICcmodavg)
  library(cAIC4)
  sample.k <-SpatialPoints(coords)
  mydps <- as.matrix(dps.out)
  id <- To.From.ID(nrow(coords))
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance')
  
  gd.bGA <- Run_gdistance(gdist.inputs = gdist.inputs, r = optim.resist)
  gd.fbm <- Run_gdistance(gdist.inputs = gdist.inputs, r = optim.fbm)
  gd.combine_surface <- Run_gdistance(gdist.inputs = gdist.inputs,  r = optim.combine_surface)
  
  
  
  #mod5
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.bGA), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod5 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  #mod6
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.fbm), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod6 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  #mod7
  df <- data.frame(g.dist = as.numeric(lower(mydps)), resist.dist = as.numeric(gd.combine_surface), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  mod7 <- mlpe_rga(formula = g.dist ~ resist.dist + (1 | pop),  data = df)
  
  
  
  
  
  library(performance)
  perf5 <- model_performance(mod5, metric = "all")
  perf6 <- model_performance(mod6, metric = "all")
  perf7 <- model_performance(mod7, metric = "all")
  AICc_mod5 <- AICcmodavg::AICc(mod5)
  AICc_mod6 <- AICcmodavg::AICc(mod6)
  AICc_mod7 <- AICcmodavg::AICc(mod7)

  
  eval_df2 <- NULL;
  final_mod5 <- cbind(perf5, AICc_mod5)
  colnames(final_mod5) <- c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df2 <- rbind(eval_df2, final_mod5)
  final_mod6 <- cbind(perf6, AICc_mod6)
  colnames(final_mod6) <- c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df2 <- rbind(eval_df2, final_mod6)
  final_mod7 <- cbind(perf7, AICc_mod7)
  colnames(final_mod7) <- c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc")
  eval_df2 <- rbind(eval_df2, final_mod7)
  print(eval_df2)
  write.table(eval_df2, paste0(out_scenario, "/df_GA_", i, ".txt"))
  
  # write.csv(eval_df, paste0(out_scenario, "MLPE_results.csv"))
  rownames(eval_df) <- c("gaus", "fbm", "combined")
  AICdf <- eval_df[order(eval_df["AIC"]),]
  topAIC <-rownames(AICdf)[1]
  BICdf <- eval_df[order(eval_df["BIC"]),]
  topBIC <-rownames(BICdf)[1]
  myR2_conditional <- eval_df[order(-eval_df["R2_conditional"]),]
  topR2c <-rownames(myR2_conditional)[1]
  myR2_marginal <- eval_df[order(-eval_df["R2_marginal"]),]
  topR2m <-rownames(myR2_marginal)[1]
  myAICc <- eval_df[order(eval_df["AICc"]),]
  topAICc <-rownames(myAICc)[1]
  resultsGF <- c(topAIC, topBIC, topR2c, topR2m, topAICc)
  
  rownames(eval_df2) <- c("gaus", "fbm", "combined")
  AICdf <- eval_df2[order(eval_df2["AIC"]),]
  topAIC <-rownames(AICdf)[1]
  BICdf <- eval_df2[order(eval_df2["BIC"]),]
  topBIC <-rownames(BICdf)[1]
  myR2_conditional <- eval_df2[order(-eval_df2["R2_conditional"]),]
  topR2c <-rownames(myR2_conditional)[1]
  myR2_marginal <- eval_df2[order(-eval_df2["R2_marginal"]),]
  topR2m <-rownames(myR2_marginal)[1]
  myAICc <- eval_df2[order(eval_df2["AICc"]),]
  topAICc <-rownames(myAICc)[1]
  resultsGA <- c(topAIC, topBIC, topR2c, topR2m, topAICc)
  
  
  mydf_top_GF <- rbind(mydf_top_GF, resultsGF)
  mydf_top_GA <- rbind(mydf_top_GA, resultsGA)
  
  # correlation genetic to resistance
  library(adegenet)
  res_m <- as.matrix(gd.gaus)
  res_ga <- as.matrix(gd.bGA)
  
  mycor_gf <- cor(lower(mydps), lower(res_m))
  mycor_ga <- cor(lower(mydps), lower(res_ga))
  
  Dgen <- lower(mydps)
  Dgeo <- dist(sample.k@coords)
  ibd <- cor(Dgen, Dgeo)
  ibd
  pdf(paste0(out_scenario, "correlation_dps_landscape_distance.pdf"))
  par(mfrow = c(1,3))
  plot(Dgen ~ Dgeo, cex = 0.2, main = paste0("IBD: ", round(ibd, 3), "; p:" ))
  abline(lm(Dgen ~ Dgeo), col = "red")
  plot(lower(mydps) ~ lower(res_ga), cex = 0.2, main = paste0("GA: ", round(mycor_ga, 3)))
  abline(lm(lower(mydps) ~ lower(res_ga)), col = "red")
  plot(lower(mydps) ~ lower(res_m), cex = 0.2, main = paste0("GF: ", round(mycor_gf, 3)))
  abline(lm(lower(mydps) ~ lower(res_m)), col = "red")
  dev.off()
  par(mfrow = c(1,1))
  
  iter_mycor <- c(ibd, mycor_gf, mycor_ga)
  cor_ld <- rbind(cor_ld, iter_mycor)
  
  # plotting
  geodist <- as.matrix(sample.k@coords)
  tr <- transition(gaus_gf, function(x) 1/mean(x) , directions = 8)
  trC <- geoCorrection(tr, "c", scl = TRUE)
  cosDist_gf <- costDistance(trC, geodist)
  path.1a <- shortestPath( trC, sample.k[3], sample.k[16], output="SpatialLines")
  path.2a <- shortestPath( trC, sample.k[12], sample.k[33], output="SpatialLines")
  path.3a <- shortestPath( trC, sample.k[1], sample.k[45], output="SpatialLines")
  path.4a <- shortestPath( trC, sample.k[6], sample.k[23], output="SpatialLines")
  path.5a <- shortestPath( trC, sample.k[11], sample.k[37], output="SpatialLines")
  path.6a <- shortestPath( trC, sample.k[4], sample.k[49], output="SpatialLines")
  path.7a <- shortestPath( trC, sample.k[28], sample.k[45], output="SpatialLines")
  
  geodist <- as.matrix(sample.k@coords)
  tr <- transition(optim.resist, function(x) 1/mean(x) , directions = 8)
  trC <- geoCorrection(tr, "c", scl = TRUE)
  cosDist_gf <- costDistance(trC, geodist)
  path.1b <- shortestPath( trC, sample.k[3], sample.k[16], output="SpatialLines")
  path.2b <- shortestPath( trC, sample.k[12], sample.k[33], output="SpatialLines")
  path.3b <- shortestPath( trC, sample.k[1], sample.k[45], output="SpatialLines")
  path.4b <- shortestPath( trC, sample.k[6], sample.k[23], output="SpatialLines")
  path.5b <- shortestPath( trC, sample.k[11], sample.k[37], output="SpatialLines")
  path.6b <- shortestPath( trC, sample.k[4], sample.k[49], output="SpatialLines")
  path.7b <- shortestPath( trC, sample.k[28], sample.k[45], output="SpatialLines")
  
  pdf(paste0(out_scenario, "surfaces_LCP.pdf"))
  par(mfrow = c(1,3))
  plot(gaus, main=(paste0("IBD:", round(ibd, 2))))
  colorplot( x1@other$xy, pca1$l1, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  
  plot(optim.resist, main=(paste0("GA:", round(mycor_ga, 2))))
  plot(sample.k, add = T)
  lines( path.1b, col="red")
  lines( path.2b, col="blue")
  lines( path.3b, col="black")
  lines( path.4b, col="red")
  lines( path.5b, col="blue")
  lines( path.6b, col="black")
  lines( path.7b, col="red")
  plot(gaus_gf, main=(paste0("GF:", round(mycor_gf, 2))))
  plot(sample.k, add = T)
  lines( path.1a, col="red")
  lines( path.2a, col="blue")
  lines( path.3a, col="black")
  lines( path.4a, col="red")
  lines( path.5a, col="blue")
  lines( path.6a, col="black")
  lines( path.7a, col="red")
  dev.off()
  par(mfrow = c(1,1))
  
  eval_test_df <- rbind(eval_test_df, final_mod0)
  eval_test_df <- rbind(eval_test_df, final_mod5)
  
  #########################################
  # investigating IBE
  # building MMR
  print("investigating IBE")
  library('ecodist')
  geoMat <- as.matrix(dist(sample.k@coords))
  gaus_IBE <- raster::extract(gaus, sample.k)
  rownames(mydps) <- rownames(geoMat)
  colnames(mydps) <- colnames(geoMat)
  # Explicit conversion to vector
  mydps_v <- as.vector(dist(mydps))
  geoMat_v <- as.vector(dist(geoMat))
  ecoMat_gaus <- as.vector(dist(gaus_IBE))
  
  MMR_eucl <- MRM(mydps_v ~ geoMat_v, nperm=999) # NULL model
  MMR_IBE <- MRM((mydps_v) ~ (geoMat_v) + (ecoMat_gaus), nperm=999)
  
  df <- data.frame(g.dist = as.numeric(lower(mydps)), e.dist = as.numeric(lower(geoMat)), resist.dist = as.numeric(ecoMat_gaus), pop = id$pop1)
  df$g.dist[df$g.dist == "Inf"] <- 0
  MLPE_eucl <- mlpe_rga(formula = g.dist ~ e.dist  + (1 | pop),  data = df)
  MLPE_eucl
  MLPE_IBE <- mlpe_rga(formula = g.dist ~ e.dist +  resist.dist + (1 | pop),  data = df)
  MLPE_IBE
  
  gf_contri <- as.data.frame(gf$overall.imp)
  gf_contri$names <- rownames(gf_contri)
  # mycontrib_list <- list('raster1', 'raster2', 'raster3', 'raster4', 'raster5', 'raster6', 'raster7', 'raster8', 'raster9', 'raster10')
  # gf_contri2 <- as.data.frame(gf_contri[rownames(gf_contri) %in% mycontrib_list, ])
  colnames(gf_contri) <- c('resGF')
  gf_contri$resGF_pc <- gf_contri$resGF / sum(gf_contri$resGF) *100
  colnames(gf_contri) <- c('resGF', 'variable', 'pc')
  gf_contri
  mygaus_gf <- gf_contri[gf_contri$variable == "raster1",]
  
  library(AICcmodavg)
  library(performance)
  perf_eucl <- model_performance(MLPE_eucl, metric = c("AIC", "BIC", "R2","RMSE"))
  perf_IBE <- model_performance(MLPE_IBE, metric = c("AIC", "BIC", "R2","RMSE"))
  AICc_mod0 <- AICc(MLPE_eucl)
  AICc_mod1 <- AICc(MLPE_IBE)
  
  final_eucl <- cbind(perf_eucl, AICc_mod0)
  final_IBE <- cbind(perf_IBE, AICc_mod1)
  
  
  MMR_df <- NULL;
  res_euc <- cbind(mygaus_gf[3], MMR_eucl$r.squared[1], MMR_eucl$r.squared[2], MMR_eucl$coef[1], MMR_eucl$coef[3], MMR_eucl$coef[2], MMR_eucl$coef[4],MMR_eucl$coef[2], MMR_eucl$coef[4], MMR_eucl$r.squared[1], MMR_eucl$r.squared[2], final_eucl)
  colnames(res_euc) <- c("GF_pc","mmr_R2", "pval","beta", "pval", "beta_geo", "pval", "beta_res", "pval", "R2", "pval", "MPLE_AIC", "BIC", "R2c", "R2m","RMSE", "AICc")
  MMR_df <- rbind(MMR_df, res_euc)
  res_IBE <- cbind(mygaus_gf[3], MMR_IBE$r.squared[1], MMR_IBE$r.squared[2], MMR_IBE$coef[1], MMR_IBE$coef[4], MMR_IBE$coef[2], MMR_IBE$coef[5],MMR_IBE$coef[3], MMR_IBE$coef[6], MMR_IBE$r.squared[1], MMR_IBE$r.squared[2], final_IBE)
  colnames(res_IBE) <- c("GF_pc","mmr_R2", "pval","beta", "pval", "beta_geo", "pval", "beta_res", "pval", "R2", "pval", "MPLE_AIC", "BIC", "R2c", "R2m","RMSE", "AICc")
  MMR_df <- rbind(MMR_df, res_IBE)
  rownames(MMR_df) <- c(paste0('eucl_', i), paste0('IBE_', i))
  write.csv(MMR_df, paste0(out_scenario, "IBE_results_", i , ".csv"))
  
  final_dfIBE <- rbind(final_dfIBE, MMR_df)
  # print(final_dfIBE)
  
  
  # time
  mytimetable <- matrix(c(time_GF1, time_GA1),ncol=2,byrow=TRUE)
  colnames(mytimetable) <- c("GF","GA")
  write.table(mytimetable, paste0(out_time, "time_table", i, ".txt"))
  # correlation
  mycor_fbm <- cor(values(fbm), values(gaus))
  mycor_combine_surface <-cor(values(combine_surface), values(gaus))
  mycor <- matrix(c(mycor_fbm, mycor_combine_surface),ncol=2,byrow=TRUE)
  colnames(mycor) <- c("fbm","combine_surface")
  final_cor <- rbind(final_cor, mycor)
  
  GFGA <- cor(values(gaus_gf), values(optim.resist))
  write.csv(GFGA, paste0(out_scenario, "GFGA.txt"))
  GA_GF_cor <- rbind(GA_GF_cor, GFGA)
  
  save.image(paste0(out_scenario0, "GF-vs-GA.R"))
  
  print(paste0("end: iter_", i))
  
  
}
colnames(mydf_top_GA) <- c("AIC", "BIC", "R2c", "R2m", "AICc")
colnames(mydf_top_GF) <- c("AIC", "BIC", "R2c", "R2m", "AICc")
colnames(mydf_top_GF0) <- c("AIC", "BIC", "R2c", "R2m", "AICc")
colnames(cor_ld) <- c("IBD","GF", "GA")
write.table(final_dfIBE, paste0(out_scenario0, "final_dfIBE", iters, ".txt"))
write.table(final_cor, paste0(out_scenario0, "final_cor", iters, ".txt"))
write.table(GA_GF_cor, paste0(out_scenario0, "GA_GF_cor_",iters,".txt"))
write.table(mydf_top_GA, paste0(out_scenario0, "mydf_top_GA_",iters,".txt"))
write.table(mydf_top_GF, paste0(out_scenario0, "mydf_top_GF_",iters,".txt"))
write.table(mydf_top_GF0, paste0(out_scenario0, "mydf_top_GF0_",iters,".txt"))
write.table(cor_ld, paste0(out_scenario0, "cor_ld",iters,".txt"))
write.table(GF1GF2_df, paste0(out_scenario0, "GF1GF2_df_",iters,".txt"))
write.table(imp_gf, paste0(out_scenario0, "imp_gf_",iters,".txt"))
write.table(eval_test_df, paste0(out_scenario0, "eval_test_df_",iters,".txt"))

