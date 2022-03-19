# indvidual_gradient
library(NLMR)
library(raster)
library(sdmpredictors)
library(leaflet)
library(raster)
library(vegan)
library(ResistanceGA)
library(spdep)
library(rmaxent)
library(rJava)
library(dismo)
library("coenocliner")
library(resGF)
library(climateStability)
library(topAIC) ## find package / topicmodels
library(gdistance)


setwd("~/Desktop/SAMARCH/Paper1/submission")
out <- '~/Desktop/SAMARCH/Paper1/submission/test/sim_single_cont'

setwd(out)
out_scenario0 <- paste0(out, "/individual_gradient_test/")
dir.create(out_scenario0, showWarnings = F)
out_time <- paste0(out_scenario0, "/out_time/")
dir.create(out_time, showWarnings = F)
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

iters <- 1
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
for (i in 1:iters) {
  library(adegenet)
  library(pegas)
  library(poppr)
  out_scenario <- paste0(out_scenario0, "/iter_", i, '/')
  dir.create(out_scenario, showWarnings = F)
  raster <- paste0(out_scenario, "/raster_", i, '/')
  dir.create(raster, showWarnings = F)
  
  ## Bernoulli distribution (occurrence)
  ## ===================================
  nb_ind <- 10000
  per_cent <- 0.01
  x <- seq(from = 0, to = 1, length = nb_ind)
  
  M <- 200 # snps
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
  t <- as.loci(yy,ploidy=1, col.loci = c(1:(dim(yy)[2])))
  x <- loci2genind(t, ploidy = 1) # no point in the column name
  
  ## get 50% of individual
  smp_size <- floor(per_cent * nrow(locs0))
  print(paste0("samples size: ", smp_size))
  train_ind <- sample(seq_len(nrow(locs0)), size = smp_size)
  train_ind <- train_ind[order(train_ind)]
  spam.train <- locs0[train_ind, ]
  # spam.train2 <- spam.train[order(spam.train$group),]
  locs1 <- data.matrix(spam.train)
  
  # hist(locs1$V1, breaks=20)
  # hist(values(gaus), breaks=50)
  # hist(values(gaus_gf), breaks=50)
  # snp
  myalleles <- x$tab
  colSums(myalleles)
  t.sub <- t[train_ind, ]
  snp.sub <- myalleles[train_ind, ]
  snp <- snp.sub[order(row.names(snp.sub)),]
  dim(snp)
  
  x <- loci2genind(t.sub, ploidy = 1) # no point in the column name
  dps <- propShared(x) # using adegenet
  mydps <- 1- dps
  dps.out <- mydps
  tre <- nj(dist(as.matrix(x)))
  plot(tre, typ="unrooted", cex=0.7)
  write.csv(dps.out, paste0(raster, "/dps.csv"))
  
  # using NLM
  library(NLMR)
  library(sf)
  library(exactextractr)
  # these are mines
  print("Building rasters")
  
  crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  e <- c(-5.9, 2.8, 47, 53.5)
  # bathy <- load_layers( layercodes = 'BO_bathymean', equalarea=FALSE, rasterstack=TRUE)
  bathy <- raster("../../../BO_bathymean_lonlat.tif")
  bathy.crop <- crop(bathy, e, snap="out")
  library(RandomFields)
  library(climateStability)
  # raster1 <- nlm_planargradient(ncol = 100, nrow = 100, resolution = res(bathy.crop)[1])
  raster1 <- nlm_gaussianfield(ncol = 100, nrow = 100, autocorr_range = 70, mag_var = 25, nug = 0, resolution = res(bathy.crop)[1])
  crs(raster1) <- crs.wgs
  extent(raster1) <- extent(bathy.crop)
  raster2 <-  nlm_fbm(ncol = 100, nrow = 100, fract_dim = 0.5, resolution= res(bathy.crop)[1])
  crs(raster2) <- crs.wgs
  extent(raster2) <- extent(bathy.crop)
  
  writeRaster(raster1, paste0(raster, "/raster1.tif"), overwrite=TRUE)
  
  clipped_stack <- stack(raster1, raster2)
  r1resampled <- projectRaster(clipped_stack, bathy.crop, method = 'bilinear', crs = crs.wgs)
  names(r1resampled) <- cbind("raster1", "raster2")
  r1resampled <- stack(r1resampled)
  mylist <- rep(list("A"),length(r1resampled@layers))
  
  gaus <- r1resampled$raster1
  names(gaus) <- 'gaus'
  gaus <- rescale0to1(gaus)
  fbm <- r1resampled$raster2
  names(fbm) <- 'fbm'
  fbm <- rescale0to1(fbm)
  clipped_stack <- stack(gaus, fbm)
  r.pts <- rasterToPoints(gaus, spatial=TRUE)
  
  # get xy value from raster
  library(rgdal)
  spts <- rasterToPoints(gaus, spatial = TRUE)
  llpts <- spTransform(spts, CRS(crs.wgs))
  x_raster <- as.data.frame(llpts)
  x_raster <- x_raster[order(x_raster$gaus), ]
  locs1 <- as.data.frame(locs1)
  
  # individuals coordinates
  ind_coord_df <- NULL;
  for(m in as.matrix(locs1)){
    temp_df <- NULL
    new_value <- x_raster[which.min(abs(m-x_raster$gaus)),]
    temp_df <- rbind(temp_df, c(m, new_value))
    ind_coord_df <- rbind(ind_coord_df, temp_df)
  }
  ind_coord_df <- as.data.frame(ind_coord_df)
  colnames(ind_coord_df) <- c("average_value", "nearest_value", "x", "y")
  xcoords <- as.numeric(ind_coord_df$x)
  ycoords <- as.numeric(ind_coord_df$y)
  ind_coords <- cbind(xcoords, ycoords)
  sample.ind <- SpatialPoints(ind_coords)
  
  
  
  plot(gaus)
  plot(sample.ind, add = T)
  x@other$xy <- ind_coords
  coords <- ind_coords
  pca1 <- dudi.pca(x, scale = FALSE, scannf = FALSE, nf = 3)
  plot(lower(mydps) ~ c(dist(ind_coords)), cex = 0.1)
  abline(lm(lower(mydps) ~ c(dist(ind_coords))), col = "red")
  sample.pts <- coords
  sample.k <- SpatialPoints(coords)
  eucl <- log(dist(coords))
  gen <- lower(mydps)/ (1-lower(mydps))
  plot(eucl, gen, ylab="Fst/(1/Fst)", xlab = "log distance")
  
  plot(gaus)
  colorplot(ind_coords, pca1$li, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  ind_coords_m <- as.data.frame(ind_coords)
  write.csv(ind_coords_m, paste0(raster, "/lat_lon.csv"))
  
  # combine surface
  full.dat <- paste0(out_scenario, '/full_sim_dat/')
  dir.create(full.dat, showWarnings = F)
  mystack <- stack(clipped_stack)
  mylist <- rep(list("A"),length(mystack@layers))
  PARM <- c(3,3,100,3,3,100, 3,3,100) #max parameter for continuous surface 2
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = gen,
                             method = 'costDistance')
  GA.inputs <- GA.prep(ASCII.dir = mystack,
                       #Results.dir = 'all.comb',
                       Results.dir = paste0(full.dat),
                       select.trans = mylist,
                       parallel = 2
  )
  combine_surface <- Combine_Surfaces(PARM = PARM,
                                      gdist.inputs = gdist.inputs,
                                      GA.inputs = GA.inputs,
                                      out = NULL,
                                      rescale = TRUE)
  combine_surface <- rescale0to1(combine_surface)
  names(combine_surface) <- 'combine_surface'
  mystack <- stack(mystack, combine_surface)
  plot(mystack)
  
  # gradient forest
  library(pegas)
  t <- as.loci(x)
  mysamples <- as.data.frame(sample.pts)
  colnames(mysamples) <- c('x','y')
  rownames(mysamples)
  lat <- cbind(ID = rownames(mysamples), mysamples)
  p1 <- lat
  p2 <- cbind(ID = rownames(p1), p1)
  p2[2] <- NULL
  tail(p2)
  sample.coord.sp <- SpatialPointsDataFrame(p1[,c('x','y')], proj4string=CRS(crs.wgs), data=p2)
  clipped_stack <- stack(gaus)
  clim.points <- raster::extract(clipped_stack, sample.coord.sp)
  clim.points <- cbind(p2, clim.points)
  sample.coord <- clim.points[,c('x', 'y')]
  id <- To.From.ID(nrow(coords))
  pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  # generates the PCNMs, you could stop here if you want all of them
  keep <- round(length(which(pcnm$value > 0))/2)
  pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
  head(pcnm.keep)
  colnames(clim.points) <- c("ID", "x", "y", "gaus")
  nbvar <- length(colnames(clim.points))
  env.gf <- cbind(clim.points[ , c(seq(4, nbvar, by=1)) ], pcnm.keep)
  library(gradientForest)
  maxLevel <- log2(0.368*nrow(env.gf)/2)
  env.gf <- as.data.frame(env.gf)
  colnames(env.gf)[1] <- "gaus"
  snp <- as.data.frame(snp)
  print(paste0("## performing gradient forest - gaus"))
  start_time <- Sys.time()
  gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  
  
  pdf(paste0(out_scenario, "Split_Density.pdf"))
  plot(gf,plot.type="Split.Density", imp.vars="gaus")
  dev.off()
  pdf(paste0(out_scenario, "importance.pdf"))
  plot(gf, plot.type = "O")
  dev.off()
  by.importance <- names(extendedForest::importance(gf))
  pdf(paste0(out_scenario, "cummulative_curves.pdf"))
  plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
  dev.off()
  
  plot(gf,plot.type="Cumulative.Importance", imp.vars="gaus")
  
  write.csv(gf$species.pos.rsq, paste0(out_scenario ,"species.pos.rsq.csv"))
  
  gaus_gf <- resGF(gf, gaus, save.image = F, results_dir = out_scenario)
  names(fbm) <- 'gaus'
  fbm_gf <- resGF(gf, fbm, save.image = F, results_dir = out_scenario)
  names(combine_surface) <- 'gaus'
  com_surface_gf <- resGF(gf, combine_surface, save.image = F, results_dir = out_scenario)
  
  
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
  
  gd.gaus <- Run_gdistance(gdist.inputs = gdist.inputs, r = gaus_gf)
  gd.fbm_gf <- Run_gdistance(gdist.inputs = gdist.inputs, r = fbm_gf)
  gd.com_surface <- Run_gdistance(gdist.inputs = gdist.inputs,  r = com_surface_gf)

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
  final_mod1 <- cbind(perf1, AICc_mod1)
  colnames(final_mod1) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc" )
  eval_df0 <- rbind(eval_df0, final_mod1)
  final_mod2 <- cbind(perf2, AICc_mod2)
  colnames(final_mod2) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc" )
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
  
  ##########################################################"
  ## second models
  gaus_gf <- resGF(gf, gaus, save.image = T, results_dir = out_scenario)
  gaus_imp <- get_imp(gf, gaus)
  
  writeRaster(gaus_gf, paste0(raster, "/gaus_gf.tif"), overwrite  = TRUE)
  end_time <- Sys.time()
  time_GF1 <- end_time - start_time
  save.image(paste0(out_scenario0, "GF-vs-GA.R"))
  # plot(gaus_gf_inv)
  # colorplot(x$other$xy, pca1$li, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  writeRaster(gaus_gf, paste0(out_scenario, "/gaus_gf.tif"), overwrite = T)
  
  # combined surface
  names(combine_surface) <- 'combine_surface'
  clipped_stack <- stack(combine_surface)
  clim.points <- raster::extract(clipped_stack, sample.coord.sp)
  clim.points <- cbind(p2, clim.points)
  sample.coord <- clim.points[,c('x', 'y')]
  id <- To.From.ID(nrow(coords))
  pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  #this generates the PCNMs, you could stop here if you want all of them
  keep <- round(length(which(pcnm$value > 0))/2)
  pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
  head(pcnm.keep)
  colnames(clim.points) <- c("ID", "x", "y", "combine_surface")
  nbvar <- length(colnames(clim.points))
  env.gf <- cbind(clim.points[ , c(seq(4, nbvar, by=1)) ], pcnm.keep)
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
  clipped_stack <- stack(fbm)
  clim.points <- raster::extract(clipped_stack, sample.coord.sp)
  clim.points <- cbind(p2, clim.points)
  # thingy to prepare
  sample.coord <- clim.points[,c('x', 'y')]
  id <- To.From.ID(nrow(coords))
  pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  #this generates the PCNMs, you could stop here if you want all of them
  keep <- round(length(which(pcnm$value > 0))/2)
  pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
  head(pcnm.keep)
  colnames(clim.points) <- c("ID", "x", "y", "fbm")
  nbvar <- length(colnames(clim.points))
  env.gf <- cbind(clim.points[ , c(seq(4, nbvar, by=1)) ], pcnm.keep)
  maxLevel <- log2(0.368*nrow(env.gf)/2)
  env.gf <- as.data.frame(env.gf)
  colnames(env.gf)[1] <- "fbm"
  snp <- as.data.frame(snp)
  print(paste0("## performing gradient forest - fbm"))
  start_time <- Sys.time()
  gf3 <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  fbm_gf <- resGF(gf3, fbm, save.image = F, results_dir = out_scenario)
  fbm_imp <- get_imp(gf3, fbm)
  
  gaus_table <- gf_importance(gf)
  write.csv(gaus_table, paste0(out_scenario, "gaus_table.csv"))
  comb_surface_table <- gf_importance(gf2)
  write.csv(comb_surface_table, paste0(out_scenario, "comb_surface_table.csv"))
  fbm_table <- gf_importance(gf3)
  write.csv(fbm_table, paste0(out_scenario, "fbm_table.csv"))
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
  id <- To.From.ID(nrow(coords))
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance')
  
  gd.gaus <- Run_gdistance(gdist.inputs = gdist.inputs, r = gaus_gf)
  gd.fbm_gf <- Run_gdistance(gdist.inputs = gdist.inputs, r = fbm_gf)
  gd.com_surface <- Run_gdistance(gdist.inputs = gdist.inputs,  r = com_surface_gf)
  mycor_gf <- cor(lower(mydps), gd.gaus)
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
  colnames(final_mod1) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc" )
  eval_df <- rbind(eval_df, final_mod1)
  final_mod2 <- cbind(perf2, AICc_mod2)
  colnames(final_mod2) <-  c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc" )
  eval_df <- rbind(eval_df, final_mod2)
  print(eval_df)
  write.table(eval_df, paste0(out_scenario, "/df_GF_", i, ".txt"))
  
  save.image(paste0(out_scenario, "GF-vs-GA.R"))
  
  
  
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
  ## Create Rga objects
  mydps <- as.data.frame(mydps)
  rownames(mydps) <- rownames(ind_coord_df)
  colnames(mydps) <- rownames(ind_coord_df)
  gdist.inputs <- gdist.prep(n.Pops = length(sample.k),
                             samples = sample.k,
                             response = lower(mydps),
                             method = 'costDistance') # commuteDistance

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
  optim.resist <-raster(paste0(out_scenario,"/resistanceGA/Results/gaus.asc"))
  optim.resist@data@min <- min(values(optim.resist))
  optim.resist@data@max <- max(values(optim.resist))
  optim.resist <- rescale0to1(optim.resist)
  plot(optim.resist)
  writeRaster(optim.resist, paste0(raster, "/basicRasterGA.tif"), overwrite = T)
  end_time <- Sys.time()
  time_GA1 <- end_time - start_time
  
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
  colnames(final_mod5) <- c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc"  )
  eval_df2 <- rbind(eval_df2, final_mod5)
  final_mod6 <- cbind(perf6, AICc_mod6)
  colnames(final_mod6) <- c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc" )
  eval_df2 <- rbind(eval_df2, final_mod6)
  final_mod7 <- cbind(perf7, AICc_mod7)
  colnames(final_mod7) <- c("AIC", "BIC", "R2_conditional", "R2_marginal", "ICC", "RMSE", "Sigma", "AICc" )
  eval_df2 <- rbind(eval_df2, final_mod7)
  print(eval_df2)
  write.table(eval_df2, paste0(out_scenario, "/df_GA_", i, ".txt"))
  
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
  colorplot(x$other$xy, pca1$l1, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  
  plot(optim.resist, main=(paste0("GA:", round(mycor_ga, 2))))
  plot(sample.k, add = T)
  plot(gaus_gf, main=(paste0("GF:", round(mycor_gf, 2))))
  plot(sample.k, add = T)
  dev.off()
  par(mfrow = c(1,1))
  
  eval_test_df <- rbind(eval_test_df, final_mod0)
  eval_test_df <- rbind(eval_test_df, final_mod5)
  
  pdf(paste0(out_scenario, "mysurfaces.pdf"))
  par(mfrow = c(1,3))
  plot(gaus, col = topo.colors(20), main=(paste0("IBD:", round(ibd, 2))))
  colorplot(x$other$xy, pca1$l1, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  plot(optim.resist, col = rev(terrain.colors(50)), main=(paste0("GA:", round(mycor_ga, 2))))
  colorplot(x$other$xy, pca1$l1, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  plot(gaus_gf, col = rev(terrain.colors(50)), main=(paste0("GF:", round(mycor_gf, 2))))
  colorplot(x$other$xy, pca1$l1, axes = 1:2, transp = F, add = TRUE, cex = 1.5)
  dev.off()
  par(mfrow = c(1,1))
  
  ## updated version
  ## nb - Iterative random forest optimization (IRF) not run in this script
  ## please refer to Pless (2021)
  # IRF_ped <- raster("TrainingTestingRfsrc_all/prediction_res_ITER3.tif")
  
  library(adegenet)
  res_m <- as.matrix(gd.gaus)
  gd.GA <- Run_gdistance(gdist.inputs = gdist.inputs, r = optim.resist)
  res_ga <- as.matrix(gd.GA)
  # gd.IRF <- Run_gdistance(gdist.inputs = gdist.inputs, r = IRFO_ped)
  # res_irf <- as.matrix(gd.IRF)
  mycor_gf <- cor(lower(mydps), lower(res_m))
  mycor_ga <- cor(lower(mydps), lower(res_ga))
  # mycor_irf <- cor(lower(mydps), lower(res_irf))
  
  Dgen <- lower(mydps)
  Dgeo <- dist(sample.k@coords)
  ibd <- cor(Dgen, Dgeo)
  ibd
  
  
  
  # raster correlation
  stack1 <- stack(gaus_gf, optim.resist)
  jnk=layerStats(stack1, 'pearson', na.rm=T)
  cor_GF_GA=jnk$'pearson correlation coefficient'
  cor_GF_GA
  
  # stack2 <- stack(gaus_gf, IRFO_ped)
  # jnk=layerStats(stack2, 'pearson', na.rm=T)
  # cor_GF_IRF=jnk$'pearson correlation coefficient'
  # cor_GF_IRF
  
  # stack3 <- stack(optim.resist, IRFO_ped)
  # jnk=layerStats(stack3, 'pearson', na.rm=T)
  # cor_GA_IRF=jnk$'pearson correlation coefficient'
  # cor_GA_IRF
  
  
  
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
  colnames(gf_contri) <- c('resGF')
  gf_contri$resGF_pc <- gf_contri$resGF / sum(gf_contri$resGF) *100
  colnames(gf_contri) <- c('resGF', 'variable', 'pc')
  gf_contri
  mygaus_gf <- gf_contri[gf_contri$variable == "gaus",]
  
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
  save.image(paste0(out_scenario, "GF-vs-GA.R"))
  print(paste0("end: iter_", i))
  
}
colnames(mydf_top_GA) <- c("AIC", "BIC", "R2c", "R2m", "AICc")
colnames(mydf_top_GF) <- c("AIC", "BIC", "R2c", "R2m", "AICc",)
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




