########### GLOBALS ---------------------------
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


########### FUNCTIONS ---------------------------
## Corrplot
Corrplot <-function(covs){
  set.seed(12)

  str <- deparse(substitute(covs)) 
  samp <- spatSample(covs, 5000)
  samp <- na.omit(samp)
  cor <- cor(samp)
  cor[is.na(cor)] <- 0
  
  setwd(outDir)
  png(height=800, width=800, file=paste0(str,".png"), type='cairo')
  corrplot(cor, method='number')
  dev.off()
  setwd(dataDir)
  cor
}


## Spatial splits for AOA
SpatialFolds <- function(df,str){
  set.seed(0)
  pt <- df
  coordinates(pt) <- ~POINT_X+POINT_Y
  crs(pt) <- CRS("+init=epsg:4326")
  pt <- SpatialPoints(pt)
  n <- 5000 #5 km apart
  print(n)  
  folds <- spatialBlock(speciesData = pt, theRange = n, k = 3,'systematic',
                        iteration = 2,progress=T)
  df$folds <- folds$foldID
  df <- df[,-10]
  df
}


## Function to map under-sampled areas and pseudo-sample "no change" and "loss" data
PseudoSamp <- function(df=usa_ref,covs=covs){
  str <- deparse(substitute(covs)) 
  str1 <- deparse(substitute(df)) 
  
  print(str)
  df <- df[, !(colnames(df) %in% c('fmask','defo'))]
  df$diffPlt <- df$AGB_T_HA_2018 - df$AGB_T_HA_2010
  exc <- c('diff_lvod','flii_envi','mgmt_mgmt', 'TC_envi','height_envi',
           'forType_envi', 'diff_jpl','bio_climate','forType_envi')
          
  # use all covs for feature space analysis
  p <- vect(df, geom=c("POINT_X", "POINT_Y"))
  print(paste(nrow(df), 'total plots'))
  df.covs <- terra::extract(covs, p)[,-1]
  df <- data.frame(df, df.covs)
  print(paste(nrow(df), 'total plots after masking'))
  df <- df[, !(colnames(df) %in% exc)]
  df <- na.omit(df)
  
  #clustering of ref data 
  df <- SpatialFolds(df, str)
  folds <- CreateSpacetimeFolds(df, spacevar='folds',k=length(unique(df$folds)))
  nms <- setdiff(names(covs),exc)
  
  unregister_dopar()
  
  model <- train(df[,which(names(df) %in% nms)], 
                 df$diffPlt, method="rf", importance=TRUE,tuneLength = 5,
                 trControl = trainControl(method="cv",index=folds$index))
  
  print(model)
  DI <- trainDI(model=model)
  cl <- makeCluster(6)
  registerDoParallel(cl) 
  AOA <- aoa(covs, model = model, trainDI = DI,cl=cl) #saving AOA not working!
  stopCluster(cl)  
    
  msk <-  AOA$AOA
  msk [msk==1] <- NA
  under_samp <- terra::mask(covs[[1]],msk)
  under_samp[!is.na(under_samp)] <- 1
  plot(under_samp, main='Under-sampled areas', col='red')
  plot(p,add=T)
  setwd(outDir)
  writeRaster(under_samp, 'under_samp.tif',overwrite=T)
  setwd(dataDir)
  under_samp
}


## FUNCTION TO INCLUDE PSEUDO LOSS SAMPLES WITHIN UNDER-SAMPLED AREAS
AddPseudo <- function(df.raw, mask='unmasked'){
  str <- deparse(substitute(df.raw))
  
  if (grepl('phl',str)==T){
    str <- 'PHL'
    type <- 'NFI'
    setwd('C:/AOA_PHL_1km')
    if (mask =='unmasked'){
      f <- list.files(getwd(), 'unmasked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }else{
      f <- list.files(getwd(), '_masked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }}
  
  if (grepl('bra',str)==T){
    str <- 'BRA'
    type <- 'LiDAR'
    setwd('C:/AOA_BRA_1km')
    if (mask =='unmasked'){
      f <- list.files(getwd(), 'unmasked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }else{
      f <- list.files(getwd(), '_masked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }}
  
  if (grepl('nld',str)==T){
    str <- 'NLD'
    type <- 'NFI'
    setwd('C:/AOA_NLD_1km')
    if (mask =='unmasked'){
      f <- list.files(getwd(), 'unmasked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }else{
      f <- list.files(getwd(), 'unmasked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }}
  
  if (grepl('swe',str)==T){
    str <- 'SWE'
    type <- 'NFI'
    setwd('C:/AOA_SWE_1km')
    if (mask =='unmasked'){
      f <- list.files(getwd(), 'unmasked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }else{
      f <- list.files(getwd(), '_masked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }}
  
  if (grepl('usa',str)==T){
    str <- 'USA'
    type <- 'LiDAR'
    setwd('C:/AOA_USA_1km')
    if (mask =='unmasked'){
      f <- list.files(getwd(), 'unmasked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }else{
      f <- list.files(getwd(), '_masked.+pseudoLoss')
      flist <- lapply(f, function(x) readRDS(x))
      df.pseudo <- ldply(flist,data.frame)
    }}
  
  df.pseudo1 <- df.pseudo %>%
    arrange(POINT_Y, -diffPlt) %>%
    filter(duplicated(POINT_Y) == FALSE)
  df.pseudo1 <- df.pseudo1 %>%
    arrange(POINT_X, -diffPlt) %>%
    filter(duplicated(POINT_X) == FALSE)
  
  n <- round(nrow(df.raw) * 0.1, 0) #10% of total should be added
  if (nrow(df.pseudo1) < n){ n <- round(nrow(df.raw) * 0.08, 0)}
  
  add_loss <- df.raw[sample(nrow(df.raw),n), ]
  df.pseudo2 <- df.pseudo1[sample(nrow(df.pseudo1),n), ]
  
  add_loss$POINT_X  <- df.pseudo2$POINT_X
  add_loss$POINT_Y  <- df.pseudo2$POINT_Y
  add_loss$AGB_T_HA_2010  <- df.pseudo2$AGB_T_HA_2010
  add_loss$AGB_T_HA_2018 <- df.pseudo2$AGB_T_HA_2018
#  add_loss$diff <- df.pseudo2$diffPlt
  setwd(dataDir)
  write.csv(rbind(df.raw,add_loss),paste0('ref',type,'_deltaAGB_',str,'_pseudo.csv'),row.names=F)
  rbind(df.raw,add_loss)
}


## Function to extract covariates value at plot locations
ExtractVal <- function(plt=val.ldr,covs=covs_phl_unmasked){
  str <- deparse(substitute(covs)) 
  rsl <- .00088889
  
  plt$diffPlt <- plt$AGB_T_HA_2018 - plt$AGB_T_HA_2010
  plt$diffSD <- sqrt(plt$AGB_T_HA_SD_2010^2) + sqrt(plt$AGB_T_HA_SD_2018^2)
  
  plt = plt[,c('POINT_X', 'POINT_Y','AGB_T_HA_2018', 'AGB_T_HA_2010','diffSD', 'diffPlt' )]
  plt <- subset(plt, !is.na(plt$AGB_T_HA_2018))
  plt <- na.omit(plt)
  p <- vect(plt, geom=c("POINT_X", "POINT_Y"))
  
  # extract values 
  plt <- data.frame(plt, terra::extract(covs,p))
  plt <- subset(plt, !is.na(plt$diff_cci))
  
  #assures that data is confined with UN-SEEA classes!
  # remove non-forestb based on CCI land cover dataset for now
  plt$lcov10_LC <- ifelse (plt$lcov10_LC < 40 | plt$lcov10_LC >130, NA, plt$lcov10_LC)
 subset(plt, !is.na(plt$lcov10_LC))
}


## Machine learning model fitting and cross validation 
XVal <- function(vt){ #runs all model
  str <- deparse(substitute(vt)) 
  print(paste('RUNNING:',str))
  
  #exclude these
  exc <- c('X','ID','diff_lvod','folds','texMean_flux','texVar_flux','texMean_cci','texVar_cci')
  data <- vt[!names(vt) %in% exc]
  data <- na.omit(data)
  head(data)
  
  #XG --------------------------------------------------
  unregister_dopar()
  n <- nrow(data)
  
  tune_grid <- expand.grid(nrounds = c(500,1000,2000),
                           max_depth = 5,
                           eta = c(0.005,0.01,0.03,0.1,0.5), ####!!! very sensitive #low eta value means model more robust to overfitting but slower to comput
                           gamma = c(0.01,0.1,0.5),
                           colsample_bytree = 0.75,
                           min_child_weight = c(1,5,10),
                           subsample = 0.5)
  
  
#  gbm_mod <- caret::train(diffPlt~ .,  data=data[1:n,-c(1:5)],method='gbm',
 #                         trControl= trainControl( method = "cv",number = 5,
  #                                                 savePredictions=T,search='random'))
                                                   
  xg_mod <- train(diffPlt ~., data = data[1:n,-c(1:5)], method = "xgbTree",
                  trControl=trainControl(method = "cv", number = 5,savePredictions = T),
                  tuneGrid = tune_grid)

  #RF --------------------------------------
  rf_mod <- caret::train(diffPlt~ .,  data=data[1:n,-c(1:5)],
                        method='ranger', importance='permutation',
                        trControl= trainControl( method = "cv",number = 5,savePredictions=T))

  ####SVM -----------------------------------------
  tune_grid <- expand.grid(C=seq(0.1,3,0.1))
  svm_mod <- caret::train(diffPlt~ .,  data=data[1:n,-c(1:5)],method='svmLinear',
                          trControl= trainControl(method = "cv",number = 5,savePredictions=T,
                                                  search='random',tuneGrid = tune_grid))
  
  ### Make sure training data are independent between base learners and meta learner
  pred2 <- svm_mod$pred
  tuned <- svm_mod$bestTune
  pred2 <- subset(pred2, pred2$C==tuned[1,1])
  pred2 <- pred2[order(pred2$rowIndex),]
  
  pred3 <- xg_mod$pred
  pred3 <- pred3[order(pred3$rowIndex),]
  
  pred4 <- rf_mod$pred
  tuned <- rf_mod$bestTune
  pred4 <- subset(pred4, pred4$mtry== tuned[1,1] & pred4$splitrule ==tuned[1,2]&
                    pred4$min.node.size ==tuned[1,3])
  pred4 <- pred4[order(pred4$rowIndex),]
    
  #### Ensemble model ---------------------------------
  #ifelse shuffling folds1-5
  #input index = cv_folds to train
  data1 <- data.frame(pred2$obs, pred2$Resample, #pred1$pred, 
                      pred2$pred, pred3$pred, pred4$pred)  
  names(data1) <- c('diffPlt','Resample','svmRadial', 'xgbTree','ranger')
  data1$Resample <- parse_number(data1$Resample)
  data2 <- transform(data1, folds= sample(Resample))
    
  cv_folds <- createFolds(data2$folds,k=5, returnTrain = T)
  en_mod <- caret::train(diffPlt~ .,  data=data2[,-c(2,6)],method='ranger',importance='permutation',
                        trControl=trainControl(method = "cv",5,index=cv_folds, savePredictions=T))
  
  unregister_dopar()
  if ("diffSD" %in% colnames(data)){ 
    cw <- 1/ (data$diffSD^2)
  }else{cw<-NA}
  
  rf_se <- ranger(diffPlt ~ ., data=data2[,-c(2,6)], quantreg = T, case.weights=cw,
                  keep.inbag = T,importance='permutation')
  setwd(dataDir)
  fin <- list(svm_mod,xg_mod,rf_mod,en_mod,rf_se,data)
  saveRDS(fin, paste0('Mods_',str,'.rds'))  
  return(fin)
}

  
=