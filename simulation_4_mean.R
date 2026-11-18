###
### This program runs the assessment of mean calibration in the small sample simulation
###

library(foreach)
library(doParallel)

source("/z_functions.R")
source("/z_functions_true_transition_probs.R")
source("/z_load_packages.R")

### Define time at which risks will be generated and  calibration will be assessed
t.eval <- round(365.25*7)

meanCal <- function(scen, n.cohort, n.pctls, n.sim){
  
  safe_calc <- function(expr) {
    tryCatch(expr, error = function(e) return(NULL))
  }
  ### Create list to store calibration of each estimates
  ### (include a list to store true calibration in each simulation iteration)
  calib.blr.mean.list <- vector("list", 3)
  calib.mlr.mean.list <- vector("list", 3)
  calib.aj.mean.list <- vector("list", 3)
  calib.true.mean.list <- vector("list", 3)
  
  ### Calculate the calibration intecepts and slopes
  for (i in 1:3){
    
    ### BLR-IPCW
    calib.blr.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.blr.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
    ### MLR-IPCW
    calib.mlr.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.mlr.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
    ### Pseudo-value
    calib.aj.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.aj.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
    ### True
    calib.true.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.true.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
  }
  
  ###
  ### Run the simulation loop
  ###
  
  #fail counts for every type of transition probability and every method
  fail_count <- 0
  fail_count_aj_m1 <- 0
  fail_count_aj_m2 <- 0
  fail_count_aj_m3 <- 0
  fail_count_blr_m1 <- 0
  fail_count_blr_m2 <- 0
  fail_count_blr_m3 <- 0
  fail_count_mlr_m1 <- 0
  fail_count_mlr_m2 <- 0
  fail_count_mlr_m3 <- 0
  
  for (sim in 1:n.sim){
    
    ### Choose patids at random
    patids.vec <- sample((1:1000000), n.cohort, replace = FALSE)
    
    ### Extract individuals
    data.mstate.reduc <- data.mstate[data.mstate$id %in% patids.vec, ]
    data.raw.reduc <- data.raw[data.raw$id %in% patids.vec, ]
    
    success <- TRUE
    ###
    ### Calculate calibration for each type of predicted risk
    ###
    for (i in 1:3){
      
      ### Extract predicted risks
      tp.pred[[i]]$id <- 1:nrow(tp.pred[[i]])
      tp.pred.reduc <- tp.pred[[i]][tp.pred[[i]]$id %in% patids.vec, ] |>
        dplyr::select(!id)
      
      if (scen == "C1"){
        
        ###
        ### Estimate mean calibration with BLR-IPCW
        calib.mean <- safe_calc(
          calib_msm(data_ms = data.mstate.reduc,
                    data_raw = data.raw.reduc,
                    j = 1, s = 0, t = t.eval,
                    tp_pred = tp.pred.reduc,
                    calib_type = "blr",
                    assess_moderate = FALSE,
                    assess_mean = TRUE)[["mean"]]
        )
        if (is.null(calib.mean)) {
          success <- FALSE
          if (i == 1) fail_count_blr_m1 <- fail_count_blr_m1 + 1
          if (i == 2) fail_count_blr_m2 <- fail_count_blr_m2 + 1
          if (i == 3) fail_count_blr_m3 <- fail_count_blr_m3 + 1
        } else {
          calib.blr.mean.list[[i]][sim, ] <- calib.mean
        }
        
        ###
        ### Estimate mean calibration with MLR-IPCW
        calib.mean <- safe_calc(
          calib_msm(data_ms = data.mstate.reduc,
                    data_raw = data.raw.reduc,
                    j = 1, s = 0, t = t.eval,
                    tp_pred = tp.pred.reduc,
                    calib_type = "mlr",
                    assess_moderate = FALSE,
                    assess_mean = TRUE)[["mean"]]
        )
        
        if (is.null(calib.mean)) {
          success <- FALSE
          if (i == 1) fail_count_mlr_m1 <- fail_count_mlr_m1 + 1
          if (i == 2) fail_count_mlr_m2 <- fail_count_mlr_m2 + 1
          if (i == 3) fail_count_mlr_m3 <- fail_count_mlr_m3 + 1
        } else {
          calib.mlr.mean.list[[i]][sim, ] <- calib.mean
        }
        
        ###
        ### Estimate mean calibration with Aalen-Johansen
        calib.mean <- safe_calc(
          calib_msm(data_ms = data.mstate.reduc,
                    data_raw = data.raw.reduc,
                    j = 1, s = 0, t = t.eval,
                    tp_pred = tp.pred.reduc,
                    calib_type = "aj",
                    assess_moderate = FALSE,
                    assess_mean = TRUE)[["mean"]]
        )
        if (is.null(calib.mean)) {
          success <- FALSE
          if (i == 1) fail_count_aj_m1 <- fail_count_aj_m1 + 1
          if (i == 2) fail_count_aj_m2 <- fail_count_aj_m2 + 1
          if (i == 3) fail_count_aj_m3 <- fail_count_aj_m3 + 1
        } else {
          calib.aj.mean.list[[i]][sim, ] <- calib.mean
        }
        
      } else if (scen %in% c("C2", "C3","C4","C5", "C6","C7")){
        
        ###
        ### Estimate mean calibration with BLR-IPCW
        calib.mean <- safe_calc(
          calib_msm(data_ms = data.mstate.reduc,
                    data_raw = data.raw.reduc,
                    j = 1, s = 0, t = t.eval,
                    tp_pred = tp.pred.reduc,
                    calib_type = "blr",
                    w_covs = c("x12","x13","x15","x24","x25","x34","x35","x45"),
                    assess_moderate = FALSE,
                    assess_mean = TRUE)[["mean"]]
        )
        if (is.null(calib.mean)) {
          success <- FALSE
          if (i == 1) fail_count_blr_m1 <- fail_count_blr_m1 + 1
          if (i == 2) fail_count_blr_m2 <- fail_count_blr_m2 + 1
          if (i == 3) fail_count_blr_m3 <- fail_count_blr_m3 + 1
        } else {
          calib.blr.mean.list[[i]][sim, ] <- calib.mean
        }
        
        ###
        ### Estimate mean calibration with MLR-IPCW
        calib.mean <- safe_calc(
          calib_msm(data_ms = data.mstate.reduc,
                    data_raw = data.raw.reduc,
                    j = 1, s = 0, t = t.eval,
                    tp_pred = tp.pred.reduc,
                    calib_type = "mlr",
                    w_covs = c("x12","x13","x15","x24","x25","x34","x35","x45"),
                    assess_moderate = FALSE,
                    assess_mean = TRUE)[["mean"]]
        )
        if (is.null(calib.mean)) {
          success <- FALSE
          if (i == 1) fail_count_mlr_m1 <- fail_count_mlr_m1 + 1
          if (i == 2) fail_count_mlr_m2 <- fail_count_mlr_m2 + 1
          if (i == 3) fail_count_mlr_m3 <- fail_count_mlr_m3 + 1
        } else {
          calib.mlr.mean.list[[i]][sim, ] <- calib.mean
        }
        
        ###
        ### Estimate mean calibration with Aalen-Johansen
        calib.mean <- safe_calc(
          calib_msm(data_ms = data.mstate.reduc,
                    data_raw = data.raw.reduc,
                    j = 1, s = 0, t = t.eval,
                    tp_pred = tp.pred.reduc,
                    calib_type = "aj",
                    pv_n_pctls = n.pctls,
                    assess_moderate = FALSE,
                    assess_mean = TRUE)[["mean"]]
        )
        if (is.null(calib.mean)) {
          success <- FALSE
          if (i == 1) fail_count_aj_m1 <- fail_count_aj_m1 + 1
          if (i == 2) fail_count_aj_m2 <- fail_count_aj_m2 + 1
          if (i == 3) fail_count_aj_m3 <- fail_count_aj_m3 + 1
        } else {
          calib.aj.mean.list[[i]][sim, ] <- calib.mean
        }
        
      }
      
      ###
      ### Save true calibration
      calib.true.mean.list[[i]][sim, ] <- colMeans(data.raw.reduc[, paste("p.true", 1:5, sep = "")] - tp.pred.reduc)
      
    }
    if (!success) {
      fail_count <- fail_count + 1
    }
    
  }
  save(
    list = c("calib.blr.mean.list", "calib.mlr.mean.list", "calib.aj.mean.list",
             "calib.true.mean.list", "scen", "n.cohort", "n.sim", "n.pctls", "fail_count","fail_count_blr_m1","fail_count_blr_m2","fail_count_blr_m3",
             "fail_count_mlr_m1","fail_count_mlr_m2","fail_count_mlr_m3",
             "fail_count_aj_m1","fail_count_aj_m2","fail_count_aj_m3"),
    file = paste0("/small_sample_mean_analysis_",
                  scen, "_n", n.cohort, "_npctls", n.pctls,".RData")
  )
  
}

# settings
ScenVec <- c("C1","C2","C3","C4","C5","C6","C7")
n.cohort <- 400  #500/700
n.pctls  <- 10
n.sim    <- 1000

#cl = parallel::makeCluster(6)
#doParallel::registerDoParallel(cl)
cl = makeCluster(7)
clusterSetRNGStream(cl, iseed = 1)
registerDoParallel(cl)

clusterEvalQ(cl, {
  
  # Source all functions ONCE per worker
  source("/z_functions.R")
  source("/z_functions_true_transition_probs.R")
  source("/z_load_packages.R")
  
  NULL
})
### Export specific objects used inside meanCal
clusterExport(cl, c("meanCal", "t.eval", "n.cohort", "n.pctls", "n.sim"))


# 2: switch do %dopar% and export necessary packages
results = foreach(i = ScenVec) %dopar% {
  load(paste0("/small_sample_prep_data_", i, ".RData"))
  meanCal(i, n.cohort, n.pctls, n.sim)
}

stopCluster(cl)
