###
### This program prepares three datasets of 1000000 subjects each
### with seven types of scenario: CS1, CS2, CS3, CS4, CS5, CS6, CS7
###

source("/z_functions.R")
source("/z_functions_true_transition_probs.R")
source("/z_load_packages.R")

### Define time at which risks will be generated and  calibration will be assessed
t.eval <- round(365.25*7)

### Defne the max follow up to be after t.eval
max.follow <- t.eval + 2

### Censoring is applied at a scale resulting in 40% censored by 7 years follow up
cens_scale <- calc.scale(0.6)

### Read in population data
### Write function to read in data and assign patids
extract_function <- function(seed){
  data <- readRDS(paste("/dataout", seed, ".rds", sep = ""))
  data$patid <- ((seed-1)*1000+1):(seed*1000)
  return(data)
}


#-------------C1-------------------------------------

scen <- "C1"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0
cens_beta_x13 <- 0
cens_beta_x15 <- 0
cens_beta_x24 <- 0
cens_beta_x25 <- 0
cens_beta_x34 <- 0
cens_beta_x35 <- 0
cens_beta_x45 <- 0


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC1/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC1/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC1/tmat.dgm1.RData", sep = ""))

print("FINISHED")


#----------------------C2----------------

scen <- "C2"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0.125
cens_beta_x13 <- 0.125
cens_beta_x15 <- 0.125
cens_beta_x24 <- 0.125
cens_beta_x25 <- 0.125
cens_beta_x34 <- 0.125
cens_beta_x35 <- 0.125
cens_beta_x45 <- 0.125


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC2/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC2/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC2/tmat.dgm1.RData", sep = ""))

print("FINISHED")


#---------------------C3-----------------

scen <- "C3"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0.25
cens_beta_x13 <- 0.25
cens_beta_x15 <- 0.25
cens_beta_x24 <- 0.25
cens_beta_x25 <- 0.25
cens_beta_x34 <- 0.25
cens_beta_x35 <- 0.25
cens_beta_x45 <- 0.25


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC3/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC3/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC3/tmat.dgm1.RData", sep = ""))

print("FINISHED")


#---------------------C4-----------------

scen <- "C4"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0.5
cens_beta_x13 <- 0.5
cens_beta_x15 <- 0.5
cens_beta_x24 <- 0.5
cens_beta_x25 <- 0.5
cens_beta_x34 <- 0.5
cens_beta_x35 <- 0.5
cens_beta_x45 <- 0.5


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC4/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC4/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC4/tmat.dgm1.RData", sep = ""))

print("FINISHED")


#---------------------C5-----------------

scen <- "C5"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0.75
cens_beta_x13 <- 0.75
cens_beta_x15 <- 0.75
cens_beta_x24 <- 0.75
cens_beta_x25 <- 0.75
cens_beta_x34 <- 0.75
cens_beta_x35 <- 0.75
cens_beta_x45 <- 0.75


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC5/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC5/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC5/tmat.dgm1.RData", sep = ""))

print("FINISHED")

#---------------------C6-----------------

scen <- "C6"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0.25
cens_beta_x13 <- 0.25
cens_beta_x15 <- 0.25
cens_beta_x24 <- 0.5
cens_beta_x25 <- 0.5
cens_beta_x34 <- 0.25
cens_beta_x35 <- 0.25
cens_beta_x45 <- 0.25


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC6/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC6/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC6/tmat.dgm1.RData", sep = ""))

print("FINISHED")

#---------------------C7-----------------

scen <- "C7"
print(paste("cens = ", scen, sep = " "))

### Set the coefficients for the censoring mechanism
cens_beta_x12 <- 0
cens_beta_x13 <- 0
cens_beta_x15 <- 0
cens_beta_x24 <- 0.75
cens_beta_x25 <- 0.75
cens_beta_x34 <- 0
cens_beta_x35 <- 0
cens_beta_x45 <- 0


### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
                                                    cens_shape = 1,
                                                    cens_scale = cens_scale,
                                                    cens_beta_x12 = cens_beta_x12,
                                                    cens_beta_x13 = cens_beta_x13,
                                                    cens_beta_x15 = cens_beta_x15,
                                                    cens_beta_x24 = cens_beta_x24,
                                                    cens_beta_x25 = cens_beta_x25,
                                                    cens_beta_x34 = cens_beta_x34,
                                                    cens_beta_x35 = cens_beta_x35,
                                                    cens_beta_x45 = cens_beta_x45)

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("/censC7/pop.mstate.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["data.raw"]], paste("/censC7/pop.raw.dgm1.", scen, ".RData", sep = ""))

saveRDS(data.mstate.obj[["tmat"]], paste("/censC7/tmat.dgm1.RData", sep = ""))

print("FINISHED")
