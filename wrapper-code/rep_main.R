### Set working directory to this folder (wrapper-code) before running, e.g.
###   setwd(".../QC-Screen/wrapper-code")

### Set any flag to FALSE to skip that block (useful when runs are long).
run_Part1_Model1a <- TRUE
run_Part1_Model1b <- TRUE
run_Part1_Model1c <- TRUE
run_Part2_Table4 <- TRUE
run_Part3_Model3 <- TRUE
run_Part4_FigureB1 <- TRUE

## source() in each script's folder so outputs (csv, rds, figures) land under that Part.
run_in_dir <- function(subdir, script) {
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(normalizePath(file.path(owd, subdir), winslash = "/", mustWork = TRUE))
  source(script, encoding = "UTF-8")
}

######################################
#### Part 1: Section 5.1 / Appx B.1 ####
######################################
if (isTRUE(run_Part1_Model1a)) {
  run_in_dir("Part 1/Model 1a", "model_1a_n_400_p_1000.R")
  run_in_dir("Part 1/Model 1a", "model_1a_n_400_p_5000.R")
}

if (isTRUE(run_Part1_Model1b)) {
  run_in_dir("Part 1/Model 1b", "model_1b_n_400_p_1000.R")
  run_in_dir("Part 1/Model 1b", "model_1b_n_400_p_5000_alpha_05.R")
  run_in_dir("Part 1/Model 1b", "model_1b_n_400_p_5000_alpha_075.R")
}

if (isTRUE(run_Part1_Model1c)) {
  run_in_dir("Part 1/Model 1c", "model_1c_n_400_p_1000_alpha_05.R")
  run_in_dir("Part 1/Model 1c", "model_1c_n_400_p_1000_alpha_075.R")
  run_in_dir("Part 1/Model 1c", "model_1c_n_400_p_5000_alpha_05.R")
  run_in_dir("Part 1/Model 1c", "model_1c_n_400_p_5000_alpha_075.R")
}

#################
#### Part 2 #####
#################
if (isTRUE(run_Part2_Table4)) {
  run_in_dir("Part 2", "Table4.R")
}

#################
#### Part 3 #####
#################
if (isTRUE(run_Part3_Model3)) {
  run_in_dir("Part 3", "model_3a_alpha_05.R")
  run_in_dir("Part 3", "model_3a_alpha_075.R")
  run_in_dir("Part 3", "model_3b_alpha_05.R")
  run_in_dir("Part 3", "model_3b_alpha_075.R")
  run_in_dir("Part 3", "model_3c_alpha_05.R")
  run_in_dir("Part 3", "model_3c_alpha_075.R")
}

#################
#### Part 4 #####
#################
if (isTRUE(run_Part4_FigureB1)) {
  run_in_dir("Part 4", "FigureB1.R")
}

