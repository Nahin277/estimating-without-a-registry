## =============================================================================
## run_all.R  --  Reproduce every number and figure in the post.
##
##   cd R
##   Rscript run_all.R
##
## Base R only. No packages. About 4 minutes.
## Output: ../figures/*.png  and  ../results/*.csv
##
## To force a refit after changing anything in 00_inputs.R, delete
## ../results/fit.rds first.
## =============================================================================

t0 <- Sys.time()
source("04_figures.R")      # which pulls in 00 -> 01 -> 02 -> 03
cat(sprintf("\nDone in %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
