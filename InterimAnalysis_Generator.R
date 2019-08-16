rm(list = ls())

dummy <- T #Please set this to FALSE in the real analysis

# Please open this file in RStudio with 32-bit version of R
# Click [Source] on the right hand side of the toolbar above.






############################################ CODE ZONE ###################################################
if (dummy)
  cat('Running in dummy mode.\nChange `dummy` to FALSE if this is not your intention.\n') else cat('Real analysis mode.\n')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('code/_includes/checkFunc.R')

cat('Checking for requirement...\n')
first_check()
cat('All check passed. Knitting... ')
suppressWarnings(rmarkdown::render('26TB_interim_analysis.Rmd', quiet = TRUE))
cat('Finished!')
system(paste('explorer', '26TB_interim_analysis.html'))
############################################ CODE ZONE ###################################################