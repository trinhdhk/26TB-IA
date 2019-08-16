first_check <- function(){
  if(Sys.getenv("R_ARCH") != "/i386") stop('Please switch to 32-bit R to continue.')
  auto.install <- function(pkgs){
    for (pkg in pkgs){
      if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    }
  }
  
  auto.install(c('knitr','Hmisc','survival','rmarkdown','plyr','tidyverse','haven','kableExtra','gdata','xtable','lubridate','readxl'))
}