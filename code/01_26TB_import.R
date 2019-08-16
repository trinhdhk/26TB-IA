# * Meta data ######
## Import Raw Interim Data to R
## Author Trinh Dong, based on SAS code by Thao Le
## Version 2019/02/14

# * Initialization ######
library(dplyr)
library(readxl)
rm(list=ls()[!ls()%in%c('project.dir', 'dummy')])
if(!exists('project.dir')) project.dir <- '..'
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Function to import data from Acces
import.Access <- function(dataPath){
  if (Sys.info()[['sysname']] == 'Windows'){
    library(RODBC)
    channel <- RODBC::odbcConnectAccess(dataPath)
    tableNames <- subset(sqlTables(channel), TABLE_TYPE == "TABLE")$TABLE_NAME
    tables <- list()
    for (i in tableNames) tables[[i]] <- sqlFetch(channel, i, as.is = TRUE)
    odbcClose(channel)
  } else {
    library(Hmisc)
    tables <- Hmisc::mdb.get(dataPath, dateformat = '%Y/%m/%d')
    tableNames <- names(tables)
  }
  
  return(list(tables = tables, tableNames = tableNames))
}

import.Excel <- function(dataPath){
  requireNamespace('readxl')
  sheets <- readxl::excel_sheets(dataPath)
  tables <- sapply(sheets, 
                   function(sheet) readxl::read_excel(dataPath, sheet = sheet), 
                   simplify = FALSE,
                   USE.NAMES = TRUE)
  return(list(tables = tables, tableNames = sheets))
}

import.Data <- function(dataPath){
  if(tools::file_ext(dataPath) %in% c('xls','xlsx')) import.Excel(dataPath) else import.Access(dataPath)
}


dataFolder <- paste0(project.dir, '/data')
dataFolder.fileList <- list.files(dataFolder)
dataName <- grep('_26TB', dataFolder.fileList, value = TRUE)[1]
dataName.27TB <- grep('_27TB', dataFolder.fileList, value = TRUE)[1]
dataPath <- paste(dataFolder, dataName, sep = '/')
dataPath.27TB <- paste(dataFolder, dataName.27TB, sep = '/')

# Import data from 26TB ######

#### Enable these 4 lines if wanna choose manually
# cat('Import 26TB file')
# dataPath <- file.choose()
# dataName <- basename(dataPath)
.tb26 <- import.Data(dataPath)
for (i in .tb26$tableNames) assign(i, .tb26$tables[[i]])

#Import data from 27TB #####

#### Enable these 3 lines if wanna choose manually
# cat('Import 27TB file')
# dataPath.27TB <- file.choose()


.tb27 <- import.Data(dataPath.27TB)
.tb27$tableNames <- c('DILI', 'DILIOC', 'BASE', 'SCR', 'STDR_STDRUG', 'TC', 'TBDrug_TBDRUG',  'OPFU')
for (i in .tb27$tableNames) assign(paste0(i,'.27TB'), .tb27$tables[[i]])


# Import data from 

#Replace variable names by column labels #####
for (table in c(.tb26$tableNames, .tb27$tableNames)){
  .table.get <- get(table)
  names(.table.get) <- ifelse(Hmisc::label(.table.get) != '', Hmisc::label(.table.get), names(.table.get))
  assign(table, .table.get)
}

# Sort data #####
BASE <- BASE %>% arrange(USUBJID, DATEADM)
AE <- AE %>% arrange(USUBJID, AESTDTC)
BASELAB <- BASELAB %>% arrange(USUBJID, DATESAMPLEBL)
BASELAB_JKT <- BASELAB_JKT %>% arrange(USUBJID, DATESAMPLEBL)
BASECSF <- BASECSF %>% arrange(USUBJID, DATESAMPLECSF)
BASECSF_JKT <- BASECSF_JKT %>% arrange(USUBJID, DATESAMPLECSF)
BASELABOTH <- BASELABOTH %>% arrange(USUBJID, DATESAMPLE)
XRAY_XRAY <- XRAY_XRAY %>% arrange(USUBJID, DSDTC)
STDR_STDRUG <- STDR_STDRUG %>% arrange(USUBJID, DATESTART)
ANTINFDRUG_ANTIFLAMDRUG <- ANTINFDRUG_ANTIFLAMDRUG %>% arrange(USUBJID, DATESTART)
TBDrug_TBDRUG <- TBDrug_TBDRUG %>% arrange(USUBJID, DATESTART)
ART_ART <- ART_ART %>% arrange(USUBJID, DATESTART)
IPFU <- IPFU %>% arrange(USUBJID, DATEFU)
LAB_STOOL <- LAB_STOOL %>% arrange(USUBJID, DATESTSAMPLE)
LAB_JKT_STOOL <- LAB_JKT_STOOL %>% arrange(USUBJID, DATESTSAMPLE)
LAB_CSF <- LAB_CSF %>% arrange(USUBJID, DATESAMPLECSF)
LAB_JKT_CSF <- LAB_JKT_CSF %>% arrange(USUBJID, DATESAMPLECSF)
TC <- TC %>% arrange(USUBJID, DATECOMPLETE)
DILI <- DILI %>% arrange(USUBJID, DATESTRATEGY)

BASE.27TB <- BASE.27TB %>% arrange(USUBJID, DATEADM)
DILI.27TB <- DILI.27TB %>% arrange(USUBJID, DATESTRATEGY)
TBDrug_TBDRUG.27TB <- TBDrug_TBDRUG.27TB %>% arrange(USUBJID, DATESTART)



#Save as Rdata #####
save(list = ls(), file = paste0(project.dir, "/Rdata/26TB_imported.Rdata"))
