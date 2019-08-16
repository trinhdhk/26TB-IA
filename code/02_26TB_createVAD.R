# * Meta data ######
## Import Raw Interim Data to R
## Author Trinh Dong, based on SAS code by Thao Le
## Version 2019/02/14

# * Initialization ######
library(plyr)
library(tidyverse)
library(lubridate)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls()[!ls()%in%c('project.dir', 'dummy')])
if(!exists('project.dir')) project.dir <- '..'
load(paste0(project.dir, '/Rdata/26TB_imported.Rdata'))

is.na.blank <- function(x){
  ifelse(is.na(x) | as.character(x) == '', TRUE, FALSE)
}

first <- function(data, by){
  by <- quo_name(enquo(by))
  do.call(rbind,
          lapply(unique(data[[by]]),
                 function(this)
                   filter(data, data[[by]] ==this) %>% head(1)
                 )
  )
}

log_gen <- function(data, check, ...){
  check <- enquo(check)
  filter(data, !!check) %>% select(...)
}

.iVectorize <- function(fun, ...){
  args <- list(...)
  if(length(names(args))) do.call(Vectorize(fun, names(args)), args)
  else do.call(Vectorize(fun), args)
}

isCharacter <- Vectorize(FUN = is_character, vectorize.args = 'x', SIMPLIFY = TRUE, USE.NAMES = FALSE)

## bl.date from st_drg;  
.DATEDRUG <-
  STDR_STDRUG %>%
  filter(STDR_STDRUG_SEQ == 1) %>%
  select('USUBJID', 'DATESTART') %>%
  mutate(BL_DATE_DR = as.Date(DATESTART))

## Year,gender and diabetes from BASE
switch_vector <- Vectorize(
  function(v){
    switch(v,
           '003' = 'HTD',
           '013' = 'PNT',
           '044' = 'Cip',
           '047' = 'Persahabatan')
  }, 'v')

.log.SCR <- SCR %>% filter(USUBJID == setdiff(SCR$USUBJID, BASE$USUBJID) & ENROLLED == 'Y') %>% select(USUBJID, ENROLLED)

.DEMO <- 
  BASE %>%
  select(USUBJID, BIRTHYR, SEX, DIABETES, DATEADM, SITEID) %>%
  mutate(BL_DATE_DEMO = as.Date(DATEADM),
         SITE = switch_vector(SITEID)) %>%
  select(-DATEADM, -SITEID)

## Take tbm grade information
.BASE <-
  SCR %>%
  filter(ENROLLED == 'Y') %>%
  select(USUBJID, ENROLLED, DATERANDOM, DATEICF, TBMGRADE, RANNO) %>%
  mutate(BL_DATE_BASE = as.Date(DATERANDOM)) %>%
  select(-ENROLLED)


with(.BASE,
     for (i in seq_along(DATERANDOM))
       if (DATERANDOM[[i]] != DATEICF[[i]])
         cat('date ICF and date random do not match, UBSUBJID =',
             USUBJID[[i]], ', DATERANDOM =', DATERANDOM[[i]], ', DATEICF =', DATEICF[[i]],
             '\n'
             )
)

.log.BASE <- log_gen(.BASE, DATERANDOM != DATEICF, USUBJID, DATERANDOM, DATEICF)
# .log.BASE <-
#   .BASE %>% filter(DATERANDOM != DATEICF) %>% select(c(USUBJID, DATERANDOM, DATEICF))

## * CD4 and glucose from baselab, same unit for jkt
.BASELAB <- 
  BASELAB %>%
  select(USUBJID, CD4, FASTGLUC) %>%
  rbind(
    BASELAB_JKT %>%
      select(USUBJID, CD4, FASTGLUC)
  )

# * - ART from base table                                                                                                                                       
# * naive ART: KNOWHIV = "NO" + ARV == "3" + NEVERART == 1                                                                                                      
# * ART <= 3 month: use RTSTART (after Joe clean the data)                                                                                                      
# * ART > 3 month

.BASEART <- 
  BASE %>%
  select(USUBJID, KNOWHIV, ART, NEVERART, ARTSTARTUKN, ARTSTART) %>%
  mutate(ART_BL = if_else(KNOWHIV == 'N' | ART == '3' | NEVERART == 1,
                          'naive',
                          if_else(ARTSTARTUKN == 1, 
                                  'unknown',
                                  'art'))) # modify this when date information is available;

# *- Merge all base datasets;  
.TMP <- 
  left_join(.BASE, .DEMO, by = 'USUBJID') %>%
  list(.DATEDRUG, .BASELAB, .BASEART) %>%
  join_all(by = 'USUBJID', type = 'left') %>%
  dplyr::mutate(BL_DATE = 
           as.Date(
             if_else(
             !is.na(BL_DATE_DR), 
             BL_DATE_DR,
             if_else(!is.na(BL_DATE_DEMO),
                    BL_DATE_DEMO,
                    if_else(
                      !is.na(BL_DATE_BASE),
                      BL_DATE_BASE,
                      as.Date(NA))
                    )
             )
           )
         ) %>%
  mutate(AGE = year(BL_DATE) - BIRTHYR) %>%
  unique() #to remove duplicated 044-202

with(.TMP,
     for (i in seq_along(BL_DATE_DR)){
       if (isTRUE(BL_DATE_DR[i] != BL_DATE_DEMO[i]))
         cat(
           'date start drug and date in demo do not match',
           'USUBJID =', USUBJID[i],
           'BL_DATE_DR =', BL_DATE_DR[i],
           'BL_DATE_DEMO =', BL_DATE_DEMO[i],
           '\n'
           )
     }
)  

.log.TMP <- log_gen(.TMP, BL_DATE_DR != BL_DATE_DEMO, USUBJID, BL_DATE_DR, BL_DATE_DEMO)

.TMP <- .TMP %>% select(USUBJID, SEX, AGE, TBMGRADE, SITE, DIABETES, ART_BL, CD4, FASTGLUC, RANNO, BL_DATE, ARTSTART)

# *- Export to vad.base
base <- .TMP %>% arrange(USUBJID)

# *-----------------------------------------------------------------------------/;                                                                              
# * vad.ae: Clinical adverse events;                                                                                                                            
# *-----------------------------------------------------------------------------/;            

# *- Take code name from dictionary;

.AENAMES <- 
  Category %>%
  select(category, submissionvalue, text) %>%
  filter(tolower(category) == 'aeterm') %>%
  dplyr::rename(AENAME = text,
         AETERM = submissionvalue) %>%
  select(AENAME, AETERM)

# *-- Remove duplicates???????
.AENAMES <- unique(.AENAMES)

.CLINEVENTS_1 <- 
  AE %>%
  mutate(
    STARTDATE1 = as.Date(AESTDTC),
    STOPDATE1 = as.Date(AEENDTC)
  ) %>%
  select(
    USUBJID, AESEQ, STARTDATE1, STOPDATE1, 
    AETERM, OTHAESPEC, NEWNEUROSIGN, AETOXGR, 
    AESER, RELATESTDR, RELATESTDR, RELATETBM,
    RELATEATT, RELATEART, USAE) %>%
  arrange(AETERM)

.CLINEVENTS_2 <- 
  left_join(.CLINEVENTS_1, .AENAMES, by = 'AETERM') %>%
  arrange(AETERM)

# *- !!!! check for missing AEterm and grade and the relateness;                                                               

with(
  .CLINEVENTS_2, 
  for (i in 1:length(AETERM)){
    if (!is.na(STARTDATE1[i]) && is.na.blank(AETERM[i])) cat(' Missing AETERM, please check for patient - ', USUBJID[i], '\n')
    if (is.na.blank(AETOXGR[i])) cat(' Missing AE GRADE, please check for patient - ', USUBJID[i], '\n')
    if (is.na.blank(RELATESTDR[i])) cat(' Missing RELATESTDR, please check for patient - ', USUBJID[i], '\n')
    if (is.na.blank(AESER[i])) cat(' Missing AESER, please check for patient - ', USUBJID[i], '\n')
    if (is.na.blank(RELATEART[i])) cat(' Missing RELATEART, please check for patient - ', USUBJID[i], '\n')
  }
)

.log.CLINEVENTS_2.missingAETERM <- 
  log_gen(.CLINEVENTS_2, !is.na(STARTDATE1) & is.na.blank(AETERM), USUBJID, STARTDATE1, AETERM, AETOXGR, RELATESTDR, AESER, RELATEART)
.log.CLINEVENTS_2.missingAETOXGR <- 
  log_gen(.CLINEVENTS_2, is.na.blank(AETOXGR), USUBJID, STARTDATE1, AETERM, AETOXGR, RELATESTDR, AESER, RELATEART)
.log.CLINEVENTS_2.missingRELATESTDR <- 
  log_gen(.CLINEVENTS_2, is.na.blank(RELATESTDR), USUBJID, STARTDATE1, AETERM, AETOXGR, RELATESTDR, AESER, RELATEART)
.log.CLINEVENTS_2.missingAESER <- 
  log_gen(.CLINEVENTS_2, is.na.blank(AESER), USUBJID, STARTDATE1, AETERM, AETOXGR, RELATESTDR, AESER, RELATEART)
.log.CLINEVENTS_2.missingRELATEART <- 
  log_gen(.CLINEVENTS_2, is.na.blank(RELATEART), USUBJID, STARTDATE1, AETERM, AETOXGR, RELATESTDR, AESER, RELATEART)


# *-- add whether AE lead to TB drug interruption/stop;                                                                                                         
# *-- TB drug interuption;                                                                                                                                      
# * check if there is missing reason to stop when date end is given; 

.MISSINGREASON <- 
  TBDrug_TBDRUG %>%
  filter(!is.na(DATEEND), is.na.blank(REASONSTOP)) %>%
  select(USUBJID, TBDRUG, DATEEND, AENO) %>%
  arrange(USUBJID, AENO)

.TMP <- 
  TC %>%
  select(USUBJID, DATEDEATH) %>%
  right_join(.MISSINGREASON, 'USUBJID')

#print data .tmp
.PRINT.TMP <- .TMP

# * stop Tb drug due to (S)AE

.stopAE1 <- 
  TBDrug_TBDRUG %>%
  filter(REASONSTOP == '1' | !is.na(AENO)) %>%
  select(USUBJID, AENO, TBDRUG, DATEEND, REASONSTOP) %>%
  arrange(USUBJID, AENO)

# *-- ART interruption/stop;                                                                                                                                    
# * check if there is missing reason to stop when date end is given;          
.MISSINGREASON_ART <- 
  ART_ART %>%
  filter(!is.na(DATEEND) & is.na.blank(REASONSTOP)) %>%
  select(USUBJID, ARTDRUG, DATEEND, AENO) %>%
  arrange(USUBJID, AENO)

.TMP1 <- 
  TC %>%
  select(USUBJID, DATEDEATH) %>%
  right_join(.MISSINGREASON_ART, 'USUBJID')

.PRINT.TMP1 <- .TMP1

#print .tmp1

# * stop ART drug due to (S)AE
.stopAE2 <- 
  ART_ART %>%
  filter(REASONSTOP == '1' | !is.na(AENO)) %>%
  select(USUBJID, AENO, ARTDRUG, DATEEND) %>%
  arrange(USUBJID, AENO) %>%
  na.omit()

.stopAE <- union_all(.stopAE1, .stopAE2) %>% arrange(USUBJID, AENO) %>% dplyr::rename('AESEQ' = 'AENO')
for (i in 1:nrow(.stopAE)){
  with(.stopAE, 
       if (is.na(AESEQ[i])) cat('missing AE number when treatment interuption is due to AE ',
                             USUBJID[i], TBDRUG[i], ARTDRUG[i],DATEEND[i],
                             '\n')
       )
}

.log.stopAE <- 
  log_gen(.stopAE, is.na(AESEQ), USUBJID, TBDRUG, ARTDRUG, DATEEND)


# * this is the data with subjecid, aeno, tb and art drug that have to stop due to AE;                                                                     
# * export this data in case people want to know which drug need to stop due to AE;    

.stopAE.tomerge <- 
  .stopAE %>%
  filter(!is.na(AESEQ)) %>%
  select(-TBDRUG, -ARTDRUG) %>%
  arrange(USUBJID, AESEQ) %>%
  distinct()

# * Check date interruption and date death for those with missing AE number;                                                                                    
# *data _stopae_missingAE;                                                                                                                                      
# *merge _stopAE (where = (aeseq = .) in = a)                                                                                                                   
# vad.endpoint;                                                                                                                                       
# *by usubjid;                                                                                                                                                  
# *if a;                                                                                                                                                        
# *datestop = datepart(dateend);                                                                                                                                
# *format datestop date9.;                                                                                                                                      
# *keep usubjid datestop death_date;                                                                                                                            
# *run; * death date is the same with date stop for AE, so the reason for stopping drug is death, no adverse number needed;        

.CLINEVENTS_2 <- .CLINEVENTS_2 %>% arrange(USUBJID, AESEQ)


# * Old code duplicate the number of ae by the number of drug that might cause AE ;                                                                             
# * change to: merge AE table with AENo that leading to drug interuption;                                                                                       
# * Create another table for AE and drug that need to stop;
.CLINEVENTS_3 <- 
  .CLINEVENTS_2 %>%
  left_join(.stopAE.tomerge, by = c('USUBJID', 'AESEQ')) %>%
  mutate(TBARVDRUGSTOP = ifelse(!is.na(DATEEND), 'Yes', 'No'))

for (i in 1:nrow(.stopAE.tomerge)){
  if (!.stopAE.tomerge$USUBJID[i] %in% .CLINEVENTS_2$USUBJID)
    warning('USUBJID =', .stopAE.tomerge$USUBJID[i], 'not found in AE TABLE. Please check! \n')
  else{
    if (!.stopAE.tomerge$AESEQ[i] %in%
        .CLINEVENTS_2$AESEQ[.CLINEVENTS_2$USUBJID == .stopAE.tomerge$USUBJID[i]]) 
      cat('WARNING: AEname leading to drug interruption not found in AE dataset -- PLEASE CORRECT!!!!',
          'USUBJID =', .stopAE.tomerge$USUBJID[i],
          'AESEQ =', .stopAE.tomerge$AESEQ[i],
          '\n')
  }
}

.log.stopAE.tomerge.noAE <-
  log_gen(.stopAE.tomerge, !USUBJID %in% .CLINEVENTS_2$USUBJID, USUBJID)

.log.stopAE.tomerge.mismatch <- 
  .CLINEVENTS_2 %>% select(USUBJID, AESEQ) %>% unique() %>% mutate(IN = TRUE) %>%
  right_join(.stopAE.tomerge, by = c('USUBJID', 'AESEQ')) %>%
  log_gen(., is.na(IN), USUBJID, AESEQ)
 
# *- Make vad.ae table ; 
ae <-
  .CLINEVENTS_3 %>%
  select(USUBJID, AETERM, AENAME, OTHAESPEC, NEWNEUROSIGN,
         AETOXGR, AESER, USAE, RELATESTDR, RELATETBM, RELATEATT,
         RELATEART, TBARVDRUGSTOP, STARTDATE1, STOPDATE1) %>%
  dplyr::rename(
    GRADE = AETOXGR,
    SAE = AESER,
    STARTDATE = STARTDATE1,
    STOPDATE = STOPDATE1
  ) %>%
  mutate(
    AENAME.FULL =
      if_else(AETERM %in% c('10',''), OTHAESPEC, AENAME)
  ) %>%
  arrange(USUBJID, STARTDATE, AETERM)

ae_missing <-
  .stopAE %>% arrange(USUBJID)


# *--------------------------------------------------------------------------------/                                                                            
#   * vad.labae: laboratory abnormalities (grading as per appendix interim analysis plan);                                                                        
# * Convert lab values for patient in JKT to VN units                                                                                                           
# *--------------------------------------------------------------------------------/                                                                            
#   
#   *- Get hae and bio dataset;

.ALL.LAB <-
  base %>%
  select(USUBJID, BL_DATE, SEX, SITE) %>%
  inner_join(
    LAB_LAB %>%
      bind_rows(LAB_JKT_LAB) %>%
      select(USUBJID, DATEBL, HGB, SODIUM, FASTGLUCOSE,
             CREAT, ALT),
    by = 'USUBJID'
    ) %>%
  left_join(
    POTASSIUM_LOG %>%
      select(USUBJID, KALI),
    by = 'USUBJID'
  ) %>%
  dplyr::mutate(LABDATE = if_else(isCharacter(DATEBL), as.Date(ymd_hms(DATEBL, quiet = TRUE)), as.Date(DATEBL))) %>%
  select(-DATEBL) %>%
  dplyr::mutate(
    FASTGLUCOSE = 
      FASTGLUCOSE * if_else(SITE == 'Cip', 0.0555, 1)
  ) %>%
  dplyr::mutate(
    CREAT = 
      CREAT * if_else(SITE == 'Cip', 88.4, 1)
  )

.CREAT.BL <-
  BASELAB %>%
  bind_rows(BASELAB_JKT) %>%
  select(USUBJID, CREAT) %>%
  mutate(
    CREAT_BL = 
      CREAT * if_else(grepl('^044', CREAT, perl = TRUE), 88.4, 1)
  ) %>%
  select(-CREAT)

.ALL.LAB1 <- 
  .ALL.LAB %>%
  left_join(.CREAT.BL, by = 'USUBJID')

# *- Grade all relevant lab parameters (set grade to 0 if not grade 3 or 4);     
.LAB.GRADES <-
  .ALL.LAB1 %>%
  mutate(
    HGB_GRADE = 
      with(.ALL.LAB1,
           if_else(HGB < 6.5,
                   4,
                   if_else(HGB < 8,
                           3,
                           0)
                   )
           ),
    LOW_NA_GRADE = 
      with(.ALL.LAB1,
           if_else(SODIUM < 120,
                   4,
                   if_else(SODIUM <= 130,
                           3,
                           0)
                   )
           ),
    LOW_K_GRADE = 
      with(.ALL.LAB1,
           if_else(KALI < 2.5,
                   4,
                   if_else(KALI < 3.0,
                           3,
                           0)
                   )
           ),
    GLUCOSE_GRADE =
      with(.ALL.LAB1,
           if_else(FASTGLUCOSE < 1.7,
                   4,
                   if_else(FASTGLUCOSE <= 2.2,
                           3,
                           0)
                   )
           ),
    HIGH_CREA_GRADE =
      with(.ALL.LAB1,
           if_else(
             CREAT/88.4 > 6*1.36 | 
               (CREAT/88.4 > 6*1.13 & (SEX == 'F' | SEX == '')),
             4,
             if_else(
               !is.na(CREAT_BL) & CREAT > 3*CREAT_BL | 
                 CREAT/88.4 >= 3*1.36 |
                 (CREAT/88.4 >= 3*1.13 & (SEX == 'F' | SEX == '')),
               3,
               0
             ))
           ),
    HIGH_ALT_GRADE = 
      with(.ALL.LAB1,
           if_else(ALT > 20*40,
                   4,
                   if_else(ALT > 5*40,
                           3,
                           0)
                   )
           )
  ) %>% arrange(USUBJID)

# ** Function for deriving abnormalities at baseline and new abnormalities;                                                     

lab_abnorm <- function(dataset, labgrade, aename){
  .tmp <-
    dataset %>%
    filter(dataset[[labgrade]] %in% c(3,4))
  if (nrow(.tmp)){
    .tmp <-
      .tmp %>%
      mutate(AENAME = aename) %>%
      dplyr::rename(GRADE = !!labgrade) %>%
      select(USUBJID, LABDATE, AENAME, GRADE) %>%
      arrange(USUBJID, LABDATE)
    
    .id.list <- unique(.tmp$USUBJID)
    out <-
      do.call(rbind, 
              lapply(.id.list,
                     function(id) 
                       .tmp %>% 
                       filter(USUBJID == id) %>% 
                       head(1)))
    return(out)
  }
  out <- data.frame(USUBJID = NULL, 
                    LABDATE = NULL, 
                    AENAME = NULL,
                    GRADE = NULL)
  
  return(out)
}


.LABAE <- list()
.LABAE$HGB_GRADE <- lab_abnorm(.LAB.GRADES, 'HGB_GRADE', 'Haemoglobin - low')
.LABAE$LOW_NA_GRADE <- lab_abnorm(.LAB.GRADES, 'LOW_NA_GRADE', 'Sodium - low')
.LABAE$LOW_K_GRADE <- lab_abnorm(.LAB.GRADES, 'LOW_K_GRADE', 'Potassium - low')
.LABAE$GLUCOSE_GRADE <- lab_abnorm(.LAB.GRADES, 'GLUCOSE_GRADE', 'BlGlucose')
.LABAE$HIGH_CREA_GRADE <- lab_abnorm(.LAB.GRADES, 'HIGH_CREA_GRADE', 'Creatinine - high')
.LABAE$HIGH_ALT_GRADE <- lab_abnorm(.LAB.GRADES, 'HIGH_ALT_GRADE', 'ALT - high')

# ** put all togeter
labae <- bind_rows(.LABAE) %>% arrange(USUBJID, LABDATE, AENAME)
meta.nLAB <- length(unique(.ALL.LAB$USUBJID))

# *-------------------------------------------------------------------------------/;                                                                           
# * vad.endpoints: all endpoints needed for the interim analysis                   ;                                                                           
# *-------------------------------------------------------------------------------/;                                                                           
# 
# ***************** Primary  endpoint: Time to death in 12 months;                                                                                             
# ***************** Add time to death, add date withdrawn ;                                                                                                    
# 
# 
# *- Create date last known to be alive (to be used if patient did not withdraw or month 8 not yet filled out);                                                
# *  Include fup visits, hema, chemistry, and csf;   

LAB_CSF <- LAB_CSF %>% arrange(USUBJID)
LAB_JKT_CSF <- LAB_JKT_CSF %>% arrange(USUBJID)
LAB_LAB <- LAB_LAB%>% arrange(USUBJID)
LAB_JKT_LAB <- LAB_JKT_LAB %>% arrange(USUBJID)
LAB_STOOL <- LAB_STOOL %>% arrange(USUBJID)
LAB_JKT_STOOL <- LAB_JKT_STOOL %>% arrange(USUBJID)
IPFU <- IPFU %>% arrange(USUBJID)
IPFU_EXTRA <- IPFU_EXTRA %>% arrange(USUBJID)

.ALLDATES <- 
  bind_rows(
    IPFU %>% select(USUBJID, DATEFU),
    IPFU_EXTRA %>% select(USUBJID, DATEFU),
    OPFU %>% select(USUBJID, DATEFU),
    LAB_STOOL %>% select(USUBJID, DATESTSAMPLE) %>% dplyr::rename(DATEFU = DATESTSAMPLE),
    LAB_JKT_STOOL %>% select(USUBJID, DATESTSAMPLE) %>% dplyr::rename(DATEFU = DATESTSAMPLE),
    LAB_CSF %>% select(USUBJID, DATESAMPLECSF) %>% dplyr::rename(DATEFU = DATESAMPLECSF),
    LAB_JKT_CSF %>% select(USUBJID, DATESAMPLECSF) %>% dplyr::rename(DATEFU = DATESAMPLECSF),
    LAB_LAB %>% select(USUBJID, DATEBL) %>% dplyr::rename(DATEFU = DATEBL),
    LAB_JKT_LAB %>% select(USUBJID, DATEBL) %>% dplyr::rename(DATEFU = DATEBL),
    TC %>% select(USUBJID, DATELASTFU) %>% dplyr::rename(DATEFU = DATELASTFU),
    TC %>% select(USUBJID, DATEWITHDRAW) %>% dplyr::rename(DATEFU = DATEWITHDRAW)
  ) %>%
  mutate(DATE = if_else(isCharacter(DATEFU), as.Date(ymd_hms(DATEFU, quiet = TRUE)), as.Date(DATEFU))) %>%
  filter(!is.na(DATE)) %>%
  select(-DATEFU) %>%
  arrange(USUBJID, DATE)
  
.LASTFUPDATE <- 
  .ALLDATES %>%
  group_by(USUBJID) %>%
  filter(DATE == max(DATE)) %>%
  distinct() %>%
  mutate(LASTFUP_DATE = as.Date(DATE)) %>%
  select(-DATE)

#missing LASTFUDATE
.CHECK <-
  base %>%
  anti_join(.LASTFUPDATE, by = 'USUBJID') %>%
  mutate(LASTFUP_DATE = rep(NA, nrow(.))) %>%
  select(USUBJID, LASTFUP_DATE, everything())
         
#print check
.PRINT.CHECK <- .CHECK

# *- Create time to death outcome;  
.TIMETODEATH <- 
  join_all(
    list(
      base %>% select(USUBJID, BL_DATE),
      TC %>% select(USUBJID, REASONLOSTFU, DATEDEATH, DATELASTFU, DATEWITHDRAW, DATELASTALIVE),
      .LASTFUPDATE %>% select(USUBJID, LASTFUP_DATE)
    ),
    by = 'USUBJID',
    type = 'left'
  ) %>%
  mutate(
    DEATH_DATE = 
      if_else(!is.na(REASONLOSTFU) & REASONLOSTFU == '8' | !is.na(DATEDEATH),
              if_else(isCharacter(DATEDEATH), as.Date(ymd_hms(DATEDEATH, quiet = TRUE)), as.Date(DATEDEATH)),
              as.Date(NA)),
    LASTFUP_DATE = 
      if_else(!is.na(REASONLOSTFU) & REASONLOSTFU == '8' | !is.na(DATEDEATH),
             as.Date(NA),
             LASTFUP_DATE),
    LASTALIVE_DATE = 
      if_else(!is.na(REASONLOSTFU) & REASONLOSTFU == '8' | !is.na(DATEDEATH),
              if_else(isCharacter(DATELASTALIVE), as.Date(ymd_hms(DATELASTALIVE, quiet = TRUE)), as.Date(DATELASTALIVE)),
              as.Date(NA)),
    TT_DEATH = 
      if_else(!is.na(REASONLOSTFU) & REASONLOSTFU == '8' | !is.na(DATEDEATH),
              if_else(!is.na(DEATH_DATE), DEATH_DATE - BL_DATE + 1, LASTALIVE_DATE - BL_DATE + 1),
              if_else(!is.na(LASTFUP_DATE), LASTFUP_DATE - BL_DATE + 1, LASTALIVE_DATE - BL_DATE + 1)
      ),
    EV_DEATH = 
      if_else(!is.na(REASONLOSTFU) & REASONLOSTFU == '8' | !is.na(DATEDEATH),
              1,
              0),
    LFU = 
      if_else(!is.na(REASONLOSTFU) & REASONLOSTFU %in% c('5','6') | !is.na(DATELASTFU),
              'Yes',
              'No')
  )

# Censor everyone at 12months (day 365)
.log.TIMETODEATH <- log_gen(.TIMETODEATH, TT_DEATH > 365  & EV_DEATH == 1, USUBJID, DEATH_DATE, TT_DEATH, EV_DEATH)

with(.TIMETODEATH,
     for (i in 1:nrow(.TIMETODEATH)){
       if (isTRUE(TT_DEATH[i] > 365)){
         if (isTRUE(EV_DEATH[i] == 1)){
           cat('WARNING: Death censored because it occurred after day 365:',
               'USUBJID =', USUBJID[i], 
               'DEATH_DATE =', DEATH_DATE[i],
               'TT_DEATH =', TT_DEATH[i],
               'EV_DEATH =', EV_DEATH[i],
               '\n')
           
           DEATH_DATE[i] <- NA
         }
         
         TT_DEATH[i] <- 365
         EV_DEATH[i] <- 0
         LASTFUP_DATE[i] <- BL_DATE[i] + 364
       }
     }
)

#PROC PRINT timetodeath
.PRINT.TIMETODEATH <-
  .TIMETODEATH %>%
  filter(EV_DEATH == 1 & is.na(TT_DEATH) & is.na(LASTFUP_DATE))

#Create table endpoint
endpoint <- 
  .TIMETODEATH %>%
  select(USUBJID, TT_DEATH, EV_DEATH, DEATH_DATE,
         LASTFUP_DATE, LFU) %>%
  arrange(USUBJID)



#****************** Open dex: OD;
.OD.DURING <-
  STDR %>%
  filter(
    REASONSTOP == '2'|
      !is.na(DATEOPEN)|
      FALLGCS != 0|
      NEUDEFICIT != 0|
      SEIZURES != 0|
      TUBERCULOMA != 0|
      HYDROCEPHALUS != 0|
      VASCULITIS != 0|
      INFARCTS != 0|
      PARADOXICAL != 0|
      NEUROIRIS != 0|
      EXTRANEUIRIS != 0|
      REASONOPENOTH != ''
    ) %>%
  mutate(DATESTART = if_else(isCharacter(DATEOPEN), as.Date(ymd_hms(DATEOPEN, quiet = TRUE)), as.Date(DATEOPEN))) %>%
  select(USUBJID, DATESTART, FALLGCS, NEUDEFICIT, SEIZURES,
         TUBERCULOMA, HYDROCEPHALUS, VASCULITIS,
         INFARCTS, PARADOXICAL, NEUROIRIS, EXTRANEUIRIS,
         REASONOPENOTH) 

OD.DURING <-
  base %>%
  select(USUBJID, BL_DATE) %>%
  right_join(.OD.DURING, by = 'USUBJID') %>%
  mutate(TTOPDEX = DATESTART - BL_DATE)

# * opend dex after treatment;                                                                                                                                  
# * tbm grade 1: 6 weeks.                                                                                                                                       
# * tbm grade 2&3: 8 weeks;                                                                                                                                     
# 
# * complete study drug or not;

.COMPLETE.STDRUG <- 
  STDR_STDRUG %>%
  filter(!is.na(DATEEND)) %>%
  left_join(base, by = 'USUBJID') %>%
  mutate(
    DATEEND_STDRUG = if_else(isCharacter(DATEEND), as.Date(ymd_hms(DATEEND, quiet = TRUE)), as.Date(DATEEND)),
    COMPLETESTDRUG = 
      if_else((TBMGRADE == 1 & DRUG == 'W6')|
                (TBMGRADE %in% c(2,3) & DRUG == 'W8'),
              1,
              0)
  ) %>%
  select(USUBJID, DRUG, DATEEND_STDRUG, TBMGRADE, COMPLETESTDRUG)

.OD.AFTER1 <- 
  ANTINFDRUG_ANTIFLAMDRUG %>%
  filter(!is.na(DATESTART)) %>%
  select(USUBJID, REASONSTART, DATESTART) %>%
  inner_join(.COMPLETE.STDRUG %>%
               filter(COMPLETESTDRUG == 1),
             by = 'USUBJID') %>%
  mutate(
    DATESTART_PRED = as.Date(DATESTART),
    DURATION = DATESTART_PRED - DATEEND_STDRUG
  ) %>%
  select(-DATESTART)

with(.OD.AFTER1,
     for (i in 1:nrow(.OD.AFTER1)){
       if (DURATION[i] <= 0) {
         cat('WARNING: Date start Dex was before date end study drug, please check',
             'USUBJID =', USUBJID[i], '\n')
       }
     }
)

.log.OD.AFTER1 <- log_gen(.OD.AFTER1, DURATION <=0 , USUBJID, DATESTART_PRED, DATEEND_STDRUG, DURATION)

.OD.AFTER1 <- 
  .OD.AFTER1 %>% 
  filter(DURATION > 0) %>%
  arrange(USUBJID)

.OD.AFTER2 <- first(.OD.AFTER1, USUBJID)

OD.AFTER <- 
  base %>%
  select(USUBJID, BL_DATE) %>%
  right_join(.OD.AFTER2, by = 'USUBJID') %>%
  mutate(TTOPDEX = DATESTART_PRED - BL_DATE)

with(OD.AFTER,
     for (i in 1:nrow(OD.AFTER))
       if (TTOPDEX[i] < 42)
         cat('WARNING: time to open dex less than 6 weeks since the start of TB treatment- please check',
             'USUBJID =', USUBJID[i], '\n')
)

.log.OD.AFTER <- 
  log_gen(OD.AFTER, TTOPDEX < 42, USUBJID, DATESTART_PRED, BL_DATE, TTOPDEX)

OD.AFTER <- OD.AFTER %>% filter(TTOPDEX >= 42)

# Create tables
odduring <- OD.DURING %>% arrange(USUBJID)
odafter <- OD.AFTER %>% arrange(USUBJID)


# *-------------------------------------------------------------------------------/;                                                                            
# * vad.dili: DILI (drug induced liver injury) study                   ;                                                                                        
# *-------------------------------------------------------------------------------/;                                                                           

.DILI <-
  DILI %>%
  filter(DILIPAR != 1 & !is.na.blank(DATESTRATEGY)) %>%
  mutate(
    ALT_DILI = ALT,
    BILI_DILI = BILDIR * if_else(BILIU == 'mgdL', 17.1, 1),
    INR_DILI = INR,
    DILI_DATE = as.Date(DATESTRATEGY)
  ) %>%
  select(USUBJID, STUDYID, STRATEGY, DILI_DATE, ALT_DILI, BILI_DILI, INR_DILI)

.DILI.27TB <- 
  DILI.27TB %>%
  filter(DILIPAR != 1 & !is.na.blank(DATESTRATEGY)) %>%
  mutate(
    ALT_DILI = ALT,
    BILI_DILI = BILDIR,
    INR_DILI = INR,
    DILI_DATE = as.Date(DATESTRATEGY)
  ) %>%
  select(USUBJID, STUDYID, STRATEGY, DILI_DATE, ALT_DILI, BILI_DILI, INR_DILI)

.DILI.PEAK <- 
  DILIOC %>%
  mutate(
    ALT_PEAK = ALT,
    BILI_PEAK = BILDIR * if_else(BILIU == 'mgdL', 17.1, 1),
    INR_PEAK = INR
  ) %>%
  select(USUBJID, STUDYID, ALT_PEAK, BILI_PEAK, INR_PEAK)

.DILI.PEAK.27TB <- 
  DILIOC.27TB %>%
  dplyr::rename(ALT_PEAK = ALT, BILI_PEAK = BILDIR, INR_PEAK = INR)

.ART.DATE <- 
  ART_ART %>%
  select(USUBJID, DATESTART) %>%
  arrange(USUBJID, DATESTART)

with (.ART.DATE,
      for (i in 1:nrow(.ART.DATE))
        if (is.na.blank(DATESTART[i]))
          cat('WARNING: missing ART date start for patient-',
              'USUBJID =', USUBJID[i], '\n')
      )

.log.ART.DATE <- 
  log_gen(.ART.DATE, is.na.blank(DATESTART), USUBJID, DATESTART)

.ART.DATE <- 
  .ART.DATE %>% 
  filter(!is.na.blank(DATESTART)) %>%
  mutate(DATESTART = gsub('^00/00/', '20', DATESTART, perl = TRUE)) %>%
  mutate(DATEART = if_else(isCharacter(DATESTART), as.Date(parse_date_time(DATESTART, orders=c('y','my', 'dmy'), quiet = TRUE)), as.Date(DATESTART))) %>%
  select(-DATESTART) %>%
  first(USUBJID)

.TB.DATE <- 
  TBDrug_TBDRUG %>%
  select(USUBJID, DATESTART) %>%
  mutate(DATE_TBDRUG = if_else(isCharacter(DATESTART), as.Date(ymd_hms(DATESTART, quiet = TRUE)), as.Date(DATESTART))) %>%
  select(-DATESTART) %>%
  first(USUBJID)

.DILI.ENDPOINT <-
  join_all(
    list(endpoint %>%
           select(USUBJID, EV_DEATH, TT_DEATH, DEATH_DATE),
         base %>%
           select(USUBJID, BL_DATE, ART_BL),
         .ART.DATE,
         .TB.DATE),
    
    by = 'USUBJID',
    type = 'left'
  ) %>%
  right_join(.DILI,
             by = 'USUBJID') %>%
  mutate(
    DEATH_DURINGDILI = 
      if_else(EV_DEATH == 1 & TT_DEATH <= 60,
              1,
              0),
    DEATHDAYS_TB = 
      if_else(!is.na(DEATH_DATE),
              as.numeric(DEATH_DATE - DATE_TBDRUG + 1),
              as.numeric(NA)),
    DEATHDAYS_ART = 
      if_else(!is.na(DEATH_DATE),
              as.numeric(DEATH_DATE - DATEART + 1),
              as.numeric(NA)),
    DEATHDAYS_DILI =
      if_else(!is.na(DEATH_DATE),
              as.numeric(DEATH_DATE - DILI_DATE + 1),
              as.numeric(NA))
  )

.DILI.ALL <- 
  join_all(
    list(
      .DILI.PEAK,
      .DILI.ENDPOINT
    ),
    type = 'inner',
    by = c('USUBJID', 'STUDYID')
  )

.DILI.ALL.27TB <-
  plyr::join_all(
    list(
      .DILI.27TB,
      TC.27TB %>% select(USUBJID, DATEDEATH, DATELASTALIVE, REASONLOSTFU),
      STDR_STDRUG.27TB %>% select(USUBJID, DATESTART) %>% first(USUBJID),
      .DILI.PEAK.27TB %>% select(-STUDYID)
    ),
    type = 'left',
    by = 'USUBJID'
  ) %>%
  mutate(
    DATE_TBDRUG = if_else(isCharacter(DATESTART), as.Date(ymd_hms(DATESTART, quiet = TRUE)), as.Date(DATESTART)),
    DATE_DEATH = if_else(REASONLOSTFU == '8' | !is.na(DATEDEATH),
                        if_else(!is.na(DATEDEATH), 
                                if_else(isCharacter(DATEDEATH), as.Date(ymd_hms(DATEDEATH, quiet = TRUE)), as.Date(DATEDEATH)),
                                if_else(isCharacter(DATELASTALIVE), as.Date(ymd_hms(DATELASTALIVE, quiet = TRUE)), as.Date(DATELASTALIVE))),
                        as.Date(NA)),
    DEATHDAYS_TB = if_else(!is.na(DATE_DEATH), as.numeric(DATE_DEATH - DATE_TBDRUG + 1), as.numeric(NA)),
    DEATHDAYS_DILI = if_else(!is.na(DATE_DEATH), as.numeric(DATE_DEATH - DILI_DATE + 1), as.numeric(NA)),
    EV_DEATH = if_else(!is.na(DATEDEATH) | REASONLOSTFU == '8', 1, 0),
    DEATH_DURINGDILI = if_else(DEATHDAYS_TB <= 60 & REASONLOSTFU == '8', 1, 0)
  )


dili <- .DILI.ALL %>% arrange(USUBJID)
dili.27TB <- .DILI.ALL.27TB %>% arrange(USUBJID)
# dili.combined <- .DILI.COMBINE %>% dplyr::arrange(STUDYID, USUBJID)

################################################################
#                                                              #
#          Add BASELINE INFORMATION for 27TB DILI              #
#                                                              #
################################################################

.BASE.27TB <- 
  inner_join(
    BASE.27TB %>% select(USUBJID, DATEADM),
    SCR.27TB %>% select(USUBJID, DATERANDOM, RANNO, LTA4HRESULT),
    by = 'USUBJID'
  ) %>%
  left_join(
    STDR_STDRUG.27TB %>% filter(STDR_STDRUG_SEQ == 1) %>% select(USUBJID, DATESTART),
    by = 'USUBJID'
  ) %>%
  mutate(
    BL_DATE = if_else(!is.na(DATESTART), DATESTART,
                      if_else(!is.na(DATEADM), DATEADM, DATERANDOM))
  ) %>%
  right_join(.DILI.27TB, by = 'USUBJID') %>% 
  select(USUBJID, RANNO, BL_DATE, LTA4HRESULT) 

base.27TB <- .BASE.27TB %>% arrange(USUBJID)

### This table is only needed in 3rd IA
opfu.27TB <- OPFU.27TB %>% select(USUBJID, DATEFU) %>% mutate(DATEFU = if_else(isCharacter(DATEFU), as.Date(ymd_hms(DATEFU)), as.Date(DATEFU)))

  
save(
  list = c('ae','ae_missing','base', 'base.27TB', 'SCR', 'dili', 'dili.27TB','endpoint', 'labae','odafter','odduring', 'meta.nLAB', 'opfu.27TB'),
  file = paste0(project.dir, '/Rdata/26TB_vad.Rdata')
)

save(
  list = ls(all.names = TRUE)[grepl('(^\\.log)|(^\\.PRINT)',ls(all.names = TRUE), perl = T)],
  file = paste0(project.dir,'/Rdata/26TB_vad_log.Rdata')
)
