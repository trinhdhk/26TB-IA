# /******************************************************************************                                                                               
#   Implement the consensus criteria for TBM diagnosis                                                                                                            
# (see folder 05TB/ConsensusDiagnosticCriteria/ for more documentation)                                                                                         
# 
# Author: Trinh Dong, based on SAS code by Marcel Wolbers; Thao Le                                                                                                                               
# 
# ******************************************************************************/                                                 
if(!exists('project.dir')) project.dir <- '..'
load(paste0(project.dir, '/Rdata/26TB_imported.Rdata'))
load(paste0(project.dir, '/Rdata/26TB_vad.Rdata'))

library(plyr)
library(tidyverse)

any2 <- function(..., na.rm = FALSE){
  pmap_lgl(list(...), any, na.rm = na.rm)
}

all2 <- function(..., na.rm = FALSE){
  pmap_lgl(list(...), all, na.rm = na.rm)
}

# Clinical criteria
clini_cirteria <- 
  inner_join(base, BASE, by = 'USUBJID') %>%
  mutate(
    #Clini_1
    CLINI_1 = 
      if_else(
        any2(
          (HEADACHE == 'Y' & HEADACHEDAY > 5),
          (IRRITABILITY == 'Y' & IRRITABILITYDAY > 5),
          (VOMIT == 'Y' & VOMITDAY > 5),
          (FEVER == 'Y' & FEVERDAY > 5),
          (NECKSTIF == 'Y' & NECKSTIFFDAY > 5),
          (SEIZURES == 'Y' & SEIZURESDAY > 5),
          (NEURODAY > 5),
          (ALTEREDCONSCIOUS  == 'Y' & CONSCIOUSDAY > 5),
          (LETHARGY == 'Y' & LETHARGYDAY > 5)
        ),
        4,
        0
      ),
    #Clini_2:
    CLINI_2 =
      if_else(
        pmap_lgl(
          list(
            COUGH == 'Y',
            WEIGHTLOSS == 'Y',
            NIGHTSWEATS == 'Y'
          ),
          any,
          na.rm = TRUE
        ), 
        2,
        if_else(
            any2(COUGH == 'N', COUGH == 'UNKNOWN') &
              any2(WEIGHTLOSS == 'N', WEIGHTLOSS == 'UNKNOWN') &
              any2(NIGHTSWEATS == 'N', NIGHTSWEATS == 'UNKNOWN'),
            0,
            as.numeric(NA)
          )
      ),
    #Clini_3
    CLINI_3 = 
      if_else(!is.na(LIVELUNGTB) & LIVELUNGTB == 'Y',
              2,
              if_else(!is.na(LIVELUNGTB) & any2(LIVELUNGTB == 'N', LIVELUNGTB == 'UNKNOWN'),
                      0,
                      as.numeric(NA))
      ),
    #Clini_4
    CLINI_4 = 
      if_else(
        any2(
          HEMIPLEGIA == 'Y',
          PARAPLEGIA == 'Y',
          TETRAPLEGIA == 'Y'
        ),
        1,
        if_else(
          all2(
            HEMIPLEGIA %in% c('N', 'UNKNOWN'),
            PARAPLEGIA %in% c('N', 'UNKNOWN'),
            TETRAPLEGIA %in% c('N', 'UNKNOWN')
          ),
          0,
          as.numeric(NA)
        )
      ),
    #Clini_5
    CLINI_5 = 
      if_else(
        CNP == 'Y',
        1,
        if_else(CNP %in% c('N', 'UNKNOWN'), 0, as.numeric(NA))
      ),
    #Clini_6
    CLINI_6 =
      if_else(
        any2(
          ALTEREDCONSCIOUS == 'Y',
          !is.na(GCS) & GCS < 15
        ),
        1,
        if_else(any2(
          ALTEREDCONSCIOUS %in% c('N', 'UNKNOWN'),
          GCS == 15),
          0, 
          as.numeric(NA))
      )
  ) %>%
  #Clini_total
  mutate(
    CLINI_SCORE = 
      select(., starts_with('CLINI_')) %>% 
      rowSums(., na.rm = TRUE) %>%
      if_else(.>6, 6, .)
  ) %>%
  select(USUBJID, starts_with('CLINI_'))
# 
# /* CSF criteria */                                                                                                                                            
#   * first merge vn and indo basecsf data;  

.BASECSF.ALL1 <- 
  bind_rows(BASECSF, BASECSF_JKT) %>% arrange(USUBJID)

.BASECSF.ALL2 <- 
  base %>% select(USUBJID, SITE) %>%
  right_join(.BASECSF.ALL1, by = 'USUBJID') %>%
  mutate(
    PROTEIN = PROTEIN * if_else(SITE == 'CIP', 0.01, 1),
    CSFGLUC = CSFGLUC * if_else(SITE == 'CIP', 0.0555, 1),
    PAIREDGLUC = PAIREDGLUC * if_else(SITE == 'CIP', 0.0555, 1)
  )

csf_criteria <- 
  .BASECSF.ALL2 %>%
  mutate(
    CSF_1 = 
      if_else(
        APPEARANCE == '1',
        1, 
        if_else(
          APPEARANCE != '',
          0,
          as.numeric(NA)
        )
      ),
    CSF_2 = if_else(CSFWBC >= 10 & CSFWBC <= 500, 1, 0),
    CSF_3 = if_else(CSFLYMLE > 50, 1, 0),
    CSF_4 = if_else(PROTEIN > 1, 1, 0),
    CSF_5 = 
      if_else(
        any2(
          CSFGLUC/PAIREDGLUC < 0.5,
          CSFGLUC < 2.2
        ),
        1,
        0
      )
  ) %>% 
  mutate(
    CSF_SCORE =
      select(., starts_with('CSF_')) %>%
      rowSums(., na.rm = TRUE) %>%
      if_else(.>4, 4, .)
    ) %>%
  select(USUBJID, starts_with('CSF_'))

#Evidence of TB elsewhere

#tube-1
.TUBE.1 <- 
  XRAYE %>%
  filter(VISITNUM == 'BASE') %>%
  mutate(
    TUBE_1 =
      if_else(
        XRAYRESULT == '2',
        2,
        if_else(
          XRAYRESULT == '3',
          4,
          if_else(
            XRAYRESULT != '',
            0,
            as.numeric(NA)
          )
        )
      )
  ) %>% select(USUBJID, TUBE_1)

#tube-2

#tube-3 AFB identified or Mycobacterium tuberculosis cultured from another source—;                                                                     
.OTHER <-
  BASELABOTH %>%
  mutate(
    TEST = 
      if_else(
        ZNSMEAR == '1' | MYCORESULT == '1',
        4,
        if_else(
          ZNSMEAR != '' | MYCORESULT != '',
          0,
          as.numeric(NA)
        )
      )
  ) %>% 
  select(USUBJID, TEST) %>%
  filter(!is.na(TEST))

.TUBE.3 <- .OTHER %>% group_by(USUBJID) %>% dplyr::summarise(TUBE_3 = max(TEST))

#Tube
tube_criteria <- 
  full_join(.TUBE.1, .TUBE.3, by = 'USUBJID') %>%
  mutate(
    TUBE_SCORE = 
      select(., TUBE_1, TUBE_3) %>% 
      rowSums(na.rm = TRUE) %>%
      if_else(.>4, 4, .))

# Criteria for definite TBM
.CONFIRMEDTBM <- 
  bind_rows(BASECSF, BASECSF_JKT) %>%
  mutate(
    CSF_AFB_POS = 
      if_else(
        ZNSMEAR == '1',
        1,
        if_else(
          ZNSMEAR != '',
          0,
          as.numeric(NA)
          )
      ),
    CSF_MTB_POS = 
      if_else(
        TBNAAT == '1',
        1,
        if_else(
          TBNAAT != '',
          0,
          as.numeric(NA)
        )
      ),
    CSF_XPERT_POS = 
      if_else(
        MYCORESULT == '1',
        1,
        if_else(
          MYCORESULT == '2',
          0,
          as.numeric(NA)
        )
      )
  ) %>%
  filter(!is.na(CSF_AFB_POS) | 
           !is.na(CSF_MTB_POS) | 
           !is.na(CSF_XPERT_POS)) %>%
  arrange(USUBJID)

confirmedTBM <- 
  .CONFIRMEDTBM %>%
  group_by(USUBJID) %>%
  dplyr::summarise(
    CSF_AFB_POS = max(CSF_AFB_POS),
    CSF_MTB_POS = max(CSF_MTB_POS),
    CSF_XPERT_POS = max(CSF_XPERT_POS)
  ) %>%
  arrange(USUBJID)

#Put all together
tbmdiagnosis <- 
  join_all(
    list(
      clini_cirteria,
      csf_criteria,
      tube_criteria,
      confirmedTBM,
      BASECSF %>% select(USUBJID, MYCORESULT, MYCOOTH),
      BASE %>% 
        select(USUBJID, HEADACHE, IRRITABILITY, VOMIT,
               FEVER, NECKSTIF, SEIZURES, ALTEREDCONSCIOUS,
               LETHARGY)
    ),
    
    type = 'left',
    by = 'USUBJID'
  ) %>%
  mutate(
    # Crude total score (treat missing scores as 0 but do NOT calculate total score if both csf_score and cere_score are missing)
    CRUDE_TOTAL_SCORE = 
      if_else(
        !is.na(CSF_SCORE),
        select(., CLINI_SCORE, CSF_SCORE, TUBE_SCORE) %>%
          rowSums(na.rm = TRUE),
        as.numeric(NA)
      ),
    #  crude final classification
    DIAGNOSTIC_SCORE = 
      if_else(
        any2(
          !is.na(CSF_AFB_POS) & CSF_AFB_POS == 1,
          !is.na(CSF_MTB_POS) & CSF_MTB_POS == 1,
          !is.na(CSF_XPERT_POS) & CSF_XPERT_POS == 1
        ),
        'definite TBM',
        if_else(
          !is.na(CRUDE_TOTAL_SCORE),
          if_else(
            all2(
              CRUDE_TOTAL_SCORE >= 10,
              CSF_SCORE >= 2,
              any2(
                HEADACHE == 'Y',
                IRRITABILITY == 'Y',
                VOMIT == 'Y',
                FEVER == 'Y',
                NECKSTIF == 'Y',
                SEIZURES == 'Y',
                ALTEREDCONSCIOUS == 'Y',
                LETHARGY == 'Y'
              )
            ),
            'probable TBM',
            if_else(
              all2(
                CRUDE_TOTAL_SCORE >= 6,
                any2(
                  HEADACHE == 'Y',
                  IRRITABILITY == 'Y',
                  VOMIT == 'Y',
                  FEVER == 'Y',
                  NECKSTIF == 'Y',
                  SEIZURES == 'Y',
                  ALTEREDCONSCIOUS == 'Y',
                  LETHARGY == 'Y'
                )
              ),
              'possible TBM',
              if_else(
                !is.na(MYCORESULT) & MYCORESULT == '3',
                'confirmed other diagnosis',
                NA_character_
              )
            )
          ),
          NA_character_
        )
      )
  ) %>%
  select(-MYCORESULT, -MYCOOTH)


# proc freq data=vad.tbmdiagnosis;                                                                                                                              
# table diagnostic_score;                                                                                                                                     
# run;

save(list = 'tbmdiagnosis', file = paste0(project.dir, '/Rdata/tbmdiagnosis.Rdata'))
