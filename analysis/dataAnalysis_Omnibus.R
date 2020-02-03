##### OMNIBUS TESTS WITH BOTH SINGLE AND DUAL DATA
#setwd('/Users/mayala/Desktop/seq data')
#setwd('/Users/mayala/Desktop/conseq data')
#setwd('/Users/mayala/Desktop/static data')
setwd('/Users/mayala/Desktop/explicit data')
#setwd('/Users/mayala/Desktop/preq data')

##### 
getOmniStats <- function(){
  
  library(dplyr)
  library(tidyr)
  library(Hmisc)
  library(ez)
  
  ################################
  ################################
  ################################
  
  # LOAD SINGLE LCs
  CW_ONLY <- read.csv("SINGLE_LCs_1.csv", header = TRUE)
  CCW_ONLY <- read.csv("SINGLE_LCs_-1.csv", header = TRUE)
  CCW_ONLY$participant <- CCW_ONLY$participant + 10 # EZ NEEDS DISTINCT ID FOR EVERY DISTINCT PARTICIPANT  
  
  DUAL_CW <- read.csv("DUAL_LCs_1.csv", header = TRUE)
  DUAL_CW$participant <- DUAL_CW$participant + 20 # EZ NEEDS DISTINCT ID FOR EVERY DISTINCT PARTICIPANT
  DUAL_CCW <- read.csv("DUAL_LCs_-1.csv", header = TRUE)
  DUAL_CCW$participant <- DUAL_CCW$participant + 20 # EZ NEEDS DISTINCT ID FOR EVERY DISTINCT PARTICIPANT
  
  CW_ONLY <- CW_ONLY %>% mutate(group = 'CW', dual = 0) %>% select(-rotation)
  CCW_ONLY <- CCW_ONLY %>% mutate(group = 'CCW', dual = 0) %>% select(-rotation)
  DUAL_CW <- DUAL_CW %>% mutate(group = 'DUAL_CW', dual = 1)
  DUAL_CCW <- DUAL_CCW %>% mutate(group ='DUAL_CCW', dual = 1)
  
  CCW_ONLY$blockmean <- CCW_ONLY$blockmean*-1 # FLIP THE SIGNS
  DUAL_CCW$blockmean <- DUAL_CCW$blockmean*-1
  
  omni <- rbind(CW_ONLY, CCW_ONLY, DUAL_CW, DUAL_CCW)
  
  omni <- omni %>% drop_na(blockmean) # NOTE FIX THIS MISSING VALUE!!!
  
  omni$block <- as.factor(omni$block)
  
  mod1 <- ezANOVA(data = omni,
                  dv = blockmean, # pv angle
                  wid = participant,
                  within = block, # Dual-CW, Dual-CCW, Single-CW, Single-CCW
                  between = .(dual), # Dual vs. Single Conditions
                  detailed = TRUE,
                  return_aov = TRUE)
  
  print(mod1)

}