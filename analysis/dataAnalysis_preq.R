########## FOR PREQ DATA ONLY ##########
setwd('/Users/mayala/Desktop/preq data')
subject_numbers <- c(1:18) 
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(0) # if this is 0, no scaling will be done in taskAnalysis()

##### function for plotting data #####
plotData <- function(){
  
  library(dplyr)
  library(ggplot2)
  
  filename <- sprintf('allTaggedData_n%d_%s.csv', length(subject_numbers), outfile_suffix)
  df <- read.csv(filename, header = TRUE)
  
  df <- df[-c(1),] # get rid of that random first row of NAs
  
  # first replace outlier-tagged trials with NA so they don't get plotted
  for (rowno in 1:nrow(df)) {
    if (df$isoutlier[rowno] == TRUE) {
      df$pathlength[rowno] <- NA
      df$pv_angle[rowno] <- NA
      df$pv_angle_n[rowno] <- NA
    }
  }
  
  tdf <- tbl_df(df) # convert to tibble for dplyr
  rotations <- sort(unique(tdf$prehome)) # use this because no-cursor trials are labeled rotation = 0
  
  for (rotationno in rotations) {
    print(rotationno)
    
    # filter by task
    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every trial -- use for plotting learning curves
      dfname<- sprintf('rotation%d_task%d_means', rotationno, taskno)
      print(dfname)
      
      ## NOTE: YOU NEED TO SUBTRACT baseline!!!! DO THIS HERE !! SUBTRACT PER ROTATION BASELINE VALUE
      taskmeans<- tdf %>% filter(task==taskno) %>% filter(prehome == rotationno) %>% group_by(participant) %>%
        group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                                      SEM_pl = SD_pl/sqrt(length(unique(participant))),
                                      Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
                                      SEM_pv = SD_pv/sqrt(length(unique(participant))))
      
      # plot each task's learning curve
      taskplot <- ggplot(data=taskmeans, aes(x=trial, y=Mean_pl)) +
        geom_line() + 
        geom_ribbon(aes(ymin=Mean_pl-SEM_pl, ymax=Mean_pl+SEM_pl),
                    alpha=0.4) +
        geom_line(data=taskmeans, aes(x=trial, y=Mean_pv), color="red") +
        geom_ribbon(aes(ymin=Mean_pv-SEM_pv, ymax=Mean_pv+SEM_pv),
                    alpha=0.4) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        ylim(-50, 50) +
        scale_colour_manual("", values = c("Mean_pl"="black", "Mean_pv"="red")) +
        ggtitle(dfname)
      
      print(taskplot)
    }
  }
}
##### function for statistical analysis
getStatistics <- function(){
  library(dplyr)
  
  filename <- sprintf('allTaggedData_n%d_%s.csv', length(subject_numbers), outfile_suffix)
  df <- read.csv(filename, header = TRUE)
  
  df <- df[-c(1),] # get rid of that random first row of NAs
  
  # first replace outlier-tagged trials with NA so they don't get analyzed
  for (rowno in 1:nrow(df)) {
    if (df$isoutlier[rowno] == TRUE) {
      df$pathlength[rowno] <- NA
      df$pv_angle[rowno] <- NA
      df$pv_angle_n[rowno] <- NA
    }
  }
  
  tdf <- tbl_df(df) # convert to tibble for dplyr
  ### NOTE: this is analyzing COLLAPSED rotations 
  
  # analyze adaptation
  baselinelast <- tdf %>% filter(task==0) %>% filter(trial %in% c(57,58,59)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  block1<- tdf %>% filter(task==3) %>% filter(trial %in% c(0,1,2)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  blocklast<- tdf %>% filter(task==7) %>% filter(trial %in% c(21,22,23)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  
  adaptdf<- rbind(block1,blocklast)
  
  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  RM_pv <- aov(pv ~ block + Error(participant/block), data=adaptdf)
  summary(RM_pv)
  RM_pl <- aov(pl ~ block + Error(participant/block), data=adaptdf)
  summary(RM_pl)
  
  #######################################################################
  ###############EDIT FOR REACH AFTEREFFECTS ANALYSIS####################
  #######################################################################
  
  # get df of just reach AEs first. baseline = -1; excludeAE = 0; includeAE = 1
  baselineAE <- tdf %>% filter(task==1) %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = -1)
  excludeAE <- tdf %>% filter(instruction=='exclude') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = 0)
  includeAE <- tdf %>% filter(instruction=='include') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = 1)
  
  AEdf <- rbind(baselineAE, excludeAE, includeAE)
  AEdf$block <- factor(AEdf$block)
  AEdf$participant <- factor(AEdf$participant)
  AEdf <-AEdf %>% mutate(group_instruction=1) # add group label for later 
  
  t.test(excludeAE$pv, baselineAE$pv, alternative ="greater" )
  t.test(excludeAE$pv-baselineAE$pv, mu=0, alternative ="greater" )
  t.test(includeAE$pv-baselineAE$pv, excludeAE$pv-baselineAE$pv, alternative ="greater" )
  
  ## visualize reach AEs
  tdf_NCs_rot <- tdf %>% filter(instruction == 'exclude' | instruction == 'include') %>% filter(trial == 0)
  ggplot(tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(prehome))) +
    geom_boxplot() +
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("PreQ Dual Group")

}
