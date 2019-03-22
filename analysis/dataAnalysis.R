
######### MAKE CHANGES HERE ############
########################################
setwd('/Users/mayala/Desktop/conseq data/selected data')

subject_numbers <- c(1:7, 9:31)
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()

######### END OF CHANGES  ##############
########################################

##### function for plotting data
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
  rotations <- sort(unique(tdf$garage_location)) # use this because no-cursor trials are labeled rotation = 0
  
  for (rotationno in rotations) {
    print(rotationno)

  # filter by task
    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every trial -- use for plotting learning curves
      dfname<- sprintf('rotation%d_task%d_means', rotationno, taskno)
      print(dfname)
      
      ## NOTE: YOU NEED TO SUBTRACT baseline!!!! DO THIS HERE !! SUBTRACT PER ROTATION BASELINE VALUE
      taskmeans<- tdf %>% filter(task==taskno) %>% filter(garage_location == rotationno) %>% group_by(participant) %>% group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),
                  SD_pl = sd(pathlength, na.rm=TRUE),
                  SEM_pl = SD_pl/sqrt(length(unique(participant))),
                  Mean_pv = mean(pv_angle, na.rm=TRUE), 
                  SD_pv = sd(pv_angle, na.rm=TRUE),
                  SEM_pv = SD_pv/sqrt(length(unique(participant))))
    
      # plot each task learning curve
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
  ### -- also still need to subtract baseline 
  
  ## #### NOTE: to do - add for-loop here for separating analysis for implicit & explicit experiments
  # analyze adaptation
  
  block1<- tdf %>% filter(task==3) %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  blocklast<- tdf %>% filter(task==5) %>% filter(trial %in% c(357,358,359)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
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
  consequence_baselineAE <- tdf %>% filter(task==1) %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = -1)
  consequence_excludeAE <- tdf %>% filter(instruction=='exclude') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = 0)
  consequence_includeAE <- tdf %>% filter(instruction=='include') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = 1)
  
  consequence_AEdf <- rbind(consequence_baselineAE, consequence_excludeAE, consequence_includeAE)
  consequence_AEdf$block <- factor(consequence_AEdf$block)
  consequence_AEdf$participant <- factor(consequence_AEdf$participant)
  consequence_AEdf <- consequence_AEdf %>% mutate(group_instruction=1) # add group label for later 
  
  #########################################################################################
  ############################# E/I COMPARISONS ACROSS GROUPS #############################
  #########################################################################################
  
  ## did any implicit learning occur?
  ## compare reach AEs across baseline NC and exclude-strategy NC b/w instruct and non-instruct groups
  both_groups_AEdf <- rbind(static_AEdf, instructed_AEdf)
  
  baseline_exclude_NCs <- both_groups_AEdf %>% filter(block == 0 | block == -1)
  implicit_pv <- aov(pv ~ block + Error(participant/block), data=baseline_exclude_NCs)
  summary(implicit_pv) # NO IMPLICIT LEARNING F(1,12) = 1.081, p = 0.319
  
  implicit_pv_plot <- ggplot(baseline_exclude_NCs, aes(group_instruction, pv, colour = factor(block))) +
    geom_boxplot() 
  print(implicit_pv_plot) ## REPLACE THIS PLOT 
  
  ## does instruction have an effect?
  ## compare reach AEs across exclude-strategy NC and include-strategy NC b/w instruct and non-instruct groups
  exclude_include_NCs <- both_groups_AEdf %>% filter(block == 0 | block == 1)
  instruction_pv <- aov(pv ~ block + Error(participant/block), data = exclude_include_NCs)
  summary(instruction_pv) # EFFECT OF INSTRUCTION F(1,12) = 14, p = 0.00281
  
  instruction_pv_plot <- ggplot(exclude_include_NCs, aes(group_instruction, pv, colour = factor(block))) +
    geom_boxplot()
  print(instruction_pv_plot)
  
  ## within each group, does implementing a strategy affect the size of reach AEs?
  # no instruction group
  planned_comparisons1 <- both_groups_AEdf %>% filter(group_instruction == 0) %>% filter(block == 0 | block == 1)
  t.test(planned_comparisons1[which(planned_comparisons1$block==0),]$pv,planned_comparisons1[which(planned_comparisons1$block==1),]$pv,paired=TRUE)
  ggplot(planned_comparisons1, aes(group_instruction, pv, colour = factor(block))) +
    geom_boxplot()
  # large significant difference in reach AEs between exclude and include strategy without instructions
  # BUT later in the plots you see they are not going in the right directions.
  # instruction group 
  planned_comparisons2 <- both_groups_AEdf %>% filter(group_instruction == 1) %>% filter(block == 0 | block == 1)
  t.test(planned_comparisons2[which(planned_comparisons2$block==0),]$pv,planned_comparisons2[which(planned_comparisons2$block==1),]$pv,paired=TRUE)
  # large significant difference in reach AEs between exclude and include strategy with instruction
  
  ## visualize reach AEs for non-instructed (static) group
  static_tdf_NCs_rot <- static_tdf %>% filter(instruction == 'exclude' | instruction == 'include') %>% filter(trial == 0)
  ggplot(static_tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(garage_location))) +
    geom_boxplot() +
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("No instruction (static) Dual Group")
  explicit_tdf_NCs_rot <- tdf %>% filter(instruction == 'exclude' | instruction == 'include') %>% filter(trial == 0)
  ggplot(explicit_tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(garage_location))) +
    geom_boxplot() +
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("Instruction (explicit) Dual Group")
  # DUAL ADAPTATION IS EXPLICIT !!!
}

}

