######### MAKE CHANGES HERE ############
########################################
# note to self this is the master script (04 05 2019)
setwd('/Users/mayala/Desktop/conseq data')

subject_numbers <- c(1:7, 9:31) # consequence experiment
#subject_numbers <- c(3:9, 11:17, 20:33) # sequence experiment
#subject_numbers <- c(1:12) #static experiment
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
      taskmeans<- tdf %>% filter(task==taskno) %>% filter(garage_location == rotationno) %>% group_by(participant) %>%
        group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                                      SEM_pl = SD_pl/sqrt(length(unique(participant))),
                                      Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
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
  library(ggplot2)
  library(gginnards)
  
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
  t.test(block1$pv, blocklast$pv, alternative ="greater" )
  
  # mod1 <- ezANOVA(data = adaptdf,
  #                 dv = pv,
  #                 wid = participant,
  #                 within=.(grg,block),
  #                 detailed = TRUE)
  
  # VISUALIZE LEARNING
  tdfsmoothCW <- tdf %>% filter(garage_location==1) %>% filter(task == 3 | task == 5 | task == 7)
  tdfsmoothCCW <- tdf %>% filter(garage_location==-1) %>% filter(task == 3 | task == 5 | task == 7) 
  
  traindfCW <- tdfsmoothCW %>% group_by(participant) %>% group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                                                                                   SEM_pl = SD_pl/sqrt(length(unique(participant))),
                                                                                   Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
                                                                                   SEM_pv = SD_pv/sqrt(length(unique(participant))))
  traindfCCW <- tdfsmoothCCW %>% group_by(participant) %>% group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                                                                                       SEM_pl = SD_pl/sqrt(length(unique(participant))),
                                                                                       Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
                                                                                       SEM_pv = SD_pv/sqrt(length(unique(participant))))
  
  
  ## huge plot with all trial datapoints 
  # Means <- ggplot(data=traindfCW, aes(x=trial, y=Mean_pv)) +
  #   geom_line() + 
  #   geom_ribbon(data=traindfCW, aes(ymin=Mean_pv-SEM_pv, ymax= Mean_pv+SEM_pv), alpha=0.4) +
  #   geom_point(data=tdfsmoothCW, aes(x=trial, y= pv_angle), colour='chartreuse1',alpha = 0.1) +
  #   geom_point(data=tdfsmoothCCW, aes(x=trial, y= pv_angle), colour='blue',alpha = 0.1) +
  #   geom_line(data=traindfCCW, aes(x=trial, y=Mean_pv)) + 
  #   geom_ribbon(data=traindfCCW, aes(ymin=Mean_pv-SEM_pv, ymax= Mean_pv+SEM_pv), alpha=0.4) +
  #   ylim(-100, 100) +
  #   geom_smooth(data = traindfCW, model=loess) +
  #   geom_smooth(data = traindfCCW, model=loess) +
  #   ggtitle('Sequence Experiment - Dual CCW training') +
  #   ylab("Angular error (Degrees)") +
  #   xlab("Trial") +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  # print(Means)
  # move_layers(Means,"GeomRibbon", position = "top")
  # move_layers(Means,"GeomLine", position = "top")
  # move_layers(Means,"GeomPoint", position = "bottom")

  ## plotting blocked training
  for (rot in sort(unique(tdf$garage_location))){
    
    dfplot1 <- tdf %>% filter(garage_location==rot) %>% filter(task == 3)
    dfplot2 <- tdf %>% filter(garage_location==rot) %>% filter(task == 5)
    dfplot3 <- tdf %>% filter(garage_location==rot) %>% filter(task == 7)
    dfplot2$trial <- dfplot2$trial + 359
    dfplot3$trial <- dfplot3$trial + 359 + 359
    dfplot <- rbind(dfplot1, dfplot2, dfplot3)
    
    traindf <- dfplot %>% group_by(participant) %>% group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                                                                                  SEM_pl = SD_pl/sqrt(length(unique(participant))), Mean_pv = mean(pv_angle, na.rm=TRUE),                                                                                SD_pv = sd(pv_angle, na.rm=TRUE),SEM_pv = SD_pv/sqrt(length(unique(participant))))
    bl1 <- traindf[1:3,] %>% mutate(block = 1)
    bl2 <- traindf[4:6,] %>% mutate(block = 2)
    bll <- traindf[(nrow(traindf)-2):nrow(traindf),] %>% mutate(block = 7)
    combobl <- rbind(bl1, bl2, bll)
    combobl <- combobl %>% group_by(block) %>% summarise(pv = mean(Mean_pv, na.rm=TRUE), sem = mean(SEM_pv,na.rm=TRUE))
    
    bltrain <- ggplot(data=combobl, aes(x=block, y=pv)) +
      geom_point() +
      geom_line(data=combobl, aes(x=block, y=pv)) +
      geom_ribbon(data=combobl, aes(ymin=pv-sem, ymax= pv+sem), alpha=0.4) +
      ylim(-50,50) +
      xlim(1,7) +
      coord_fixed(ratio = 1/7) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      ggtitle('Consequence experiment - Blocked training')
    print(bltrain)
  }
  #######################################################################
  ###############EDIT FOR REACH AFTEREFFECTS ANALYSIS####################
  #######################################################################
  
  # get df of just reach AEs first. baseline = -1; excludeAE = 0; includeAE = 1
  sequence_baselineAE <- tdf %>% filter(task==1) %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = -1)
  sequence_excludeAE <- tdf %>% filter(instruction=='exclude') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = 0)
  sequence_includeAE <- tdf %>% filter(instruction=='include') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), block = 1)
  
  sequence_AEdf <- rbind(sequence_baselineAE, sequence_excludeAE, sequence_includeAE)
  sequence_AEdf$block <- factor(sequence_AEdf$block)
  sequence_AEdf$participant <- factor(sequence_AEdf$participant)
  sequence_AEdf <- sequence_AEdf %>% mutate(group_instruction=1) # add group label for later 
  
  t.test(sequence_excludeAE$pv, sequence_baselineAE$pv, alternative ="greater" )
  t.test(sequence_excludeAE$pv-sequence_baselineAE$pv, mu=0, alternative ="greater" )
  t.test(sequence_includeAE$pv-sequence_baselineAE$pv, sequence_excludeAE$pv-sequence_baselineAE$pv, alternative ="greater" )
  
  ##### CROSS EXPERIMENT COMPARISONS ##############################################################
  #################################################################################################
  ####### E/I COMPARISONS ACROSS INSTRUCTED AND NON-INSTRUCTED EXPERIMENTS ########################
  #################################################################################################
  
  ## did any implicit learning occur?
  ## compare reach AEs across baseline NC and exclude-strategy NC b/w instruct and non-instruct groups
  both_groups_AEdf <- rbind(static_AEdf, instructed_AEdf)
  
  baseline_exclude_NCs <- both_groups_AEdf %>% filter(block == 0 | block == -1)
  implicit_pv <- aov(pv ~ block + Error(participant/block), data=baseline_exclude_NCs)
  summary(implicit_pv) 
  # F(1,11) = 2.21 p=0.165
  
  implicit_pv_plot <- ggplot(baseline_exclude_NCs, aes(group_instruction, pv, colour = factor(block))) +
    geom_boxplot() 
  print(implicit_pv_plot) ## REPLACE THIS PLOT 
  
  ## does instruction have an effect?
  ## compare reach AEs across exclude-strategy NC and include-strategy NC b/w instruct and non-instruct groups
  exclude_include_NCs <- both_groups_AEdf %>% filter(block == 0 | block == 1)
  instruction_pv <- aov(pv ~ block + Error(participant/block), data = exclude_include_NCs)
  summary(instruction_pv) # EFFECT OF INSTRUCTION F(1,11) = 10.72, p = 0.00742
  
  instruction_pv_plot <- ggplot(exclude_include_NCs, aes(group_instruction, pv, colour = factor(block))) +
    geom_boxplot()
  print(instruction_pv_plot)
  
  ## within each group, does implementing a strategy affect the size of reach AEs?
  # no instruction group
  planned_comparisons1 <- both_groups_AEdf %>% filter(group_instruction == 0) %>% filter(block == 0 | block == 1)
  t.test(planned_comparisons1[which(planned_comparisons1$block==0),]$pv,planned_comparisons1[which(planned_comparisons1$block==1),]$pv,paired=TRUE)
  ggplot(planned_comparisons1, aes(group_instruction, pv, colour = factor(block))) +
    geom_boxplot()
  #t = -2.6557, df = 10, p-value = 0.02408

  # instruction group 
  planned_comparisons2 <- both_groups_AEdf %>% filter(group_instruction == 1) %>% filter(block == 0 | block == 1)
  t.test(planned_comparisons2[which(planned_comparisons2$block==0),]$pv,planned_comparisons2[which(planned_comparisons2$block==1),]$pv,paired=TRUE)
  #t = -3.5606, df = 10, p-value = 0.005176
  
  ########### visualize reach AEs for non-instructed (static) group
  static_tdf_NCs_rot <- tdf %>% filter(instruction == 'exclude' | instruction == 'include') %>% filter(trial == 0)
  ggplot(static_tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(garage_location))) +
    geom_boxplot() +
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("No instruction (static) Dual Group")
  
## REACH AE Visualizations
##  get SEMs for errorbar ##
SEMs <- NA
 for (garage in sort(unique(tdf$garage_location))) {
   for (instruct in sort(unique(tdf$instruction))) {
     x <- tdf %>% filter(garage_location == garage ) %>% filter(instruction == instruct) %>% filter(trial == 0) %>% group_by(participant) %>%
     group_by(trial) %>% summarise(Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
                                   SEM_pv = SD_pv/sqrt(length(unique(participant))), 
                                   instruction = instruct, garage_location = garage, lowerSEM = Mean_pv-SEM_pv, upperSEM = Mean_pv + SEM_pv)
     
     if (is.data.frame(SEMs) == TRUE ) {
       SEMs <- rbind(SEMs, x)
     } else {
       SEMs <- x
     }
     }
   }
  ## bar plot I/E Reach aftereffects ##
IEbars<- ggplot(data=SEMs, aes(x=instruction, y=Mean_pv, fill=as.factor(garage_location))) +
          geom_bar(stat="identity", position ="dodge") +
          geom_errorbar(data=SEMs, mapping=aes(x=instruction, y=Mean_pv, ymin=SEMs$lowerSEM, ymax=SEMs$upperSEM),
                        width=0.1, size=1, color="grey", position = position_dodge(width = 0.9)) +
          ylim(-50, 50) +
          ylab("Angular error (Degrees)") +
          ggtitle("No instruction (Seq) Dual Group")
print(IEbars)

}




