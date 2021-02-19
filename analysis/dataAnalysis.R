##### MASTER SCRIPT FOR ANALYSIS AND PLOTS FOR SEQ, CONSEQ, STATIC EXPERIMENTS
setwd('~/science/repos/sequence')
setwd('conseq data')
rm(list=ls())

subject_numbers <- c(1:7, 9:31) # CONSEQUENCE EXPERIMENT
#subject_numbers <- c(3:9, 11:17, 20:33) # SEQUENCE EXPERIMENT
#subject_numbers <- c(1:12) # STATIC EXPERIMENT
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # IF THIS IS 0, NO SCALING WILL BE DONE IN TASKANALYSIS()

library(dplyr)
library(ggplot2)
library(svglite)
library(gginnards)
library(Hmisc)
library(ggbeeswarm)
library(tidyr)
library(lsr)
library(OneR)
library(TOSTER)
library(pwr)

##### FUNCTION FOR PLOTTING BLOCKED RAW DATA ONLY
plotData <- function(){
  
  filename <- sprintf('md_analysis/allTaggedData_n%d_%s.csv',
                      length(subject_numbers), outfile_suffix)
  
  df <- read.csv(filename, header = TRUE)

  df <- df[-c(1),] # get rid of that random first row of NAs
  # note for static - due to experiment error p6 did double the baseline NCs 
  
  df <- df %>%
    group_by(participant) %>%
    group_by(task) %>%
    mutate(binno = bin(trial,  # make bins for plotting
                       nbins = (max(trial)+1)/3,
                       labels = c(1:((max(trial)+1)/3))))
  
  # first replace outlier-tagged trials with NA so they don't get plotted
  for (rowno in 1:nrow(df)) {
    if (df$isoutlier[rowno] == TRUE) {
      df$pathlength[rowno] <- NA
      df$pv_angle[rowno] <- NA
      df$pv_angle_n[rowno] <- NA
    }
  }
  
  tdf <- tbl_df(df) 
  
  rotations <- sort(unique(tdf$garage_location)) # use this because no-cursor trials
  # are labeled rotation = 0
  
  wholePI <- NA
  
  for (rotationno in rotations) {
    
    print(rotationno)
    
    # PLOT BLOCKED LEARNING CURVES
    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every trial -- use for plotting learning curves
      dfname <- sprintf('rotation%d_task%d_means', rotationno, taskno)
      print(dfname)
      
      taskmeans <- tdf %>%
        filter(task == taskno) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        group_by(binno) %>%
        summarise(Mean_pl = mean(pathlength, na.rm=TRUE),
                  SD_pl = sd(pathlength, na.rm=TRUE),
                  SEM_pl = SD_pl/sqrt(length(unique(participant))),
                  Mean_pv = mean(pv_angle, na.rm=TRUE),
                  SD_pv = sd(pv_angle, na.rm=TRUE),
                  SEM_pv = SD_pv/sqrt(length(unique(participant))))
      
      # plot each task's learning curve
      taskplot <- ggplot(data = taskmeans, aes(x = binno, y = Mean_pv, group = 1)) +
        geom_line() + 
        geom_ribbon(aes(ymin=Mean_pv-SEM_pv, ymax=Mean_pv+SEM_pv),
                    alpha=0.4) +
        coord_fixed(ratio = 2) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        ylim(-50, 50) +
        ylab("Angular error (Degrees)") +
        xlab("Block") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        ggtitle(dfname)
      
      #print(taskplot)
      plotfilename <- sprintf('md_analysis/Full_LCs_%s.svg', dfname)
      # ggsave(file = plotfilename,
      #       plot = taskplot,
      #       height = 10, dpi = 96, units = "cm")
       
    }
    
    # GET PI FOR EACH PARTICIPANT, STORE, AND PLOT
    if (rotationno == 1){
      
      PI_inblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == min(binno_no) | binno_no == min(binno_no)+1) %>%
        summarise(Mean_pv_inblock = mean(pv_angle, na.rm=TRUE))
      
      PI_midblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
        summarise(Mean_pv_midblock = mean(pv_angle, na.rm=TRUE))
      
      PI_finblock <- tdf %>%
        filter(task == 7) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
        summarise(Mean_pv_finblock = mean(pv_angle, na.rm=TRUE))
      
      PI_df <- PI_inblock %>%
        inner_join(PI_finblock) %>%
        inner_join(PI_midblock) %>%
        mutate(rot = rotationno,
               PI_1 = ((Mean_pv_inblock - Mean_pv_midblock)/Mean_pv_inblock)*100,
               PI_2 = ((Mean_pv_inblock - Mean_pv_finblock)/Mean_pv_inblock)*100,
               PI_r30_1 = ((30 - Mean_pv_finblock)/30)*100,
               PI_r30_2 = ((30 - Mean_pv_midblock)/30)*100)
      
    } else { 
      
      tdf$pv_angle_neg <- tdf$pv_angle*-1
      
      PI_inblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == min(binno_no) | binno_no == min(binno_no)+1) %>%
        summarise(Mean_pv_inblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_midblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
        summarise(Mean_pv_midblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_finblock <- tdf %>%
        filter(task == 7) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
        summarise(Mean_pv_finblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_df <- PI_inblock %>%
        inner_join(PI_finblock) %>%
        inner_join(PI_midblock) %>%
        mutate(rot = rotationno,
               PI_1 = ((Mean_pv_inblock - Mean_pv_midblock)/Mean_pv_inblock)*100,
               PI_2 = ((Mean_pv_inblock - Mean_pv_finblock)/Mean_pv_inblock)*100,
               PI_r30_1 = ((30 - Mean_pv_finblock)/(30))*100,
               PI_r30_2 = ((30 - Mean_pv_midblock)/(30))*100)
      
    }
    
    if (is.data.frame(wholePI) == TRUE) {
      
      wholePI <- rbind(PI_df, wholePI)
      
    } else {
      
      wholePI <- PI_df
      
    }
    
    
  }
  
# PLOT PI PLOT FOR WHOLE CONDITION
 wholePI$isoutlier <- FALSE # remove any intense outliers
 wholePI$isoutlier[which(wholePI$PI_2 %in% boxplot(wholePI$PI_2)$out)] <- TRUE
 
 #write.csv(wholePI, "md_analysis/PI_inblock.csv", row.names = FALSE)
 
 wholePI <- wholePI %>% 
   filter(isoutlier == FALSE) %>%
   group_by(rot) %>%
   mutate(PImean = mean(PI_2, na.rm = TRUE),
          PIsd = sd(PI_2, na.rm = TRUE),
          SEM = PIsd/sqrt(length(unique(participant))))
 
  PIplot <- ggplot(data = wholePI, aes(x = rot, y = PI_2, fill = rot)) +
    stat_summary(fun.y = mean, geom = "bar", na.rm = TRUE) +
    geom_errorbar(data = wholePI,
                  mapping = aes(x = rot, y = PI_2,
                                ymin =PImean - SEM , ymax = PImean + SEM),
                  width = 0.1, size = 0.5, color = "black",
                  position = position_dodge(width = 0.9)) +
    geom_beeswarm(data = wholePI, aes(x = rot, y = PI_2),
                  alpha = 1/7,
                  dodge.width = 2, cex = 3,
                  stroke = 0.3) +
    ylab("Percentage Improvement") +
    ggtitle("PI relative to initial error") +
    coord_fixed(ratio = 1/30) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_blank(),
          legend.position = "none") +
    scale_y_continuous(breaks = seq(-150, +150, 50), limits = c(-150, 150))
  
  print(PIplot)
  plotfilename <- sprintf('md_analysis/PI_v2.svg')
  ggsave(file = plotfilename,
         plot = PIplot,
         height = 10, dpi = 96, units = "cm")
  
}

##### FUNCTION FOR GETTING STATS AND CLEAN PLOTS
getStatistics <- function(){

  filename <- sprintf('md_analysis/allTaggedData_n%d_%s.csv', length(subject_numbers), outfile_suffix)
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
  
  tdf <- tbl_df(df) 
  ### NOTE: this is analyzing COLLAPSED rotations 

  ################################
  ################################
  ################################
  
  # ANALYZE LEARNING
  baselinelast <- tdf %>%
    filter(task == 0) %>%
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2,max(trial)-1,max(trial))) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), 
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
  
  block1 <- tdf %>%
    filter(task == 3) %>%
    group_by(participant) %>%
    filter(trial %in% c(min(trial), min(trial)+1, min(trial)+2)) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), 
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
  
  blocklast <- tdf %>% 
    filter(task == 7) %>%
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2,max(trial)-1,max(trial))) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
 
  adaptdf <- rbind(block1, blocklast)

  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  # FOR STATIC - save file for easy import for cross experiment comparisons
  #write.csv(block1, "md_analysis/block1_noninstructed.csv", row.names = FALSE)
  
  TOSTpaired(n = length(unique(subject_numbers)),
                         m1 = mean(baselinelast$pv), m2 = mean(blocklast$pv),
                         sd1 = sd(baselinelast$pv), sd2 = sd(blocklast$pv),
                         low_eqbound_dz = -0.2, high_eqbound_dz = 0.2,
                         r12 = cor(baselinelast$pv, blocklast$pv))
  
  
  
  # error proxy 1 - peak velocity angle
  shapiro.test(adaptdf$pv) # if significant - do wilcoxon rank test
  t.test(block1$pv, blocklast$pv,
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS
  wilcox.test(block1$pv, blocklast$pv, # if you fail normality tests
              alternative = "greater",
              paired = TRUE)
  cohensD(block1$pv, blocklast$pv,
          method = "paired")
  
  pwr.t.test(n = length(unique(tdf$participant)), 
             d = 0.5,
             #sig.level = 0.02406,
             type = "paired",
             alternative = "greater")
  
  # error proxy 2 - pathlength
  boxplot(adaptdf$pl)
  shapiro.test(adaptdf$pl) # if significant - log to get a more normal dist'n
  t.test(log(block1$pl), log(blocklast$pl),
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS
  
  # error proxy - maximum deviation angle
  boxplot(adaptdf$md)
  shapiro.test(adaptdf$md) # if significant - log to get a more normal dist'n
  t.test(log(block1$md), log(blocklast$md),
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS

  # VISUALIZE LEARNING
  tdfsmoothCW <- tdf %>%
    filter(garage_location == 1) %>%
    filter(task == 3 | task == 5 | task == 7)
  
  tdfsmoothCCW <- tdf %>%
    filter(garage_location==-1) %>%
    filter(task == 3 | task == 5 | task == 7) 
  
  traindfCW <- tdfsmoothCW %>%
    group_by(participant) %>%
    group_by(trial) %>% 
    summarise(Mean_pl = mean(pathlength, na.rm = TRUE),
              SD_pl = sd(pathlength, na.rm = TRUE),
              SEM_pl = SD_pl/sqrt(length(unique(participant))),
              Mean_pv = mean(pv_angle, na.rm = TRUE), SD_pv = sd(pv_angle, na.rm = TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))))
  
  traindfCCW <- tdfsmoothCCW %>%
    group_by(participant) %>%
    group_by(trial) %>%
    summarise(Mean_pl = mean(pathlength, na.rm = TRUE),
              SD_pl = sd(pathlength, na.rm = TRUE),
              SEM_pl = SD_pl/sqrt(length(unique(participant))),
              Mean_pv = mean(pv_angle, na.rm = TRUE),
              SD_pv = sd(pv_angle, na.rm = TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))))
  
  # PLOT TRUNCATED LCs
  for (rot in sort(unique(tdf$garage_location))){
    
    dfplot1 <- tdf %>% filter(garage_location==rot) %>% filter(task == 3)
    dfplot2 <- tdf %>% filter(garage_location==rot) %>% filter(task == 5)
    dfplot3 <- tdf %>% filter(garage_location==rot) %>% filter(task == 7)
    dfplot2$trial <- dfplot2$trial + 359
    dfplot3$trial <- dfplot3$trial + 359 + 359
    dfplot <- rbind(dfplot1, dfplot2, dfplot3)
    
   
    ppdf1 <- dfplot %>%
      filter(trial %in% c(0,1,2)) %>%
      group_by(participant) %>%
      mutate(block = 1, blockmean = mean(pv_angle, na.rm = TRUE)) %>%
      select(participant, block, blockmean) %>%
      distinct(participant, .keep_all = TRUE)
    ppdf2 <- dfplot %>%
      filter(trial %in% c(3,4,5)) %>%
      group_by(participant) %>%
      mutate(block = 2, blockmean = mean(pv_angle, na.rm = TRUE)) %>%
      select(participant, block, blockmean) %>%
      distinct(participant, .keep_all = TRUE)
    
    # comparing learning before and after participants are made aware that the cursors are weird
    # ppdf_final720 <- dfplot %>%
    #   filter(trial %in% c((359-2):359)) %>%
    #   group_by(participant) %>%
    #   mutate(block = 6, blockmean = mean(pv_angle, na.rm = TRUE)) %>%
    #   select(participant, block, blockmean) %>%
    #   distinct(participant, .keep_all = TRUE)
    # ppdf7_24 <- dfplot %>%
    #   group_by(participant) %>%
    #   mutate(block = 8, blockmean = mean(pv_angle, na.rm = TRUE)) %>%
    #   select(participant, block, blockmean) %>%
    #   distinct(participant, .keep_all = TRUE)
    # ytemp<- t.test(ppdf_final720$blockmean, ppdf7_24$blockmean, paired = TRUE)
    # print(ytemp)
    
    ppdf7 <- dfplot %>%
      group_by(participant) %>%
      filter(trial %in% c(max(trial), max(trial) - 1, max(trial) - 2)) %>%
      mutate(block = 7, blockmean = mean(pv_angle, na.rm = TRUE)) %>%
      select(participant, block, blockmean) %>%
      distinct(participant, .keep_all = TRUE)
    ppdf_full <- rbind(ppdf1, ppdf2, ppdf7)
    
    traindf <- ppdf_full %>%
      group_by(block) %>%
      mutate(pv = mean(blockmean, na.rm = TRUE),
             sdpv = sd(blockmean, na.rm = TRUE),
             sem = sdpv/sqrt(length(unique(participant)))) 
    
    #outfile_name = sprintf('DUAL_LCs_%s.csv', rot)
    #write.csv(traindf, file = outfile_name, row.names = FALSE)  
    
    # block training plots
    bltrain <- ggplot(data=traindf, aes(x=block, y=pv)) +
      geom_point() +
      geom_line(data=traindf, aes(x=block, y=pv)) +
      geom_line(data=ppdf_full, aes(x=block, y=blockmean, colour=as.factor(participant)), alpha = 0.3) + 
      geom_ribbon(data=traindf, aes(ymin=pv-sem, ymax= pv+sem), alpha=0.4) +
      coord_fixed(ratio = 1/7) +
      ylim(-50,50) +
      scale_x_continuous(labels = c(0,1,2,3,4,5,6,7), breaks = c(0,1,2,3,4,5,6,7))+
      #xlim(1,7) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_blank(), legend.position = "none") +
      ggtitle('Dual Static') 
    
    print(bltrain)
    move_layers(bltrain, "GeomRibbon", position = "top")
    move_layers(bltrain, "GeomPoint", position = "top")
    
    out_name <- sprintf('DUAL_LCs_%s.csv', rot)
    write.csv(traindf, out_name, row.names = FALSE)
  }

  ################################
  ################################
  ################################
  
  # ANALYZE REACH AFTEREFFECTS
  
  # COLLAPSED ROTATIONS FIRST
  sequence_baselineAE <- tdf %>%
    filter(task == 1) %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              instruction = 'baseline')
  
  sequence_excludeAE <- tdf %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              instruction = 'exclude')
  
  sequence_includeAE <- tdf %>%
    filter(instruction =='include') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              instruction = 'include')
  
  # AE ANALYSIS BUT BASELINE SEPARATE
  t.test(sequence_excludeAE$pv, sequence_baselineAE$pv,
         alternative = "greater",
         paired = TRUE) 
  
  t.test(sequence_includeAE$pv - sequence_baselineAE$pv,
         sequence_excludeAE$pv - sequence_baselineAE$pv,
         alternative = "greater",
         paired = TRUE )
  
  # SUBTRACT BASELINE FOR FIGURES & MORE AE ANALYSIS
  sequence_excludeAE$pv <- sequence_excludeAE$pv - sequence_baselineAE$pv
  sequence_includeAE$pv <- sequence_includeAE$pv - sequence_baselineAE$pv
  
  # AE ANALYSIS, BASELINE SUBTRACTED
  t.test(sequence_excludeAE$pv,
         mu=0,
         alternative ="greater") # IMPLICIT LEARNING MEASURE
  
  ################################
  ################################
  ################################
  
  # REACH AE VISUALIZATIONS - BASELINE DEDUCTED - SIGNED ERRORS
  
  sequence_baselineAE_CCW <- tdf %>%
    filter(garage_location == -1) %>%
    filter(task == 1) %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>% # BASELINE HAS BOTH GARAGES & NEED 3 TRIALS FOR BASELINE BLOCK
    summarise(pvB = mean(pv_angle, na.rm = TRUE))
  
  sequence_baselineAE_CW <- tdf %>%
    filter(garage_location == 1) %>%
    filter(task == 1) %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
    summarise(pvB = mean(pv_angle, na.rm = TRUE))
  
  sequence_excludeAE_CCW <- tdf %>%
    filter(garage_location == -1) %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              garage_location = unique(garage_location),
              instruction = 'exclude')
  sequence_excludeAE_CCW$pv <- sequence_excludeAE_CCW$pv - sequence_baselineAE_CCW$pv
  
  sequence_excludeAE_CW <- tdf %>%
    filter(garage_location == 1) %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              garage_location = unique(garage_location),
              instruction = 'exclude')
  sequence_excludeAE_CW$pv <- sequence_excludeAE_CW$pv - sequence_baselineAE_CW$pv
  
  sequence_includeAE_CCW <- tdf %>%
    filter(garage_location == -1) %>%
    filter(instruction == 'include') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              garage_location = unique(garage_location),
              instruction = 'include')
  sequence_includeAE_CCW$pv <- sequence_includeAE_CCW$pv - sequence_baselineAE_CCW$pv

  sequence_includeAE_CW <- tdf %>%
    filter(garage_location == 1) %>%
    filter(instruction == 'include') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm=TRUE),
              garage_location = unique(garage_location),
              instruction = 'include')
  sequence_includeAE_CW$pv <- sequence_includeAE_CW$pv - sequence_baselineAE_CW$pv
  
  swarms <- rbind(sequence_excludeAE_CCW,
                  sequence_excludeAE_CW,
                  sequence_includeAE_CCW,
                  sequence_includeAE_CW)
  
  ## COLLECT BASELINE-SUBTRACTED REACH AEs
  
  sequence_excludeAE_CCW.summary <- sequence_excludeAE_CCW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "exclude",
              garage_location = "-1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)

  sequence_excludeAE_CW.summary <- sequence_excludeAE_CW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "exclude",
              garage_location = "1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  sequence_includeAE_CCW.summary <- sequence_includeAE_CCW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "include",
              garage_location = "-1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  sequence_includeAE_CW.summary <- sequence_includeAE_CW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "include",
              garage_location = "1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  SEMs <- rbind(sequence_excludeAE_CCW.summary, 
                sequence_excludeAE_CW.summary,
                sequence_includeAE_CCW.summary,
                sequence_includeAE_CW.summary)

  ## bar plot I/E Reach aftereffects ##
IEbars <- ggplot(data = SEMs,
                aes(x = instruction, y = Mean_pv, fill = as.factor(garage_location))) +
            geom_bar(stat = "identity", position = "dodge") +
            geom_errorbar(data = SEMs, mapping = aes(x = instruction, y = Mean_pv, 
                                                     ymin = SEMs$lowerSEM, ymax = SEMs$upperSEM),
                                                     width = 0.2, size = 0.5, color = "black",
                                                     position = position_dodge(width = 0.9)) +
            geom_beeswarm(data = swarms, aes(x = instruction, y = pv),
                alpha = 1/7,
                dodge.width = .9, cex = 3,
                stroke = 0.3) +
            # geom_point(data = sequence_excludeAE_CCW, size = 1, stroke = 0,
            #            aes(x = instruction, y = pv), alpha = 1/20,
            #            position = position_dodge(width = 0.5, preserve = "single")) +
            # geom_point(data = sequence_excludeAE_CW, size = 1, stroke = 0,
            #            aes(x = instruction, y = pv), alpha = 1/20,
            #            position = position_dodge(width = -0.5, preserve = "single")) +
            # geom_point(data = sequence_includeAE_CCW, size = 1, stroke = 0,
            #            aes(x = instruction, y = pv), alpha = 1/20,
            #            position = position_dodge(width = 0.5, preserve = "single")) +
            # geom_point(data = sequence_includeAE_CW, size = 1, stroke = 0,
            #            aes(x = instruction, y = pv), alpha = 1/20,
            #            position = position_dodge(width = -0.5, preserve = "single")) +
            ylab("Angular error (Degrees)") +
            ggtitle("Dual Sequence") +
            coord_fixed(ratio = 1/13) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  legend.title = element_blank(), legend.position = "none") +
            scale_y_continuous(breaks = seq(-30, +30, 10), limits = c(-30, 30))
#move_layers(IEbars, "GeomPoint", position = "bottom")

print(IEbars)

## ANALYZE AE BASELINE DEDUCTED:
# GETS SIGNED ERRORS FIRST,
# DEDUCTS BASELINE BASED ON GARAGE/PREHOME,
# FLIPS THE SIGN ON THE NEGATIVE ROTATION.

# EXCLUDE-STRATEGY REACH AEs - PV 
YY <- tdf %>%
  filter(instruction == 'exclude') %>%
  drop_na(pv_angle_n) %>% 
  group_by(participant) %>%
  filter(task == min(task)) %>%
  filter(trial == min(trial)) %>%
  select(participant, garage_location, instruction, pv_angle, pv_angle_n)


YY <- YY %>% distinct(participant, .keep_all = T)
YY <- inner_join(YY, sequence_baselineAE_CW, by = "participant")
YY <- inner_join(YY, sequence_baselineAE_CCW, by = "participant")

YY$pv_angle_nR <- NA
for (rowno in 1:nrow(YY)) {
  if (YY$garage_location[rowno] == -1){
    YY$pv_angle_nR[rowno] <- (YY$pv_angle[rowno] - YY$pvB.x[rowno])
  } else {
    YY$pv_angle_nR[rowno] <- (YY$pv_angle[rowno] + YY$pvB.y[rowno])*-1
    
  }
}


t.test(YY$pv_angle_nR,
       mu = 0,
       alternative = "greater") # IMPLICIT LEARNING MEASURE
cohensD(x = YY$pv_angle_nR)

# NOW DO SAME BUT DO ENDPOINT ANGLE - first, get baseline endpoint angle
sequence_baselineAE_CCW <- tdf %>%
  filter(garage_location == -1) %>%
  filter(task == 1) %>%
  drop_na(endpointang) %>% 
  group_by(participant) %>%
  filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>% # BASELINE HAS BOTH GARAGES & NEED 3 TRIALS FOR BASELINE BLOCK
  summarise(eaB = mean(endpointang, na.rm = TRUE))

sequence_baselineAE_CW <- tdf %>%
  filter(garage_location == 1) %>%
  filter(task == 1) %>%
  drop_na(endpointang) %>% 
  group_by(participant) %>%
  filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
  summarise(eaB = mean(endpointang, na.rm = TRUE))

YY <- tdf %>%
  filter(instruction == 'exclude') %>%
  drop_na(endpointang) %>%
  group_by(participant) %>%
  filter(task == min(task)) %>%
  filter(trial == min(trial)) %>%
  select(participant, garage_location, instruction, endpointang) %>%
  distinct(participant, .keep_all = T)

YY <- inner_join(YY, sequence_baselineAE_CW, by = "participant")
YY <- inner_join(YY, sequence_baselineAE_CCW, by = "participant")

YY$endpointang_nR <- NA
for (rowno in 1:nrow(YY)) {
  if (YY$garage_location[rowno] == -1){
    YY$endpointang_nR[rowno] <- (YY$endpointang[rowno] - YY$eaB.x[rowno])
  } else {
    YY$endpointang_nR[rowno] <- (YY$endpointang[rowno] + YY$eaB.y[rowno])*-1
    
  }
}

t.test(YY$endpointang_nR,
       mu = 0,
       alternative = "greater") # IMPLICIT LEARNING MEASURE
cohensD(x = YY$endpointang_nR)

# INCLUDE-STRATEGY REACH AEs
ZZ <- tdf %>%
  filter(instruction == 'include') %>%
  drop_na(pv_angle_n) %>% 
  group_by(participant) %>%
  filter(trial == min(trial)) %>%
  filter(task == min(task)) %>%
  select(participant, garage_location, instruction, pv_angle, pv_angle_n)

ZZ <- ZZ %>% distinct(participant, .keep_all = T)

ZZ <- inner_join(ZZ, sequence_baselineAE_CW, by = "participant")
ZZ <- inner_join(ZZ, sequence_baselineAE_CCW, by = "participant")

ZZ$pv_angle_nR <- NA
for (rowno in 1:nrow(ZZ)) {
  if (ZZ$garage_location[rowno] == -1){
    ZZ$pv_angle_nR[rowno] <- (ZZ$pv_angle[rowno] - ZZ$pvB.x[rowno])
  } else {
    ZZ$pv_angle_nR[rowno] <- (ZZ$pv_angle[rowno] + ZZ$pvB.y[rowno])*-1
    
  }
}

t.test(ZZ$pv_angle_nR,
       YY$pv_angle_nR,
       paired = TRUE,
       alternative = "greater") # EXPLICIT LEARNING MEASURE
cohensD(x = ZZ$pv_angle_nR,
        y = YY$pv_angle_nR,
        method = "paired")

# now do the same but with endpoint angle
# INCLUDE-STRATEGY REACH AEs
ZZ <- tdf %>%
  filter(instruction == 'include') %>%
  drop_na(endpointang) %>% 
  group_by(participant) %>%
  filter(trial == min(trial)) %>%
  filter(task == min(task)) %>%
  select(participant, garage_location, instruction, endpointang)

ZZ <- ZZ %>% distinct(participant, .keep_all = T)

ZZ <- inner_join(ZZ, sequence_baselineAE_CW, by = "participant")
ZZ <- inner_join(ZZ, sequence_baselineAE_CCW, by = "participant")

ZZ$endpointang_nR <- NA
for (rowno in 1:nrow(ZZ)) {
  if (ZZ$garage_location[rowno] == -1){
    ZZ$endpointang_nR[rowno] <- (ZZ$endpointang[rowno] - ZZ$eaB.x[rowno])
  } else {
    ZZ$endpointang_nR[rowno] <- (ZZ$endpointang[rowno] + ZZ$eaB.y[rowno])*-1
    
  }
}

t.test(ZZ$endpointang_nR,
       YY$endpointang_nR,
       paired = TRUE,
       alternative = "greater") # EXPLICIT LEARNING MEASURE

################################
################################
################################


}

##### CROSS EXPERIMENT COMPARISONS
crossComparisons <- function(){
################################
################################
################################

# CROSS EXPERIMENT COMPARISONS 

instructed_AEdf_pv <- read.csv("md_analysis/explicit_AEdf_pv.csv", header = TRUE)
static_AEdf_pv <- read.csv("md_analysis/static_AEdf_pv.csv", header = TRUE)
static_AEdf_pv$participant <- static_AEdf_pv$participant + 20

instructed_AEdf_ea <- read.csv("md_analysis/explicit_AEdf_ea.csv", header = TRUE)
static_AEdf_ea <- read.csv("md_analysis/static_AEdf_ea.csv", header = TRUE)
static_AEdf_ea$participant <- static_AEdf_ea$participant + 20

# E/I COMPARISONS ACROSS INSTRUCTED AND NON-INSTRUCTED EXPERIMENTS 

## did any implicit learning occur?
## compare reach AEs across baseline NC and exclude-strategy NC b/w instruct and non-instruct groups
both_groups_AEdf <- rbind(static_AEdf_pv, instructed_AEdf_pv) # CHANGE PV TO EA 
#both_groups_AEdf <- rbind(static_AEdf_ea, instructed_AEdf_ea) 

#baseline_exclude_NCs <- both_groups_AEdf %>% filter(block == 0 | block == -1)
baseline_exclude_NCs <- both_groups_AEdf %>%
  filter(instruction == "baseline" | instruction == "exclude")
#implicit_pv <- aov(pv ~ block + Error(participant/block), data=baseline_exclude_NCs)
mod1 <- ezANOVA(data = baseline_exclude_NCs,
                dv = pv, # pv angle OR ea
                wid = participant,
                within = instruction, 
                between = .(group), 
                detailed = TRUE,
                return_aov = TRUE)
print(mod1)

implicit_pv_plot <- ggplot(baseline_exclude_NCs, aes(group,ea, colour = factor(instruction))) +
  geom_boxplot() 
print(implicit_pv_plot) ## REPLACE THIS PLOT 

## does instruction have an effect?
## compare reach AEs across exclude-strategy NC and include-strategy NC b/w instruct and non-instruct groups
#exclude_include_NCs <- both_groups_AEdf %>% filter(block == 0 | block == 1)
exclude_include_NCs <- both_groups_AEdf %>%
  filter(instruction == "exclude" | instruction == "include")
#instruction_pv <- aov(pv ~ block + Error(participant/block), data = exclude_include_NCs)
#summary(instruction_pv) # EFFECT OF INSTRUCTION F(1,11) = 10.72, p = 0.00742
mod2 <- ezANOVA(data = exclude_include_NCs,
                dv = pv, # pv angle
                wid = participant,
                within = instruction, # Dual-CW, Dual-CCW, Single-CW, Single-CCW
                between = .(group), # Dual vs. Single Conditions
                detailed = TRUE,
                return_aov = TRUE)
print(mod2)

pwr.anova.test(k = 2, n= 12, f = 0.06, sig.level =  3.625135e-04)
pwr.anova.test(k = 2,  f = 0.06, sig.level =  0.05, power = 0.8)

instruction_pv_plot <- ggplot(exclude_include_NCs,
                              aes(group_instruction, pv, colour = factor(block))) +
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

}

##### PASSIVE VS. ACTIVE COMPARISONS
passiveVactive <- function() {
  
  setwd('~/science/repos/sequence')
  setwd('preq data')
  PIpreq <- read.csv("md_analysis/PI_inblock.csv", header = TRUE)
  
  setwd('~/science/repos/sequence')
  setwd('seq data')
  PIseq <- read.csv("md_analysis/PI_inblock.csv", header = TRUE)
  
  setwd('~/science/repos/sequence')
  setwd('conseq data')
  PIconseq <- read.csv("md_analysis/PI_inblock.csv", header = TRUE)
  
  PIconseq <- PIconseq %>%
    filter(isoutlier == FALSE) %>%
    group_by(participant) %>% 
    summarise(meanPI = mean(PI_2)) %>%
    mutate(group = "conseq", active = "no")
  
  PIpreq <- PIpreq %>%
    filter(isoutlier == FALSE) %>%
    group_by(participant) %>% 
    summarise(meanPI = mean(PI_2)) %>%
    mutate(group = "preq", active = "yes") 
  
  PIseq <- PIseq %>%
    filter(isoutlier == FALSE) %>%
    group_by(participant) %>% 
    summarise(meanPI = mean(PI_2)) %>%
    mutate(group = "seq", active = "yes") 
  
  PIdf <- rbind(PIconseq, PIpreq, PIseq)
  
  shapiro.test(PIdf$meanPI) # significant! do nonparam. test
  
  wilcox.test(PIdf$meanPI[which(PIdf$active == "yes")], PIdf$meanPI[which(PIdf$active == "no")], # if you fail normality tests
              alternative = "greater",
              paired = F)
  
}