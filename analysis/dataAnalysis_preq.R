##### FOR DUAL PREQUENCE DATA ANALYSIS
# NOTE: REVERSED LABEL TAGS 
setwd('~/science/repos/sequence')
setwd('preq data')
rm(list=ls())

subject_numbers <- c(1:31)  # PARTICIPANT 1 EXCLUDED DUE TO HUGE BASELINE BIASES
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(0) # IF THIS IS 0, NO SCALING WILL BE DONE IN TASKANALYSIS()

library(dplyr)
library(ggplot2)
library(svglite)
library(gginnards)
library(Hmisc)
library(ggbeeswarm)
library(tidyr)

##### FUNCTION FOR PLOTTING BINNED RAW DATA & PI
plotData <- function(){
  
  filename <- sprintf('md_analysis/allTaggedData_n%d_%s.csv',
                      length(subject_numbers), outfile_suffix)
  
  df <- read.csv(filename, header = TRUE)
  
  df <- df[-c(1),] # get rid of that random first row of NAs
  
  df <- df %>%
    filter(participant != 4) %>% # outlier participants excluded 
    filter(participant != 6) %>%
    filter(participant != 24) %>%
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
  
  rotations <- sort(unique(tdf$prehome)) # use this because no-cursor trials 
  # are labeled rotation = 0
  
  wholePI <- NA
  
  for (rotationno in rotations) {
    
    print(rotationno)
    
    # PLOT BLOCKED LEARNING CURVES
    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every bin -- use for plotting LCs
      dfname <- sprintf('rotation%d_task%d_means', rotationno, taskno)
      print(dfname)
      
      taskmeans <- tdf %>% 
        filter(task == taskno) %>%
        filter(prehome == rotationno) %>% 
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
        geom_ribbon(aes(ymin = Mean_pv-SEM_pv, ymax = Mean_pv+SEM_pv),
                    alpha = 0.4) +
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
      
       print(taskplot)
      plotfilename <- sprintf('md_analysis/Full_LCs_%s.svg', dfname)
      ggsave(file = plotfilename,
             plot = taskplot,
             height = 10, dpi = 96, units = "cm")
      
    }
    
  # GET PI FOR EACH PARTICIPANT, STORE, AND PLOT
  if (rotationno == 1){ # flipped for preq
    
    
    tdf$pv_angle_neg <- tdf$pv_angle*-1
    
    PI_inblock <- tdf %>%
      filter(task == 3) %>%
      filter(prehome == rotationno) %>%
      group_by(participant) %>%
      mutate(binno_no = as.double(binno)) %>%
      filter(binno_no == min(binno_no) | binno_no == min(binno_no)+1) %>%
      summarise(Mean_pv_inblock = mean(pv_angle_neg, na.rm=TRUE))
    
    PI_midblock <- tdf %>%
      filter(task == 3) %>%
      filter(prehome == rotationno) %>%
      group_by(participant) %>%
      mutate(binno_no = as.double(binno)) %>%
      filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
      summarise(Mean_pv_midblock = mean(pv_angle_neg, na.rm=TRUE))
    
    PI_finblock <- tdf %>%
      filter(task == 7) %>%
      filter(prehome == rotationno) %>%
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
  } else { 
    
    PI_inblock <- tdf %>%
      filter(task == 3) %>%
      filter(prehome == rotationno) %>%
      group_by(participant) %>%
      mutate(binno_no = as.double(binno)) %>%
      filter(binno_no == min(binno_no) | binno_no == min(binno_no)+1) %>%
      summarise(Mean_pv_inblock = mean(pv_angle, na.rm=TRUE))
    
    PI_midblock <- tdf %>%
      filter(task == 3) %>%
      filter(prehome == rotationno) %>%
      group_by(participant) %>%
      mutate(binno_no = as.double(binno)) %>%
      filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
      summarise(Mean_pv_midblock = mean(pv_angle, na.rm=TRUE))
    
    PI_finblock <- tdf %>%
      filter(task == 7) %>%
      filter(prehome == rotationno) %>%
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
  
  filename <- sprintf('md_analysis/allTaggedData_n%d_%s.csv',
                      length(subject_numbers), outfile_suffix)
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
  tdf <- tdf %>% filter(participant != 1)  # PARTICIPANT 1 EXCLUDED DUE TO HUGE BASELINE BIASES
  
  # NOW ANALYZING COLLAPSED ROTATIONS
  
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
  
  # IS DUAL ADAPTATION EQUIVALENT TO BASELINE?
  TOSTpaired(n = length(unique(subject_numbers)),
             m1 = mean(baselinelast$pv), m2 = mean(blocklast$pv),
             sd1 = sd(baselinelast$pv), sd2 = sd(blocklast$pv),
             low_eqbound_dz = -0.2, high_eqbound_dz = 0.2,
             r12 = cor(baselinelast$pv, blocklast$pv))
  
  adaptdf<- rbind(block1, blocklast)
  
  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)

  boxplot(adaptdf$pv) # remove the extreme outliers
  
  adaptdf2 <- adaptdf %>%
    filter(participant != "4") %>%
    filter(participant != "6") %>%
    filter(participant != "24")    # adaptdf[-c(35,37,55),] # outliers
  
  t.test(adaptdf2$pv[which(adaptdf2$block == "3")],
         adaptdf2$pv[which(adaptdf2$block == "7")],
         alternative ="greater",
         paired = TRUE )
  cohensD(x = adaptdf2$pv[which(adaptdf2$block == "3")],
          y= adaptdf2$pv[which(adaptdf2$block == "7")],
          method = "paired")
  
  pwr.t.test(n = length(unique(tdf$participant)), 
             d = 0.5,
             sig.level = 0.0009797,
             type = "paired",
             alternative = "greater") # report this
  
  pwr.t.test(d = 0.5,
             sig.level = 0.0009797,
             type = "paired",
             alternative = "greater",
             power = 0.8)
  
  # error proxy 2 - pathlength
  boxplot(adaptdf2$pl) # couple of outliers but not too intense
  shapiro.test(adaptdf2$pl) # significant - log to get a more normal dist'n
  t.test(log(adaptdf2$pl[which(adaptdf2$block == "3")]),
         log(adaptdf2$pl[which(adaptdf2$block == "7")]),
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS
  cohensD(log(adaptdf2$pl[which(adaptdf2$block == "3")]),
          log(adaptdf2$pl[which(adaptdf2$block == "7")]),
          method = "paired")
  
  # error proxy - maximum deviation angle
  boxplot(adaptdf$md) # some very intense outliers, remove them
  
  block1 <- tdf %>%
    filter(task == 3) %>%
    filter(isoutlier_md == FALSE) %>%
    filter(selection_1 == 1) %>%
    group_by(participant) %>%
    filter(trial %in% c(min(trial), min(trial)+1, min(trial)+2)) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), 
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
  
  blocklast <- tdf %>%
    filter(task == 7) %>%
    filter(isoutlier_md == FALSE) %>%
    filter(selection_1 == 1) %>%
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2,max(trial)-1,max(trial))) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
  
  adaptdfmd<- rbind(block1, blocklast)
  
  shapiro.test(adaptdfmd$md) # significant - log to get a more normal dist'n

  t.test(log(adaptdfmd$md[which(adaptdfmd$block == "3")]),
         log(adaptdfmd$md[which(adaptdfmd$block == "7")]),
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS
  cohensD(log(adaptdf$md[which(adaptdf$block == "3")]),
          log(adaptdf$md[which(adaptdf$block == "7")]),
          method = "paired")
  
  #### visualize blocked learning
  for (rot in sort(unique(tdf$prehome))){
    
    dfplot1 <- tdf %>%
      filter(prehome==rot) %>%
      filter(task == 3)  %>%
      filter(participant != "4") %>%
      filter(participant != "6") %>%
      filter(participant != "24")
    
    dfplot2 <- tdf %>%
      filter(prehome==rot) %>%
      filter(task == 5)  %>%
      filter(participant != "4") %>%
      filter(participant != "6") %>%
      filter(participant != "24")
    
    dfplot3 <- tdf %>%
      filter(prehome==rot) %>%
      filter(task == 7) %>%
      filter(participant != "4") %>%
      filter(participant != "6") %>%
      filter(participant != "24")
    
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
    
    outfile_name = sprintf('DUAL_LCs_%s.csv', rot)
    #write.csv(traindf, file = outfile_name, row.names = FALSE)  
    
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
      ggtitle('Dual Lead-in') 
    
    print(bltrain)
    move_layers(bltrain, "GeomRibbon", position = "top")
    move_layers(bltrain, "GeomPoint", position = "top")
  }
  

  ################################
  ################################
  ################################
  
  # ANALYZE REACH AFTEREFFECTS
  
  # COLLAPSED ROTATIONS FIRST
  sequence_baselineAE <- tdf %>%
    filter(task == 1) %>%
    drop_na(pv_angle_n) %>% 
    filter(trial %in% c(21,22,23)) %>%
    group_by(participant) %>%
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

  ## visualize reach AEs
  tdf_NCs_rot <- tdf %>% filter(instruction == 'exclude' | instruction == 'include') %>% filter(trial == 0) 
  ggplot(tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(prehome))) +
    geom_boxplot() +
    geom_point()+
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("PreQ Dual Group")  # note : participant 6 has one outlier and the remaining also an outlier in boxplots
  
  ## AEs LOOK WEIRD - DOUBLE CHECK THAT THEY'RE GOING THE EXPECTED DIRECTIONS
  
  tdf_NCs_baseline.summary1 <- tdf %>%
    filter(prehome=="1") %>%
    filter(task=="1") %>%
    filter(trial %in% c(21,22,23)) %>%
    group_by(participant) %>%
    summarise(meanaebase = mean(pv_angle, na.rm=T))
  
  tdf_NCs_baseline.summary2 <- tdf %>%
    filter(prehome=="-1") %>%
    filter(task=="1") %>%
    filter(trial %in% c(21,22,23)) %>%
    group_by(participant) %>%
    summarise(meanaebase = mean(pv_angle, na.rm=T))
  
  tdf_NCs_rot.summary1 <- tdf_NCs_rot %>%   # NOTE - ONLY FIRST TRIAL
    filter(prehome=="1") %>%
    filter(instruction == "include") %>%
    group_by(participant) %>%
    summarise(meanae = mean(pv_angle, na.rm=T))
  
  tdf_NCs_rot.summary2 <- tdf_NCs_rot %>%
    filter(prehome=="1") %>%
    filter(instruction == "exclude") %>%
    group_by(participant) %>%
    summarise(meanae = mean(pv_angle, na.rm=T))
  
  tdf_NCs_rot.summary3 <- tdf_NCs_rot %>%
    filter(prehome=="-1") %>%
    filter(instruction == "include") %>%
    group_by(participant) %>%
    summarise(meanae = mean(pv_angle, na.rm=T))
  
  tdf_NCs_rot.summary4 <- tdf_NCs_rot %>%
    filter(participant != 14) %>%
    filter(prehome=="-1") %>%
    filter(instruction == "exclude") %>%
    group_by(participant) %>%
    summarise(meanae = mean(pv_angle, na.rm=T))
  # double check that AEs are indeed not significant from baseline 
  #yy <- left_join(tdf_NCs_rot.summary4,tdf_NCs_baseline.summary2 )
  #t.test(yy$meanae,yy$meanaebase,paired = TRUE,alternative = "greater")
  

  # }   # note that due to the counterbalance mistake for participants 1-9, group by participants (so they only get one score in the mean)

  sequence_baselineAE_CCW <- tdf %>%
    filter(prehome == -1) %>%
    filter(task == 1) %>%
    drop_na(pv_angle_n) %>% 
    filter(trial %in% c(18:23)) %>% # BASELINE HAS BOTH GARAGES & NEED 3 TRIALS FOR BASELINE BLOCK
    group_by(participant) %>%
    summarise(pvB = mean(pv_angle, na.rm = TRUE))
  
  sequence_baselineAE_CW <- tdf %>%
    filter(prehome == 1) %>%
    filter(task == 1) %>%
    drop_na(pv_angle_n) %>% 
    filter(trial %in% c(18:23)) %>%
    group_by(participant) %>%
    summarise(pvB = mean(pv_angle, na.rm = TRUE))
  
  sequence_excludeAE_CCW <- tdf %>% # NOTE HERE A COUPLE PARTICIPANTS DO NOT HAVE ALL CONDITIONS
    filter(prehome == -1) %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              prehome = unique(prehome),
              instruction = 'exclude')
  sequence_excludeAE_CCW <- inner_join(sequence_baselineAE_CCW,
                                       sequence_excludeAE_CCW,
                                       by = "participant")
  sequence_excludeAE_CCW$pv <- sequence_excludeAE_CCW$pv - sequence_excludeAE_CCW$pvB
  
  sequence_excludeAE_CW <- tdf %>%
    filter(prehome == 1) %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              prehome = unique(prehome),
              instruction = 'exclude')
  sequence_excludeAE_CW <- inner_join(sequence_baselineAE_CW,
                                       sequence_excludeAE_CW,
                                       by = "participant")
  sequence_excludeAE_CW$pv <- sequence_excludeAE_CW$pv - sequence_excludeAE_CW$pvB
  
  sequence_includeAE_CCW <- tdf %>%
    filter(prehome == -1) %>%
    filter(instruction == 'include') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              prehome = unique(prehome),
              instruction = 'include')
  sequence_includeAE_CCW <- inner_join(sequence_baselineAE_CCW,
                                      sequence_includeAE_CCW,
                                      by = "participant")
  sequence_includeAE_CCW$pv <- sequence_includeAE_CCW$pv - sequence_includeAE_CCW$pvB
  
  sequence_includeAE_CW <- tdf %>%
    filter(prehome == 1) %>%
    filter(instruction == 'include') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm=TRUE),
              prehome = unique(prehome),
              instruction = 'include')
  sequence_includeAE_CW <- inner_join(sequence_baselineAE_CW,
                                       sequence_includeAE_CW,
                                       by = "participant")
  sequence_includeAE_CW$pv <- sequence_includeAE_CW$pv - sequence_includeAE_CW$pvB
  
  swarms <- rbind(sequence_excludeAE_CCW,
                  sequence_excludeAE_CW,
                  sequence_includeAE_CCW,
                  sequence_includeAE_CW)
    
  ## ANALYSE REACH AE FOR SEPARATED ROTATIONS
  t.test(sequence_excludeAE_CW$pv,
         mu = 0,
         alternative = "greater")
  
  t.test(sequence_excludeAE_CCW$pv,
         mu = 0,
         alternative = "less")
  
  ## COLLECT BASELINE-SUBTRACTED REACH AEs
  sequence_excludeAE_CCW.summary <- sequence_excludeAE_CCW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "exclude",
              prehome = "-1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  sequence_excludeAE_CW.summary <- sequence_excludeAE_CW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "exclude",
              prehome = "1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  sequence_includeAE_CCW.summary <- sequence_includeAE_CCW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "include",
              prehome = "-1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  sequence_includeAE_CW.summary <- sequence_includeAE_CW %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "include",
              prehome = "1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  SEMs <- rbind(sequence_excludeAE_CCW.summary, 
                sequence_excludeAE_CW.summary,
                sequence_includeAE_CCW.summary,
                sequence_includeAE_CW.summary)

  
  ## bar plot I/E Reach aftereffects ##
  IEbars <- ggplot(data = SEMs,
                   aes(x = instruction, y = Mean_pv, fill = as.factor(prehome))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(data = SEMs, 
                  mapping = aes(x = instruction,
                                y = Mean_pv,
                                ymin = SEMs$lowerSEM,
                                ymax = SEMs$upperSEM),
                  width = 0.2,
                  size = 0.5,
                  color = "black",
                  position = position_dodge(width = 0.9)) +
    geom_beeswarm(data = swarms, aes(x = instruction, y = pv),
                  alpha = 1/7,
                  dodge.width = .9, cex = 3,
                  stroke = 0.3) +
  geom_point(data = sequence_excludeAE_CCW, size = 1, stroke = 0,
             aes(x = instruction, y = pv), alpha = 1/20,
             position = position_dodge(width = 0.5, preserve = "single")) +
  geom_point(data = sequence_excludeAE_CW, size = 1, stroke = 0,
             aes(x = instruction, y = pv), alpha = 1/20,
             position = position_dodge(width = -0.5, preserve = "single")) +
  geom_point(data = sequence_includeAE_CCW, size = 1, stroke = 0,
             aes(x = instruction, y = pv), alpha = 1/20,
             position = position_dodge(width = 0.5, preserve = "single")) +
  geom_point(data = sequence_includeAE_CW, size = 1, stroke = 0,
             aes(x = instruction, y = pv), alpha = 1/20,
             position = position_dodge(width = -0.5, preserve = "single")) +
  ylab("Angular error (Degrees)") +
  ggtitle("Dual Prequence") +
  coord_fixed(ratio = 1/13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.position = "none") +
  scale_y_continuous(breaks = seq(-30, +30, 10), limits = c(-30, 30))
  move_layers(IEbars, "GeomPoint", position = "bottom")

  print(IEbars)
  
## ANALYZE AE BASELINE DEDUCTED:
  # GETS SIGNED ERRORS FIRST,
  # DEDUCTS BASELINE BASED ON PREHOME,
  # FLIPS THE SIGN ON THE NEGATIVE ROTATION.
  
  YY <- tdf %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    filter(task == min(task)) %>%
    select(participant, prehome, instruction, pv_angle, pv_angle_n)
  
  YY <- inner_join(YY, sequence_baselineAE_CW, by = "participant")
  YY <- inner_join(YY, sequence_baselineAE_CCW, by = "participant")
  
  YY$pv_angle_nR <- NA
  for (rowno in 1:nrow(YY)) {
    if (YY$prehome[rowno] == -1){
      YY$pv_angle_nR[rowno] <- (YY$pv_angle[rowno] - YY$pvB.y[rowno])*-1 # NOTE - REVERSED LABELS
    } else {
      YY$pv_angle_nR[rowno] <- YY$pv_angle[rowno] + YY$pvB.x[rowno]
      
    }
  }
  
  t.test(YY$pv_angle_nR,
         mu = 0,
         alternative = "greater") # IMPLICIT LEARNING MEASURE
  cohensD(x = YY$pv_angle_nR)
  
  # NOW DO SAME BUT DO ENDPOINT ANGLE - first, get baseline endpoint angle
  sequence_baselineAE_CCW <- tdf %>%
    filter(prehome == -1) %>%
    filter(task == 1) %>%
    drop_na(endpointang) %>% 
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>% # BASELINE HAS BOTH GARAGES & NEED 3 TRIALS FOR BASELINE BLOCK
    summarise(eaB = mean(endpointang, na.rm = TRUE))
  
  sequence_baselineAE_CW <- tdf %>%
    filter(prehome == 1) %>% 
    filter(task == 1) %>%
    drop_na(endpointang) %>% 
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
    summarise(eaB = mean(endpointang, na.rm = TRUE))
  
  YY <- tdf %>%
    filter(instruction == 'exclude') %>%
    drop_na(endpointang) %>%
    group_by(participant) %>%
    filter(task == min(task)) %>% # the first exclude-strategy task
    filter(trial == min(trial)) %>% # the first exclude-strategy trial
    select(participant, task, prehome, instruction, endpointang) %>%
    distinct(participant, .keep_all = T)
  
  YY <- inner_join(YY, sequence_baselineAE_CW, by = "participant")
  YY <- inner_join(YY, sequence_baselineAE_CCW, by = "participant")
  
  YY$endpointang_nR <- NA
  for (rowno in 1:nrow(YY)) {
    if (YY$prehome[rowno] == 1){ #CW for this condition
      YY$endpointang_nR[rowno] <- (YY$endpointang[rowno] - YY$eaB.x[rowno])*-1
    } else { # it was a CCW-associated trial
      YY$endpointang_nR[rowno] <- (YY$endpointang[rowno] + YY$eaB.y[rowno])
      
    }
  }
  shapiro.test(abs(YY$endpointang_nR)) 
  
  t.test(log(YY$endpointang_nR + 100), # make normal a distn with negative values
         mu = log(100),
         alternative = "greater") # IMPLICIT LEARNING MEASURE

  t.test(YY$endpointang_nR, mu = 0, alternative = "greater")
  # INCLUDE-STRATEGY REACH AEs
  ZZ <- tdf %>%
    filter(instruction == 'include') %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    filter(task == min(task)) %>%
    select(participant, prehome, instruction, pv_angle, pv_angle_n)
  
  ZZ <- inner_join(ZZ, sequence_baselineAE_CW, by = "participant")
  ZZ <- inner_join(ZZ, sequence_baselineAE_CCW, by = "participant")
  
  ZZ$pv_angle_nR <- NA
  for (rowno in 1:nrow(ZZ)) {
    if (ZZ$prehome[rowno] == -1){
      ZZ$pv_angle_nR[rowno] <- (ZZ$pv_angle[rowno] - ZZ$pvB.x[rowno])*-1 # NOTE - REVERSED LABELS
    } else {
      ZZ$pv_angle_nR[rowno] <- (ZZ$pv_angle[rowno] + ZZ$pvB.y[rowno])
      
    }
  }
  
  t.test(ZZ$pv_angle_nR,
         YY$pv_angle_nR,
         paired = TRUE,
         alternative = "greater") # EXPLICIT LEARNING MEASURE
  
  # now do the same but with endpoint angle
  ZZ <- tdf %>%
    filter(instruction == 'include') %>%
    drop_na(endpointang) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    filter(task == min(task)) %>%
    select(participant, prehome, instruction, endpointang)
  
  ZZ <- ZZ %>% distinct(participant, .keep_all = T)
  
  ZZ <- inner_join(ZZ, sequence_baselineAE_CW, by = "participant")
  ZZ <- inner_join(ZZ, sequence_baselineAE_CCW, by = "participant")
  
  ZZ$endpointang_nR <- NA
  for (rowno in 1:nrow(ZZ)) {
    if (ZZ$prehome[rowno] == -1){
      ZZ$endpointang_nR[rowno] <- (ZZ$endpointang[rowno] - ZZ$eaB.x[rowno])
    } else {
      ZZ$endpointang_nR[rowno] <- (ZZ$endpointang[rowno] + ZZ$eaB.y[rowno])*-1
    }
  }
  
  
  t.test(ZZ$endpointang_nR,
         YY$endpointang_nR,
         paired = TRUE,
         alternative = "greater") # EXPLICIT LEARNING MEASURE
  # WITHIN-STRATEGY REACH AEs
  
  sequence_withinAE <- tdf %>%
    filter(task == 6) %>%
    drop_na(pv_angle_n) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    select(participant, prehome, instruction, pv_angle, pv_angle_n)
  
  sequence_withinAE <- inner_join(sequence_withinAE, sequence_baselineAE_CW, by = "participant")
  sequence_withinAE <- inner_join(sequence_withinAE, sequence_baselineAE_CCW, by = "participant")
  
  sequence_withinAE$pv_angle_nR <- NA
  for (rowno in 1:nrow(ZZ)) {
    if (sequence_withinAE$prehome[rowno] == -1){
      sequence_withinAE$pv_angle_nR[rowno] <- (sequence_withinAE$pv_angle[rowno] - sequence_withinAE$pvB.x[rowno])*-1
    } else {
      sequence_withinAE$pv_angle_nR[rowno] <- (sequence_withinAE$pv_angle[rowno] + sequence_withinAE$pvB.y[rowno])
      
    }
  }
  
  t.test(sequence_withinAE$pv_angle_nR,
         YY$pv_angle_nR,
         paired = TRUE,
         alternative = "greater")

}
