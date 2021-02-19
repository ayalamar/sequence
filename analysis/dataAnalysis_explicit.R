##### FOR DUAL EXPLICIT DATA ANALYSIS
setwd('~/science/repos/sequence')
setwd('explicit data')
rm(list=ls())

subject_numbers <- c(1:7,9:12,16) # DUAL EXPLICIT PARTICIPANTS
tasks <- c(0, 1, 2, 3, 4, 5, 6, 7, 8) 
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

##### FUNCTION FOR PLOTTING BINNED RAW DATA ONLY
plotData <- function(){
  
  filename <- sprintf('md_analysis/allTaggedData_n%d_%s.csv', 
                      length(subject_numbers), outfile_suffix)
  
  df <- read.csv(filename, header = TRUE)
  
  df <- df[-c(1),] # get rid of that random first row of NAs
  
  df <- df %>% # make bins for plotting
    group_by(participant) %>%
    group_by(task) %>%
    mutate(binno = bin(trial,
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
  rotations <- sort(unique(tdf$garage_location)) # use this because no-cursor 
  # trials are labeled rotation = 0
  
  wholePI <- NA
  
  for (rotationno in rotations) {
    
    print(rotationno)
    
    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every bin -- use for plotting 
      #learning curves
      dfname <- sprintf('rotation%d_task%d_means', rotationno, taskno)
      print(dfname)

      taskmeans <- tdf %>%
        filter(task == taskno) %>% 
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        group_by(binno) %>%
        summarise(Mean_pl = mean(pathlength, na.rm = TRUE),
                  SD_pl = sd(pathlength, na.rm = TRUE),
                  SEM_pl = SD_pl/sqrt(length(unique(participant))),
                  Mean_pv = mean(pv_angle, na.rm = TRUE),
                  SD_pv = sd(pv_angle, na.rm = TRUE),
                  SEM_pv = SD_pv/sqrt(length(unique(participant))))
      
      # plot each task learning curve
      taskplot <- ggplot(data = taskmeans,
                         aes(x = binno, y = Mean_pv, group = 1)) +
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
      
      # print(taskplot)
      plotfilename <- sprintf('md_analysis/Full_LCs_%s.svg', dfname)
      # ggsave(file = plotfilename,
      #        plot = taskplot,
      #        height = 10, dpi = 96, units = "cm")
      # 
    }
    
    # GET PI FOR EACH PARTICIPANT, STORE, AND PLOT
    if (rotationno == 1){
      
      PI_inblock <- tdf %>%
        filter(task == 2) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == min(binno_no) | binno_no == min(binno_no)+1) %>%
        summarise(Mean_pv_inblock = mean(pv_angle, na.rm=TRUE))
      
      PI_midblock <- tdf %>%
        filter(task == 2) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
        summarise(Mean_pv_midblock = mean(pv_angle, na.rm=TRUE))
      
      PI_finblock <- tdf %>%
        filter(task == 6) %>%
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
               PI_r30_2 = ((30 - Mean_pv_midblock)/30)*100,
               oldPI_2 = ((Mean_pv_inblock - Mean_pv_finblock)/30)*100)
      
    } else { 
      
      tdf$pv_angle_neg <- tdf$pv_angle*-1
      
      PI_inblock <- tdf %>%
        filter(task == 2) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == min(binno_no) | binno_no == min(binno_no)+1) %>%
        summarise(Mean_pv_inblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_midblock <- tdf %>%
        filter(task == 2) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)-1 | binno_no == max(binno_no)) %>%
        summarise(Mean_pv_midblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_finblock <- tdf %>%
        filter(task == 6) %>%
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
               PI_r30_2 = ((30 - Mean_pv_midblock)/(30))*100,
               oldPI_2 = ((Mean_pv_inblock - Mean_pv_finblock)/30)*100)
      
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
  
  write.csv(wholePI, "md_analysis/PI_inblock.csv", row.names = FALSE)
  
  wholePI <- wholePI %>% 
    #filter(isoutlier == FALSE) %>%
    group_by(rot) %>%
    mutate(PImean = mean(oldPI_2, na.rm = TRUE),
           PIsd = sd(oldPI_2, na.rm = TRUE),
           SEM = PIsd/sqrt(length(unique(participant))))
  
  PIplot <- ggplot(data = wholePI, aes(x = rot, y = oldPI_2, fill = rot)) +
    stat_summary(fun.y = mean, geom = "bar", na.rm = TRUE) +
    geom_errorbar(data = wholePI,
                  mapping = aes(x = rot, y = oldPI_2,
                                ymin = PImean - SEM , ymax = PImean + SEM),
                  width = 0.1, size = 0.5, color = "black",
                  position = position_dodge(width = 0.9)) +
    geom_beeswarm(data = wholePI, aes(x = rot, y = oldPI_2),
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
  
  tdf <- tbl_df(df) # convert to tibble for dplyr
  ### NOTE: this is analyzing COLLAPSED rotations 
  
  #### ANALYZE LEARNING ####
  baselinelast <- tdf %>%
    filter(task == 0) %>%
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), 
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
  
  block1 <- tdf %>%
    filter(task == 2) %>%
    group_by(participant) %>%
    filter(trial %in% c(min(trial), min(trial)+1, min(trial)+2)) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), 
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))
  
  blocklast <- tdf %>%
    filter(task == 6) %>%
    group_by(participant) %>%
    filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              md = mean(maxdev, na.rm = TRUE),
              block = mean(task))

  adaptdf <- rbind(block1, blocklast)

  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  # error proxy 1 - peak velocity angle
  shapiro.test(adaptdf$pv) # if sig., wicoxon rank test
  t.test(block1$pv, blocklast$pv,
         alternative = "greater", paired = TRUE )
  
  # error proxy 2 - pathlength
  boxplot(adaptdf$pl)
  shapiro.test(adaptdf$pl) # if significant - log to get a more normal dist'n
  t.test(log(block1$pl), log(blocklast$pl),
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS
  cohensD(x = log(block1$pl), y = log(blocklast$pl),
          method = "paired")
  
  # error proxy 3 - maximum deviation angle
  boxplot(adaptdf$md)
  shapiro.test(adaptdf$md) #  if sig., wicoxon rank test
  t.test(log(block1$md), log(blocklast$md),
         alternative = "greater",
         paired = TRUE) # COLLAPSED ROTATIONS
  wilcox.test(log(block1$md), log(blocklast$md), 
              alternative = "greater", 
              paired = TRUE)
  
  # save file for easy import for cross experiment comparisons
  # write.csv(block1, "md_analysis/block1_instructed.csv", row.names = FALSE)
  #### #### #### #### #### #### ####
  
  #### COMPARING BETWEEN GROUPS ####
  # compare block 1 of instructed vs. noninstructed
  block1_noninstructed <- read.csv("md_analysis/block1_noninstructed.csv", header = TRUE)
  
  # error proxy 1 - pv
  t.test(block1$pv, block1_noninstructed$pv, paired = FALSE)
  cohensD(x = block1$pv, y = block1_noninstructed$pv, method = "corrected")
  # error proxy 2 - pathlength
  t.test(log(block1$pl), log(block1_noninstructed$pl), paired = FALSE)
  # error proxy 3 - maxdev
  t.test(block1$md, block1_noninstructed$md, paired = FALSE)
  wilcox.test(block1$md, block1_noninstructed$md, paired = FALSE)
  #### #### #### #### #### #### ####
  
  #### VISUALIZE LEARNING ####
    tdfsmooth <- tdf %>%
      filter(garage_location == 1) %>%
      filter(task == 2)
    
    tdfsmooth2 <- tdf %>%
      filter(garage_location == 1) %>%
      filter(task == 3)
    
    tdfsmooth2$trial <- tdfsmooth2$trial + 359
    tdfsmooth <- rbind(tdfsmooth, tdfsmooth2)
    
    traindf <- tdfsmooth %>%
      group_by(participant) %>%
      group_by(trial) %>%
      summarise(Mean_pl = mean(pathlength, na.rm=TRUE), 
                SD_pl = sd(pathlength, na.rm=TRUE),
                SEM_pl = SD_pl/sqrt(length(unique(participant))),
                Mean_pv = mean(pv_angle, na.rm=TRUE),
                SD_pv = sd(pv_angle, na.rm=TRUE),
                SEM_pv = SD_pv/sqrt(length(unique(participant))))
    
    ## summarized plot (block 1, 2, last(topup))
    for (rot in sort(unique(tdf$garage_location))){
      
      dfplot1 <- tdf %>% filter(garage_location == rot) %>% filter(task == 2)
      dfplot2 <- tdf %>% filter(garage_location == rot) %>% filter(task == 3)
      dfplot3 <- tdf %>% filter(garage_location == rot) %>% filter(task == 6)
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
      write.csv(traindf, file = outfile_name, row.names = FALSE)  
    
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
        ggtitle('Explicit experiment') 

      
      print(bltrain)
      move_layers(bltrain, "GeomRibbon", position = "top")
      move_layers(bltrain, "GeomPoint", position = "top")
    }

    # ANALYZE PERCENT IMPROVEMENT
    PI <- c()
    for (ppno in sort(unique(dfplot$participant))) {
      dfplot1 <- tdf  %>% filter(task == 2)
      dfplot2 <- tdf  %>% filter(task == 3)
      dfplot3 <- tdf  %>% filter(task == 6)
      dfplot2$trial <- dfplot2$trial + 359
      dfplot3$trial <- dfplot3$trial + 359 + 359
      dfplot <- rbind(dfplot1, dfplot2, dfplot3)
      
      inblock <- dfplot %>% filter(participant == ppno) %>% filter(trial==0|trial==1|trial==2)
      inblock <- mean(inblock$pv_angle_n, na.rm = TRUE)
      finblock <- dfplot %>% filter(participant == ppno) %>% filter(trial == 740|trial==741|trial==742) 
      finblock <- mean(finblock$pv_angle_n, na.rm = TRUE)
      
      y <- ((inblock - finblock)/30)*100
      #y <- ((inblock - finblock)/inblock)*100
      
      if (is.null(PI) == TRUE ) {
        PI <- y
      } else {
        PI <- c(PI, y)
      }
    }
    
    # PP4 ONLY UP TO 12 
    PI4_1 <- tdf %>%
      filter(participant == 4) %>%
      filter(task==2) %>%
      filter(trial %in% c(0,1,2)) %>%
      summarise(meanpi=mean(pv_angle_n, na.rm = TRUE))
    PI4_2 <- tdf %>%
      filter(participant == 4) %>%
      filter(task==6) %>%
      filter(trial %in% c(10,11,12)) %>%
      summarise(meanpi=mean(pv_angle_n, na.rm = TRUE))
    y <- ((PI4_1 - PI4_2)/PI4_1)*100
    
    PI[4] <- y$meanpi[1]
    boxplot(PI) # NO OUTLIERS
    
    meanPI <- mean(PI, na.rm=TRUE)
    semPI <- sd(PI, na.rm=TRUE)/sqrt(length(PI))
    
    t.test(PI, mu = 0,
           alternative = "greater") # IMPROVEMENT ACROSS TRAINING
    
    PI <- tbl_df(PI) # ADD LABELS TO PLOT
    PI <- PI %>% mutate(participant = 1:n(),
                        group = "explicit",
                        lowerSEM = meanPI - semPI,
                        upperSEM = meanPI + semPI)
    
    PIplot <- ggplot(data = PI, aes(x = group, y = value)) +
      stat_summary(fun.y = mean, geom = "bar", na.rm = TRUE) +
      geom_errorbar(data = PI, mapping = aes(x = group, y = value, 
                                             ymin = PI$lowerSEM, ymax = PI$upperSEM),
                    width = 0.1, size = 0.5, color = "black",
                    position = position_dodge(width = 0.9)) +
      geom_beeswarm(data = PI, aes(x = group, y = value),
                    alpha = 1/7,
                    dodge.width = 2, cex = 3,
                    stroke = 0.3) +
      #geom_point(data = PI, aes(x = group, y = value), size = 1, alpha = 1/20) +
      ylab("Percentage Improvement") +
      ggtitle("Dual Explicit") +
      coord_fixed(ratio = 1/30) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_blank(), legend.position = "none") +
    scale_y_continuous(breaks = seq(-150, +150, 50), limits = c(-150,150))
    
    print(PIplot)

    ################################
    ################################
    ################################
    
    # ANALYZE REACH AFTEREFFECTS
    
    # COLLAPSED ROTATIONS FIRST

    baselineAE <- tdf %>%
      filter(task == 1) %>%
      drop_na(pv_angle_n) %>% 
      filter(trial %in% c(21,22,23)) %>%
      group_by(participant) %>%
      summarise(pv = mean(pv_angle_n, na.rm = TRUE),
                instruction = 'baseline')
    
    excludeAE <- tdf %>%
      filter(instruction == 'exclude') %>%
      drop_na(pv_angle_n) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      summarise(pv = mean(pv_angle_n, na.rm = TRUE),
                instruction = 'exclude')
    
    includeAE <- tdf %>%
      filter(instruction == 'include') %>%
      drop_na(pv_angle_n) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      summarise(pv = mean(pv_angle_n, na.rm = TRUE),
                instruction = 'include')
    
    # AE ANALYSIS BUT BASELINE SEPARATE
    t.test(excludeAE$pv, baselineAE$pv,
           alternative = "greater",
           paired = TRUE) 
    
    t.test(includeAE$pv - baselineAE$pv,
           excludeAE$pv - baselineAE$pv,
           alternative = "greater",
           paired = TRUE )

    # SUBTRACT BASELINE FOR FIGURES & MORE AE ANALYSIS
    excludeAE$pv <- excludeAE$pv - baselineAE$pv
    includeAE$pv <- includeAE$pv - baselineAE$pv
    
    # AE ANALYSIS, BASELINE SUBTRACTED
    t.test(excludeAE$pv,
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
      filter(trial %in% c(18:23)) %>% # BASELINE HAS BOTH GARAGES & NEED 3 TRIALS FOR BASELINE BLOCK
      group_by(participant) %>%
      summarise(pvB = mean(pv_angle, na.rm = TRUE))
    
    sequence_baselineAE_CW <- tdf %>%
      filter(garage_location == 1) %>%
      filter(task == 1) %>%
      drop_na(pv_angle_n) %>% 
      filter(trial %in% c(18:23)) %>%
      group_by(participant) %>%
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
      ggtitle("Dual Explicit") +
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
    
    # EXCLUDE-STRATEGY REACH AEs
    YY <- tdf %>%
      filter(instruction == 'exclude') %>%
      drop_na(pv_angle_n) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      filter(task == min(task)) %>%
      select(participant, garage_location, instruction, pv_angle, pv_angle_n)
    
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
    
    # INCLUDE-STRATEGY REACH AEs
    ZZ <- tdf %>%
      filter(instruction == 'include') %>%
      drop_na(pv_angle_n) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      filter(task == min(task)) %>%
      select(participant, garage_location, instruction, pv_angle, pv_angle_n)
    
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
    
    #### NOW DO SAME BUT DO ENDPOINT ANGLE - first, get baseline endpoint angle  ####
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
    
    
    ZZ <- tdf %>%
      filter(instruction == 'include') %>%
      drop_na(endpointang) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      filter(task == min(task)) %>%
      select(participant, garage_location, instruction, endpointang)
    
    ZZ <- ZZ %>%
      distinct(participant, .keep_all = T)
    
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
    cohensD(ZZ$endpointang_nR,
            YY$endpointang_nR,
            method = "paired")
    
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
    
    #### prep for cross experiment comparisons for PDP no-cursors ####
    ### PV FILE
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
    
    explicit_AEdf <- rbind(sequence_baselineAE,
                         sequence_excludeAE,
                         sequence_includeAE)
    
    explicit_AEdf <- explicit_AEdf %>% 
      mutate(group = "instructed")
    
    write.csv(explicit_AEdf,
              "md_analysis/explicit_AEdf_pv.csv",
              row.names = FALSE)
    
    ### ENDPOINT ANGLE FILE
    tdf$endpointang_n <- abs(tdf$endpointang)
    
    sequence_baselineAE <- tdf %>%
      filter(task == 1) %>%
      drop_na(endpointang_n) %>% 
      group_by(participant) %>%
      filter(trial %in% c(max(trial)-2, max(trial)-1, max(trial))) %>%
      summarise(pv = mean(endpointang_n, na.rm = TRUE),
                instruction = 'baseline')
    
    sequence_excludeAE <- tdf %>%
      filter(instruction == 'exclude') %>%
      drop_na(endpointang_n) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      summarise(pv = mean(endpointang_n, na.rm = TRUE),
                instruction = 'exclude')
    
    sequence_includeAE <- tdf %>%
      filter(instruction =='include') %>%
      drop_na(endpointang_n) %>% 
      group_by(participant) %>%
      filter(trial == min(trial)) %>%
      summarise(pv = mean(endpointang_n, na.rm = TRUE),
                instruction = 'include')
    
    explicit_AEdf <- rbind(sequence_baselineAE,
                         sequence_excludeAE,
                         sequence_includeAE)
    
    explicit_AEdf <- explicit_AEdf %>% 
      mutate(group = "instructed")
    
    write.csv(explicit_AEdf,
              "md_analysis/explicit_AEdf_ea.csv",
              row.names = FALSE)
}



