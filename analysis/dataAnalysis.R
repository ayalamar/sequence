##### MASTER SCRIPT FOR ANALYSIS AND PLOTS FOR SEQ, CONSEQ, STATIC EXPERIMENTS
setwd('/Users/mayala/Desktop/conseq data')

subject_numbers <- c(1:7, 9:31) # CONSEQUENCE EXPERIMENT
#subject_numbers <- c(3:9, 11:17, 20:33) # SEQUENCE EXPERIMENT
#subject_numbers <- c(1:12) # STATIC EXPERIMENT
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # IF THIS IS 0, NO SCALING WILL BE DONE IN TASKANALYSIS()

##### FUNCTION FOR PLOTTING RAW DATA ONLY
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
  
  tdf <- tbl_df(df) 
  rotations <- sort(unique(tdf$garage_location)) # use this because no-cursor trials are labeled rotation = 0
  
  for (rotationno in rotations) {
    print(rotationno)
    
    # filter by task
    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every trial -- use for plotting learning curves
      dfname<- sprintf('rotation%d_task%d_means', rotationno, taskno)
      print(dfname)
      
      taskmeans<- tdf %>%
        filter(task == taskno) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        group_by(trial) %>%
        summarise(Mean_pl = mean(pathlength, na.rm=TRUE),
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

##### FUNCTION FOR GETTING STATS AND CLEAN PLOTS
getStatistics <- function(){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gginnards)
  library(Hmisc)
  library(ggbeeswarm)

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
  
  tdf <- tbl_df(df) 
  ### NOTE: this is analyzing COLLAPSED rotations 

  ################################
  ################################
  ################################
  
  # ANALYZE LEARNING
  baselinelast <- tdf %>%
    filter(task == 0) %>%
    filter(trial %in% c(57,58,59)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), pl = mean(pathlength, na.rm = TRUE), block = mean(task))
  
  block1 <- tdf %>%
    filter(task == 3) %>%
    filter(trial %in% c(0,1,2)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), pl = mean(pathlength, na.rm = TRUE), block = mean(task))
  
  blocklast <- tdf %>% 
    filter(task == 7) %>%
    filter(trial %in% c(21,22,23)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), pl = mean(pathlength, na.rm = TRUE), block = mean(task))
 
  adaptdf <- rbind(block1, blocklast)

  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  # RM_pv <- aov(pv ~ block + Error(participant/block), data=adaptdf)
  # summary(RM_pv)
  # RM_pl <- aov(pl ~ block + Error(participant/block), data=adaptdf)
  # summary(RM_pl)
   t.test(block1$pv, blocklast$pv, alternative = "greater", paired = TRUE ) # COLLAPSED ROTATIONS
  
  # mod1 <- ezANOVA(data = adaptdf,
  #                 dv = pv,
  #                 wid = participant,
  #                 within=.(block),
  #                 detailed = TRUE)
  
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
    
    # traindf <- dfplot %>%
    #   group_by(participant) %>%
    #   group_by(trial) %>%
    #   summarise(Mean_pl = mean(pathlength, na.rm=TRUE),
    #             SD_pl = sd(pathlength, na.rm=TRUE),
    #             SEM_pl = SD_pl/sqrt(length(unique(participant))),
    #             Mean_pv = mean(pv_angle, na.rm=TRUE),
    #             SD_pv = sd(pv_angle, na.rm=TRUE),
    #             SEM_pv = SD_pv/sqrt(length(unique(participant))))
    # 
    # bl1 <- traindf[1:3,] %>% mutate(block = 1)
    # bl2 <- traindf[4:6,] %>% mutate(block = 2)
    # bll <- traindf[(nrow(traindf)-2):nrow(traindf),] %>% mutate(block = 7)
    # 
    # combobl <- rbind(bl1, bl2, bll)
    # combobl <- combobl %>%
    #   group_by(block) %>% 
    #   summarise(pv = mean(Mean_pv, na.rm=TRUE),sem = mean(SEM_pv,na.rm=TRUE))
    
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
      geom_line(data=ppdf_full, aes(x=block, y=blockmean, colour=as.factor(participant)), alpha = 0.1) + 
      geom_ribbon(data=traindf, aes(ymin=pv-sem, ymax= pv+sem), alpha=0.4) +
      coord_fixed(ratio = 1/7) +
      ylim(-50,50) +
      xlim(1,7) +
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
  
  # ANALYZE PERCENT IMPROVEMENT
  
  PI <- c()
  for (ppno in sort(unique(dfplot$participant))) {
    
    dfplot1 <- tdf  %>% filter(task == 3)
    dfplot2 <- tdf  %>% filter(task == 5)
    dfplot3 <- tdf  %>% filter(task == 7)
    dfplot2$trial <- dfplot2$trial + 359
    dfplot3$trial <- dfplot3$trial + 359 + 359
    dfplot <- rbind(dfplot1, dfplot2, dfplot3)
    
    inblock <- dfplot %>%
      filter(participant == ppno) %>%
      filter(trial == 0|trial == 1|trial == 2)
    
    inblock <- mean(inblock$pv_angle_n, na.rm = TRUE)
    
    finblock <- dfplot %>%
      filter(participant == ppno) %>%
      filter(trial == 740|trial == 741|trial == 742) 
    
    finblock <- mean(finblock$pv_angle_n, na.rm = TRUE)
    y <- ((inblock - finblock)/30)*100
    #y <- ((inblock - finblock)/inblock)*100
    
    if (is.null(PI) == TRUE ) {
      PI <- y
    } else {
      PI <- c(PI, y)
    }
  }
  
  boxplot(PI) # CHECK FOR ANY OUTLIERS

  meanPI <- mean(PI, na.rm = TRUE)
  semPI <- sd(PI, na.rm = TRUE)/sqrt(length(PI))
  
  t.test(PI, alternative = "greater") # IMPROVEMENT ACROSS TRAINING
  
  PI <- tbl_df(PI) # ADD LABELS TO PLOT
  PI <- PI %>% mutate(participant = 1:n(),
                      group = "sequence",
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
    ggtitle("Dual Static - PI relative to 30") +
    coord_fixed(ratio = 1/30) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank(), legend.position = "none") +
  scale_y_continuous(breaks = seq(-150, +150, 50), limits = c(-150,150))
  
  #move_layers(PIplot,"GeomPoint", position = "bottom")
  
  print(PIplot)
  
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

  
# SEMs <- NA
#  for (garage in sort(unique(tdf$garage_location))) {
#    for (instruct in sort(unique(tdf$instruction))) {
#      x <- tdf %>% filter(garage_location == garage ) %>% filter(instruction == instruct) %>% filter(trial == 0) %>% group_by(participant) %>%
#      group_by(trial) %>% summarise(Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
#                                    SEM_pv = SD_pv/sqrt(length(unique(participant))), 
#                                    instruction = instruct, garage_location = garage, lowerSEM = Mean_pv-SEM_pv, upperSEM = Mean_pv + SEM_pv)
#      
#      if (is.data.frame(SEMs) == TRUE ) {
#        SEMs <- rbind(SEMs, x)
#      } else {
#        SEMs <- x
#      }
#      }
#    }

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


################################
################################
################################

# ANALYZE WITHIN-AE VS. EXCLUDE-AE

# WITHIN-STRATEGY REACH AEs

sequence_withinAE <- tdf %>%
  filter(task == 6) %>%
  drop_na(pv_angle_n) %>% 
  group_by(participant) %>%
  filter(trial == min(trial)) %>%
  select(participant, garage_location, instruction, pv_angle, pv_angle_n)

sequence_withinAE <- inner_join(sequence_withinAE, sequence_baselineAE_CW, by = "participant")
sequence_withinAE <- inner_join(sequence_withinAE, sequence_baselineAE_CCW, by = "participant")

sequence_withinAE$pv_angle_nR <- NA
for (rowno in 1:nrow(ZZ)) {
  if (sequence_withinAE$garage_location[rowno] == -1){
    sequence_withinAE$pv_angle_nR[rowno] <- (sequence_withinAE$pv_angle[rowno] - sequence_withinAE$pvB.x[rowno])
  } else {
    sequence_withinAE$pv_angle_nR[rowno] <- (sequence_withinAE$pv_angle[rowno] + sequence_withinAE$pvB.y[rowno])*-1
    
  }
}

t.test(sequence_withinAE$pv_angle_nR,
       YY$pv_angle_nR,
       paired = TRUE,
       alternative = "greater")
}

##### CROSS EXPERIMENT COMPARISONS
crossComparisons <- function(){
################################
################################
################################

# CROSS EXPERIMENT COMPARISONS 

# E/I COMPARISONS ACROSS INSTRUCTED AND NON-INSTRUCTED EXPERIMENTS 

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

}