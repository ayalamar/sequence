##### FOR DUAL PREQUENCE DATA ANALYSIS
# NOTE: REVERSED LABEL TAGS 
setwd('/Users/mayala/Desktop/preq data')
subject_numbers <- c(1:31)  # PARTICIPANT 1 EXCLUDED DUE TO HUGE BASELINE BIASES
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(0) # IF THIS IS 0, NO SCALING WILL BE DONE IN TASKANALYSIS()

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

##### FUNCTION FOR GETTING STATS AND CLEAN PLOTS
getStatistics <- function(){
  
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gginnards)
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
  tdf <- tdf %>% filter(participant != 1)  # PARTICIPANT 1 EXCLUDED DUE TO HUGE BASELINE BIASES
  # NOW ANALYZING COLLAPSED ROTATIONS
  
  ################################
  ################################
  ################################
  
  # ANALYZE LEARNING
  baselinelast <- tdf %>%
    filter(task == 0) %>%
    filter(trial %in% c(57, 58, 59)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              block = mean(task))
  
  block1 <- tdf %>%
    filter(task == 3) %>%
    filter(trial %in% c(0, 1, 2)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              block = mean(task))
  
  blocklast <- tdf %>%
    filter(task == 7) %>%
    filter(trial %in% c(21, 22, 23)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              block = mean(task))
  
  adaptdf<- rbind(block1, blocklast)
  
  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  # RM_pv <- aov(pv ~ block + Error(participant/block), data=adaptdf)
  # summary(RM_pv)
  # RM_pl <- aov(pl ~ block + Error(participant/block), data=adaptdf)
  # summary(RM_pl)
  boxplot(adaptdf$pv) # remove the intense outliers
  
  adaptdf2 <- adaptdf %>%
    filter(participant != "4") %>%
    filter(participant != "6") %>%
    filter(participant != "24")    # adaptdf[-c(35,37,55),] # outliers
  
  t.test(adaptdf2$pv[which(adaptdf2$block == "3")],
         adaptdf2$pv[which(adaptdf2$block == "7")],
         alternative ="greater", paired = TRUE )
  
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
    # traindf <- dfplot %>% group_by(participant) %>% group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
    #                                                                               SEM_pl = SD_pl/sqrt(length(unique(participant))), Mean_pv = mean(pv_angle, na.rm=TRUE),                                                                                SD_pv = sd(pv_angle, na.rm=TRUE),SEM_pv = SD_pv/sqrt(length(unique(participant))))
    # bl1 <- traindf[1:3,] %>% mutate(block = 1)
    # bl2 <- traindf[4:6,] %>% mutate(block = 2)
    # bll <- traindf[(nrow(traindf)-2):nrow(traindf),] %>% mutate(block = 7)
    # combobl <- rbind(bl1, bl2, bll)
    # combobl <- combobl %>% group_by(block) %>% summarise(pv = mean(Mean_pv, na.rm=TRUE), sem = mean(SEM_pv,na.rm=TRUE))
    # 
    # bltrain <- ggplot(data=combobl, aes(x=block, y=pv)) +
    #   geom_point() +
    #   geom_line(data=combobl, aes(x=block, y=pv)) +
    #   geom_ribbon(data=combobl, aes(ymin=pv-sem, ymax= pv+sem), alpha=0.4) +
    #   ylim(-50,50) +
    #   xlim(1,7) +
    #   coord_fixed(ratio = 1/7) +
    #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    #   ggtitle('Dual Lead-in - Blocked training')
    
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
      ggtitle('Dual Lead-in') 
    
    print(bltrain)
    move_layers(bltrain, "GeomRibbon", position = "top")
    move_layers(bltrain, "GeomPoint", position = "top")
  }
  
  ################################
  ################################
  ################################
  
  # ANALYZE PERCENT IMPROVEMENT
  
  dfplot1 <- tdf  %>%
    filter(task == 3) %>%
    filter(participant != "4") %>%
    filter(participant != "6") %>%
    filter(participant != "24")
  
  dfplot2 <- tdf  %>%
    filter(task == 5)%>%
    filter(participant != "4") %>%
    filter(participant != "6") %>%
    filter(participant != "24")
  
  dfplot3 <- tdf  %>%
    filter(task == 7) %>%
    filter(participant != "4") %>%
    filter(participant != "6") %>%
    filter(participant != "24")
  
  dfplot2$trial <- dfplot2$trial + 359
  dfplot3$trial <- dfplot3$trial + 359 + 359
  dfplot <- rbind(dfplot1, dfplot2, dfplot3)
  
  PI <- c()
  
  for (ppno in sort(unique(dfplot$participant))) {
    inblock <- dfplot %>% filter(participant == ppno) %>% filter(trial==0|trial==1|trial==2)
    inblock <- mean(inblock$pv_angle_n, na.rm=TRUE)
    finblock <- dfplot %>% filter(participant == ppno) %>% filter(trial == 740|trial==741|trial==742) 
    finblock <- mean(finblock$pv_angle_n, na.rm=TRUE)
    y <- ((inblock - finblock)/inblock)*100
    
    if (is.null(PI) == TRUE ) {
      PI <- y
    } else {
      PI <- c(PI, y)
    }
  }
  # outliers
  boxplot(PI) # one outlier here but we already removed outliers above- don't wanna overdo it
  
  meanPI <- mean(PI, na.rm=TRUE)
  semPI <- sd(PI, na.rm=TRUE)/sqrt(length(PI))
  
  t.test(PI, alternative = "greater") # IMPROVEMENT ACROSS TRAINING
  
  PI <- tbl_df(PI) # ADD LABELS TO PLOT
  PI <- PI %>% mutate(participant = 1:n(),
                      group = "prequence",
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
    ggtitle("Dual Prequence") +
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
  
  ## more visualizations of AE
  # SEMs <- NA
  # for (ph in sort(unique(tdf$prehome))) {
  #   for (instruct in sort(unique(tdf$instruction))) {
  #     x <- tdf %>% filter(prehome == ph ) %>%
  #       filter(instruction == instruct) %>%
  #       filter(trial == 0) %>%
  #       group_by(participant) %>% 
  #       group_by(trial) %>%
  #       summarise(Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
  #                                     SEM_pv = SD_pv/sqrt(length(unique(participant))), 
  #                                     instruction = instruct, prehome = ph, lowerSEM = Mean_pv-SEM_pv, upperSEM = Mean_pv + SEM_pv)
  #     
  #     if (is.data.frame(SEMs) == TRUE ) {
  #       SEMs <- rbind(SEMs, x)
  #     } else {
  #       SEMs <- x
  #     }
  #   }
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
  
  
  # SEMs <- NA
  #  for (garage in sort(unique(tdf$prehome))) {
  #    for (instruct in sort(unique(tdf$instruction))) {
  #      x <- tdf %>% filter(prehome == garage ) %>% filter(instruction == instruct) %>% filter(trial == 0) %>% group_by(participant) %>%
  #      group_by(trial) %>% summarise(Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
  #                                    SEM_pv = SD_pv/sqrt(length(unique(participant))), 
  #                                    instruction = instruct, prehome = garage, lowerSEM = Mean_pv-SEM_pv, upperSEM = Mean_pv + SEM_pv)
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
                   aes(x = instruction, y = Mean_pv, fill = as.factor(prehome))) +
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
    ggtitle("Dual Prequence") +
    coord_fixed(ratio = 1/13) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank(), legend.position = "none") +
    scale_y_continuous(breaks = seq(-30, +30, 10), limits = c(-30, 30))
  # move_layers(IEbars, "GeomPoint", position = "bottom")
  
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

return()
}
