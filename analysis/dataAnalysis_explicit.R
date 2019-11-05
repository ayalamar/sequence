##### FOR DUAL EXPLICIT DATA ANALYSIS
setwd('/Users/mayala/Desktop/explicit data')

subject_numbers <- c(1:7,9:12,16) # DUAL EXPLICIT PARTICIPANTS
tasks <- c(0, 1, 2, 3, 4, 5, 6, 7, 8) 
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

##### FUNCTION FOR GETTING STATS AND CLEAN PLOTS
getStatistics <- function(){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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
  
  tdf <- tbl_df(df) # convert to tibble for dplyr
  ### NOTE: this is analyzing COLLAPSED rotations 
  ### -- also still need to subtract baseline 
  
  
  ################################
  ################################
  ################################
  
  # ANALYZE LEARNING
  baselinelast <- tdf %>%
    filter(task == 0) %>%
    filter(trial %in% c(57, 58, 59)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm = TRUE), block = mean(task))
  
  block1 <- tdf %>%
    filter(task == 2) %>%
    filter(trial %in% c(0, 1, 2)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), pl = mean(pathlength, na.rm = TRUE), block = mean(task))
 #block1<- tdf %>% filter(task==2) %>% filter(trial %in% c(3,4,5)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  
  blocklast <- tdf %>%
    filter(task == 6) %>%
    filter(trial %in% c(21, 22, 23)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm = TRUE), pl = mean(pathlength, na.rm = TRUE), block = mean(task))
  
  # participant 4 - file corrupted - problem with task 6 (only trials 0-12 got copied)
  blocklast4 <- tdf %>% 
    filter(participant == 4) %>%
    filter(task == 6) %>%
    filter(trial %in% c(10, 11, 12)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  
  adaptdf <- rbind(block1, blocklast, blocklast4)
  blocklast <- rbind(blocklast, blocklast4)
  
  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  # RM_pv <- aov(pv ~ block + Error(participant/block), data=adaptdf)
  # summary(RM_pv)
  # RM_pl <- aov(pl ~ block + Error(participant/block), data=adaptdf)
  # summary(RM_pl)
  t.test(block1$pv, blocklast$pv, alternative = "greater", paired = TRUE )
  
  # VISUALIZE LEARNING
    
    tdfsmooth <- tdf %>%
      filter(garage_location == 1) %>%
      filter(task == 2)
    
    tdfsmooth2 <- tdf %>%
      filter(garage_location == 1) %>%
      filter(task == 3)
    
    tdfsmooth2$trial <- tdfsmooth2$trial + 359
    tdfsmooth <- rbind(tdfsmooth,tdfsmooth2)
    
    traindf <- tdfsmooth %>%
      group_by(participant) %>%
      group_by(trial) %>%
      summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                SEM_pl = SD_pl/sqrt(length(unique(participant))), Mean_pv = mean(pv_angle, na.rm=TRUE),
                SD_pv = sd(pv_angle, na.rm=TRUE),SEM_pv = SD_pv/sqrt(length(unique(participant))))
    
    # ## giant plot that include all trial datapoints
    # Means <- ggplot(data=traindf, aes(x=trial, y=Mean_pv)) +
    #   geom_line() + 
    #   geom_ribbon(data=traindf, aes(ymin=Mean_pv-SEM_pv, ymax= Mean_pv+SEM_pv), alpha=0.4) +
    #   geom_point(data=tdfsmooth, aes(x=trial, y= pv_angle), colour='steelblue',alpha = 0.1) +
    #   #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, linetype="dashed", color="blue", aes(outfit=fit<<-..y..),n=359) +
    #   #geom_line(aes(x=trial, y=exp(decay$fitted.values)), color = "red")
    #   #stat_smooth(method = "nls", formula = y ~ a * exp(-S * x), 
    #   # method.args = list(start = list(a = 78, S = 0.02)), se = FALSE, #starting values obtained from fit above
    #   # color = "dark red", linetype ="dashed")+
    #   ylim(-50, 50) +
    #   geom_smooth(model=loess) +
    #   ggtitle('Explicit Experiment - Dual CW training') +
    #   ylab("Angular error (Degrees)") +
    #   xlab("Trial") +
    #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    # print(Means)
    # move_layers(Means,"GeomRibbon", position = "top")
    # move_layers(Means,"GeomLine", position = "top")
    # move_layers(Means,"GeomPoint", position = "bottom")
    
    ## summarized plot (block 1, 2, last(topup))
    for (rot in sort(unique(tdf$garage_location))){
      
      dfplot1 <- tdf %>% filter(garage_location == rot) %>% filter(task == 2)
      dfplot2 <- tdf %>% filter(garage_location == rot) %>% filter(task == 3)
      dfplot3 <- tdf %>% filter(garage_location == rot) %>% filter(task == 6)
      dfplot2$trial <- dfplot2$trial + 359
      dfplot3$trial <- dfplot3$trial + 359 + 359
      dfplot <- rbind(dfplot1, dfplot2, dfplot3)
      
      # traindf <- dfplot %>% # GROUP MEANS
      #   group_by(participant) %>%
      #   group_by(trial) %>%
      #   summarise(Mean_pl = mean(pathlength, na.rm = TRUE), 
      #             SD_pl = sd(pathlength, na.rm = TRUE),
      #             SEM_pl = SD_pl/sqrt(length(unique(participant))),
      #             Mean_pv = mean(pv_angle, na.rm = TRUE),
      #             SD_pv = sd(pv_angle, na.rm = TRUE),
      #             SEM_pv = SD_pv/sqrt(length(unique(participant))))
      
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
      
      # bl1 <- traindf[1:3,] %>% mutate(block = 1)
      # bl2 <- traindf[4:6,] %>% mutate(block = 2)
      # bll <- traindf[(nrow(traindf)-2):nrow(traindf),] %>% mutate(block = 7)
      # combobl <- rbind(bl1, bl2, bll)
      # combobl <- combobl %>% group_by(block) %>% summarise(pv = mean(Mean_pv), sem = mean(SEM_pv))
    
    
      
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
        ggtitle('Explicit experiment') 
        
        
      
      print(bltrain)
      move_layers(bltrain, "GeomRibbon", position = "top")
      move_layers(bltrain, "GeomPoint", position = "top")
    }

    ################################
    ################################
    ################################
    
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
      y <- ((inblock - finblock)/inblock)*100
      
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
    
    #move_layers(PIplot,"GeomPoint", position = "bottom")
    
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
    

  # SEMs <- NA
  # for (garage in sort(unique(tdf$garage_location))) {
  #   for (instruct in sort(unique(tdf$instruction))) {
  #     x <- tdf %>% filter(garage_location == garage ) %>% filter(instruction == instruct) %>% filter(trial == 0) %>% group_by(participant) %>%
  #       group_by(trial) %>% summarise(Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
  #                                     SEM_pv = SD_pv/sqrt(length(unique(participant))), 
  #                                     instruction = instruct, garage_location = garage, lowerSEM = Mean_pv-SEM_pv, upperSEM = Mean_pv + SEM_pv)
  #     
  #     if (is.data.frame(SEMs) == TRUE ) {
  #       SEMs <- rbind(SEMs, x)
  #     } else {
  #       SEMs <- x
  #     }
  #   }
  # }
  
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

