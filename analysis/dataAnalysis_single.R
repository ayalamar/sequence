##### FOR SINGLE CW AND CCW CONTROL GROUPS
setwd('~/science/repos/sequence')
setwd('single CW data')
rm(list=ls())

subject_numbers <- c(1:10) # SAME FOR SINGLE CW & CCW GROUPS
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9) 
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
library(OneR)

##### FUNCTION FOR PLOTTING BLOCKED RAW DATA & PI
plotData <- function(){
  
  filename <- sprintf('md_analysis/allTaggedData_n%d_%s.csv',
                      length(subject_numbers), outfile_suffix)
  
  df <- read.csv(filename, header = TRUE)
  
  df <- df[-c(1),] # get rid of that random first row of NAs
  
  df <- df %>%
    group_by(participant) %>%
    group_by(task) %>%
    mutate(binno = bin(trial,  # make bins for plotting
                       nbins = (max(trial)+1)/3,
                       labels = c(1:((max(trial)+1)/3))))
  
  # first replace outlier-tagged trials with NA so they don't get plotted
  # for (rowno in 1:nrow(df)) {
  #   if (df$isoutlier[rowno] == TRUE) {
  #     df$pathlength[rowno] <- NA
  #     df$pv_angle[rowno] <- NA
  #   }
  # }
  
  tdf <- tbl_df(df) # convert to tibble for dplyr
  rotations <- sort(unique(tdf$garage_location)) # use this because no-cursor trials are labeled rotation = 0
  # tdf <- tdf %>%
  #   filter(participant != 6) %>% # super outlying participants for CCW
  #   filter(participant != 7) %>% # super outlying participants for CCW
  #   filter(participant != 7) %>% # super outlying participants for CW
  #   filter(participant != 2) %>% # super outlying participants for CW
  
  for (rotationno in rotations) {
    
    print(rotationno)

    for (taskno in sort(unique(tdf$task))) {
      # this creates a column of MEANS for every bin -- use for plotting learning curves
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
      
      # plot each task learning curve
      taskplot <- ggplot(data = taskmeans, aes(x = binno, y = Mean_pv, group = 1)) +
        geom_line() + 
        geom_ribbon(aes(ymin = Mean_pv-SEM_pv, ymax = Mean_pv+SEM_pv),
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
      
      # print(taskplot)
      # plotfilename <- sprintf('md_analysis/Full_LCs_%s.svg', dfname)
      # ggsave(file = plotfilename,
      #        plot = taskplot,
      #        height = 10, dpi = 96, units = "cm")
            
    }
    
    if (rotationno == 1){
      
      PI_inblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == min(binno_no)) %>%
        summarise(Mean_pv_inblock = mean(pv_angle, na.rm=TRUE))
      
      PI_midblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)) %>%
        summarise(Mean_pv_midblock = mean(pv_angle, na.rm=TRUE))
      
      PI_finblock <- tdf %>%
        filter(task == 7) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)) %>%
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
        filter(binno_no == min(binno_no)) %>%
        summarise(Mean_pv_inblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_midblock <- tdf %>%
        filter(task == 3) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)) %>%
        summarise(Mean_pv_midblock = mean(pv_angle_neg, na.rm=TRUE))
      
      PI_finblock <- tdf %>%
        filter(task == 7) %>%
        filter(garage_location == rotationno) %>%
        group_by(participant) %>%
        mutate(binno_no = as.double(binno)) %>%
        filter(binno_no == max(binno_no)) %>%
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
    

    wholePI <- PI_df
    wholePI$isoutlier <- FALSE # remove any intense outliers
    wholePI$isoutlier[which(wholePI$PI_2 %in% boxplot(wholePI$PI_2)$out)] <- TRUE
    
    wholePI <- wholePI %>% 
      filter(isoutlier == FALSE) %>%
      group_by(rot) %>%
      mutate(PImean = mean(PI_2, na.rm = TRUE),
             PIsd = sd(PI_2, na.rm = TRUE),
             SEM = PIsd/sqrt(length(unique(participant))))
    
    PIplot <- ggplot(data = wholePI, aes(x = rot, y = PImean, fill = rot)) +
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
      scale_y_continuous(breaks = seq(-150, +150, 50), limits = c(-150, 200))
    
    print(PIplot)
    plotfilename <- sprintf('md_analysis/PI_v2.svg')
    ggsave(file = plotfilename,
           plot = PIplot,
           height = 10, dpi = 96, units = "cm")
    
  }
  
}

  
  

##### FUNCTION FOR GETTING STATS AND CLEAN PLOTS
getStatistics <- function(){
  
  filename <- sprintf('check/allTaggedData_n%d_%s.csv', length(subject_numbers), outfile_suffix)
  df <- read.csv(filename, header = TRUE)
  
  df <- df[-c(1),] # get rid of that random first row of NAs
  
  # first replace outlier-tagged trials with NA so they don't get analyzed
  for (rowno in 1:nrow(df)) {
    if (df$isoutlier[rowno] == TRUE) {
      df$pathlength[rowno] <- NA
      df$pv_angle[rowno] <- NA
    }
  }
  
  tdf <- tbl_df(df) 
  
  ################################
  ################################
  ################################
  
  # ANALYZE LEARNING
  # note: Errors are signed since groups are analyzed separately
  
  baselinelast <- tdf %>%
    filter(task==0) %>%
    filter(trial %in% c(57,58,59)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle, na.rm=TRUE),
              pl = mean(pathlength, na.rm=TRUE),
              block = mean(task))
  
  block1 <- tdf %>%
    filter(task == 3) %>%
    filter(trial %in% c(0,1,2)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              pl = mean(pathlength, na.rm = TRUE),
              block = mean(task))
  
  blocklast <- tdf %>%
    filter(task == 7) %>%
    filter(trial %in% c(10,11,12)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle, na.rm=TRUE),
              pl = mean(pathlength, na.rm=TRUE),
              block = mean(task))
  
  adaptdf <- rbind(block1, blocklast)
  
  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  RM_pv <- aov(pv ~ block + Error(participant/block), data=adaptdf)
  summary(RM_pv)

  
  # VISUALIZE LEARNING

  tdfsmooth <- tdf %>% filter(task == 3)
  tdfsmooth2 <- tdf %>% filter(task == 5)
  tdfsmooth2$trial <- tdfsmooth2$trial + 179
  tdfsmooth <- rbind(tdfsmooth,tdfsmooth2)

  traindf <- tdfsmooth %>%
    group_by(participant) %>%
    group_by(trial) %>%
    summarise(Mean_pl = mean(pathlength, na.rm=TRUE),
              SD_pl = sd(pathlength, na.rm=TRUE),
              SEM_pl = SD_pl/sqrt(length(unique(participant))),
              Mean_pv = mean(pv_angle, na.rm=TRUE),
              SD_pv = sd(pv_angle, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))))

 ##plotting blocked training
  for (rot in sort(unique(tdf$garage_location))){
    
    dfplot1 <- tdf %>% filter(garage_location==rot) %>% filter(task == 3)
    dfplot2 <- tdf %>% filter(garage_location==rot) %>% filter(task == 5)
    dfplot3 <- tdf %>% filter(garage_location==rot) %>% filter(task == 7)
    dfplot2$trial <- dfplot2$trial + 179
    dfplot3$trial <- dfplot3$trial + 179 + 179
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
      ggtitle('CW ONLY') 
    
    print(bltrain)
    move_layers(bltrain, "GeomRibbon", position = "top")
    move_layers(bltrain, "GeomPoint", position = "top")
    
    traindf <- traindf %>% # NEED THIS FOR OMNIBUS TESTS IN DUAL DATA ANALYSIS LATER
      mutate(rotation = rot)
    #outfile_name <- sprintf('SINGLE_LCs_%s.csv', rot)
    #write.csv(traindf, file = outfile_name, row.names = FALSE)  
  }

  
  ################################
  ################################
  ################################
  
  # ANALYZE REACH AFTEREFFECTS
  # baseline = -1; excludeAE = 0; includeAE = 1
  # note: Errors are signed since groups are analyzed separately
  
  baselineAE <- tdf %>%
    filter(task == 1) %>%
    filter(trial %in% c(10,11,12)) %>%
    group_by(participant) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE), 
              garage_location = unique(garage_location),
              instruction = 'baseline')
  
  excludeAE <- tdf %>%
    filter(instruction == 'exclude') %>%
    drop_na(pv_angle) %>% 
    group_by(participant) %>% 
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              garage_location = unique(garage_location),
              instruction = 'exclude')
  
  includeAE <- tdf %>%
    filter(instruction == 'include') %>%
    drop_na(pv_angle) %>% 
    group_by(participant) %>%
    filter(trial == min(trial)) %>%
    summarise(pv = mean(pv_angle, na.rm = TRUE),
              garage_location = unique(garage_location),
              instruction = 'include')
  
  #subtract baseline
  excludeAE$pv <- excludeAE$pv - baselineAE$pv
  includeAE$pv <- includeAE$pv - baselineAE$pv
  
  t.test(excludeAE$pv, mu=0, alternative ="greater") # is there implicit learning?
  t.test(includeAE$pv, excludeAE$pv, alternative ="greater", paired = TRUE) # is there explicit learning?
  
  swarms <- rbind(excludeAE,
                  includeAE)
  
  ## VISUALIZE RAW I/E NO-CURSORS
  tdf_NCs_rot <- tdf %>%
    filter(instruction == 'exclude' | instruction == 'include') %>%
    filter(trial == 0)
  
  ggplot(tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(garage_location))) +
    geom_boxplot() +
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("Single rotation group")
  #boxplot(includeAE$pv) #CW - A COUPLE OF OUTLIERS
  
  ## COLLECT BASELINE-SUBTRACTED REACH AEs
  includeAE.summary <- includeAE %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "include",
              garage_location = "1",
              lowerSEM = Mean_pv-SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  excludeAE.summary <- excludeAE %>%
    summarise(Mean_pv = mean(pv, na.rm=TRUE),
              SD_pv = sd(pv, na.rm=TRUE),
              SEM_pv = SD_pv/sqrt(length(unique(participant))),
              instruction = "exclude",
              garage_location = "1",
              lowerSEM = Mean_pv - SEM_pv,
              upperSEM = Mean_pv + SEM_pv)
  
  SEMs2 <- rbind(includeAE.summary, excludeAE.summary)

  ## bar plot I/E Reach aftereffects ##
  IEbars<- ggplot(data = SEMs2, aes(x = instruction, y = Mean_pv,
                                  fill = as.factor(garage_location))) +
            geom_bar(stat = "identity", position = "dodge") +
            geom_errorbar(data = SEMs2, mapping = aes(x = instruction, y = Mean_pv,
                                                  ymin = SEMs2$lowerSEM, ymax = SEMs2$upperSEM),
                                                  width = 0.2, size = 0.5, color = "black",
                                                  position = position_dodge(width = 0.9)) +
            geom_beeswarm(data = swarms, aes(x = instruction, y = pv),
                  alpha = 1/7,
                  dodge.width = .9, cex = 3,
                  stroke = 0.3) +
            # geom_point(data = includeAE, aes(x = instruction, y = pv),
            #            size = 1, stroke = 0,
            #            alpha = 1/20) +
            # geom_point(data = excludeAE, aes(x = instruction, y = pv),
            #            size = 1, stroke = 0,
            #            alpha = 1/20) +
            ylab("Angular error (Degrees)") +
            ggtitle("Single CW") +
            coord_fixed(ratio = 1/13) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  legend.title = element_blank(), legend.position = "none") +
            scale_y_continuous(breaks = seq(-30, +30, 10), limits = c(-30, 30))
  print(IEbars)
  
}
