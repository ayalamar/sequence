# for single CW and CCW
setwd('/Users/mayala/Desktop/single CCW data')

subject_numbers <- c(1:10) # same for CW & CCW
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()

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
    }
  }
  
  tdf <- tbl_df(df) # convert to tibble for dplyr
  
  # ANALYZE LEARNING
  # note: Errors are signed since groups are analyzed separately
  baselinelast <- tdf %>% filter(task==0) %>% filter(trial %in% c(57,58,59)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  block1<- tdf %>% filter(task==3) %>% filter(trial %in% c(0,1,2)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  blocklast<- tdf %>% filter(task==7) %>% filter(trial %in% c(10,11,12)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  
  adaptdf<- rbind(block1,blocklast)
  
  adaptdf$block <- factor(adaptdf$block)
  adaptdf$participant <- factor(adaptdf$participant)
  
  RM_pv <- aov(pv ~ block + Error(participant/block), data=adaptdf)
  summary(RM_pv)
  RM_pl <- aov(pl ~ block + Error(participant/block), data=adaptdf)
  summary(RM_pl)
  
  # VISUALIZE LEARNING

  tdfsmooth <- tdf %>% filter(task == 3)
  tdfsmooth2 <- tdf %>% filter(task == 5)
  tdfsmooth2$trial <- tdfsmooth2$trial + 179
  tdfsmooth <- rbind(tdfsmooth,tdfsmooth2)

  traindf <- tdfsmooth %>% group_by(participant) %>% group_by(trial) %>% summarise(Mean_pl = mean(pathlength, na.rm=TRUE),SD_pl = sd(pathlength, na.rm=TRUE),
                                  SEM_pl = SD_pl/sqrt(length(unique(participant))),
                                  Mean_pv = mean(pv_angle, na.rm=TRUE), SD_pv = sd(pv_angle, na.rm=TRUE),
                                  SEM_pv = SD_pv/sqrt(length(unique(participant))))

# #   exp.model <-lm(pv_angle ~ exp(trial), tdfsmooth)
# #   decay <- lm(log(pv_angle) ~ trial, data=tdfsmooth) 
#   
#   Means <- ggplot(data=traindf, aes(x=trial, y=Mean_pv)) +
#     geom_line() + 
#     geom_ribbon(data=traindf, aes(ymin=Mean_pv-SEM_pv, ymax= Mean_pv+SEM_pv), alpha=0.4) +
#     geom_point(data=tdfsmooth, aes(x=trial, y= pv_angle), colour='chartreuse1',alpha = 0.1) +
#     #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, linetype="dashed", color="blue", aes(outfit=fit<<-..y..),n=359) +
#     #geom_line(aes(x=trial, y=exp(decay$fitted.values)), color = "red")
#     #stat_smooth(method = "nls", formula = y ~ a * exp(-S * x), 
#                 # method.args = list(start = list(a = 78, S = 0.02)), se = FALSE, #starting values obtained from fit above
#                 # color = "dark red", linetype ="dashed")+
#     ylim(-50, 50) +
#     ggtitle('Single CW training') +
#     ylab("Angular error (Degrees)") +
#     xlab("Trial") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#   print(Means)
#   move_layers(Means,"GeomRibbon", position = "top")
#   move_layers(Means,"GeomLine", position = "top")
#   move_layers(Means,"GeomPoint", position = "bottom")
  
  ##plotting blocked training
  for (rot in sort(unique(tdf$garage_location))){
    
    dfplot1 <- tdf %>% filter(garage_location==rot) %>% filter(task == 3)
    dfplot2 <- tdf %>% filter(garage_location==rot) %>% filter(task == 5)
    dfplot3 <- tdf %>% filter(garage_location==rot) %>% filter(task == 7)
    dfplot2$trial <- dfplot2$trial + 179
    dfplot3$trial <- dfplot3$trial + 179 + 179
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
  
  #smooth_vals <- predict(loess(pv_angle~trial,tdfsmooth), tdfsmooth$trial)

  # ANALYZE REACH AFTEREFFECTS
  # get df of just reach AEs first. baseline = -1; excludeAE = 0; includeAE = 1
  # note: Errors are signed since groups are analyzed separately
  baselineAE <- tdf %>% filter(task==1) %>% filter(trial %in% c(10,11,12)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle, na.rm=TRUE), block = -1)
  excludeAE <- tdf %>% filter(instruction=='exclude') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle, na.rm=TRUE), block = 0)
  includeAE <- tdf %>% filter(instruction=='include') %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle, na.rm=TRUE), block = 1)
  
  AEdf <- rbind(baselineAE, excludeAE, includeAE)
  AEdf$block <- factor(AEdf$block)
  AEdf$participant <- factor(AEdf$participant)
  AEdf <- AEdf %>% mutate(group_instruction=1) # add group label for later 
  
  t.test(excludeAE$pv, baselineAE$pv, alternative ="greater") # is there implicit learning?
  t.test(excludeAE$pv - baselineAE$pv, mu=0, alternative ="greater") 
  t.test(includeAE$pv - baselineAE$pv, excludeAE$pv - baselineAE$pv, alternative ="greater") 
  
  # VISUALIZE REACH AFTEREFFECTS (to see if in the expected directions)
  tdf_NCs_rot <- tdf %>% filter(instruction == 'exclude' | instruction == 'include') %>% filter(trial == 0)
  ggplot(tdf_NCs_rot, aes(instruction, pv_angle, colour = factor(garage_location))) +
    geom_boxplot() +
    ylim(-50, 50) +
    theme_classic() +
    ggtitle("Single rotation group")
  
}