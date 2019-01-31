
######### MAKE CHANGES HERE ############
########################################
setwd('/Users/mayala/Desktop/DATS')

subject_numbers <- c(1,2)
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) ## THIS NEEDS TO CHANGE FOR EXPLICIT VERSION OF EXP (9 TASKS)
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()

######### END OF CHANGES  ##############
########################################

# TASK COMBINE SEQUENCE . M
# AFTER having been selected, this combines all tasks per participant into one huge
# participant file that includes ALL samples for every trial and block

taskCombine <- function() {
  for (ppno in 1:length(subject_numbers)) {
    
    subject_id <- subject_numbers[ppno]
  
    ppdf <- NA
  
    
    for (taskno in tasks) {
  
      task <- tasks[taskno]
      
      filename <- sprintf('p%03d-%d_selected.txt', subject_id, taskno)
      print(filename)
      
      taskdf <- read.table(filename, header = FALSE)
      
      colnames(taskdf) <- c("trial","time","cursorx","cursory","handx","handy","rotation",
                               "step","target_angle","targetx","targety","garage_location","garagex","garagey","homex","homey",
                               "group","trial_type","selection_1","selection_2","selection_3","selection_peakvel","selection_unsure")
  
      taskdf$participant <- subject_id
      taskdf$task <- taskno
      taskdf$pathlength <- NA
      taskdf$instruction <- NA
      
      ## calculate path length per trial and store it 
      for (trialno in sort(unique(taskdf$trial))) {

        # select only the one current trial for analysis
        trialsamples <- taskdf[ which(taskdf$trial == trialno), ]

        # create empty column to store pathlength values

        trial_pathlength <- sum(sqrt(diff(trialsamples$cursorx)^2 + diff(trialsamples$cursory)^2))
        taskdf$pathlength[which(taskdf$trial == trialno)] <- trial_pathlength
      
        # separate out labelling for implicit experiments and explicit experiment
        if (max(tasks)==12) { ## THIS IS AN IMPLICIT EXPERIMENT
          if (trialsamples$participant%%2 == 1) { # this is an ODD-numbered participant
            # print include or exclude for no-cursor blocks 8,9 11,12
            # if subject number is ODD, store exclude-include-include-exclude instruction sequence
            # else, store include-exclude-include-exclude instruction sequence
            taskdf$instruction[which(taskdf$task == 8)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 9)] <- 'include'
            taskdf$instruction[which(taskdf$task == 11)] <- 'include'
            taskdf$instruction[which(taskdf$task == 12)] <- 'exclude'
            
          } else { #this is an even-numbered participant
            taskdf$instruction[which(taskdf$task == 8)] <- 'include'
            taskdf$instruction[which(taskdf$task == 9)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 11)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 12)] <- 'include'
          }
          
        } else { ## THIS IS AN EXPLICIT EXPERIMENT
          if (trialsamples$participant%%2 == 1) { # this is an ODD-numbered participant
            # for block 8, print 'exclude'
            taskdf$instruction[which(taskdf$task == 4)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 5)] <- 'include'
            taskdf$instruction[which(taskdf$task == 7)] <- 'include'
            taskdf$instruction[which(taskdf$task == 8)] <- 'exclude'
            
          } else { #this is an even-numbered participant
            taskdf$instruction[which(taskdf$task == 4)] <- 'include'
            taskdf$instruction[which(taskdf$task == 5)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 7)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 8)] <- 'include'
          }
          
        }

      }
      
      if (is.data.frame(ppdf)==TRUE) {
        
        ppdf <- rbind(ppdf, taskdf)
      
      } else {
        
        ppdf <- taskdf
        
      }
  }
    
  outfile_name = sprintf('combined_p0%02d_%s.csv', subject_id, outfile_suffix)
  write.csv(ppdf, file = outfile_name, row.names = FALSE)  
  
  }
  
}

# QUICK ANALYSIS SEQUENCE . M
# After having been combined, this function subsets the variable(s) of interest and 
# creates an output w/ a single sample per trial

# load the csv files you just made above ^
taskPreprocess <- function() {
for (ppno in 1:length(subject_numbers)) {
  
  subject_id <- subject_numbers[ppno]

  filename <- sprintf('combined_p0%02d_%s.csv', subject_id, outfile_suffix)
  print(filename)
  
  taskdf <- read.csv(filename, header = TRUE)
  
  # create subset where you only select samples that occur at peak velocity
  pvsamples <- taskdf[ which(taskdf$selection_peakvel == 1), ]
  
  if (homey != 0) {
    
    # calculate new target Y because (0,0) is not the origin/home & store it
    pvsamples$new_targety <- pvsamples$targety + homey
    pvsamples$new_target_angle <- (atan2(pvsamples$new_targety, pvsamples$targetx))/(pi/180)
    
    # calculate new cursor Y because (0,0) is not origin/home & store it
    
    pvsamples$relative_cursory <- pvsamples$cursory + homey
    pvsamples$pv_angle_OG <- (atan2(pvsamples$relative_cursory, pvsamples$cursorx))/(pi/180)
    
    # calculate angle at peak velocity relative to target location & store it
    
    pvsamples$pv_angle <- pvsamples$target_angle - pvsamples$pv_angle_OG # positive angles to the right of target; negative to left
    
  } else {
    
    # no need to adjust Y, just calculate angle at peak velocity relative to target location & store it
    pvsamples$pv_angle_OG <- (atan2(pvsamples$cursory, pvsamples$cursorx))/(pi/180)
    pvsamples$pv_angle <- pvsamples$target_angle - pvsamples$pv_angle_OG # positive angles to the right of target; negative to left
    
  }

  outfile_name = sprintf('trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
  
  write.csv(pvsamples, file = outfile_name, row.names = FALSE)  
} 
  
}

# outlierRemove . R 
tagOutliers <- function() {
  
  for (ppno in 1:length(subject_numbers)) {
    
    subject_id <- subject_numbers[ppno]
  
    filename <- sprintf('trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    print(filename)
    
    ppdf <- read.csv(filename, header = TRUE)
    
    ppdf$pv_angle_n <- abs(ppdf$pv_angle)
    boxplot(pv_angle_n~rotation, ylab='angle at peak velocity', xlab='rotation',data = ppdf)
    outlier_values <- boxplot.stats(ppdf$pv_angle_n)$out
    
    ppdf$isoutlier <- FALSE
    ppdf$isoutlier[which(ppdf$pv_angle_n %in% outlier_values)] <- TRUE
    ppdf$pv_angle_n[which(ppdf$isoutlier == TRUE)] <- NA
    
    outfile_name = sprintf('tagged_trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    
    write.csv(ppdf, file = outfile_name, row.names = FALSE) 
    
  }
}

# combine ALL ppdf with tagged outliers into one big file that includes ALL participants (use for plots)
combineTagged <- function() {
  groupdf <- NA
  
  for (ppno in 1:length(subject_numbers)) {
    
    subject_id <- subject_numbers[ppno]
    
    filename <- sprintf('tagged_trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    print(filename)
    
    ppdf <- read.csv(filename, header = TRUE)
    
    if (is.data.frame(ppdf)==TRUE) {
      
      groupdf <- rbind(groupdf, ppdf)
      
    } else {
      
      groupdf <- ppdf
      
    }
    
  }
  
  outfile_name = sprintf('allTaggedData_n%d_%s.csv', length(subject_numbers), outfile_suffix)

  write.csv(groupdf, file = outfile_name, row.names = FALSE) 
  # note : this has a first row of NAs
}

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
  
  # analyze reach aftereffects
  baselineAE<- tdf %>% filter(task==1) %>% filter(trial %in% c(0)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  block1AE<- tdf %>% filter(task==) %>% filter(trial %in% c(357,358,359)) %>% group_by(participant) %>% summarise(pv = mean(pv_angle_n, na.rm=TRUE), pl = mean(pathlength, na.rm=TRUE), block = mean(task))
  adaptdf<- rbind(block1,blocklast)

  
  # analyze effect of instruction
  
  tdf %>% filter(instruction == include) ...
  
  
  
}

