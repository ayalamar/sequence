# For prepping single CW and CCW data for analysis; run these functions in order
setwd('~/science/repos/sequence')
setwd('single CW data')
rm(list=ls())

subject_numbers <- c(1:10) # same for CW & CCW
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()

library(dplyr)

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
      taskdf$angdev <- NA
      
      ## calculate error proxies (path length & ang deviation) per trial and store it 
      for (trialno in sort(unique(taskdf$trial))) {
        
        # select only the one current trial for analysis
        trialsamples <- taskdf[ which(taskdf$trial == trialno), ]
        
        # calculate pathlength and store it in column "pathlength"
        trial_pathlength <- sum(sqrt(diff(trialsamples$cursorx)^2 + diff(trialsamples$cursory)^2))
        taskdf$pathlength[which(taskdf$trial == trialno)] <- trial_pathlength
        
        # calculate maximum angular deviation and store it in column "maxdev"
        if (homey != 0) {
          # calculate new target Y because (0,0) is not the origin/home
          trialtemp <- trialsamples
          trialtemp$newtargety <- trialtemp$targety + homey
          trialtemp$newtargetangle <- (atan2(trialtemp$newtargety, trialtemp$targetx))/(pi/180)
          # calculate new cursor Y because (0,0) is not origin/home
          trialtemp$relativecursory <- trialtemp$cursory + homey
          trialtemp$angdev_OG <- (atan2(trialtemp$relativecursory, trialtemp$cursorx))/(pi/180) 
          # calculate angular deviation relative to target
          trialtemp$angdev <- trialtemp$target_angle - trialtemp$angdev_OG 
          trialsamples$angdev <- trialtemp$angdev
        } else {
          # no need to adjust y, calculate maxdev relative to target location and store it
          trialtemp <- trialsamples
          trialtemp$angdev_OG <- (atan2(trialtemp$cursory, trialtemp$cursorx))/(pi/180)
          trialtemp$angdev <- trialtemp$target_angle - trialtemp$angdev_OG
          trialsamples$angdev <- trialtemp$angdev
        }
        
        # put back the angular deviations and ...
        taskdf$angdev[which(taskdf$trial == trialno)] <- trialsamples$angdev
        
        if (trialsamples$participant%%2 == 1) { # this is an ODD-numbered participant
            # print include or exclude for no-cursor tasks 8 & 9
            # if subject number is ODD, store exclude-include instruction sequence
            # else, store include-exclude instruction sequence
            taskdf$instruction[which(taskdf$task == 8)] <- 'exclude'
            taskdf$instruction[which(taskdf$task == 9)] <- 'include'
            
        } else { #this is an even-numbered participant
          
            taskdf$instruction[which(taskdf$task == 8)] <- 'include'
            taskdf$instruction[which(taskdf$task == 9)] <- 'exclude'
            
        }
        
      }
      
      if (is.data.frame(ppdf)==TRUE) {
        
        ppdf <- rbind(ppdf, taskdf)
        
      } else {
        
        ppdf <- taskdf
        
      }
    }
    
    outfile_name = sprintf('md_analysis/combined_p0%02d_%s.csv', subject_id, outfile_suffix)
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
    
    filename <- sprintf('md_analysis/combined_p0%02d_%s.csv', subject_id, outfile_suffix)
    print(filename)
    
    taskdf <- read.csv(filename, header = TRUE)
    taskdf <- taskdf %>%
      group_by(task,trial) %>%
      filter(step %in% c(2,3)) %>%
      mutate(maxdev = max(abs(angdev), na.rm = TRUE))
    
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
    
    outfile_name = sprintf('md_analysis/trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    
    write.csv(pvsamples, file = outfile_name, row.names = FALSE)  
  } 
  
}

# outlierRemove . R 
tagOutliers <- function() {
  
  for (ppno in 1:length(subject_numbers)) {
    
    subject_id <- subject_numbers[ppno]
    
    filename <- sprintf('md_analysis/trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    print(filename)
    
    ppdf <- read.csv(filename, header = TRUE)
    
    boxplot(pv_angle~rotation, ylab='angle at peak velocity', xlab='rotation',data = ppdf)
    outlier_values <- boxplot.stats(ppdf$pv_angle)$out
    
    ppdf$isoutlier <- FALSE
    ppdf$isoutlier[which(ppdf$pv_angle %in% outlier_values)] <- TRUE
    #ppdf$pv_angle[which(ppdf$isoutlier == TRUE)] <- NA # NOTE : USING SIGNED ERRORS
    ppdf$pv_angle[which(ppdf$selection_1 != 1)] <- NA
    
    outfile_name = sprintf('md_analysis/tagged_trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    
    write.csv(ppdf, file = outfile_name, row.names = FALSE) 
    
  }
}

# combine ALL ppdf with tagged outliers into one big file that includes ALL participants (use for plots)
combineTagged <- function() {
  groupdf <- NA
  
  for (ppno in 1:length(subject_numbers)) {
    
    subject_id <- subject_numbers[ppno]
    
    filename <- sprintf('md_analysis/tagged_trialdata_p0%02d_%s.csv', subject_id, outfile_suffix)
    print(filename)
    
    ppdf <- read.csv(filename, header = TRUE)
    
    if (is.data.frame(ppdf)==TRUE) {
      
      groupdf <- rbind(groupdf, ppdf)
      
    } else {
      
      groupdf <- ppdf
      
    }
    
  }
  
  outfile_name = sprintf('md_analysis/allTaggedData_n%d_%s.csv', length(subject_numbers), outfile_suffix)
  
  write.csv(groupdf, file = outfile_name, row.names = FALSE) 
  # note : this has a first row of NAs
}
