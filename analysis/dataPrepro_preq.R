########## FOR PREQ DATA SPECIFICALLY ##########
setwd('/Users/mayala/Desktop/preq data')
subject_numbers <- c(1:2) 
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()

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
      
      colnames(taskdf) <- c("cursorx","cursory","homex","homey","handx","handy","participant",
                            "prehome", "prehomex", "prehomey", "rotation", "step", "target_angle",
                            "targetx", "targety", "time", "trial", "selection_1","selection_2",
                            "selection_3","selection_peakvel","selection_unsure")
      
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
        
        ## PARTICIPANTS BEFORE SUBJECT # 19 NEED TO HAVE CB PRINTED IN OUTPUT
        if (unique(taskdf$participant) < 18.5) { 
          
          if (unique(taskdf$participant) < 9.5) {
            
            orderIndex <- unique(taskdf$participant) %% 4
            instrorder <- (list(c(0,1), c(1,0), c(0,1), c(1,0)))[orderIndex + 1] #python has it 0-indexed
            
            if (instrorder[[1]][1] == 0) {
              
              taskdf$instruction[which(taskdf$task == 8)] <- 'include'
              taskdf$instruction[which(taskdf$task == 9)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 11)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 12)] <- 'include'
              
            } else {
              
              taskdf$instruction[which(taskdf$task == 8)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 9)] <- 'include'
              taskdf$instruction[which(taskdf$task == 11)] <- 'include'
              taskdf$instruction[which(taskdf$task == 12)] <- 'exclude'
              
            }
            
          } else {
            
            orderIndex <- unique(taskdf$participant) %% 8
            instrorder <- (list(c(0,1,0,1), c(1,0,1,0), c(0,1,0,1), c(1,0,1,0),
                                c(0,1,1,0), c(1,0,0,1), c(0,1,1,0), c(1,0,0,1)))[orderIndex + 1]
            instrorder1 <- list(c(0,1,0,1))
            instrorder2 <- list(c(1,0,1,0))
            instrorder3 <- list(c(0,1,1,0))
            instrorder4 <- list(c(1,0,0,1))
            
            if (all.equal(instrorder, instrorder1) == TRUE) {
              
              taskdf$instruction[which(taskdf$task == 8)] <- 'include'
              taskdf$instruction[which(taskdf$task == 9)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 11)] <- 'include'
              taskdf$instruction[which(taskdf$task == 12)] <- 'exclude'
              
            } else if (all.equal(instrorder, instrorder2) == TRUE) {
              
              taskdf$instruction[which(taskdf$task == 8)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 9)] <- 'include'
              taskdf$instruction[which(taskdf$task == 11)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 12)] <- 'include'
              
            } else if (all.equal(instrorder, instrorder3) == TRUE) {
              
              taskdf$instruction[which(taskdf$task == 8)] <- 'include'
              taskdf$instruction[which(taskdf$task == 9)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 11)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 12)] <- 'include'
              
            } else { # instrorder4
              
              taskdf$instruction[which(taskdf$task == 8)] <- 'exclude'
              taskdf$instruction[which(taskdf$task == 9)] <- 'include'
              taskdf$instruction[which(taskdf$task == 11)] <- 'include'
              taskdf$instruction[which(taskdf$task == 12)] <- 'exclude'
              
            }
            
          }
          
        } else { 
          print('this participant was run after fix')
          ## THIS PARTICIPANT WAS RUN AFTER THE EXPERIMENT WAS FIXED - NO ADJUSTMENTS NEEDED
          
          
          
        }
        
      }
      
  if (is.data.frame(ppdf)==TRUE) {
        
    ppdf <- rbind(ppdf, taskdf)
        
  } else {
        
    ppdf <- taskdf
        
  }
    
    outfile_name = sprintf('combined_p0%02d_%s.csv', subject_id, outfile_suffix)
    write.csv(ppdf, file = outfile_name, row.names = FALSE)
    
  }
  
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
    ppdf$pv_angle_n[which(ppdf$selection_1 != 1)] <- NA
    
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
