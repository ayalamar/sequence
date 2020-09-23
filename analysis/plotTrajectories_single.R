# for single CW and CCW - plot trajectories
setwd('/Users/mayala/Documents/single CW data')

subject_numbers <- c(1:10) # same for CW & CCW
tasks <- c(0, 1, 3, 4, 5, 6, 7, 8, 9) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()
npoints <- 50 # how many points to interpolate for reach trajectory plots
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gginnards)

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
        
        trial_pathlength <- sum(sqrt(diff(trialsamples$cursorx)^2 + diff(trialsamples$cursory)^2))
        taskdf$pathlength[which(taskdf$trial == trialno)] <- trial_pathlength

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
    
    outfile_name = sprintf('combined_p0%02d_%s.csv', subject_id, outfile_suffix)
    write.csv(ppdf, file = outfile_name, row.names = FALSE)  
    
  }
  
}

# QUICK ANALYSIS SEQUENCE . M
# After having been combined, this function subsets the variable(s) of interest and 
# creates an output w/ a single sample per trial
# also gets splined trajectories

# load the csv files you just made above ^
taskPreprocess <- function() {
  
  allpaths <- NA
  
  for (ppno in 1:length(subject_numbers)) {
    
    avgpaths <- NA
    
    subject_id <- subject_numbers[ppno]
    
    filename <- sprintf('combined_p0%02d_%s.csv', subject_id, outfile_suffix)
    print(filename)
    
    taskdf <- read.csv(filename, header = TRUE)
    
    #get some excluded trials out
    taskdf <- taskdf %>% 
      filter(step %in% c(1,2)) %>% # get outward movement only
      filter(selection_1 == 1) %>% # no weird trials
      filter(selection_unsure == 0) %>% # no weird trials
      filter(pathlength < boxplot(taskdf$pathlength)$stats[4]) # get rid of very wonky reaches that weren't selected out for some reason
    
    ## calculate new target Y because (0,0) is not the origin/home & store it
    taskdf$new_targety <- taskdf$targety + homey
    taskdf$new_target_angle <- (atan2(taskdf$new_targety, taskdf$targetx))/(pi/180)
    taskdf$target_angle <- taskdf$new_target_angle
    ##calculate new cursor Y because (0,0) is not origin/home & store it
    taskdf$relative_cursory <- taskdf$cursory + homey
    taskdf$cursory <- taskdf$relative_cursory
    
    ## Get some exemplary hand paths for each task
    
      for (taskno in tasks) {
        print(tasks[taskno])
 
        # ####FIRST BLOCKS - EARLY ADAPTATION
        # pathdf1 <- taskdf %>%
        #   filter(garage_location == unique(garage_location)) %>%
        #   filter(task == taskno) %>%
        #   filter(trial == unique(trial)[1])
        # pathdf2 <- taskdf %>%
        #   filter(garage_location == unique(garage_location)) %>%
        #   filter(task == taskno) %>%
        #   filter(trial == unique(trial)[2])
        # pathdf3 <- taskdf %>%
        #   filter(garage_location == unique(garage_location)) %>%
        #   filter(task == taskno) %>%
        #   filter(trial == unique(trial)[3])

        ####LAST BLOCKS - LATE ADAPTATION
        pathdf1 <- taskdf %>%
          filter(garage_location == unique(garage_location)) %>%
          filter(task == taskno) %>%
          filter(trial == unique(trial)[length(unique(trial))]) # last trial
        pathdf2 <- taskdf %>%
          filter(garage_location == unique(garage_location)) %>%
          filter(task == taskno) %>%
          filter(trial == unique(trial)[length(unique(trial))-1]) # 2nd last trial
        pathdf3 <- taskdf %>%
          filter(garage_location == unique(garage_location)) %>%
          filter(task == taskno) %>%
          filter(trial == unique(trial)[length(unique(trial))-2]) # 3rd last trial
        
        
        ## spline interpolate new X and Y
        # path t1
        XX1 <- pathdf1$cursorx
        YY1 <- pathdf1$cursory
        TT1 <- pathdf1$time
        angle1 <- unique(pathdf1$target_angle)
      
        trajectory1 <- rotateTrajectory(XX1, YY1, (-1*angle1))
        XX1 <- trajectory1[,1]
        YY1 <- trajectory1[,2]
        
        newXX1 <- spline(TT1, XX1, n = npoints)$y
        newYY1 <- spline(TT1, YY1, n = npoints)$y
        
        newpathdf1 <- tibble(newXX1, newYY1)
        
        # path t2
        XX2 <- pathdf2$cursorx
        YY2 <- pathdf2$cursory
        TT2 <- pathdf2$time
        angle2 <- unique(pathdf2$target_angle)
        
        trajectory2 <- rotateTrajectory(XX2, YY2, (-1*angle2))
        XX2 <- trajectory2[,1]
        YY2 <- trajectory2[,2]
        
        newXX2 <- spline(TT2, XX2, n = npoints)$y
        newYY2 <- spline(TT2, YY2, n = npoints)$y
        
        newpathdf2 <- tibble(newXX2, newYY2)
        
        # path t3
        XX3 <- pathdf3$cursorx
        YY3 <- pathdf3$cursory
        TT3 <- pathdf3$time
        angle3 <- unique(pathdf3$target_angle)
        
        trajectory3 <- rotateTrajectory(XX3, YY3, (-1*angle3))
        XX3 <- trajectory3[,1]
        YY3 <- trajectory3[,2]
        
        newXX3 <- spline(TT3, XX3, n = npoints)$y
        newYY3 <- spline(TT3, YY3, n = npoints)$y
        
        newpathdf3 <- tibble(newXX3, newYY3)
        
        # get the average path across 3 trials -- block, and plot
        avgpathdfX <- tibble(newXX1, newXX2, newXX3)
        avgpathdfY <- tibble(newYY1, newYY2, newYY3)
        
        avgpathdfX <- avgpathdfX %>%
          mutate(meanX = rowMeans(avgpathdfX))
        
        avgpathdfY <- avgpathdfY %>%
          mutate(meanY = rowMeans(avgpathdfY))
        
        pppath <- tibble(avgpathdfX, avgpathdfY)
        
        pppath <- pppath %>%
          select(meanX, meanY) %>%
          mutate(task = taskno,
                 npoint = 1:n(),
                 participant = ppno)
        
        # combine into one df per pp
        if (is.data.frame(avgpaths) == TRUE){
          avgpaths <- rbind(avgpaths, pppath)
        } else {
          avgpaths <- pppath 
        }
        
        pathfilename <- sprintf('splinedPaths_%d.csv', ppno)
        #write.csv(avgpaths, file = pathfilename, row.names = FALSE)
      
      }
    
      # combine into one df across all pp
      if (is.data.frame(allpaths) == TRUE){
        allpaths <- rbind(allpaths, avgpaths)
      } else {
        allpaths <- avgpaths 
      }
    }

    # plot averaged paths
    for (tasknum in tasks){
      library(gginnards)
      taskpaths <- allpaths %>%
        filter(task == tasknum) %>%
        group_by(npoint) %>%
        summarise(meansampX = mean(meanX, na.rm = TRUE),
                  meansampY = mean(meanY, na.rm = TRUE))
      
      indivpaths <- allpaths %>%
        filter(task == tasknum)
      
      targetloc <- data.frame(
        x0 = rep(12, 20),
        y0 = rep(0, 20))
      
      homeloc <- data.frame(
        xh = rep(0, 20),
        yh = rep(0, 20))
      
      pathplot <- ggplot(data = taskpaths, 
             aes(x = meansampX, y = meansampY)) +
          geom_line(size = 1) +
          geom_point(data = indivpaths, 
                    aes(x = meanX, y = meanY, colour = as.factor(participant)),
                    alpha = 0.3) +
          geom_point(data = targetloc,
                   aes(x = x0, y= y0), size = 5, shape = 21, colour = "black", fill = "white", stroke = 2, alpha = 0.2) +
        geom_point(data = homeloc,
                   aes(x = xh, y= yh), size = 5, shape = 21, colour = "black", fill = "white", stroke = 2, alpha = 0.2) +
          xlim(-15, 15) +
          ylim(-15, 15) +
          ylab("Y (cm)") +
          xlab("X (cm)") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
                #legend.title = element_blank(), legend.position = "none") 
      
      print(pathplot)
      #move_layers(Means,"GeomLine", position = "top")
      outfile_suffix = "singleCW"
      plotfilename <- sprintf('paths_%s_task%d_lateadaptation.svg', outfile_suffix,tasknum)
      ggsave(file=plotfilename, plot=pathplot, width=20, height=20, units = "cm")
      
    }
    
    
    ####
    
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
    
    #write.csv(pvsamples, file = outfile_name, row.names = FALSE)  
  } 
  


