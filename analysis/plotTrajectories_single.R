# for single CW and CCW - plot trajectories
setwd('~/science/repos/sequence')
setwd('single CCW data')

subject_numbers <- c(1:10) # same for CW & CCW
tasks <- c(0, 1, 3, 4, 5, 6, 7) 
outfile_suffix <- sprintf('ALL')
homex <- c(0)
homey <- c(7.3438) # if this is 0, no scaling will be done in taskAnalysis()
npoints <- 20 # how many points to interpolate for reach trajectory plots

library(dplyr)
library(ggplot2)
library(tidyverse)
library(gginnards)
  
periods <- c('early', 'late')
allpaths <- NA
  
  for (periodno in periods){
    
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
   
          if (periodno == 'early'){
            
            ##FIRST BLOCKS - EARLY ADAPTATION
            pathdf1 <- taskdf %>%
              filter(garage_location == unique(taskdf$garage_location)) %>%
              filter(task == taskno) %>%
              filter(trial == unique(trial)[1])
            pathdf2 <- taskdf %>%
              filter(garage_location == unique(taskdf$garage_location)) %>%
              filter(task == taskno) %>%
              filter(trial == unique(trial)[2])
            pathdf3 <- taskdf %>%
              filter(garage_location == unique(taskdf$garage_location)) %>%
              filter(task == taskno) %>%
              filter(trial == unique(trial)[3])
            
          } else {
            
            #LAST BLOCKS - LATE ADAPTATION
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
          }
          
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
        
        }
      
        # combine into one df across all pp
        if (is.data.frame(allpaths) == TRUE){
          allpaths <- rbind(allpaths, avgpaths)
        } else {
          allpaths <- avgpaths 
        }
    }

    for (tasknum in tasks){
      
      dev.off()
      
      taskpaths <- allpaths %>%
        filter(task == tasknum) %>%
        group_by(npoint) %>%
        summarise(meansampX = mean(meanX, na.rm = TRUE),
                  meansampY = mean(meanY, na.rm = TRUE))
      
      indivpaths <- allpaths %>%
        filter(task == tasknum)
      
      targetloc <- data.frame(
        x0 = rep(12, 1),
        y0 = rep(0, 1))
      
      homeloc <- data.frame(
        xh = rep(0, 1),
        yh = rep(0, 1))
      
      pathplot <- ggplot(data = taskpaths, 
             aes(x = meansampX, y = meansampY)) +
          geom_line(size = 1) +
          geom_point(data = indivpaths,
                    aes(x = meanX, y = meanY, colour = as.factor(participant)),
                    alpha = 0.3, size = 2, stroke = 0) +
          # geom_line(data = indivpaths,
          #           aes(x = meanX, y = meanY,
          #               colour = as.factor(participant)),
          #               linetype = "dotted", size = 1) +
          geom_point(data = targetloc,
                   aes(x = x0, y= y0), size = 7, shape = 21, colour = "black", fill = "black", stroke = 0, alpha = 0.3) +
          geom_point(data = homeloc,
                   aes(x = xh, y= yh), size = 7, shape = 21, colour = "black", fill = "black", stroke = 0, alpha = 0.3) +
          xlim(-15, 15) +
          ylim(-15, 15) +
          ylab("Y (cm)") +
          xlab("X (cm)") +
        theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             legend.title = element_blank(),
             legend.position = "none",
             axis.line=element_blank(),
             axis.text.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks=element_blank(),
             axis.title.x=element_blank(),
             axis.title.y=element_blank())
      
      
      move_layers(pathplot,"GeomLine", position = "top")
      print(pathplot)
      
      plotfilename <- sprintf('_paths/paths_task%d_%sadaptation.svg', tasknum, periodno)
      ggsave(file=plotfilename, plot=pathplot, width=9, height=9, units = "cm",dpi=20)
      
    }
    
  }  

rotateTrajectory <- function(X, Y, angle){
  
  # create rotation matrix to rotate X,Y coords
  th <- (angle/180) * pi
  R <- t(matrix(data=c(cos(th),
                       sin(th), 
                       -sin(th),
                       cos(th)),
                nrow = 2, ncol =2))
  
  # put coords in a matrix 
  coords <- matrix(data = c(X,Y),
                   ncol = 2)
  
  # rotate coords
  Rcoords <- coords %*% R
  
  return(Rcoords)
}
