#### make runPreQuence output compatible with MATLAB ComputeData GUI
#### convert .csv experiment output to .txt 
#### remove the headers 
library(dplyr)
subject_numbers <- c(1:18) # edit participant #s

# output txt file with no header
for (subject in 1:length(subject_numbers)) {

  for (task in 0:12){
    
    infile_name <- sprintf('p0%02d-%d.csv', subject, task)
    data <- read.csv(infile_name, header = TRUE)

    ## RESTART TIME STAMP FOR EVERY TRIAL
    for (trialno in sort(unique(data$trial))){
      #sample1<- data[which(data$trial == trialno),] 
      data$time_ms[which(data$trial == trialno)] <- data$time_ms[which(data$trial == trialno)] - data$time_ms[which(data$trial == trialno)][1] # need to fix
    }
    
    ## SELECT ONLY STEPS 1 & 2 & 6 (outward movement)
    data <-  data %>% filter(step == c(1,2,6))
    
   ## only do for rotated trials & subjects run in the old version of the experiment that was
   ## missing the cursorxy stamps
      if (data$participant < 9.5){
        
        if (data$rotation[1] != 0 ) {
          
          ## convert mouse data to rotated cursor data (not initially recorded)
          for (sample in 1:nrow(data)){
            ref <- c(data$homex_px[1],data$homey_px[1])
            
            t <- data$rotation[sample]*pi/180 # convert to radians
            R <- c(cos(t), -1*sin(t), sin(t), cos(t)) # rotation matrix
            R <- matrix(data = R, nrow =2,ncol = 2, byrow = TRUE)
            data$cursorx_px[sample] <-  data$cursorx_px[sample] - ref[1]
            data$cursory_px[sample] <-  data$cursory_px[sample] - ref[2]
            
            newpoint <- R %*% c(data$cursorx_px[sample],data$cursory_px[sample])
            data$cursorx_px[sample] <- newpoint[1] + ref[1]
            data$cursory_px[sample] <- newpoint[2] + ref[2]
          }
          
        }
      }
    ## remove header for GUI compatibility
    outfile_name = sprintf('p0%02d-%d.txt', subject, task)
    write.table(data, file = outfile_name, sep = "\t",
                row.names = FALSE, col.names = FALSE)
  }
  
}