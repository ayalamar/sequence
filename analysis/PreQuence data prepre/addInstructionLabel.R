
if (participant < 18.5) { # subjects without instruction counterbalance printed in output
  
  if (participant < 9.5) {
    
    orderIndex <- participant %% 4
    instrorder <- (list(c(0,1), c(1,0), c(0,1),c(1,0)))[orderIndex + 1] #python has it 0-indexed
      
      if (instrorder[[1]][1] == 0){
        
        df$instruction[which(taskdf$task == 8)] <- 'include'
        df$instruction[which(taskdf$task == 9)] <- 'exclude'
        df$instruction[which(taskdf$task == 11)] <- 'exclude'
        df$instruction[which(taskdf$task == 12)] <- 'include'
        
      } else {
        
        df$instruction[which(taskdf$task == 8)] <- 'exclude'
        df$instruction[which(taskdf$task == 9)] <- 'include'
        df$instruction[which(taskdf$task == 11)] <- 'include'
        df$instruction[which(taskdf$task == 12)] <- 'exclude'
        
      }
      
  } else {
    
    orderIndex <- participant %% 8
    instrorder <- (list(c(0,1,0,1), c(1,0,1,0), c(0,1,0,1),c(1,0,1,0),
                        c(0,1,1,0),c(1,0,0,1),c(0,1,1,0),c(1,0,0,1)))[orderIndex + 1]
    instrorder1 <- list(c(0,1,0,1))
    instrorder2 <- list(c(1,0,1,0))
    instrorder3 <- list(c(0,1,1,0))
    instrorder4 <- list(c(1,0,0,1))
    
    if (all.equal(instrorder, instrorder1) == TRUE) {
      
      df$instruction[which(taskdf$task == 8)] <- 'include'
      df$instruction[which(taskdf$task == 9)] <- 'exclude'
      df$instruction[which(taskdf$task == 11)] <- 'include'
      df$instruction[which(taskdf$task == 12)] <- 'exclude'
      
    } else if (all.equal(instrorder, instrorder2) == TRUE) {
      
      df$instruction[which(taskdf$task == 8)] <- 'exclude'
      df$instruction[which(taskdf$task == 9)] <- 'include'
      df$instruction[which(taskdf$task == 11)] <- 'exclude'
      df$instruction[which(taskdf$task == 12)] <- 'include'
      
    } else if (all.equal(instrorder, instrorder3) == TRUE) {
      
      df$instruction[which(taskdf$task == 8)] <- 'include'
      df$instruction[which(taskdf$task == 9)] <- 'exclude'
      df$instruction[which(taskdf$task == 11)] <- 'exclude'
      df$instruction[which(taskdf$task == 12)] <- 'include'
    
    } else { # instrorder4
      
      df$instruction[which(taskdf$task == 8)] <- 'exclude'
      df$instruction[which(taskdf$task == 9)] <- 'include'
      df$instruction[which(taskdf$task == 11)] <- 'include'
      df$instruction[which(taskdf$task == 12)] <- 'exclude'
      
    }
    
  }
  
}
