## plotting PIs across exps
library(ggplot2)
meanPIs <- c(CCWmeanPI, CWmeanPI,seqmeanPI,conseqmeanPI,preqmeanPI,expmeanPI,statmeanPI)
exps <- c('CCW', 'CWPI','seqPI','conseqPI','preqPI','expPI','statPI')
semPIs <- c(CCWsemPI, CWsemPI,seqsemPI,conseqsemPI,preqsemPI,expsemPI,statsemPI)

PI <- data.frame(exps,semPIs,meanPIs)
meanPIs <- c(CCWmeanPI, CWmeanPI,seqmeanPI,conseqmeanPI,preqmeanPI,expmeanPI,statmeanPI)
g<- ggplot(data=PI, aes(x=exps, y=meanPIs)) +
  geom_bar(stat="identity") +
  geom_errorbar(data=PI, mapping=aes(x=exps, y=meanPIs, ymin=meanPIs-semPIs, ymax=meanPIs+semPIs),
  width=0.1, size=1, color="grey", position = position_dodge(width = 0.9)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(-100,100)