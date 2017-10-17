#### nick's ROC demo code
library(dplyr)
library(ggplot2)
library(pROC)

set.seed(1)
example <- data.frame(type = c(rep('Best',100), rep('Expected',100), rep('Random',100)),
                      Severity = factor(rep(c(rep('Mild',50),rep('Severe',50)),3)),
                      votes = c(c(sample(c(0:40), 50, replace = T)/100,
                                  sample(c(60:100), 50, replace = T)/100),
                                c(sample(c(0:65), 50, replace = T)/100,
                                  sample(c(35:100), 50, replace = T)/100),
                                c(sample(c(0:100), 50, replace = T)/100,
                                  sample(c(0:100), 50, replace = T)/100)))

ggplot(example, aes(color = Severity, y = votes, x=1)) + 
  geom_jitter(width = 0.5, size = 3) + theme_bw() + 
  labs(x=NULL, y='Voted Classification Value') + 
  facet_grid(.~type) + xlim(0,2) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=20))

roc_best <- roc(example[example$type=='Best','Severity'], 
                example[example$type=='Best','votes'])
roc_exp <- roc(example[example$type=='Expected','Severity'], 
               example[example$type=='Expected','votes'])
roc_rand <- roc(example[example$type=='Random','Severity'], 
                example[example$type=='Random','votes'])

par(mar=c(1,4,4,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), 
     xaxs='i', yaxs='i', ylab='', xlab='', axes = F)
axis(2)
axis(3)
box()
plot(roc_best, col='purple4', lwd=3, add=T, lty=1)
plot(roc_exp, col='green3', lwd=3, add=T, lty=1)
plot(roc_rand, col='red', lwd=3, add=T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=3, text="Specificity", line=2.5)
legend('bottomright', legend=c(sprintf('Best'), #(AUC: %.3g)',roc_best$auc),
                               sprintf('Expected'), #(AUC: %.3g)',roc_exp$auc),
                               sprintf('Random')), # (AUC: %.3g)',roc_rand$auc)),
       lty=c(1,1,1), lwd=3, col=c('purple4','green3','red'), bty='n')

ggplot(example, aes(color = Severity, y = votes, x=1)) + 
  geom_jitter(width = 0.5) + theme_bw() + labs(x=NULL, y='Voted Classification Value') + 
  facet_grid(.~type) + xlim(0,2) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

roc_best <- roc(example[example$type=='Best','Severity'], 
                example[example$type=='Best','votes'])
roc_exp <- roc(example[example$type=='Expected','Severity'], 
               example[example$type=='Expected','votes'])
roc_rand <- roc(example[example$type=='Random','Severity'], 
                example[example$type=='Random','votes'])

par(mar=c(1,4,4,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), 
     xaxs='i', yaxs='i', ylab='', xlab='', axes = F)
axis(2);axis(3);box()
plot(roc_exp, col='green3', lwd=3, add=T, lty=1)
plot(roc_best, col = 'black', lwd = 3, add = T, lty=2)
plot(roc_rand, col = 'red', lwd=3, add = T, lty=3)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=3, text="Specificity", line=2.5)

example %>% 
  filter(type == 'Expected') %>% 
  ggplot(aes(color = Severity, y = votes, x=1)) + xlim(0.5,1.5) + 
  geom_jitter(width = 0.5, size=3) + theme_bw() + 
  labs(x=NULL, y='Voted Classification Value') + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=20))