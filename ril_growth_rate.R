# Analysis of rate data
library(ggplot2)
# Read in data
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
BP14<-read.csv("ril_best.fit_logistic_estimates_height.csv")

dap<-unique(BP14$dap_i) 

# Get data for AGR into QTL format
BP14_qtl<-BP14[,c(1:2)]
for(d in 1:length(dap)) {
  treat.day<-BP14[BP14$dap_i == dap[d],c(1,2,7)]
  t.names<-colnames(treat.day[3])
  #t.names[1]<-c("height")
  t.names<-paste(t.names, dap[d], sep="_")
  colnames(treat.day)[3]<-t.names
  BP14_qtl<-merge(BP14_qtl, treat.day, by=c('ril', 'treatment'), all=T)
}

# Remove duplicate rows
BP14_qtl<-unique(BP14_qtl)
BP14_qtl<-BP14_qtl[,c(1:28)]
#colnames(BP14_qtl)[28]<-c('AGR_33')

BP14_qtl$Obs<-c(1:nrow(BP14_qtl))
BP14_qtl$experiment<-rep('DP14', nrow(BP14_qtl))
BP14_qtl$year<-rep('2014', nrow(BP14_qtl))
BP14_qtl$plot<-rep('unknown', nrow(BP14_qtl))
BP14_qtl$plot_id<-rep('unknown', nrow(BP14_qtl))
BP14_qtl$measurement<-rep('phenotyper', nrow(BP14_qtl))

# reorder cols
BP14_AGR_qtl<-BP14_qtl[,c(29:31,2,32,33,1,34,3:28)]
colnames(BP14_AGR_qtl)[7]<-c("id")
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
write.csv(BP14_AGR_qtl, file="BP14_AGR_qtl.csv", quote=F, row.names=F)

########################
# Now lets get maximum AGR and Day at which growth rate is maximized in each treatment
########################


dap<-unique(BP14$dap_i)
rils<-unique(BP14$ril)
treatments<-unique(BP14$treatment)


BP14_max.rate_dap<-c()
for(r in 1:length(rils)) {
  g<-rils[r]
  temp1<-BP14[BP14$ril == g,]
  for(t in 1:length(treatments)) {
    temp2<-temp1[temp1$treatment == treatments[t],]
    if(length(temp2$dap_i > 0)) {
      #mean.AGR<-mean(temp2$AGR)
      max.AGR<-max(temp2$AGR)
      max.AGR.dap<-temp2[temp2$AGR == max.AGR, 'dap_i']
      out<-c(as.character(g),as.character(treatments[t]), max.AGR,max.AGR.dap)
      out<-t(as.data.frame(out))
      BP14_max.rate_dap<-rbind(BP14_max.rate_dap, out)
    } else (next)
  }
}

BP14_max.rate_dap<-as.data.frame(BP14_max.rate_dap)
rownames(BP14_max.rate_dap)<-c(1:nrow(BP14_max.rate_dap))
colnames(BP14_max.rate_dap)<-c('id', 'treatment', 'maxAGR', 'day_of_maxAGR')



BP14_max.rate_dap$Obs<-c(1:nrow(BP14_max.rate_dap))
BP14_max.rate_dap$experiment<-rep('BP14', nrow(BP14_max.rate_dap))
BP14_max.rate_dap$year<-rep('2014', nrow(BP14_max.rate_dap))
BP14_max.rate_dap$plot<-rep('unknown', nrow(BP14_max.rate_dap))
BP14_max.rate_dap$plot_id<-rep('unknown', nrow(BP14_max.rate_dap))
BP14_max.rate_dap$measurement<-rep('phenotyper', nrow(BP14_max.rate_dap))

BP14_max.rate_dap<-BP14_max.rate_dap[,c(5,6,7,2,8,9,1,10,3:4)]
rownames(BP14_max.rate_dap)<-c(1:nrow(BP14_max.rate_dap))
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
write.csv(BP14_max.rate_dap, file='BP14_max.rate_dap_qtl.csv', row.names=F, quote=F)


# Lets make some plots of average rate, max rate, and day of max rate

BP14_rate_dap<-c()
for(r in 1:length(rils)) {
  g<-rils[r]
  temp1<-BP14[BP14$ril == g,]
  for(t in 1:length(treatments)) {
    temp2<-temp1[temp1$treatment == treatments[t],]
    if(length(temp2$dap_i > 0)) {
      mean.AGR<-mean(temp2$AGR)
      max.AGR<-max(temp2$AGR)
      max.AGR.dap<-temp2[temp2$AGR == max.AGR, 'dap_i']
      out<-c(as.character(g),as.character(treatments[t]), mean.AGR, max.AGR,max.AGR.dap)
      out<-t(as.data.frame(out))
      BP14_rate_dap<-rbind(BP14_rate_dap, out)
    } else (next)
  }
}

BP14_rate_dap<-as.data.frame(BP14_rate_dap)
rownames(BP14_rate_dap)<-c(1:nrow(BP14_rate_dap))
colnames(BP14_rate_dap)<-c('id', 'treatment','meanAGR', 'maxAGR', 'day_of_maxAGR')

BP14_rate_dap$meanAGR<-as.numeric(as.character(BP14_rate_dap$meanAGR))
BP14_rate_dap$maxAGR<-as.numeric(as.character(BP14_rate_dap$maxAGR))
BP14_rate_dap$day_of_maxAGR<-as.numeric(as.character(BP14_rate_dap$day_of_maxAGR)) 

# Make the plots
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")

# Average growth rate per day
pdf("Figure_8a.pdf")
p<-ggplot(BP14, aes(dap_i, AGR, col=treatment)) + geom_smooth(size=1.5) +  scale_colour_manual(name="Treatment", labels=c('Dry', 'Wet'), values=c('orange', 'blue')) + ylab("Growth rate (cm/day)") + xlab("Days after planting") + theme_bw()
print(p)
dev.off()

# Max growth rate boxplot
pdf("Figure_8b.pdf")
p<-ggplot(BP14_rate_dap, aes(treatment, maxAGR, fill=factor(treatment))) + geom_boxplot() + scale_fill_manual(name="Treatment", labels=c('Dry', 'Wet'), values=c('orange', 'blue')) + theme_bw() +xlab("Treatment")+ylab("Maximum growth rate (cm/day)")
print(p)
dev.off()

# Day of maximum growth rate
pdf("Figure_8c.pdf")
p<-ggplot(BP14_rate_dap, aes(day_of_maxAGR, fill=treatment)) + geom_density(alpha=0.2) + scale_fill_manual(name="Treatment", labels=c('Dry', 'Wet'), values=c('orange', 'blue')) + theme_bw() +xlab("Day of maximum growth rate") + ylab("Density")
print(p)
dev.off()

