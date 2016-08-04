library(ggplot2)
library(lme4)
library(nlme)

########################################################################
# Plant height
########################################################################

# Load in illinois data
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")
DR13<-read.csv('DR13_qtl.csv')
DR14<-read.csv('DR14_qtl.csv')
DN13<-read.csv('DN13_qtl.csv')
DN14<-read.csv('DN14_qtl.csv')

# Lets start with DR13
# Get a list of all genotypes and treatments
genotypes<-sort(unique(DR13$id))
treatments<-sort(unique(DR13$treatment))

DR13.l<-c()
for(g in 1:length(genotypes)) {
  temp1<-DR13[DR13$id == genotypes[g],]
    # Lets extract all the timepoints from 
    day_list<-strsplit(colnames(DR13)[9:16], split="_")
    days<-c()
    for(d in 1:length(day_list)) {
      days<-c(days, day_list[[d]][2])
    }
    for(d in 1:length(days)) {
      colnames(temp1)[8]<-c("dap_i")
      temp2<-temp1[,c(1:8,8+d)]
      temp2$dap_i<-rep(days[d], nrow(temp1))
      colnames(temp2)[9]<-c('height')
      DR13.l<-rbind(DR13.l, temp2)
    }
}

days<-unique(DR13.l$dap_i)
DR13.blup<-c()
for(d in 1:length(days)){
  temp<-DR13.l[DR13.l$dap_i == days[d],]
  model<-lmer(height ~ (1|id) + treatment + plot %in% treatment + (1|id:treatment), data=temp)
  temp<-temp[complete.cases(temp),]
  temp$predicted<-predict(model)
  temp$residual<-residuals(model)
  DR13.blup<-rbind(DR13.blup, temp)
}

boxplot(DR13.blup$predicted~DR13.blup$treatment, col=c("orange", "blue"))
# Make plot of residuals
plot(DR13.blup$residual~DR13.blup$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DR13.blup$predicted~DR13.blup$height, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Measured')
abline(a=0, b=1, col='red')

# Lets take the average
DR13.l.ag<-aggregate(DR13.l[9], by=list(DR13.l$experiment, DR13.l$year, DR13.l$treatment, DR13.l$dap_i, DR13.l$id), mean, na.action="na.pass", na.rm=TRUE)
DR13.blup.ag<-aggregate(DR13.blup[,9:11], by=list(DR13.blup$experiment, DR13.blup$year, DR13.blup$treatment, DR13.blup$dap_i, DR13.blup$id), mean, na.action="na.pass", na.rm=TRUE)
plot(DR13.blup.ag$residual~DR13.blup.ag$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DR13.blup.ag$predicted~DR13.blup.ag$height, pch=20, cex=0.3, col="blue", ylab="Predicted", xlab="Measured")
abline(a=0, b=1, col='red')


colnames(DR13.blup.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'id', 'height', 'predicted', 'residual')
DR13_treatment_time.ag<-aggregate(DR13.blup.ag[,c(6:8)], by=list(DR13.blup.ag$experiment, DR13.blup.ag$year,DR13.blup.ag$treatment, DR13.blup.ag$dap_i), mean, na.action="na.pass", na.rm=TRUE)
colnames(DR13_treatment_time.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'height', 'predicted', 'residual')

# Plot raw  and BLUP average
plot(DR13_treatment_time.ag$height~DR13_treatment_time.ag$dap_i, col=c("blue", "orange"), pch=20, ylab="Height (cm)", xlab="Days after planting")
plot(DR13_treatment_time.ag$predicted~DR13_treatment_time.ag$dap_i, col=c("green", "brown"), pch=20, ylab="Height (cm)", xlab="Days after planting")

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# Lets make a ggplot version
p<-ggplot(DR13.blup, aes(x=dap_i, y=predicted, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height BLUP (cm)") + xlab("Days after planting")
pdf("DR13_height_blup_over_time.pdf")
print(p)
dev.off()

p<-ggplot(DR13.blup, aes(x=dap_i, y=height, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height (cm)") + xlab("Days after planting")
pdf("DR13_height_over_time.pdf")
print(p)
dev.off()


days<-sort(unique(DR13.blup$dap_i))
DR13.blup_qtl<-c()
for(d in 1:length(days)) {
  temp<-DR13.blup[DR13.blup$dap_i == days[d],]
  colnames(temp)[10]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DR13.blup_qtl<-temp[,-c(1,8,9,11)]
  }
  if (d > 1 ) {
    DR13.blup_qtl<-merge(DR13.blup_qtl, temp[,-c(1,8,9,11)], by=c('experiment', 'year', 'treatment', 'plot', 'subplot_id', 'id'), all=T)
  }
}

# Okay now lets finish this other crap...
DR13.blup_qtl$Obs<-c(1:nrow(DR13.blup_qtl))
DR13.blup_qtl$sampling<-rep('unknown', nrow(DR13.blup_qtl))
DR13.blup_qtl<-DR13.blup_qtl[,c(15,1:6,16,7:14)]

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# write out 
write.csv(DR13.blup_qtl, file="DR13.height.blup_qtl.csv", row.names=F, quote=F)


# Now DR14
# Get a list of all genotypes and treatments
genotypes<-sort(unique(DR14$id))
treatments<-sort(unique(DR14$treatment))

DR14.l<-c()
for(g in 1:length(genotypes)) {
  temp1<-DR14[DR14$id == genotypes[g],]
  # Lets extract all the timepoints from 
  day_list<-strsplit(colnames(DR14)[9:13], split="_")
  days<-c()
  for(d in 1:length(day_list)) {
    days<-c(days, day_list[[d]][2])
  }
  for(d in 1:length(days)) {
    colnames(temp1)[8]<-c("dap_i")
    temp2<-temp1[,c(1:8,8+d)]
    temp2$dap_i<-rep(days[d], nrow(temp1))
    colnames(temp2)[9]<-c('height')
    DR14.l<-rbind(DR14.l, temp2)
  }
}

days<-unique(DR14.l$dap_i)
DR14.blup<-c()
for(d in 1:length(days)){
  temp<-DR14.l[DR14.l$dap_i == days[d],]
  model<-lmer(height ~ (1|id) + treatment + plot %in% treatment + (1|id:treatment), data=temp)
  temp<-temp[complete.cases(temp),]
  temp$predicted<-predict(model)
  temp$residual<-residuals(model)
  DR14.blup<-rbind(DR14.blup, temp)
}

boxplot(DR14.blup$predicted~DR14.blup$treatment, col=c("orange", "blue"))
# Make plot of residuals
plot(DR14.blup$residual~DR14.blup$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
pdf("DR14_BLUP.v.actual.pdf")
plot(DR14.blup$predicted~DR14.blup$height, pch=20, cex=0.3, col='blue', xlab='Measured height (cm)', ylab='Predicted height (cm)')
abline(a=0, b=1, col='red')
dev.off()
# Lets take the average
DR14.l.ag<-aggregate(DR14.l[9], by=list(DR14.l$experiment, DR14.l$year, DR14.l$treatment, DR14.l$dap_i, DR14.l$id), mean, na.action="na.pass", na.rm=TRUE)
DR14.blup.ag<-aggregate(DR14.blup[9:11], by=list(DR14.blup$experiment, DR14.blup$year, DR14.blup$treatment, DR14.blup$dap_i, DR14.blup$id), mean, na.action="na.pass", na.rm=TRUE)
plot(DR14.blup.ag$residual~DR14.blup.ag$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DR14.blup.ag$predicted~DR14.blup.ag$height, pch=20, cex=0.3, col="blue", ylab="Predicted", xlab="Measured")
abline(a=0, b=1, col='red')

colnames(DR14.blup.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'id', 'height', 'predicted', 'residual')
DR14_treatment_time.ag<-aggregate(DR14.blup.ag[,c(6:8)], by=list(DR14.blup.ag$experiment, DR14.blup.ag$year,DR14.blup.ag$treatment, DR14.blup.ag$dap_i), mean, na.action="na.pass", na.rm=TRUE)
colnames(DR14_treatment_time.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'height', 'predicted', 'residual')

# Plot raw  and BLUP average
plot(DR14_treatment_time.ag$height~DR14_treatment_time.ag$dap_i, col=c("orange", "blue"), pch=20, ylab="Height (cm)", xlab="Days after planting")
plot(DR14_treatment_time.ag$predicted~DR14_treatment_time.ag$dap_i, col=c("orange", "blue"), pch=20, ylab="Height (cm)", xlab="Days after planting")

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# Lets make a ggplot version
p<-ggplot(DR14.blup, aes(x=dap_i, y=predicted, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height BLUP (cm)") + xlab("Days after planting")
pdf("DR14_height_blup_over_time.pdf")
print(p)
dev.off()

p<-ggplot(DR14.blup, aes(x=dap_i, y=height, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height (cm)") + xlab("Days after planting")
pdf("DR14_height_over_time.pdf")
print(p)
dev.off()

days<-sort(unique(DR14.blup$dap_i))
DR14.blup_qtl<-c()
for(d in 1:length(days)) {
  temp<-DR14.blup[DR14.blup$dap_i == days[d],]
  colnames(temp)[10]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DR14.blup_qtl<-temp[,-c(1,8,9,11)]
  }
  if (d > 1 ) {
    DR14.blup_qtl<-merge(DR14.blup_qtl, temp[,-c(1,8,9,11)], by=c('experiment', 'year', 'treatment', 'plot', 'subplot_id', 'id'), all=T)
  }
}

# Okay now lets finish this other crap...
DR14.blup_qtl$Obs<-c(1:nrow(DR14.blup_qtl))
DR14.blup_qtl$sampling<-rep('unknown', nrow(DR14.blup_qtl))
DR14.blup_qtl<-DR14.blup_qtl[,c(12,1:6,13,7:11)]

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# write out 
write.csv(DR14.blup_qtl, file="DR14.height.blup_qtl.csv", row.names=F, quote=F)


# Lets now go with DN13
# Get a list of all genotypes and treatments
genotypes<-sort(unique(DN13$id))
treatments<-sort(unique(DN13$treatment))

DN13.l<-c()
for(g in 1:length(genotypes)) {
  temp1<-DN13[DN13$id == genotypes[g],]
  # Lets extract all the timepoints from 
  day_list<-strsplit(colnames(DN13)[9:16], split="_")
  days<-c()
  for(d in 1:length(day_list)) {
    days<-c(days, day_list[[d]][2])
  }
  for(d in 1:length(days)) {
    colnames(temp1)[8]<-c("dap_i")
    temp2<-temp1[,c(1:8,8+d)]
    temp2$dap_i<-rep(days[d], nrow(temp1))
    colnames(temp2)[9]<-c('height')
    DN13.l<-rbind(DN13.l, temp2)
  }
}

days<-unique(DN13.l$dap_i)
DN13.blup<-c()
for(d in 1:length(days)){
  temp<-DN13.l[DN13.l$dap_i == days[d],]
  model<-lmer(height ~ (1|id) + treatment + plot %in% treatment + (1|id:treatment), data=temp)
  temp<-temp[complete.cases(temp),]
  temp$predicted<-predict(model)
  temp$residual<-residuals(model)
  DN13.blup<-rbind(DN13.blup, temp)
}

boxplot(DN13.blup$predicted~DN13.blup$treatment, col=c("orange", "blue"))
# Make plot of residuals
plot(DN13.blup$residual~DN13.blup$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DN13.blup$predicted~DN13.blup$height, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Measured')
abline(a=0, b=1, col='red')


# Lets take the average
DN13.l.ag<-aggregate(DN13.l[9], by=list(DN13.l$experiment, DN13.l$year, DN13.l$treatment, DN13.l$dap_i, DN13.l$id), mean, na.action="na.pass", na.rm=TRUE)
DN13.blup.ag<-aggregate(DN13.blup[,9:11], by=list(DN13.blup$experiment, DN13.blup$year, DN13.blup$treatment, DN13.blup$dap_i, DN13.blup$id), mean, na.action="na.pass", na.rm=TRUE)
plot(DN13.blup.ag$residual~DN13.blup.ag$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DN13.blup.ag$predicted~DN13.blup.ag$height, pch=20, cex=0.3, col="blue", ylab="Predicted", xlab="Measured")
abline(a=0, b=1, col='red')

colnames(DN13.blup.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'id', 'height', 'predicted', 'residual')
DN13_treatment_time.ag<-aggregate(DN13.blup.ag[,c(6:8)], by=list(DN13.blup.ag$experiment, DN13.blup.ag$year,DN13.blup.ag$treatment, DN13.blup.ag$dap_i), mean, na.action="na.pass", na.rm=TRUE)
colnames(DN13_treatment_time.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'height', 'predicted', 'residual')

# Plot raw average
plot(DN13_treatment_time.ag$height~DN13_treatment_time.ag$dap_i, col=c("orange", "blue"), pch=20, ylab="Height (cm)", xlab="Days after planting")
plot(DN13_treatment_time.ag$predicted~DN13_treatment_time.ag$dap_i, col=c("orange", "blue"), pch=20, ylab="Height (cm)", xlab="Days after planting")

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# Lets make a ggplot version
p<-ggplot(DN13.blup, aes(x=dap_i, y=predicted, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height BLUP (cm)") + xlab("Days after planting")
pdf("DN13_height_blup_over_time.pdf")
print(p)
dev.off()

p<-ggplot(DN13.blup, aes(x=dap_i, y=height, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height (cm)") + xlab("Days after planting")
pdf("DN13_height_over_time.pdf")
print(p)
dev.off()

days<-sort(unique(DN13.blup$dap_i))
DN13.blup_qtl<-c()
for(d in 1:length(days)) {
  temp<-DN13.blup[DN13.blup$dap_i == days[d],]
  colnames(temp)[10]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DN13.blup_qtl<-temp[,-c(1,8,9,11)]
  }
  if (d > 1 ) {
    DN13.blup_qtl<-merge(DN13.blup_qtl, temp[,-c(1,8,9,11)], by=c('experiment', 'year', 'treatment', 'plot', 'subplot_id', 'id'), all=T)
  }
}

# Okay now lets finish this other crap...
DN13.blup_qtl$Obs<-c(1:nrow(DN13.blup_qtl))
DN13.blup_qtl$sampling<-rep('unknown', nrow(DN13.blup_qtl))
DN13.blup_qtl<-DN13.blup_qtl[,c(15,1:6,16,7:14)]

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# write out 
write.csv(DN13.blup_qtl, file="DN13.height.blup_qtl.csv", row.names=F, quote=F)

# Lets start with DN14
# Get a list of all genotypes and treatments
genotypes<-sort(unique(DN14$id))
treatments<-sort(unique(DN14$treatment))

DN14.l<-c()
for(g in 1:length(genotypes)) {
  temp1<-DN14[DN14$id == genotypes[g],]
  # Lets extract all the timepoints from 
  day_list<-strsplit(colnames(DN14)[9:11], split="_")
  days<-c()
  for(d in 1:length(day_list)) {
    days<-c(days, day_list[[d]][2])
  }
  for(d in 1:length(days)) {
    colnames(temp1)[8]<-c("dap_i")
    temp2<-temp1[,c(1:8,8+d)]
    temp2$dap_i<-rep(days[d], nrow(temp1))
    colnames(temp2)[9]<-c('height')
    DN14.l<-rbind(DN14.l, temp2)
  }
}

days<-unique(DN14.l$dap_i)
DN14.blup<-c()
for(d in 1:length(days)){
  temp<-DN14.l[DN14.l$dap_i == days[d],]
  model<-lmer(height ~ (1|id) + treatment + plot %in% treatment + (1|id:treatment), data=temp)
  temp<-temp[complete.cases(temp),]
  temp$predicted<-predict(model)
  temp$residual<-residuals(model)
  DN14.blup<-rbind(DN14.blup, temp)
}

boxplot(DN14.blup$predicted~DN14.blup$treatment, col=c("orange", "blue"))
# Make plot of residuals
plot(DN14.blup$residual~DN14.blup$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DN14.blup$predicted~DN14.blup$height, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Measured')
abline(a=0, b=1, col='red')


# Lets take the average
DN14.l.ag<-aggregate(DN14.l[9], by=list(DN14.l$experiment, DN14.l$year, DN14.l$treatment, DN14.l$dap_i, DN14.l$id), mean, na.action="na.pass", na.rm=TRUE)
DN14.blup.ag<-aggregate(DN14.blup[10], by=list(DN14.blup$experiment, DN14.blup$year, DN14.blup$treatment, DN14.blup$dap_i, DN14.blup$id), mean, na.action="na.pass", na.rm=TRUE)

# Lets take the average
DN14.l.ag<-aggregate(DN14.l[9], by=list(DN14.l$experiment, DN14.l$year, DN14.l$treatment, DN14.l$dap_i, DN14.l$id), mean, na.action="na.pass", na.rm=TRUE)
DN14.blup.ag<-aggregate(DN14.blup[,9:11], by=list(DN14.blup$experiment, DN14.blup$year, DN14.blup$treatment, DN14.blup$dap_i, DN14.blup$id), mean, na.action="na.pass", na.rm=TRUE)
plot(DN14.blup.ag$residual~DN14.blup.ag$predicted, pch=20, cex=0.3, col='blue', xlab='Predicted', ylab='Residual')
abline(h=0, col='red')
plot(DN14.blup.ag$predicted~DN14.blup.ag$height, pch=20, cex=0.3, col="blue", ylab="Predicted", xlab="Measured")
abline(a=0, b=1, col='red')

colnames(DN14.blup.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'id', 'height', 'predicted', 'residual')
DN14_treatment_time.ag<-aggregate(DN14.blup.ag[,c(6:8)], by=list(DN14.blup.ag$experiment, DN14.blup.ag$year,DN14.blup.ag$treatment, DN14.blup.ag$dap_i), mean, na.action="na.pass", na.rm=TRUE)
colnames(DN14_treatment_time.ag)<-c('experiment', 'year', 'treatment', 'dap_i', 'height', 'predicted', 'residual')
DN14_treatment_time.ag[DN14_treatment_time.ag$dap_i == 'culm_height','dap_i']<-67

# Plot raw average
plot(DN14_treatment_time.ag$height~DN14_treatment_time.ag$dap_i, col=c("orange", "blue"), pch=20, ylab="Height (cm)", xlab="Days after planting")
plot(DN14_treatment_time.ag$predicted~DN14_treatment_time.ag$dap_i, col=c("orange", "blue"), pch=20, ylab="Height (cm)", xlab="Days after planting")

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# Lets make a ggplot version
p<-ggplot(DN14.blup, aes(x=dap_i, y=predicted, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height BLUP (cm)") + xlab("Days after planting")
pdf("DN14_height_blup_over_time.pdf")
print(p)
dev.off()

p<-ggplot(DN14.blup, aes(x=dap_i, y=height, group=treatment, colour = treatment))  + geom_smooth(method="loess") + theme_bw() + ylab("Height (cm)") + xlab("Days after planting")
pdf("DN14_height_over_time.pdf")
print(p)
dev.off()

days<-sort(unique(DN14.blup$dap_i))
DN14.blup_qtl<-c()
for(d in 1:length(days)) {
  temp<-DN14.blup[DN14.blup$dap_i == days[d],]
  colnames(temp)[10]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DN14.blup_qtl<-temp[,-c(1,8,9,11)]
  }
  if (d > 1 ) {
    DN14.blup_qtl<-merge(DN14.blup_qtl, temp[,-c(1,8,9,11)], by=c('experiment', 'year', 'treatment', 'plot', 'subplot_id', 'id'), all=T)
  }
}

# Okay now lets finish this other crap...
DN14.blup_qtl$Obs<-c(1:nrow(DN14.blup_qtl))
DN14.blup_qtl$sampling<-rep('unknown', nrow(DN14.blup_qtl))
DN14.blup_qtl<-DN14.blup_qtl[,c(10,1:6,11,7:9)]

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
# write out 
write.csv(DN14.blup_qtl, file="DN14.height.blup_qtl.csv", row.names=F, quote=F)


########################################################################
# Flowering time
########################################################################


setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")
DR13_ft<-read.csv("13DR_panicle_emergence_step_2.csv")
DR13_ft[DR13_ft$treatment == "wet",'treatment']<-c("wet")
DR13_ft[DR13_ft$treatment == "dry",'treatment']<-c("dry")

DR14_ft<-read.csv("14DR_panicle_emergence_step_2.csv")
DR14_ft[DR14_ft$treatment == "wet",'treatment']<-c("wet")
DR14_ft[DR14_ft$treatment == "dry",'treatment']<-c("dry")

DN13_ft<-read.csv("13DN_panicle_emergence_step_2.csv", colClasses=c("treatment"="character"))
DN13_ft[DN13_ft$treatment == "high_density","treatment"]<-c("dense")
DN13_ft[DN13_ft$treatment == "low_density","treatment"]<-c("sparse")

DN14_ft<-read.csv("14DN_panicle_emergence_step_2.csv", colClasses=c("treatment"="character"), na.strings = ".")
DN14_ft[DN14_ft$treatment == "high_density","treatment"]<-c("dense")
DN14_ft[DN14_ft$treatment == "low_density","treatment"]<-c("sparse")


######## REFORMAT AND CALCULATE BLUP

##### DR13
DR13_ft<-DR13_ft[,c(8,9,10,13)]
colnames(DR13_ft)[4]<-c("panicle_emergence")

DR13_ft_model<-lmer(panicle_emergence ~ (1|genotype) + treatment + plot %in% treatment + (1|genotype:treatment), data=DR13_ft)
DR13_ft$predicted<-predict(DR13_ft_model)
DR13_ft$residual<-residuals(DR13_ft_model)

plot(DR13_ft$predicted~DR13_ft$panicle_emergence, pch=20, cex=0.3, col=c("blue"), ylab=c("Predicted DAS"), xlab=c("DAS"))
plot(DR13_ft$residual~DR13_ft$predicted, pch=20, cex=0.3, col=c("blue"))


treatments<-unique(DR13_ft$treatment)
DR13_ft.blup<-c()
for(i in 1:length(treatments)){
  temp<-DR13_ft[DR13_ft$treatment == treatments[i],]
  temp.ag<-aggregate(temp[,5], by=list(temp$genotype), FUN=mean, na.rm=TRUE, na.action=NULL )
  colnames(temp.ag)[1:2]<-c("genotype", paste("DR13", "panicle_emergence", treatments[i], sep="_"))
  if (i == 1){
    DR13_ft.blup<-temp.ag
  }
  if (i > 1){
    DR13_ft.blup<-merge(DR13_ft.blup, temp.ag, by=c("genotype"))
  }
}
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_flowering/illinois_field")
write.csv(DR13_ft.blup, file="DR13.flower.blup.csv", row.names=F, quote=F)

##### DR14

DR14_ft<-DR14_ft[,c(8,9,10,13)]
colnames(DR14_ft)[4]<-c("panicle_emergence")

DR14_ft_model<-lmer(panicle_emergence ~ (1|genotype) + treatment + plot %in% treatment + (1|genotype:treatment), data=DR14_ft)
DR14_ft$predicted<-predict(DR14_ft_model)
DR14_ft$residual<-residuals(DR14_ft_model)

plot(DR14_ft$predicted~DR14_ft$panicle_emergence, pch=20, cex=0.3, col=c("blue"), ylab=c("Predicted DAS"), xlab=c("DAS"))
plot(DR14_ft$residual~DR14_ft$predicted, pch=20, cex=0.3, col=c("blue"))


treatments<-unique(DR14_ft$treatment)
DR14_ft.blup<-c()
for(i in 1:length(treatments)){
  temp<-DR14_ft[DR14_ft$treatment == treatments[i],]
  temp.ag<-aggregate(temp[,5], by=list(temp$genotype), FUN=mean, na.rm=TRUE, na.action=NULL )
  colnames(temp.ag)[1:2]<-c("genotype", paste("DR14", "panicle_emergence", treatments[i], sep="_"))
  if (i == 1){
    DR14_ft.blup<-temp.ag
  }
  if (i > 1){
    DR14_ft.blup<-merge(DR14_ft.blup, temp.ag, by=c("genotype"))
  }
}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_flowering/illinois_field")
write.csv(DR14_ft.blup, file="DR14.flower.blup.csv", row.names=F, quote=F)

##### DN13

DN13_ft<-DN13_ft[,c(8,9,10,17)]
colnames(DN13_ft)[4]<-c("panicle_emergence")

DN13_ft_model<-lmer(panicle_emergence ~ (1|genotype) + treatment + plot %in% treatment + (1|genotype:treatment), data=DN13_ft)
DN13_ft$predicted<-predict(DN13_ft_model)
DN13_ft$residual<-residuals(DN13_ft_model)

plot(DN13_ft$predicted~DN13_ft$panicle_emergence, pch=20, cex=0.3, col=c("blue"), ylab=c("Predicted DAS"), xlab=c("DAS"))
plot(DN13_ft$residual~DN13_ft$predicted, pch=20, cex=0.3, col=c("blue"))

treatments<-unique(DN13_ft$treatment)
DN13_ft.blup<-c()
for(i in 1:length(treatments)){
  temp<-DN13_ft[DN13_ft$treatment == treatments[i],]
  temp.ag<-aggregate(temp[,5], by=list(temp$genotype), FUN=mean, na.rm=TRUE, na.action=NULL )
  colnames(temp.ag)[1:2]<-c("genotype", paste("DN13", "panicle_emergence", treatments[i], sep="_"))
  if (i == 1){
    DN13_ft.blup<-temp.ag
  }
  if (i > 1){
    DN13_ft.blup<-merge(DN13_ft.blup, temp.ag, by=c("genotype"))
  }
}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_flowering/illinois_field")
write.csv(DN13_ft.blup, file="DN13.flower.blup.csv", row.names=F, quote=F)

##### DN14

DN14_ft<-DN14_ft[,c(8,9,10,14)]
colnames(DN14_ft)[4]<-c("panicle_emergence")
DN14_ft$plot<-as.factor(DN14_ft$plot)
DN14_ft<-DN14_ft[complete.cases(DN14_ft),]
DN14_ft_model<-lmer(panicle_emergence ~ (1|genotype) + treatment + plot %in% treatment + (1|genotype:treatment), data=DN14_ft)
DN14_ft$predicted<-predict(DN14_ft_model)
DN14_ft$residual<-residuals(DN14_ft_model)

plot(DN14_ft$predicted~DN14_ft$panicle_emergence, pch=20, cex=0.3, col=c("blue"), ylab=c("Predicted DAS"), xlab=c("DAS"))
plot(DN14_ft$residual~DN14_ft$predicted, pch=20, cex=0.3, col=c("blue"))

treatments<-unique(DN14_ft$treatment)
DN14_ft.blup<-c()
for(i in 1:length(treatments)){
  temp<-DN14_ft[DN14_ft$treatment == treatments[i],]
  temp.ag<-aggregate(temp[,5], by=list(temp$genotype), FUN=mean, na.rm=TRUE, na.action=NULL )
  colnames(temp.ag)[1:2]<-c("genotype", paste("DN14", "panicle_emergence", treatments[i], sep="_"))
  if (i == 1){
    DN14_ft.blup<-temp.ag
  }
  if (i > 1){
    DN14_ft.blup<-merge(DN14_ft.blup, temp.ag, by=c("genotype"))
  }
}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_flowering/illinois_field")
write.csv(DN14_ft.blup, file="DN14.flower.blup.csv", row.names=F, quote=F)








##########################################################################
# functions
##########################################################################

get.loess.fit<-function(fit, times, geno, cond) {
  return_df<-c()
  predict.vals<-predict(fit, data.frame(times), se=TRUE)
  accession<-rep(geno, length(times))
  condition<-rep(as.character(cond), length(times))
  M<-predict.vals$fit
  M.lo<-M - predict.vals$se.fit
  M.hi<-M + predict.vals$se.fit
  slope<-c(0)
  for(s in 2:length(times)) {
    s.temp<-(M[s] - M[s-1]) / 2
    slope<-c(slope, s.temp)
  }
  return_df<-cbind(accession, condition,times, M, M.lo, M.hi, slope)
  return(return_df)
}


get.mean.slope<-function(df, r.days, geno, cond) {
  return_df<-c()
  accession<-rep(geno, length(r.days))
  condition<-rep(as.character(cond), length(r.days))
  M<-df$height_above_bound
  slope<-c(0)
  for(s in 2:length(r.days)) {
    s.temp<-(M[s] - M[s-1]) / 2
    slope<-c(slope, s.temp)
  }
  return_df<-cbind(accession, condition,r.days, M, slope)
  return(return_df)
}
