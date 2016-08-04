########################################################
# Correlation
########################################################
# Load in Illinois BLUPs
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
DR13.blup<-read.csv("DR13.height.blup_qtl.csv")
DR14.blup<-read.csv("DR14.height.blup_qtl.csv")
DN13.blup<-read.csv("DN13.height.blup_qtl.csv")
DN14.blup<-read.csv("DN14.height.blup_qtl.csv")

# Peformed on BLUP values
treatments<-unique(DR13.blup$treatment)
DR13.blup.ag<-aggregate(DR13.blup[,c(9:ncol(DR13.blup))], by=list(DR13.blup$id, DR13.blup$treatment, DR13.blup$experiment), FUN=mean, na.rm=TRUE, na.action=NULL)
colnames(DR13.blup.ag)<-c("genotype", "treatment", "experiment", colnames(DR13.blup)[9:ncol(DR13.blup)])
DR13.blup.ag.ld<-DR13.blup.ag[,c(1:3,10)]
DR13.final_estimate<-c()
for(i in 1:length(treatments)) {
  temp<-DR13.blup.ag.ld[DR13.blup.ag.ld$treatment == treatments[i],]
  colnames(temp)[4]<-paste(unique(temp$experiment), treatments[i], 'height', sep="_")
  if (i == 1){
    DR13.final_estimate<-temp[,c(1,4)]
  }
  if (i > 1) {
    DR13.final_estimate<-merge(DR13.final_estimate, temp[,c(1,4)], by=c("genotype"))
  }
}


DR14.blup.ag<-aggregate(DR14.blup[,c(9:ncol(DR14.blup))], by=list(DR14.blup$id, DR14.blup$treatment, DR14.blup$experiment), mean, na.rm=TRUE, na.action=NULL)
colnames(DR14.blup.ag)<-c("genotype", "treatment", "experiment", colnames(DR14.blup)[9:ncol(DR14.blup)])
DR14.blup.ag.ld<-DR14.blup.ag[,c(1:3,7)]
DR14.final_estimate<-c()
for(i in 1:length(treatments)) {
  temp<-DR14.blup.ag.ld[DR14.blup.ag.ld$treatment == treatments[i],]
  colnames(temp)[4]<-paste(unique(temp$experiment), treatments[i], 'height', sep="_")
  if (i == 1){
    DR14.final_estimate<-temp[,c(1,4)]
  }
  if (i > 1) {
    DR14.final_estimate<-merge(DR14.final_estimate, temp[,c(1,4)], by=c("genotype"))
  }
}

BP14.est<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound/ril_best.fit_logistic_estimates_height.csv")
BP14.est.ld<-BP14.est[BP14.est$dap_i == 30,]
BP14.est.ld$experiment<-rep("BP14", nrow(BP14.est.ld))
BP14.est.ld<-BP14.est.ld[,c(1,2,16,4)]
colnames(BP14.est.ld)[1]<-c("genotype")

BP14.final_estimate<-c()
for(i in 1:length(treatments)) {
  temp<-BP14.est.ld[BP14.est.ld$treatment == treatments[i],]
  colnames(temp)[4]<-paste(unique(temp$experiment), treatments[i], 'height', sep="_")
  if (i == 1){
    BP14.final_estimate<-temp[,c(1,4)]
  }
  if (i > 1) {
    BP14.final_estimate<-merge(BP14.final_estimate, temp[,c(1,4)], by=c("genotype"))
  }
}

# Need to read int DL height
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/dinneny")
dl.height<-read.csv('dinneny_qtl.csv')

#dl.height.ag<-aggregate(dl.height[,9], by=list(dl.height$id, dl.height$treatment, dl.height$experiment), FUN=mean, na.rm=TRUE, na.action=NULL)
dl.height.ag<-aggregate(dl.height[,10], by=list(dl.height$id, dl.height$treatment, dl.height$experiment), FUN=mean, na.rm=TRUE, na.action=NULL)

colnames(dl.height.ag)<-c("genotype","treatment","experiment","height")

treatments<-unique(dl.height.ag$treatment)
DL13.final_estimate<-c()
for(i in 1:length(treatments)) {
  temp<-dl.height.ag[dl.height.ag$treatment == treatments[i],]
  colnames(temp)[4]<-paste(unique(temp$experiment), treatments[i], 'height', sep="_")
  if (i == 1){
    DL13.final_estimate<-temp[,c(1,4)]
  }
  if (i > 1) {
    DL13.final_estimate<-merge(DL13.final_estimate, temp[,c(1,4)], by=c("genotype"))
  }
}


####################### DENSITY
treatments<-unique(DN13.blup$treatment)

DN13.blup.ag<-aggregate(DN13.blup[,c(9:ncol(DN13.blup))], by=list(DN13.blup$id, DN13.blup$treatment, DN13.blup$experiment), FUN=mean, na.rm=TRUE, na.action=NULL)
colnames(DN13.blup.ag)<-c("genotype", "treatment", "experiment",colnames(DN13.blup)[9:ncol(DN13.blup)])
DN13.blup.ag.ld<-DN13.blup.ag[,c(1:3,11)]

DN13.final_estimate<-c()
for(i in 1:length(treatments)) {
  temp<-DN13.blup.ag.ld[DN13.blup.ag.ld$treatment == treatments[i],]
  colnames(temp)[4]<-paste(unique(temp$experiment), treatments[i], 'height', sep="_")
  if (i == 1){
    DN13.final_estimate<-temp[,c(1,4)]
  }
  if (i > 1) {
    DN13.final_estimate<-merge(DN13.final_estimate, temp[,c(1,4)], by=c("genotype"))
  }
}

DN14.blup.ag<-aggregate(DN14.blup[,c(9:ncol(DN14.blup))], by=list(DN14.blup$id, DN14.blup$treatment, DN14.blup$experiment), FUN=mean, na.rm=TRUE, na.action=NULL)
colnames(DN14.blup.ag)<-c("genotype", "treatment", "experiment", colnames(DN14.blup)[9:ncol(DN14.blup)])
DN14.blup.ag.ld<-DN14.blup.ag[,c(1:3,5)]


DN14.final_estimate<-c()
for(i in 1:length(treatments)) {
  temp<-DN14.blup.ag.ld[DN14.blup.ag.ld$treatment == treatments[i],]
  colnames(temp)[4]<-paste(unique(temp$experiment), treatments[i], 'height', sep="_")
  if (i == 1){
    DN14.final_estimate<-temp[,c(1,4)]
  }
  if (i > 1) {
    DN14.final_estimate<-merge(DN14.final_estimate, temp[,c(1,4)], by=c("genotype"))
  }
}

DL13.final_estimate[1:5,]
DR13.final_estimate[1:5,]
DR14.final_estimate[1:5,]
BP14.final_estimate[1:5,]
DN13.final_estimate[1:5,]
DN14.final_estimate[1:5,]

final_day_heights<-merge(DL13.final_estimate, DR13.final_estimate, by=c("genotype"))
final_day_heights<-merge(final_day_heights, DR14.final_estimate, by=c("genotype"))
final_day_heights<-merge(final_day_heights, BP14.final_estimate, by=c("genotype"))
final_day_heights<-merge(final_day_heights, DN13.final_estimate, by=c("genotype"))
final_day_heights<-merge(final_day_heights, DN14.final_estimate, by=c("genotype"))

final_day_height_rank<-c()
for (i in 2:ncol(final_day_heights)) {
  temp<-final_day_heights[,c(1,i)]
  temp[,2]<-rank(temp[,2], ties.method=c("first"))
  if (i == 2) {
    final_day_height_rank<-temp
  }
  
  if (i > 2){
    final_day_height_rank<-merge(final_day_height_rank, temp, by=c("genotype"))
  }
}

rownames(final_day_height_rank)<-final_day_height_rank$genotype
final_day_height_rank<-final_day_height_rank[,2:ncol(final_day_height_rank)]
final_day_height_rank<-as.matrix(final_day_height_rank)
colnames(final_day_height_rank)<-sub("_height", "", colnames(final_day_height_rank))
my_palette<-colorRampPalette(c("blue", "white", "red"))(n = 1000)
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/")
pdf("heatmap_final_day_height_rank.pdf")
heatmap(final_day_height_rank,Colv=F, scale='none', col=my_palette)
dev.off()


setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
my_palette<-colorRampPalette(c("blue", "white", "red"))(n = 1000)
pdf("Figure_S4.pdf")
heatmap(final_day_height_rank,Colv=F, scale='none', col=my_palette)
dev.off()

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")

cor.final_day_height_rank<-cor(final_day_height_rank)
my_palette<-colorRampPalette(c("white", "cyan3"))(n = 1000)
pdf("Figure_2.pdf")
heatmap(cor.final_day_height_rank, Colv=T, Rowv=T, col=my_palette, symm=T)
dev.off()



setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
ril<-read.csv("ril_best.fit_logistic_estimates_height.csv")
ril.anova<-get.field_anova.regression_format(ril)


########################################################################################
########################################################################################
########################################################################################

# Load in flowering time BLUPs
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_flowering/illinois_field")
DR13_ft.blup<-read.csv("DR13.flower.blup.csv")
DR14_ft.blup<-read.csv("DR14.flower.blup.csv")
DN13_ft.blup<-read.csv("DN13.flower.blup.csv")
DN14_ft.blup<-read.csv("DN14.flower.blup.csv")

# Put all values in common data frame
height_v_ft<-merge(DR13_ft.blup, DR14_ft.blup, by=c("genotype"))
height_v_ft<-merge(height_v_ft, DN13_ft.blup, by=c("genotype"))
height_v_ft<-merge(height_v_ft, DN14_ft.blup, by=c("genotype"))

# Subset to only field experiments
height_v_ft<-merge(final_day_heights, height_v_ft, by=c("genotype"))
ht_v_ft_field<-height_v_ft[,c(4:7,10:21)]
ht_v_ft_field<-ht_v_ft_field[complete.cases(ht_v_ft_field),]
cor.height_v_ft<-cor(ht_v_ft_field)
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/")
write.csv(cor.height_v_ft, file="height_v_flowering_time_field_ril.csv", quote=F, row.names=T)

# Lets make a plot of this
heights<-c(height_v_ft[,4], height_v_ft[,5], height_v_ft[,6], height_v_ft[,7], height_v_ft[,12], height_v_ft[,13])
flowering_time<-c(height_v_ft[,14], height_v_ft[,15], height_v_ft[,17], height_v_ft[,16], height_v_ft[,20], height_v_ft[,21])
plot(heights~flowering_time, pch=20, cex=0.5, col="blue", ylab="Height (cm)", xlab="Days to panicle emergence")


heights<-c(height_v_ft$DR13_dry_height, height_v_ft$DR13_wet_height, height_v_ft$DR14_dry_height, height_v_ft$DR14_wet_height, height_v_ft$DN13_dense_height, height_v_ft$DN13_sparse_height, height_v_ft$DN14_dense_height, height_v_ft$DN14_sparse_height)
colors<-c(rep("orange", 74), rep("blue", 74), rep("orange", 74), rep("blue", 74), rep("darkgray", 74), rep("green", 74), rep("darkgray", 74), rep("green", 74))

flowering<-c(height_v_ft$DR13_panicle_emergence_dry, height_v_ft$DR13_panicle_emergence_wet, height_v_ft$DR14_panicle_emergence_dry, height_v_ft$DR14_panicle_emergence_wet, height_v_ft$DN13_panicle_emergence_dense, height_v_ft$DN13_panicle_emergence_sparse, height_v_ft$DN14_panicle_emergence_dense, height_v_ft$DN14_panicle_emergence_sparse)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/")
pdf("height_v_flowering_time.pdf")
plot(heights~flowering, pch=20, col=colors, cex=0.5, ylab=c("Height (cm)"), xlab=c("Days to panicle emeregence"))
abline(a=0, b=1, col="red")
dev.off()

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril")
save.image("correlation_between_exp_and_flowering.Rdata")
