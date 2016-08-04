# Summary of growouts for table 1
library(ggplot2)
library(lme4)
library(nlme)
library(car)


# Lets calculate some of the metrics in Table1
####################################################
# Illinois field

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field")
DR13.blup<-read.csv('DR13.height.blup_qtl.csv')
DR14.blup<-read.csv('DR14.height.blup_qtl.csv')
DN13.blup<-read.csv('DN13.height.blup_qtl.csv')
DN14.blup<-read.csv('DN14.height.blup_qtl.csv')

DR13.w<-DR13.blup[DR13.blup$treatment == 'wet','height_67']
DR13.d<-DR13.blup[DR13.blup$treatment == 'dry','height_67']
#dry
quantile(DR13.d, c(0.05,0.5, 0.95), na.rm=T)
mean(DR13.d, na.rm=T)
var(DR13.d, na.rm=T)
sqrt(var(DR13.d, na.rm=T))/mean(DR13.d, na.rm=T)

#wet
quantile(DR13.w, c(0.05,0.5, 0.95), na.rm=T)
mean(DR13.w, na.rm=T)
var(DR13.w, na.rm=T)
sqrt(var(DR13.w, na.rm=T))/mean(DR13.w, na.rm=T)


DR14.w<-DR14.blup[DR14.blup$treatment == 'wet',"height_47"]
DR14.d<-DR14.blup[DR14.blup$treatment == 'dry',"height_47"]

#dry
quantile(DR14.d, c(0.05,0.5, 0.95), na.rm=T)
mean(DR14.d, na.rm=T)
var(DR14.d, na.rm=T)
sqrt(var(DR14.d, na.rm=T))/mean(DR14.d, na.rm=T)

#wet
quantile(DR14.w, c(0.05,0.5, 0.95), na.rm=T)
mean(DR14.w, na.rm=T)
var(DR14.w, na.rm=T)
sqrt(var(DR14.w, na.rm=T))/mean(DR14.w, na.rm=T)


DN13.s<-DN13.blup[DN13.blup$treatment == 'sparse',"height_67"]
DN13.d<-DN13.blup[DN13.blup$treatment == 'dense',"height_67"]

#dense
quantile(DN13.d, c(0.05,0.5, 0.95), na.rm=T)
mean(DN13.d, na.rm=T)
var(DN13.d, na.rm=T)
sqrt(var(DN13.d, na.rm=T))/mean(DN13.d, na.rm=T)

#sparse
quantile(DN13.s, c(0.05,0.5, 0.95), na.rm=T)
mean(DN13.s, na.rm=T)
var(DN13.s, na.rm=T)
sqrt(var(DN13.s, na.rm=T))/mean(DN13.s, na.rm=T)


DN14.s<-DN14.blup[DN14.blup$treatment == 'sparse',"height_67"]
DN14.d<-DN14.blup[DN14.blup$treatment == 'dense',"height_67"]

#dense
quantile(DN14.d, c(0.05,0.5, 0.95), na.rm=T)
mean(DN14.d, na.rm=T)
var(DN14.d, na.rm=T)
sqrt(var(DN14.d, na.rm=T))/mean(DN14.d, na.rm=T)

#sparse
quantile(DN14.s, c(0.05,0.5, 0.95), na.rm=T)
mean(DN14.s, na.rm=T)
var(DN14.s, na.rm=T)
sqrt(var(DN14.s, na.rm=T))/mean(DN14.s, na.rm=T)

####################################################
# Dinneny
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/dinneny")
dl<-read.csv('dinneny_qtl.csv')
dl.d<-dl[dl$treatment == 'dry', 'plant_height']
dl.w<-dl[dl$treatment == 'wet', 'plant_height']

quantile(dl.d, c(0.05,0.5, 0.95), na.rm=T)
mean(dl.d, na.rm=T)
var(dl.d, na.rm=T)
sqrt(var(dl.d, na.rm=T))/mean(dl.d, na.rm=T)


#wet
quantile(dl.w, c(0.05,0.5, 0.95), na.rm=T)
mean(dl.w, na.rm=T)
var(dl.w, na.rm=T)
sqrt(var(dl.w, na.rm=T))/mean(dl.w, na.rm=T)

####################################################
# Bellweather RIL
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
ril<-read.csv("ril_best.fit_logistic_estimates_height.csv")
ril.30<-ril[ril$dap_i == 30,]
ril.30.d<-ril.30[ril.30$treatment == 'dry', 'M']
ril.30.w<-ril.30[ril.30$treatment == 'wet', 'M']

#dry
quantile(ril.30.d, c(0.05,0.5, 0.95), na.rm=T)
mean(ril.30.d, na.rm=T)
var(ril.30.d, na.rm=T)
sqrt(var(ril.30.d, na.rm=T))/mean(ril.30.d, na.rm=T)

#wet
quantile(ril.30.w, c(0.05,0.5, 0.95), na.rm=T)
mean(ril.30.w, na.rm=T)
var(ril.30.w, na.rm=T)
sqrt(var(ril.30.w, na.rm=T))/mean(ril.30.w, na.rm=T)


####################################################
# Lets calculate average replication within experiment

# Illinois
# Load in illinois data
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")
DR13<-read.csv('DR13_qtl.csv')
DR14<-read.csv('DR14_qtl.csv')
DN13<-read.csv('DN13_qtl.csv')
DN14<-read.csv('DN14_qtl.csv')

DR13.w<-DR13[DR13$treatment == 'wet',]
DR13.d<-DR13[DR13$treatment == 'dry',]

mean(table(DR13.d$genotype)[3:length(table(DR13.d$id))])
mean(table(DR13.w$id)[3:length(table(DR13.w$id))])

DR14.w<-DR14[DR14$treatment == 'wet',]
DR14.d<-DR14[DR14$treatment == 'dry',]

mean(table(DR14.d$id)[3:length(table(DR14.d$id))])
mean(table(DR14.w$id)[3:length(table(DR14.w$id))])

DN13.s<-DN13[DN13$treatment == 'sparse',]
DN13.d<-DN13[DN13$treatment == 'dense',]

mean(table(DN13.d$id)[3:length(table(DN13.d$id))])
mean(table(DN13.s$id)[3:length(table(DN13.s$id))])


DN14.s<-DN14[DN14$treatment == 'sparse',]
DN14.d<-DN14[DN14$treatment == 'dense',]

mean(table(DN14.d$id)[3:length(table(DN14.d$id))])
mean(table(DN14.s$id)[3:length(table(DN14.s$id))])

#############
# Bellweather RIL
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril/")
ril.raw<-read.csv("ril_height_qtl_raw.csv")

ril.w<-ril.raw[ril.raw$treatment == 'wet',]
ril.d<-ril.raw[ril.raw$treatment == 'dry',]

mean(table(ril.d$id)[3:length(table(ril.d$id))])
mean(table(ril.w$id)[3:length(table(ril.w$id))])

######## How much variability throughout time and across experiments?

# Lets calculate CV (coefficient of variation)
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/table_data/cv_data")
# DR13
t.points<-c(9:ncol(DR13.blup))
DR13.cv<-getCV.qtl_format(DR13.blup, t.points)
DR13.cv<-as.data.frame(DR13.cv)
DR13.cv$location<-rep('Illinois', nrow(DR13.cv))
DR13.cv$year<-rep('2013', nrow(DR13.cv))
DR13.cv<-DR13.cv[,c(10,11,1:9)]
write.csv(DR13.cv,"DR13_cv.csv", row.names=F, quote=F)

#plot(y=DR13.cv[1,4:ncol(DR13.cv)], x=c(4:ncol(DR13.cv)), type="l", col="orange", ylim=c(0,0.5))
#points(DR13.cv[2,4:ncol(DR13.cv)]~c(4:ncol(DR13.cv)), type="l", col="blue")

# DN13
t.points<-c(9:ncol(DN13.blup))
DN13.cv<-getCV.qtl_format(DN13.blup, t.points)
DN13.cv<-as.data.frame(DN13.cv)
DN13.cv$location<-rep('Illinois', nrow(DN13.cv))
DN13.cv$year<-rep('2013', nrow(DN13.cv))
DN13.cv<-DN13.cv[,c(10,11,1:9)]
write.csv(DN13.cv,"DN13_cv.csv", row.names=F, quote=F)

#plot(DN13.cv[1,4:ncol(DN13.cv)]~c(4:ncol(DN13.cv)), type="l", col="grey", ylim=c(0,0.5))
#points(DN13.cv[2,2:ncol(DN13.cv)]~c(2:ncol(DN13.cv)), type="l", col="green")

# DR14
t.points<-c(9:ncol(DR14.blup))
DR14.cv<-getCV.qtl_format(DR14.blup, t.points)
DR14.cv<-as.data.frame(DR14.cv)
DR14.cv$location<-rep('Illinois', nrow(DR14.cv))
DR14.cv$year<-rep('2014', nrow(DR14.cv))
DR14.cv<-DR14.cv[,c(7,8,1:6)]
write.csv(DR14.cv,"DR14_cv.csv", row.names=F, quote=F)

#plot(DR14.cv[1,2:ncol(DR14.cv)]~c(2:ncol(DR14.cv)), type="l", col="orange", ylim=c(0,0.5))
#points(DR14.cv[2,2:ncol(DR14.cv)]~c(2:ncol(DR14.cv)), type="l", col="blue")

# DN14
t.points<-c(9:ncol(DN14.blup))
DN14.cv<-getCV.qtl_format(DN14.blup, t.points)
DN14.cv<-as.data.frame(DN14.cv)
DN14.cv$location<-rep('Illinois', nrow(DN14.cv))
DN14.cv$year<-rep('2014', nrow(DN14.cv))
DN14.cv<-DN14.cv[,c(5,6,1:4)]
write.csv(DN14.cv,"DN14_cv.csv", row.names=F, quote=F)

#plot(DN14.cv[1,2:ncol(DN14.cv)]~c(2:ncol(DN14.cv)), type="l", col="grey", ylim=c(0,0.5))
#points(DN14.cv[2,2:ncol(DN14.cv)]~c(2:ncol(DN14.cv)), type="l", col="green")

# Lets look at dinneny lab 
dl.height<-dl[,c(1:8,10)]
#dl.height$empty1<-rep(NA, nrow(dl.height))
#dl.height$empty2<-rep(NA, nrow(dl.height))
#dl.height<-dl.height[,c(8,9,1:7)]
t.points<-c(9:ncol(dl.height))
dl.height.cv<-getCV.qtl_format(dl.height, t.points)
dl.height.cv<-as.data.frame(dl.height.cv)
dl.height.cv$location<-rep('Carnagie', nrow(dl.height.cv))
dl.height.cv$year<-rep('2013', nrow(dl.height.cv))
dl.height.cv<-dl.height.cv[,c(3,4,1:2)]
write.csv(dl.height.cv,"dl.height_cv.csv", row.names=F, quote=F)

# Now lets reformat bellweather to qtl format data to calculate CV
ril_height_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound/ril_height_above_bound_qtl.csv")
#ril_height_qtl<-c()
#days<-sort(unique(ril.mean$dap_i))
#days_e<-days[seq(1,length(days),by=2)]

#for(d in 1:length(days_e)) {
#  day<-days_e[d]
#  data<-ril.mean[ril.mean$dap_i == day | ril.mean$dap == (day + 1), ]
#  colnames(data)[6]<-paste('height', day, sep="_")
#  data<-data[,c(1,3,4,6)]
#  if (d == 1) {
#    ril_height_qtl<-data
#  }
#  
#  if (d > 1) {
#    ril_height_qtl<-merge(ril_height_qtl, data, by=c('plant_id','genotype', 'treatment'), all=T)
#  }
#}

#ril_height_qtl$Obs<-c(1:nrow(ril_height_qtl))
#ril_height_qtl$experiment<-rep('BP14', nrow(ril_height_qtl))
#ril_height_qtl$year<-rep('2014', nrow(ril_height_qtl))
#ril_height_qtl$plot_id<-rep('phenotyper', nrow(ril_height_qtl))
#ril_height_qtl$subplot_id<-rep('phenotyper', nrow(ril_height_qtl))

#ril_height_qtl<-ril_height_qtl[,c(17,18,19,3,20,21,2,1,4:16)]
#colnames(ril_height_qtl)[7]<-c("id")
#setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril")
#write.csv(ril_height_qtl, 'ril_height_qtl_raw.csv', quote=F, row.names=F)
#setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")

# Get CV for phenotyper
t.points<-c(9:ncol(ril_height_qtl))
bp14.cv<-getCV.qtl_format(ril_height_qtl, t.points)
bp14.cv<-as.data.frame(bp14.cv)
bp14.cv$location<-rep('DDPSC', nrow(bp14.cv))
bp14.cv$year<-rep('2014', nrow(bp14.cv))
bp14.cv<-bp14.cv[,c(28,29,1:27)]
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/table_data/cv_data")
write.csv(bp14.cv,"BP14_cv.csv", row.names=F, quote=F)

# Plot it
#plot(bp14.cv[1,2:ncol(bp14.cv)]~c(2:ncol(bp14.cv)), type="l", col="orange", ylim=c(0,0.5))
#points(bp14.cv[2,2:ncol(bp14.cv)]~c(2:ncol(bp14.cv)), type="l", col="blue")

#################################
# Make boxplots for Figure 1a
#################################

# Load data from the top of the script.

dn13.b<-DN13.blup[,c(2,7,4,16)]
colnames(dn13.b)<-c("experiment", "id", "treatment", "height")
dn13.b<-dn13.b[complete.cases(dn13.b),]
dn13.b$height<-as.numeric(as.character(dn13.b$height))

dn14.b<-DN14.blup[,c(2,7,4,11)]
colnames(dn14.b)<-c("experiment", "id", "treatment", "height")
dn14.b<-dn14.b[complete.cases(dn14.b),]
dn14.b$height<-as.numeric(as.character(dn14.b$height))

dr13.b<-DR13.blup[,c(2,7,4,15)]
colnames(dr13.b)<-c("experiment", "id", "treatment", "height")
dr13.b<-dr13.b[complete.cases(dr13.b),]
dr13.b$height<-as.numeric(as.character(dr13.b$height))

dr14.b<-DR14.blup[,c(2,7,4,13)]
colnames(dr14.b)<-c("experiment", "id", "treatment", "height")
dr14.b<-dr14.b[complete.cases(dr14.b),]
dr14.b$height<-as.numeric(as.character(dr14.b$height))

dl.h<-dl[,c(2,7,4,10)]
colnames(dl.h)<-c("experiment", "id", "treatment", "height")
dl.h<-dl.h[complete.cases(dl.h),]
dl.h$height<-as.numeric(as.character(dl.h$height))

ril$experiment<-rep("BP14", nrow(ril))
ril.fh<-ril[ril$dap_i==30,c(16,1,2,4)]
colnames(ril.fh)<-c("experiment", "id", "treatment", "height")
ril.fh<-ril.fh[complete.cases(ril.fh),]
ril.fh$height<-as.numeric(as.character(ril.fh$height))

final_heights<-rbind(dn13.b, dn14.b, dr13.b, dr14.b, dl.h, ril.fh)
final_heights<-final_heights[complete.cases(final_heights),]
final_heights$group<-paste(final_heights$experiment, final_heights$treatment, sep="_")
final_heights$group<-factor(final_heights$group, levels= c("DN13_sparse","DN13_dense","DN14_sparse","DN14_dense","DR13_wet","DR13_dry","DR14_wet","DR14_dry","DL13_wet","DL13_dry","BP14_wet","BP14_dry"), ordered=TRUE)

groups<-unique(final_heights$group)
final_heights$a10_mean<-rep('NA', nrow(final_heights))
final_heights$b100_mean<-rep('NA', nrow(final_heights))

for(g in 1:length(groups)){
  a10<-c()
  a10<-final_heights[final_heights$group == groups[g] & final_heights$id == 'A10', 'height']
  if (length(a10) > 0) {
    final_heights[final_heights$group == groups[g],'a10_mean']<-mean(a10)
  } 
  if (length(a10) == 0){
    final_heights[final_heights$group == groups[g],'a10_mean']<-c('NA')
  }
  b100<-final_heights[final_heights$group == groups[g] & final_heights$id == 'B100', 'height']
  if (length(b100) > 0) {
    final_heights[final_heights$group == groups[g],'b100_mean']<-mean(b100)
  } 
  if (length(b100) == 0){
    final_heights[final_heights$group == groups[g],'b100_mean']<-c('NA')
  }
}

# A10 Final height missing from DN14 sparse
a10_dn14_sparse_day46.mean<-mean(DN14.blup[DN14.blup$id == 'A10' & DN14.blup$treatment == 'sparse', 'height_46'])
final_heights[final_heights$group == 'DN14_sparse', 'a10_mean']<-a10_dn14_sparse_day46.mean

#setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table/")
#pdf("Figure_1a.pdf")
p<-ggplot(final_heights, aes(factor(group), height)) + geom_boxplot(aes(fill=factor(treatment))) + scale_fill_manual('Treatments', values=c('grey', 'green', 'orange','blue'))+ theme_bw() + xlab("") + ylab("Height (cm)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  geom_point(aes(factor(group), as.numeric(a10_mean)), colour = "red", size = 3, na.rm=TRUE) +  geom_point(aes(factor(group), as.numeric(b100_mean)), colour = "purple", size = 3, na.rm=T)   
p
#dev.off()

#################################
# Make boxplots for Figure 1b
#################################
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")
DR13_height_total_variance.l<-read.csv("DR13_height_total_variance_np.l.csv")
DR14_height_total_variance.l<-read.csv("DR14_height_total_variance_np.l.csv")
DN13_height_total_variance.l<-read.csv("DN13_height_total_variance_np.l.csv")
DN14_height_total_variance.l<-read.csv("DN14_height_total_variance_np.l.csv")
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_bellweather")
BP14_height_total_variance.l<-read.csv("BP14_height_total_variance.l.csv")

all_variance<-rbind(DN13_height_total_variance.l[29:32,], DN14_height_total_variance.l[9:12,], DR13_height_total_variance.l[25:28,], DR14_height_total_variance.l[17:20,], BP14_height_total_variance.l[45:48,])
all_variance$Trait<-c(rep("DN13", 4), rep("DN14", 4), rep("DR13", 4), rep("DR14", 4), rep("BP14", 4))
all_variance$Trait<-factor(all_variance$Trait, levels = c('DN13', 'DN14', 'DR13', 'DR14', 'BP14'), ordered=TRUE)


p = ggplot(data=all_variance, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
p=p + geom_bar(width = 1, stat="identity")
p=p+facet_grid(facets=. ~ Trait)
p=p+xlab("Trait") + theme(legend.title=element_blank())


setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril")
pdf("all_ril_variance_at_final_timepoint.pdf")
p
dev.off()

#setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")
#pdf("Figure_1b.pdf")
#p
#dev.off()








#save.image(file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/ril_correlation_H2_and_total_variance.Rdata")
#rm(list=ls())











getCV.qtl_format<-function(dataset, times) {
  pheno<-dataset
  t.points<-times
  i.treat<-unique(pheno$treatment)
  mean.trait<-c()
  stdev.trait<-c()
  coef.var<-c()
  
  for (t in 1:length(i.treat)) {
    treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
    m2<-c()
    stdev2<-c()
    cv2<-c()
    for(d in t.points) {
      m<-mean(treatment.pheno[,d], na.rm =T)
      stdev<-sqrt(var(treatment.pheno[,d], na.rm =T))
      cv<-stdev/m
      m2<-c(m2,m)
      stdev2<-c(stdev2, stdev)
      cv2<-c(cv2,cv)
    }
    m2<-c(as.character(i.treat[t]), m2)
    stdev2<-c(as.character(i.treat[t]), stdev2)
    cv2<-c(as.character(i.treat[t]), cv2)
    mean.trait<-rbind(mean.trait, m2)
    stdev.trait<-rbind(stdev.trait, stdev2)
    coef.var<-rbind(coef.var, cv2)
    rownames(coef.var)<-c()
    colnames(coef.var)<-c('treatment', colnames(pheno)[t.points])
  }
  return(coef.var)
}  


