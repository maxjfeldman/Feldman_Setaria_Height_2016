##### Lets start with height_above_bound

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
height.l<-read.csv('ril_best.fit_logistic_estimates_height.csv')
height.l<-height.l[,c(1:4)]

days<-sort(unique(height.l$dap_i))

ril_height_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  data<-height.l[height.l$dap_i == day, ]
  colnames(data)[4]<-paste('height', day, sep="_")
  data<-data[,c(1,2,4)]
  if (d == 1) {
    ril_height_qtl<-data
  }
  
  if (d > 1) {
    ril_height_qtl<-merge(ril_height_qtl, data, by=c('ril', 'treatment'), all=T)
  }
}

# Add misc columns
ril_height_qtl$Obs<-c(1:nrow(ril_height_qtl))
ril_height_qtl$experiment<-rep("BP14", nrow(ril_height_qtl))
ril_height_qtl$year<-rep("2014", nrow(ril_height_qtl))
ril_height_qtl$plot<-rep("bellweater", nrow(ril_height_qtl))
ril_height_qtl$plot_id<-rep("bellweater", nrow(ril_height_qtl))
ril_height_qtl$sampling<-rep("bellweater", nrow(ril_height_qtl))

ril_height_qtl<-ril_height_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_height_qtl)[7]<-c("id")

write.csv(ril_height_qtl, file="ril_height_above_bound_qtl.csv", quote=F, row.names=F)

##### Now centroid_y

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/centroid_y")
height.l<-read.csv('ril_best.fit_logistic_estimates_height.csv')
height.l<-height.l[,c(1:4)]

days<-sort(unique(height.l$dap_i))

ril_height_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  data<-height.l[height.l$dap_i == day, ]
  colnames(data)[4]<-paste('height', day, sep="_")
  data<-data[,c(1,2,4)]
  if (d == 1) {
    ril_height_qtl<-data
  }
  
  if (d > 1) {
    ril_height_qtl<-merge(ril_height_qtl, data, by=c('ril', 'treatment'), all=T)
  }
}

# Add misc columns
ril_height_qtl$Obs<-c(1:nrow(ril_height_qtl))
ril_height_qtl$experiment<-rep("BP14", nrow(ril_height_qtl))
ril_height_qtl$year<-rep("2014", nrow(ril_height_qtl))
ril_height_qtl$plot<-rep("bellweater", nrow(ril_height_qtl))
ril_height_qtl$plot_id<-rep("bellweater", nrow(ril_height_qtl))
ril_height_qtl$sampling<-rep("bellweater", nrow(ril_height_qtl))

ril_height_qtl<-ril_height_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_height_qtl)[7]<-c("id")

write.csv(ril_height_qtl, file="ril_centroid_y_qtl.csv", quote=F, row.names=F)


##### Now canopy_height

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/canopy_height")
height.l<-read.csv('ril_loess_estimates_canopy_height.csv')
height.l<-height.l[,c(1:4)]

days<-sort(unique(height.l$dap_i))

ril_canopy_height_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  data<-height.l[height.l$dap_i == day, ]
  colnames(data)[4]<-paste('canopy_height', day, sep="_")
  data<-data[,c(1,2,4)]
  if (d == 1) {
    ril_canopy_height_qtl<-data
  }
  
  if (d > 1) {
    ril_canopy_height_qtl<-merge(ril_canopy_height_qtl, data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_canopy_height_qtl$Obs<-c(1:nrow(ril_canopy_height_qtl))
ril_canopy_height_qtl$experiment<-rep("BP14", nrow(ril_canopy_height_qtl))
ril_canopy_height_qtl$year<-rep("2014", nrow(ril_canopy_height_qtl))
ril_canopy_height_qtl$plot<-rep("bellweater", nrow(ril_canopy_height_qtl))
ril_canopy_height_qtl$plot_id<-rep("bellweater", nrow(ril_canopy_height_qtl))
ril_canopy_height_qtl$sampling<-rep("bellweater", nrow(ril_canopy_height_qtl))

ril_canopy_height_qtl<-ril_canopy_height_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_canopy_height_qtl)[7]<-c("id")

write.csv(ril_canopy_height_qtl, file="ril_loess_estimates_canopy_height_qtl.csv", quote=F, row.names=F)



##### Now canopy_width

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/canopy_height")
height.l<-read.csv('ril_loess_estimates_canopy_width.csv')
height.l<-height.l[,c(1:4)]

days<-sort(unique(height.l$dap_i))

ril_canopy_width_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  data<-height.l[height.l$dap_i == day, ]
  colnames(data)[4]<-paste('canopy_width', day, sep="_")
  data<-data[,c(1,2,4)]
  if (d == 1) {
    ril_canopy_width_qtl<-data
  }
  
  if (d > 1) {
    ril_canopy_width_qtl<-merge(ril_canopy_width_qtl, data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_canopy_width_qtl$Obs<-c(1:nrow(ril_canopy_width_qtl))
ril_canopy_width_qtl$experiment<-rep("BP14", nrow(ril_canopy_width_qtl))
ril_canopy_width_qtl$year<-rep("2014", nrow(ril_canopy_width_qtl))
ril_canopy_width_qtl$plot<-rep("bellweater", nrow(ril_canopy_width_qtl))
ril_canopy_width_qtl$plot_id<-rep("bellweater", nrow(ril_canopy_width_qtl))
ril_canopy_width_qtl$sampling<-rep("bellweater", nrow(ril_canopy_width_qtl))

ril_canopy_width_qtl<-ril_canopy_width_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_canopy_width_qtl)[7]<-c("id")

write.csv(ril_canopy_width_qtl, file="ril_loess_estimates_canopy_width_qtl.csv", quote=F, row.names=F)


##### Now canopy_ratio

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/canopy_height")
height.l<-read.csv('ril_loess_estimates_canopy_ratio.csv')
height.l<-height.l[,c(1:4)]

days<-sort(unique(height.l$dap_i))

ril_canopy_ratio_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  data<-height.l[height.l$dap_i == day, ]
  colnames(data)[4]<-paste('canopy_ratio', day, sep="_")
  data<-data[,c(1,2,4)]
  if (d == 1) {
    ril_canopy_ratio_qtl<-data
  }
  
  if (d > 1) {
    ril_canopy_ratio_qtl<-merge(ril_canopy_ratio_qtl, data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_canopy_ratio_qtl$Obs<-c(1:nrow(ril_canopy_ratio_qtl))
ril_canopy_ratio_qtl$experiment<-rep("BP14", nrow(ril_canopy_ratio_qtl))
ril_canopy_ratio_qtl$year<-rep("2014", nrow(ril_canopy_ratio_qtl))
ril_canopy_ratio_qtl$plot<-rep("bellweater", nrow(ril_canopy_ratio_qtl))
ril_canopy_ratio_qtl$plot_id<-rep("bellweater", nrow(ril_canopy_ratio_qtl))
ril_canopy_ratio_qtl$sampling<-rep("bellweater", nrow(ril_canopy_ratio_qtl))

ril_canopy_ratio_qtl<-ril_canopy_ratio_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_canopy_ratio_qtl)[7]<-c("id")

write.csv(ril_canopy_ratio_qtl, file="ril_loess_estimates_canopy_ratio_qtl.csv", quote=F, row.names=F)






