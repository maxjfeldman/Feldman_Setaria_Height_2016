
library(lme4)
library(ggplot2)
library(lattice)

setwd('/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril')
# Read in  RIL data
ril<-read.csv('sv.vis_data_ril.csv')
ril<-ril[ril$in_bounds == "True",]

# Lets remove all columns that arn't to do with height
ril<-subset(ril, select=-c(camera, imgtype, area, hull.area, solidity, perimeter, width, height, longest_axis, center.of.mass.x, hull_vertices, ellipse_center_x, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity, y.position, height_below_bound,above_bound_area,percent_above_bound_area, below_bound_area,percent_below_bound_area,date,dap,extent_x))

ril<-ril[,c(4,2,13,5,12,3,12,1,11,15,6,8,9,10,7)]

# Calculate correlation matrix
ril.cor<-cor(ril[,9:13])
write.csv(ril.cor, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.csv")
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.pdf")
levelplot(ril.cor, scales=list(x=list(rot=90)), ylab="", xlab="", at=seq(0,1,0.01))
dev.off()


# Calculate median, mean, max, min
ril.median<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), median)
ril.mean<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), mean)
ril.max<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), max)
ril.min<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), min)
ril.var<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), var)

colnames(ril.median)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')
colnames(ril.mean)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')
colnames(ril.max)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')
colnames(ril.min)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')
colnames(ril.var)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')

# Lets calculate correlation between height types, write them out and make some plots
ril.cor.mean<-cor(ril.mean[,6:10])
write.csv(ril.cor.mean, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.mean.csv")
# write.csv(ril.cor.mean, file="/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table/Table_S1.csv", quote=F)
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.mean.pdf")
levelplot(ril.cor.mean, scales=list(x=list(rot=90)), ylab="", xlab="", at=seq(0,1,0.01))
dev.off()
ril.cor.median<-cor(ril.median[,6:10])
write.csv(ril.cor.median, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.median.csv")
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.median.pdf")
levelplot(ril.cor.median, scales=list(x=list(rot=90)), ylab="", xlab="", at=seq(0,1,0.01))
dev.off()
ril.cor.max<-cor(ril.max[,6:10])
write.csv(ril.cor.max, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.max.csv")
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.max.pdf")
levelplot(ril.cor.max, scales=list(x=list(rot=90)), ylab="", xlab="", at=seq(0,1,0.01))
dev.off()
ril.cor.min<-cor(ril.min[,6:10])
write.csv(ril.cor.min, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.min.csv")
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/trait_correlation/ril/ril.cor.min.pdf")
levelplot(ril.cor.min, scales=list(x=list(rot=90)), ylab="", xlab="", at=seq(0,1,0.01))
dev.off()

# Lets make a data frame with at each name
ril_day<-sort(unique(ril.mean$dap_i))
ril_day<-ril_day[seq(1,length(ril_day),by=2)]
height_type<-colnames(ril.mean)[6:10]
ril_mean_day<-c()
temp3<-c()
for(h in 6:ncol(ril.mean)) {
  temp<-ril.mean[,c(1,2,3,4,5,h)]
  for (i in 1:(length(ril_day))) {
    day<-ril_day[i]
    temp1<-temp[temp$dap_i == day | temp$dap_i == (day + 1),]
    temp2<-temp1[,c(1:4,6)]
    colnames(temp2)[5]<-day
    if (i == 1) {
      temp3<-temp2
    }
    if (i > 1){
      temp3<-merge(temp3, temp2, by=c('plant_id', 'cartag', 'genotype', 'treatment'), all=T)
    }
  }
  type_col<-rep(colnames(temp1)[6], nrow(temp3))
  temp3$height_type<-type_col
  ril_mean_day<-rbind(ril_mean_day, temp3)
}

ril_mean_day<-ril_mean_day[,c(1:4,18,5:17)]



############################ MEAN ############################ 

#### Now partition just genetic variance in each treatment
# Calculate heritability based on MEAN
ril_mean_H2<-get_h2_bellweather_ril(ril_mean_day)
write.csv(ril_mean_H2, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/Table_S2.csv", quote=F, row.names=T)

# Checked mean... looks like center of mass y is the most heritable measure
# *** This is also known as Figure_S2a ***
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_type_mean_H2.pdf")
plot(ril_mean_H2[1,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="H2 (mean)", ylab="H2", xlab="Days After Planting")
points(ril_mean_H2[3,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_mean_H2[5,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_mean_H2[7,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_mean_H2[9,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Calculate total variance based on MEAN
ril_mean_tot_var<-get_total_var_bellweather_ril(ril_mean_day)
write.csv(ril_mean_tot_var, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/Table_S3.csv", quote=F, row.names=T)

# Lets plot these and write out to file
# *** This is also known as Figure_S2b ***
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_mean.pdf")
plot(ril_mean_tot_var[2,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment (mean)", ylab="% of Variance", xlab="Days After Planting")
points(ril_mean_tot_var[6,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_mean_tot_var[10,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_mean_tot_var[14,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_mean_tot_var[18,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# *** This is also known as Figure_S2c ***
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_mean.pdf")
plot(ril_mean_tot_var[3,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment (mean)", ylab="% of Variance", xlab="Days After Planting")
points(ril_mean_tot_var[7,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_mean_tot_var[11,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_mean_tot_var[15,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_mean_tot_var[19,2:ncol(ril_mean_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

############################ MEDIAN ############################ 
# Lets make a data frame with at each name
ril_day<-sort(unique(ril.median$dap_i))
ril_day<-ril_day[seq(1,length(ril_day),by=2)]
height_type<-colnames(ril.median)[6:10]
ril_median_day<-c()
temp3<-c()
for(h in 6:ncol(ril.median)) {
  temp<-ril.median[,c(1,2,3,4,5,h)]
  for (i in 1:(length(ril_day))) {
    day<-ril_day[i]
    temp1<-temp[temp$dap_i == day | temp$dap_i == (day + 1),]
    temp2<-temp1[,c(1:4,6)]
    colnames(temp2)[5]<-day
    if (i == 1) {
      temp3<-temp2
    }
    if (i > 1){
      temp3<-merge(temp3, temp2, by=c('plant_id', 'cartag', 'genotype', 'treatment'), all=T)
    }
  }
  type_col<-rep(colnames(temp1)[6], nrow(temp3))
  temp3$height_type<-type_col
  ril_median_day<-rbind(ril_median_day, temp3)
}

ril_median_day<-ril_median_day[,c(1:4,18,5:17)]

#### Now partition just genetic variance in each treatment
# Calculate heritability based on median
ril_median_H2<-get_h2_bellweather_ril(ril_median_day)

# Checked median... looks like center of mass y is the most heritable measure
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_type_median_H2.pdf")
plot(ril_median_H2[1,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="H2 (median)", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[3,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_median_H2[5,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_median_H2[7,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_median_H2[9,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Calculate total variance based on median
ril_median_tot_var<-get_total_var_bellweather_ril(ril_median_day)

# Lets plot these and write out to file
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_median.pdf")
plot(ril_median_tot_var[2,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment (median)", ylab="% of Variance", xlab="Days After Planting")
points(ril_median_tot_var[6,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_median_tot_var[10,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_median_tot_var[14,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_median_tot_var[18,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_median.pdf")
plot(ril_median_tot_var[3,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment (median)", ylab="% of Variance", xlab="Days After Planting")
points(ril_median_tot_var[7,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_median_tot_var[11,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_median_tot_var[15,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_median_tot_var[19,2:ncol(ril_median_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


############################ MAX ############################ 
# Lets make a data frame with at each name
ril_day<-sort(unique(ril.max$dap_i))
ril_day<-ril_day[seq(1,length(ril_day),by=2)]
height_type<-colnames(ril.max)[6:10]
ril_max_day<-c()
temp3<-c()
for(h in 6:ncol(ril.max)) {
  temp<-ril.max[,c(1,2,3,4,5,h)]
  for (i in 1:(length(ril_day))) {
    day<-ril_day[i]
    temp1<-temp[temp$dap_i == day | temp$dap_i == (day + 1),]
    temp2<-temp1[,c(1:4,6)]
    colnames(temp2)[5]<-day
    if (i == 1) {
      temp3<-temp2
    }
    if (i > 1){
      temp3<-merge(temp3, temp2, by=c('plant_id', 'cartag', 'genotype', 'treatment'), all=T)
    }
  }
  type_col<-rep(colnames(temp1)[6], nrow(temp3))
  temp3$height_type<-type_col
  ril_max_day<-rbind(ril_max_day, temp3)
}

ril_max_day<-ril_max_day[,c(1:4,18,5:17)]

#### Now partition just genetic variance in each treatment
# Calculate heritability based on max
ril_max_H2<-get_h2_bellweather_ril(ril_max_day)

# Checked max... looks like center of mass y is the most heritable measure
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_type_max_H2.pdf")
plot(ril_max_H2[1,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="H2 (max)", ylab="H2", xlab="Days After Planting")
points(ril_max_H2[3,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_max_H2[5,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_max_H2[7,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_max_H2[9,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Calculate total variance based on max
ril_max_tot_var<-get_total_var_bellweather_ril(ril_max_day)

# Lets plot these and write out to file
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_max.pdf")
plot(ril_max_tot_var[2,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment (max)", ylab="% of Variance", xlab="Days After Planting")
points(ril_max_tot_var[6,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_max_tot_var[10,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_max_tot_var[14,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_max_tot_var[18,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_max.pdf")
plot(ril_max_tot_var[3,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment (max)", ylab="% of Variance", xlab="Days After Planting")
points(ril_max_tot_var[7,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_max_tot_var[11,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_max_tot_var[15,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_max_tot_var[19,2:ncol(ril_max_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


############################ MIN ############################ 
# Lets make a data frame with at each name
ril_day<-sort(unique(ril.min$dap_i))
ril_day<-ril_day[seq(1,length(ril_day),by=2)]
height_type<-colnames(ril.min)[6:10]
ril_min_day<-c()
temp3<-c()
for(h in 6:ncol(ril.min)) {
  temp<-ril.min[,c(1,2,3,4,5,h)]
  for (i in 1:(length(ril_day))) {
    day<-ril_day[i]
    temp1<-temp[temp$dap_i == day | temp$dap_i == (day + 1),]
    temp2<-temp1[,c(1:4,6)]
    colnames(temp2)[5]<-day
    if (i == 1) {
      temp3<-temp2
    }
    if (i > 1){
      temp3<-merge(temp3, temp2, by=c('plant_id', 'cartag', 'genotype', 'treatment'), all=T)
    }
  }
  type_col<-rep(colnames(temp1)[6], nrow(temp3))
  temp3$height_type<-type_col
  ril_min_day<-rbind(ril_min_day, temp3)
}

ril_min_day<-ril_min_day[,c(1:4,18,5:17)]

#### Now partition just genetic variance in each treatment
# Calculate heritability based on min
ril_min_H2<-get_h2_bellweather_ril(ril_min_day)

# Checked min... looks like center of mass y is the most heritable measure
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_type_min_H2.pdf")
plot(ril_min_H2[1,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="H2 (min)", ylab="H2", xlab="Days After Planting")
points(ril_min_H2[3,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_min_H2[5,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_min_H2[7,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_min_H2[9,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Calculate total variance based on min
ril_min_tot_var<-get_total_var_bellweather_ril(ril_min_day)

# Lets plot these and write out to file
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_min.pdf")
plot(ril_min_tot_var[2,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment (min)", ylab="% of Variance", xlab="Days After Planting")
points(ril_min_tot_var[6,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_min_tot_var[10,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_min_tot_var[14,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_min_tot_var[18,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_min.pdf")
plot(ril_min_tot_var[3,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment (min)", ylab="% of Variance", xlab="Days After Planting")
points(ril_min_tot_var[7,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_min_tot_var[11,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_min_tot_var[15,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_min_tot_var[19,2:ncol(ril_min_tot_var)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Lets see which mean, median, max, min is the most heritable

########## Height above bound
# *** This is also known as Figure_S3 ***
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table/Figure_S3.pdf")
plot(ril_mean_H2[1,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="height_above_bound", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[1,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[1,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[1,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

########## extent_y

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/ril/height_types/ril_extent_y_summary_stat_H2.pdf")
plot(ril_mean_H2[3,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="extent_y", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[3,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[3,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[3,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

########## centroid_y

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/ril/height_types/ril_centroid_y_summary_stat_H2.pdf")
plot(ril_mean_H2[5,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="centroid_y", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[5,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[5,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[5,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

########## ellipse_center_y

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/ril/height_types/ril_ellipse_center_y_summary_stat_H2.pdf")
plot(ril_mean_H2[7,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="ellipse_center_y", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[1,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[7,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[7,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

########## canopy_height

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/ril/height_types/ril_canopy_height_summary_stat_H2.pdf")
plot(ril_mean_H2[9,2:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="canopy_height", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[9,2:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[9,2:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[9,2:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()




############################################################################################
# Calculate a pooled average for every 2 day block, and write it out in QTL format
############################################################################################

# Now lets reformat bellweather to qtl format data to calculate CV
ril_height_qtl<-c()
days<-sort(unique(ril.mean$dap_i))
days_e<-days[seq(1,length(days),by=2)]

for(d in 1:length(days_e)) {
  day<-days_e[d]
  data<-ril.mean[ril.mean$dap_i == day | ril.mean$dap == (day + 1), ]
  colnames(data)[6]<-paste('height', day, sep="_")
  data<-data[,c(1,3,4,6)]
  if (d == 1) {
    ril_height_qtl<-data
  }
  
  if (d > 1) {
    ril_height_qtl<-merge(ril_height_qtl, data, by=c('plant_id','genotype', 'treatment'), all=T)
  }
}

ril_height_qtl$Obs<-c(1:nrow(ril_height_qtl))
ril_height_qtl$experiment<-rep('BP14', nrow(ril_height_qtl))
ril_height_qtl$year<-rep('2014', nrow(ril_height_qtl))
ril_height_qtl$plot_id<-rep('phenotyper', nrow(ril_height_qtl))
ril_height_qtl$subplot_id<-rep('phenotyper', nrow(ril_height_qtl))

ril_height_qtl<-ril_height_qtl[,c(17,18,19,3,20,21,2,1,4:16)]
colnames(ril_height_qtl)[7]<-c("id")
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril")
write.csv(ril_height_qtl, 'ril_height_qtl_raw.csv', quote=F, row.names=F)


############################################################################################
# Define fxns
############################################################################################

get_h2_bellweather_ril<-function(data){
  h.types<-unique(data$height_type)
  i.treat<-unique(as.character(data$treatment))
  variance.type<-c()
  for(j in 1:length(h.types)){
    H2<-c()
    e2<-c()
    type.pheno<-data[data$height_type == h.types[j],]
    
    for (i in 6:length(colnames(type.pheno))){
      # Get complete cases
      cc.pheno<-type.pheno[complete.cases(type.pheno[,i]),c(1:5,i)]
    
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.pheno[,6]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.pheno)
    
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
    
      # Extract individual components (order will remain the same)
      gxt.var<-re[1]
      geno.var<-re[2]
      treat.var<-re[3]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
    
      reps.t1<-table(cc.pheno[cc.pheno$treatment == i.treat[1], 'genotype'])
      reps.t2<-table(cc.pheno[cc.pheno$treatment == i.treat[2], 'genotype'])
      reps.treatment<-c(reps.t1, reps.t2)
    
      reps.t1<-as.character(unique(cc.pheno[cc.pheno$treatment == i.treat[1], 'genotype']))
      reps.t2<-as.character(unique(cc.pheno[cc.pheno$treatment == i.treat[2], 'genotype']))
      unique.combined <- c(as.character(reps.t1), as.character(reps.t2))
    
      freq.unique.combined <- table(unique.combined)
    
      # Calculate the harmonic mean replication within treatment blocks
      hm_treatment<-harmonic.mean(freq.unique.combined)$harmean
    
      # Now get a count of total genotypic replication
      reps.total<-table(cc.pheno[,'genotype'])
      # Get the harmonic mean of this quantity
      hm_total<-harmonic.mean(reps.total)$harmean
    
      # Calculate heritability as described by AEL 
      # H2 = geno.var/(geno.var + (gxt.var/harmonic mean of treatment block replication) + (residual.var/harmonic mean of total genotype replication) )
      h2<-((geno.var)/(geno.var + (gxt.var/hm_treatment) + (res/hm_total)))
      e<-1-h2
      # This is the heritability
      H2<-c(H2,h2)
      e2<-c(e2, e)
    
    }
    names(H2)<-colnames(type.pheno)[6:ncol(type.pheno)]
    names(e2)<-colnames(type.pheno)[6:ncol(type.pheno)]
    variance.t<-rbind(H2, e2)
    v.type<-rep(as.character(h.types[j]), 2)
    #v.treat<-rep(as.character(i.treat[t]), 2)
    #colnames(treatment.pheno)[6:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[6:length(treatment.pheno)], sep="_")
    variance.t<-cbind(v.type, variance.t)
    colnames(variance.t)<-c('height_type', colnames(type.pheno)[6:length(type.pheno)])
    rownames(variance.t)<-c('Genotype', 'Error')
    variance.type<-rbind(variance.type, variance.t)
  }
  return(variance.type)
}

  
get_total_var_bellweather_ril<-function(data){
  h.types<-unique(data$height_type)
  i.treat<-unique(as.character(data$treatment))
  variance.type<-c()
  for(j in 1:length(h.types)){
    H2<-c()
    t2<-c()
    gxt2<-c()
    e2<-c()
    type.pheno<-data[data$height_type == h.types[j],]
    
    for (i in 6:length(colnames(type.pheno))){
      # Get complete cases
      cc.pheno<-type.pheno[complete.cases(type.pheno[,i]),c(1:5,i)]
      
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.pheno[,6]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.pheno)
      
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      
      # Extract individual components (order will remain the same)
      gxt.var<-re[1]
      geno.var<-re[2]
      treat.var<-re[3]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      
      # Get proportion of variance for all factors
      h<-geno.var/tot.var
      t<-treat.var/tot.var
      e<-res/tot.var
      gxt<-gxt.var/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      t2<-c(t2,t)
      e2<-c(e2,e)
      gxt2<-c(gxt2, gxt)
      
    }
    names(H2)<-colnames(type.pheno)[6:ncol(type.pheno)]
    names(t2)<-colnames(type.pheno)[6:ncol(type.pheno)]
    names(gxt2)<-colnames(type.pheno)[6:ncol(type.pheno)]
    names(e2)<-colnames(type.pheno)[6:ncol(type.pheno)]
    variance.t<-rbind(H2,t2,gxt2,e2)
    v.type<-rep(as.character(h.types[j]), 4)
    variance.t<-cbind(v.type, variance.t)
    colnames(variance.t)<-c('height_type', colnames(type.pheno)[6:length(type.pheno)])
    rownames(variance.t)<-c('Genotype', 'Treatment', "G X Treatment", 'Error')
    variance.type<-rbind(variance.type, variance.t)
  }
  return(variance.type)
}

  
  