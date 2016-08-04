# Script to analyze heritability and total variance from LemnaTec Setaria RIL dataset

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

# Get vector of individual treatments and height types
i.treat<-unique(ril_mean_day$treatment)
h.types<-unique(ril_mean_day$height_type)
variance.treat<-c()
ril_mean_H2<-c()

for (j in 1:length(h.types)) {
  #for (j in 2:2) {
  type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
  #H2<-c()
  #e2<-c()
  variance.treat<-c()
  for (t in 1:length(i.treat)) {
    # Create variables to store values
    treatment.pheno<-type.pheno[type.pheno$treatment == i.treat[t],]
    variance.t<-c()
    H2<-c()
    e2<-c()
    
    # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
    for(d in 6:length(colnames(treatment.pheno))){
      # Use only RILs with all measurements for each treatment.phenotype
      cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,d]),c(1:3,d)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.treatment.pheno[,4]~(1|genotype), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
    
    variance.t<-rbind(H2, e2)
    v.type<-rep(as.character(h.types[j]), 2)
    v.treat<-rep(as.character(i.treat[t]), 2)
    #colnames(treatment.pheno)[6:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[6:length(treatment.pheno)], sep="_")
    variance.t<-cbind(v.type, v.treat, variance.t)
    colnames(variance.t)<-c('height_type', 'treatment', colnames(treatment.pheno)[6:length(treatment.pheno)])
    rownames(variance.t)<-c('Genotype', 'Error')
    variance.treat<-rbind(variance.treat, variance.t)
  }
  ril_mean_H2<-rbind(ril_mean_H2, variance.treat)
}  

# Checked mean... looks like center of mass y is the most heritable measure
# Wet
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_height_type_mean_H2.pdf")
plot(ril_mean_H2[1,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="well_watered_mean", ylab="H2", xlab="Days After Planting")
points(ril_mean_H2[5,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_mean_H2[9,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_mean_H2[13,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_mean_H2[17,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)

dev.off()

# Dry
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_height_type_mean_H2.pdf")
plot(ril_mean_H2[3,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="drought_mean", ylab="H2", xlab="Days After Planting")
points(ril_mean_H2[7,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_mean_H2[11,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_mean_H2[15,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_mean_H2[19,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

write.csv(ril_mean_H2, "/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_mean_H2.csv", quote=F)

##########################################################################

# Get vector of individual treatments and height types
i.treat<-unique(ril_mean_day$treatment)
h.types<-unique(ril_mean_day$height_type)

variance_final<-c()
variance_final_ne<-c()
for (j in 1:length(h.types)) {
  type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
  variance<-c()
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  # Now one with all environmental variance subtracted
  H2_ne<-c()
  e2_ne<-c()
  for(d in 6:length(colnames(type.pheno))){
    cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(3:4,d)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.type.pheno[,3]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
    
    # Now only genetic and error
    # Now get proportion of variance - environment fx
    h_ne<-geno.var/(tot.var - treat.var - gxt.var)
    e_ne<-res/(tot.var - treat.var - gxt.var)
    # Append variables to a vector of variables
    H2_ne<-c(H2_ne,h_ne)
    e2_ne<-c(e2_ne,e_ne)
  }
  
  v.type<-rep(as.character(h.types[j]), 4)
  variance<-rbind(H2, t2, gxt2, e2)
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  variance.h<-cbind(v.type, variance)
  colnames(variance.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final<-rbind(variance_final, variance.h)
  v.type_ne<-rep(as.character(h.types[j]), 2)
  variance_ne<-rbind(H2_ne, e2_ne)
  rownames(variance_ne)<-c('Genotype', 'Error')
  variance_ne.h<-cbind(v.type_ne, variance_ne)
  colnames(variance_ne.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final_ne<-rbind(variance_final_ne, variance_ne.h)
}


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/H2_fx_height_type_mean.pdf")
plot(variance_final_ne[1,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="Heritability", ylab="H2", xlab="Days After Planting")
points(variance_final_ne[3,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final_ne[5,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final_ne[7,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final_ne[9,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_mean.pdf")
plot(variance_final[2,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[6,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[10,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[14,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[18,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_mean.pdf")
plot(variance_final[3,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[7,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[11,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[15,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[19,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


##### Lets also write these out to manuscripts file:
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table/Figure_S2a.pdf")
plot(variance_final_ne[1,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="Heritability", ylab="H2", xlab="Days After Planting")
points(variance_final_ne[3,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final_ne[5,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final_ne[7,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final_ne[9,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table/Figure_S2b.pdf")
plot(variance_final[2,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.25), col=c('green'), main="Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[6,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[10,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[14,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[18,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('topright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table/Figure_S2c.pdf")
plot(variance_final[3,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[7,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[11,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[15,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[19,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('topright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

variance_final_mean<-variance_final
varience_final_mean_ne<-variance_final_ne

# Looks like center of mass y is the most heritable measure
# Lets look at how this compares with median, max, min etc...

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
# Calculate heritability

# Get vector of individual treatments and height types
i.treat<-unique(ril_median_day$treatment)
h.types<-unique(ril_median_day$height_type)
variance.treat<-c()
ril_median_H2<-c()

for (j in 1:length(h.types)) {
  #for (j in 2:2) {
  type.pheno<-ril_median_day[ril_median_day$height_type == h.types[j],]
  #H2<-c()
  #e2<-c()
  variance.treat<-c()
  for (t in 1:length(i.treat)) {
    # Create variables to store values
    treatment.pheno<-type.pheno[type.pheno$treatment == i.treat[t],]
    variance.t<-c()
    H2<-c()
    e2<-c()
    
    # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
    for(d in 6:length(colnames(treatment.pheno))){
      # Use only RILs with all measurements for each treatment.phenotype
      cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,d]),c(1:3,d)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.treatment.pheno[,4]~(1|genotype), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
    
    variance.t<-rbind(H2, e2)
    v.type<-rep(as.character(h.types[j]), 2)
    v.treat<-rep(as.character(i.treat[t]), 2)
    #colnames(treatment.pheno)[6:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[6:length(treatment.pheno)], sep="_")
    variance.t<-cbind(v.type, v.treat, variance.t)
    colnames(variance.t)<-c('height_type', 'treatment', colnames(treatment.pheno)[6:length(treatment.pheno)])
    rownames(variance.t)<-c('Genotype', 'Error')
    variance.treat<-rbind(variance.treat, variance.t)
  }
  ril_median_H2<-rbind(ril_median_H2, variance.treat)
}  

# Checked median... looks like center of mass y is the most heritable measure
# Checked mean... looks like center of mass y is the most heritable measure
# Wet
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_height_type_median_H2.pdf")
plot(ril_median_H2[1,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="well_watered_median", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[5,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_median_H2[9,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_median_H2[13,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_median_H2[17,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Dry
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_height_type_median_H2.pdf")
plot(ril_median_H2[3,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="drought_median", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[7,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_median_H2[11,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_median_H2[15,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_median_H2[19,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

write.csv(ril_median_H2, "/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_median_H2.csv", quote=F)

##########################################################################

# Get vector of individual treatments and height types
i.treat<-unique(ril_median_day$treatment)
h.types<-unique(ril_median_day$height_type)

variance_final<-c()
variance_final_ne<-c()
for (j in 1:length(h.types)) {
  type.pheno<-ril_median_day[ril_median_day$height_type == h.types[j],]
  variance<-c()
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  # Now one with all environmental variance subtracted
  H2_ne<-c()
  e2_ne<-c()
  for(d in 6:length(colnames(type.pheno))){
    cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(3:4,d)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.type.pheno[,3]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
    
    # Now only genetic and error
    # Now get proportion of variance - environment fx
    h_ne<-geno.var/(tot.var - treat.var - gxt.var)
    e_ne<-res/(tot.var - treat.var - gxt.var)
    # Append variables to a vector of variables
    H2_ne<-c(H2_ne,h_ne)
    e2_ne<-c(e2_ne,e_ne)
  }
  
  v.type<-rep(as.character(h.types[j]), 4)
  variance<-rbind(H2, t2, gxt2, e2)
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  variance.h<-cbind(v.type, variance)
  colnames(variance.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final<-rbind(variance_final, variance.h)
  v.type_ne<-rep(as.character(h.types[j]), 2)
  variance_ne<-rbind(H2_ne, e2_ne)
  rownames(variance_ne)<-c('Genotype', 'Error')
  variance_ne.h<-cbind(v.type_ne, variance_ne)
  colnames(variance_ne.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final_ne<-rbind(variance_final_ne, variance_ne.h)
}


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/H2_fx_height_type_median.pdf")
plot(variance_final_ne[1,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="Heritability", ylab="H2", xlab="Days After Planting")
points(variance_final_ne[3,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final_ne[5,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final_ne[7,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final_ne[9,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_median.pdf")
plot(variance_final[2,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[6,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[10,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[14,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[18,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_median.pdf")
plot(variance_final[3,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[7,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[11,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[15,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[19,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

variance_final_median<-variance_final
varience_final_median_ne<-variance_final_ne

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
# Calculate heritability

# Get vector of individual treatments and height types
i.treat<-unique(ril_max_day$treatment)
h.types<-unique(ril_max_day$height_type)
variance.treat<-c()
ril_max_H2<-c()

for (j in 1:length(h.types)) {
  #for (j in 2:2) {
  type.pheno<-ril_max_day[ril_max_day$height_type == h.types[j],]
  #H2<-c()
  #e2<-c()
  variance.treat<-c()
  for (t in 1:length(i.treat)) {
    # Create variables to store values
    treatment.pheno<-type.pheno[type.pheno$treatment == i.treat[t],]
    variance.t<-c()
    H2<-c()
    e2<-c()
    
    # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
    for(d in 6:length(colnames(treatment.pheno))){
      # Use only RILs with all measurements for each treatment.phenotype
      cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,d]),c(1:3,d)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.treatment.pheno[,4]~(1|genotype), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
    
    variance.t<-rbind(H2, e2)
    v.type<-rep(as.character(h.types[j]), 2)
    v.treat<-rep(as.character(i.treat[t]), 2)
    #colnames(treatment.pheno)[6:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[6:length(treatment.pheno)], sep="_")
    variance.t<-cbind(v.type, v.treat, variance.t)
    colnames(variance.t)<-c('height_type', 'treatment', colnames(treatment.pheno)[6:length(treatment.pheno)])
    rownames(variance.t)<-c('Genotype', 'Error')
    variance.treat<-rbind(variance.treat, variance.t)
  }
  ril_max_H2<-rbind(ril_max_H2, variance.treat)
}  

# Checked max... looks like center of mass y is the most heritable measure
# Wet
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_height_type_max_H2.pdf")
plot(ril_max_H2[1,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="well_watered_max", ylab="H2", xlab="Days After Planting")
points(ril_max_H2[5,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_max_H2[9,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_max_H2[13,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_max_H2[17,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Dry
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_height_type_max_H2.pdf")
plot(ril_max_H2[3,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="drought_max", ylab="H2", xlab="Days After Planting")
points(ril_max_H2[7,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_max_H2[11,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_max_H2[15,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_max_H2[19,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

write.csv(ril_max_H2, "/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_max_H2.csv", quote=F)

##########################################################################

# Get vector of individual treatments and height types
i.treat<-unique(ril_max_day$treatment)
h.types<-unique(ril_max_day$height_type)

variance_final<-c()
variance_final_ne<-c()
for (j in 1:length(h.types)) {
  type.pheno<-ril_max_day[ril_max_day$height_type == h.types[j],]
  variance<-c()
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  # Now one with all environmental variance subtracted
  H2_ne<-c()
  e2_ne<-c()
  for(d in 6:length(colnames(type.pheno))){
    cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(3:4,d)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.type.pheno[,3]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
    
    # Now only genetic and error
    # Now get proportion of variance - environment fx
    h_ne<-geno.var/(tot.var - treat.var - gxt.var)
    e_ne<-res/(tot.var - treat.var - gxt.var)
    # Append variables to a vector of variables
    H2_ne<-c(H2_ne,h_ne)
    e2_ne<-c(e2_ne,e_ne)
  }
  
  v.type<-rep(as.character(h.types[j]), 4)
  variance<-rbind(H2, t2, gxt2, e2)
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  variance.h<-cbind(v.type, variance)
  colnames(variance.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final<-rbind(variance_final, variance.h)
  v.type_ne<-rep(as.character(h.types[j]), 2)
  variance_ne<-rbind(H2_ne, e2_ne)
  rownames(variance_ne)<-c('Genotype', 'Error')
  variance_ne.h<-cbind(v.type_ne, variance_ne)
  colnames(variance_ne.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final_ne<-rbind(variance_final_ne, variance_ne.h)
}


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/H2_fx_height_type_max.pdf")
plot(variance_final_ne[1,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="Heritability", ylab="H2", xlab="Days After Planting")
points(variance_final_ne[3,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final_ne[5,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final_ne[7,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final_ne[9,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_max.pdf")
plot(variance_final[2,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[6,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[10,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[14,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[18,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_max.pdf")
plot(variance_final[3,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[7,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[11,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[15,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[19,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


variance_final_max<-variance_final
varience_final_max_ne<-variance_final_ne
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
# Calculate heritability

# Get vector of individual treatments and height types
i.treat<-unique(ril_min_day$treatment)
h.types<-unique(ril_min_day$height_type)
variance.treat<-c()
ril_min_H2<-c()

for (j in 1:length(h.types)) {
  #for (j in 2:2) {
  type.pheno<-ril_min_day[ril_min_day$height_type == h.types[j],]
  #H2<-c()
  #e2<-c()
  variance.treat<-c()
  for (t in 1:length(i.treat)) {
    # Create variables to store values
    treatment.pheno<-type.pheno[type.pheno$treatment == i.treat[t],]
    variance.t<-c()
    H2<-c()
    e2<-c()
    
    # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
    for(d in 6:length(colnames(treatment.pheno))){
      # Use only RILs with all measurements for each treatment.phenotype
      cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,d]),c(1:3,d)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.treatment.pheno[,4]~(1|genotype), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
    
    variance.t<-rbind(H2, e2)
    v.type<-rep(as.character(h.types[j]), 2)
    v.treat<-rep(as.character(i.treat[t]), 2)
    #colnames(treatment.pheno)[6:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[6:length(treatment.pheno)], sep="_")
    variance.t<-cbind(v.type, v.treat, variance.t)
    colnames(variance.t)<-c('height_type', 'treatment', colnames(treatment.pheno)[6:length(treatment.pheno)])
    rownames(variance.t)<-c('Genotype', 'Error')
    variance.treat<-rbind(variance.treat, variance.t)
  }
  ril_min_H2<-rbind(ril_min_H2, variance.treat)
}  

# Checked min... looks like center of mass y is the most heritable measure
# Wet
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_height_type_min_H2.pdf")
plot(ril_min_H2[1,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="well_watered_min", ylab="H2", xlab="Days After Planting")
points(ril_min_H2[5,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_min_H2[9,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_min_H2[13,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_min_H2[17,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='purple')
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

# Dry
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_height_type_min_H2.pdf")
plot(ril_min_H2[3,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="drought_min", ylab="H2", xlab="Days After Planting")
points(ril_min_H2[7,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(ril_min_H2[11,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('red'))
points(ril_min_H2[15,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(ril_min_H2[19,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

write.csv(ril_min_H2, "/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_min_H2.csv", quote=F)

##########################################################################

# Get vector of individual treatments and height types
i.treat<-unique(ril_min_day$treatment)
h.types<-unique(ril_min_day$height_type)

variance_final<-c()
variance_final_ne<-c()
for (j in 1:length(h.types)) {
  type.pheno<-ril_min_day[ril_min_day$height_type == h.types[j],]
  variance<-c()
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  # Now one with all environmental variance subtracted
  H2_ne<-c()
  e2_ne<-c()
  for(d in 6:length(colnames(type.pheno))){
    cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(3:4,d)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.type.pheno[,3]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
    
    # Now only genetic and error
    # Now get proportion of variance - environment fx
    h_ne<-geno.var/(tot.var - treat.var - gxt.var)
    e_ne<-res/(tot.var - treat.var - gxt.var)
    # Append variables to a vector of variables
    H2_ne<-c(H2_ne,h_ne)
    e2_ne<-c(e2_ne,e_ne)
  }
  
  v.type<-rep(as.character(h.types[j]), 4)
  variance<-rbind(H2, t2, gxt2, e2)
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  variance.h<-cbind(v.type, variance)
  colnames(variance.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final<-rbind(variance_final, variance.h)
  v.type_ne<-rep(as.character(h.types[j]), 2)
  variance_ne<-rbind(H2_ne, e2_ne)
  rownames(variance_ne)<-c('Genotype', 'Error')
  variance_ne.h<-cbind(v.type_ne, variance_ne)
  colnames(variance_ne.h)<-c('height_type',colnames(type.pheno)[6:length(type.pheno)] )
  variance_final_ne<-rbind(variance_final_ne, variance_ne.h)
}


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/H2_fx_height_type_min.pdf")
plot(variance_final_ne[1,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col=c('green'), main="Heritability", ylab="H2", xlab="Days After Planting")
points(variance_final_ne[3,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final_ne[5,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final_ne[7,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final_ne[9,2:ncol(variance_final_ne)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/treatment_fx_height_type_min.pdf")
plot(variance_final[2,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.3), col=c('green'), main="Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[6,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[10,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[14,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[18,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()


pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/gxt_fx_height_type_min.pdf")
plot(variance_final[3,2:ncol(variance_final)]~c(seq(8,32,2)), type="l",ylim=c(0,.15), col=c('green'), main="G X Treatment", ylab="% of Variance", xlab="Days After Planting")
points(variance_final[7,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('blue'))
points(variance_final[11,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('red'))
points(variance_final[15,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('orange'))
points(variance_final[19,2:ncol(variance_final)]~c(seq(8,32,2)), type="l", col=c('purple'))
legend('bottomright', c('Height_Above_Bound', 'Extent_Y', 'Centroid_Y', 'Ellipse_Center_Y', 'Canopy_Height'), col=c('green','blue','red','orange','purple'), lty=1, bty='n', cex=.75)
dev.off()

variance_final_min<-variance_final
varience_final_min_ne<-variance_final_ne
######################################################## 
# Lets see which mean, median, max, min is the most heritable

########## Height above bound

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table/Figure_S3.pdf")
plot(varience_final_mean_ne[1,2:ncol(varience_final_mean_ne)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="height_above_bound", ylab="H2", xlab="Days After Planting")
points(varience_final_median_ne[1,2:ncol(varience_final_median_ne)]~c(seq(8,32,2)), type="l", col='blue')
points(varience_final_max_ne[1,2:ncol(varience_final_max_ne)]~c(seq(8,32,2)), type="l", col='green')
points(varience_final_min_ne[1,2:ncol(varience_final_min_ne)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

























######################################################## 
# Looks like centroid_y is the most heritable measurment of plant height in both conditions and in (mean, median, max, min)
# Lets see which mean, median, max, min is the most heritable

########## CENTROID_Y
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_centroid_y_summary_stat_H2.pdf")
plot(ril_mean_H2[9,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="well_watered_centroid_y_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[9,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[9,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[9,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_centroid_y_summary_stat_H2.pdf")
plot(ril_mean_H2[11,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="drought_centroid_y_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[11,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[11,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[11,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()


########## HEIGHT ABOVE BOUND
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_height_above_bound_summary_stat_H2.pdf")
plot(ril_mean_H2[1,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="well_watered_height_above_bound_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[1,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[1,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[1,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_height_above_bound_summary_stat_H2.pdf")
plot(ril_mean_H2[3,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="drought_height_above_bound_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[3,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[3,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[3,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

########## CANOPY HEIGHT
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_canopy_height_summary_stat_H2.pdf")
plot(ril_mean_H2[17,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="well_watered_canopy_height_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[17,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[17,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[17,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_canopy_height_summary_stat_H2.pdf")
plot(ril_mean_H2[19,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="drought_canopy_height_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[19,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[19,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[19,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

########## Extent_Y
pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_wet_extent_y_summary_stat_H2.pdf")
plot(ril_mean_H2[5,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="well_watered_extent_y_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[5,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[5,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[5,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

pdf("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_dry_extent_y_summary_stat_H2.pdf")
plot(ril_mean_H2[7,3:ncol(ril_mean_H2)]~c(seq(8,32,2)), type="l",ylim=c(0,1), col='red', main="drought_extent_y_summary_stats", ylab="H2", xlab="Days After Planting")
points(ril_median_H2[7,3:ncol(ril_median_H2)]~c(seq(8,32,2)), type="l", col='blue')
points(ril_max_H2[7,3:ncol(ril_max_H2)]~c(seq(8,32,2)), type="l", col='green')
points(ril_min_H2[7,3:ncol(ril_min_H2)]~c(seq(8,32,2)), type="l", col='orange')
legend('bottomright', c('Mean', 'Median', 'Max', 'Min'), col=c('red','blue','green','orange'), lty=1, bty='n', cex=.75)
dev.off()

# NO DIFFERENCE

#################################################
# TOTAL VARIANCE
#################################################


##### Start with height_above_bound_mean

# Get vector of individual treatments and height types
h.types<-unique(ril_mean_day$height_type)
#ril_mean_var<-c()

#for (j in 1:length(h.types)) {
  #for (j in 2:2) {
  # Lets focus on height_above_bound
  j<-1
  type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
    # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
    for(d in 6:length(colnames(type.pheno))){
      # Use only RILs with all measurements for each type.phenotype
      cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(1,3,4,d)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.type.pheno[,4]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      gxt.var<-re[1]
      geno.var<-re[2]
      treat.var<-re[3]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
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
  variance<-rbind(H2, t2, gxt2, e2)
  #height_type<-rep(h.types[j], nrow(variance))
  #variance<-cbind(height_type,variance)
  colnames(variance)<-colnames(type.pheno)[6:length(type.pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  #her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
  #write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
  #DR13_height_var<-variance
#}  

  ril_height_above_bound_day_var<-variance
  write.csv(ril_height_above_bound_day_var, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_above_bound_variance.csv", quote=F)


  traits<-colnames(ril_height_above_bound_day_var)
  types<-rownames(ril_height_above_bound_day_var)
  
  variance.l<-c()
  t<-c()
  for(i in 1:length(traits)) {
    t<-rep(traits[i], length(types))
    rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
    variance.l<-rbind(variance.l, rf)
  }
  colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
  rownames(variance.l)<-c(1:nrow(variance.l))
  
ril_height_above_bound_day_var.l<-variance.l

pdf('/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_above_bound_variance.pdf')
p<-ggplot(ril_height_above_bound_day_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1)
p
dev.off()

write.csv(ril_height_above_bound_day_var.l, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_height_above_bound_variance.l.csv", quote=F, row.names=F)


##### Now extent_y

# Get vector of individual treatments and height types
h.types<-unique(ril_mean_day$height_type)
#ril_mean_var<-c()

#for (j in 1:length(h.types)) {
#for (j in 2:2) {
# Lets focus on height_above_bound
j<-2
type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
H2<-c()
t2<-c()
e2<-c()
gxt2<-c()
# For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
for(d in 6:length(colnames(type.pheno))){
  # Use only RILs with all measurements for each type.phenotype
  cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(1,3,4,d)]
  # Build linear model each cofactor is a random effect
  model<-lmer(cc.type.pheno[,4]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
  # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
  re<-as.numeric(VarCorr(model))
  res<-attr(VarCorr(model), "sc")^2
  # Extract individual components (order will remain the same)
  gxt.var<-re[1]
  geno.var<-re[2]
  treat.var<-re[3]
  # Total variance is sum of all variances
  tot.var<-sum(re, res)
  # Get proportion of variance
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
variance<-rbind(H2, t2, gxt2, e2)
#height_type<-rep(h.types[j], nrow(variance))
#variance<-cbind(height_type,variance)
colnames(variance)<-colnames(type.pheno)[6:length(type.pheno)]
rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
#her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
#write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
#DR13_height_var<-variance
#}  

ril_extent_y_day_var<-variance

write.csv(ril_extent_y_day_var, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_extent_y_variance.csv", quote=F)

traits<-colnames(ril_extent_y_day_var)
types<-rownames(ril_extent_y_day_var)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}
colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

ril_extent_y_day_var.l<-variance.l

pdf('/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_extent_y_variance.pdf')
p<-ggplot(ril_extent_y_day_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1)
p
dev.off()

write.csv(ril_extent_y_day_var.l, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_extent_y_variance.l.csv", quote=F, row.names=F)



##### Now centroid_y_mean

# Get vector of individual treatments and height types
h.types<-unique(ril_mean_day$height_type)
#ril_mean_var<-c()

#for (j in 1:length(h.types)) {
#for (j in 2:2) {
# Lets focus on height_above_bound
j<-3
type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
H2<-c()
t2<-c()
e2<-c()
gxt2<-c()
# For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
for(d in 6:length(colnames(type.pheno))){
  # Use only RILs with all measurements for each type.phenotype
  cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(1,3,4,d)]
  # Build linear model each cofactor is a random effect
  model<-lmer(cc.type.pheno[,4]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
  # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
  re<-as.numeric(VarCorr(model))
  res<-attr(VarCorr(model), "sc")^2
  # Extract individual components (order will remain the same)
  gxt.var<-re[1]
  geno.var<-re[2]
  treat.var<-re[3]
  # Total variance is sum of all variances
  tot.var<-sum(re, res)
  # Get proportion of variance
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
variance<-rbind(H2, t2, gxt2, e2)
#height_type<-rep(h.types[j], nrow(variance))
#variance<-cbind(height_type,variance)
colnames(variance)<-colnames(type.pheno)[6:length(type.pheno)]
rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
#her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
#write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
#DR13_height_var<-variance
#}  

ril_centroid_y_day_var<-variance
write.csv(ril_centroid_y_day_var, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_centroid_y_variance.csv", quote=F)

traits<-colnames(ril_centroid_y_day_var)
types<-rownames(ril_centroid_y_day_var)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}
colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

ril_centroid_y_day_var.l<-variance.l

pdf('/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_centroid_y_variance.pdf')
p<-ggplot(ril_centroid_y_day_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1)
p
dev.off()

write.csv(ril_centroid_y_day_var.l, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_centroid_y_variance.l.csv", quote=F, row.names=F)



##### Now ellipse_y_mean

# Get vector of individual treatments and height types
h.types<-unique(ril_mean_day$height_type)
#ril_mean_var<-c()

#for (j in 1:length(h.types)) {
#for (j in 2:2) {
# Lets focus on height_above_bound
j<-4
type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
H2<-c()
t2<-c()
e2<-c()
gxt2<-c()
# For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
for(d in 6:length(colnames(type.pheno))){
  # Use only RILs with all measurements for each type.phenotype
  cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(1,3,4,d)]
  # Build linear model each cofactor is a random effect
  model<-lmer(cc.type.pheno[,4]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
  # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
  re<-as.numeric(VarCorr(model))
  res<-attr(VarCorr(model), "sc")^2
  # Extract individual components (order will remain the same)
  gxt.var<-re[1]
  geno.var<-re[2]
  treat.var<-re[3]
  # Total variance is sum of all variances
  tot.var<-sum(re, res)
  # Get proportion of variance
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
variance<-rbind(H2, t2, gxt2, e2)
#height_type<-rep(h.types[j], nrow(variance))
#variance<-cbind(height_type,variance)
colnames(variance)<-colnames(type.pheno)[6:length(type.pheno)]
rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
#her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
#write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
#DR13_height_var<-variance
#}  

ril_ellipse_y_day_var<-variance
write.csv(ril_ellipse_y_day_var, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_ellipse_y_variance.csv", quote=F)


traits<-colnames(ril_ellipse_y_day_var)
types<-rownames(ril_ellipse_y_day_var)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}
colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

ril_ellipse_y_day_var.l<-variance.l

pdf('/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_ellipse_y_variance.pdf')
p<-ggplot(ril_ellipse_y_day_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1)
p
dev.off()
write.csv(ril_ellipse_y_day_var.l, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_ellipse_y_variance.l.csv", quote=F, row.names=F)



##### Start with canopy_height_mean

# Get vector of individual treatments and height types
h.types<-unique(ril_mean_day$height_type)
#ril_mean_var<-c()

#for (j in 1:length(h.types)) {
#for (j in 2:2) {
# Lets focus on height_above_bound
j<-5
type.pheno<-ril_mean_day[ril_mean_day$height_type == h.types[j],]
H2<-c()
t2<-c()
e2<-c()
gxt2<-c()
# For each treatment phenotype(date), and height_type calculate the proportion of variance associated with the factor.
for(d in 6:length(colnames(type.pheno))){
  # Use only RILs with all measurements for each type.phenotype
  cc.type.pheno<-type.pheno[complete.cases(type.pheno[,d]),c(1,3,4,d)]
  # Build linear model each cofactor is a random effect
  model<-lmer(cc.type.pheno[,4]~(1|genotype)+(1|treatment)+(1|genotype:treatment), data=cc.type.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
  # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
  re<-as.numeric(VarCorr(model))
  res<-attr(VarCorr(model), "sc")^2
  # Extract individual components (order will remain the same)
  gxt.var<-re[1]
  geno.var<-re[2]
  treat.var<-re[3]
  # Total variance is sum of all variances
  tot.var<-sum(re, res)
  # Get proportion of variance
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
variance<-rbind(H2, t2, gxt2, e2)
#height_type<-rep(h.types[j], nrow(variance))
#variance<-cbind(height_type,variance)
colnames(variance)<-colnames(type.pheno)[6:length(type.pheno)]
rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
#her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
#write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
#DR13_height_var<-variance
#}  

ril_canopy_height_day_var<-variance
write.csv(ril_canopy_height_day_var, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_canopy_height_variance.csv", quote=F)

traits<-colnames(ril_canopy_height_day_var)
types<-rownames(ril_canopy_height_day_var)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}
colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

ril_canopy_height_day_var.l<-variance.l

pdf('/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_canopy_height_variance.pdf')
p<-ggplot(ril_canopy_height_day_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1)
p
dev.off()

write.csv(ril_canopy_height_day_var.l, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/height_types/ril_canopy_height_variance.l.csv", quote=F, row.names=F)


############ Lets make a data.frame of Treatment effect by method
ril_treatment_fx_height_type<-rbind(ril_height_above_bound_day_var[2,],ril_extent_y_day_var[2,],ril_centroid_y_day_var[2,],ril_ellipse_y_day_var[2,],ril_canopy_height_day_var[2,])
row.names(ril_treatment_fx_height_type)<-c("height_above_bound", "extent_y", "centroid_y", "ellipse_center_y", "canopy_height")
write.csv(ril_treatment_fx_height_type, file="/Users/mfeldman/Dropbox/setaria_height_paper/results/height_types/ril_treatment_fx_height_type.csv", quote=F, row.names=T)





save.image(file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/ril_correlation_H2_total_variance.Rdata")




