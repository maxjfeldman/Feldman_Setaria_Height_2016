# Analysis of variance 
library(ggplot2)
library(lme4)
library(nlme)
library(car)

#############################################################################
# Load in data and format
#############################################################################

# Load in Illinois raw data
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")
DR13<-read.csv('DR13_qtl.csv')
DR14<-read.csv('DR14_qtl.csv')
DN13<-read.csv('DN13_qtl.csv')
DN14<-read.csv('DN14_qtl.csv')

# Load in Dinneny Lab data
dl.height<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/data/dinneny/dinneny_qtl.csv")
dl.height<-dl.height[,c(2:7,10)]
colnames(dl.height)[c(5,7)]<-c('subplot_id','height')
dl.height<-dl.height[complete.cases(dl.height),]

# Load in Bellweather data
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril")
ril<-read.csv("sv.vis_data_ril.csv")
# Redo first few steps from variance scripts
# Need to remove plants out of bound
ril<-ril[ril$in_bound == 'True',]

# Lets remove all columns that arn't to do with height
ril<-subset(ril, select=-c(camera, imgtype, area, hull.area, solidity, perimeter, width, height, longest_axis, center.of.mass.x, hull_vertices, ellipse_center_x, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity, y.position, height_below_bound,above_bound_area,percent_above_bound_area, below_bound_area,percent_below_bound_area,date,dap,extent_x))

# This order makes more sense
ril<-ril[,c(4,2,13,5,12,3,12,1,11,15,6,8,9,10,7)]

# From the previous section we believe that the mean of 'height_above_bound' is the best measure
# Perform this operation
ril.mean<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), mean)
colnames(ril.mean)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')

# Lets do our best to remove empty cars...
plot(ril.mean$height_above_bound~ril.mean$dap_i, pch=20, cex=0.3, ylab="height", xlab="day")

# Need to zoom in 17 - 32 dap_i
close_up<-ril.mean[ril.mean$dap_i > 16,]
plot(close_up$height_above_bound~close_up$dap_i, ylim=c(0,25), pch=20, cex=0.3)
filter_1<-close_up[close_up$height_above_bound > 5,]
plot(filter_1$height_above_bound~filter_1$dap_i, ylim=c(0,25), pch=20, cex=0.3)
filter_2<-filter_1[(filter_1$height_above_bound > 12 & (filter_1$dap_i == 28 | filter_1$dap_i == 29)),'plant_id']
filter_2<-filter_1[filter_1$plant_id %in% filter_2, ]
plot(filter_2$height_above_bound~filter_2$dap_i, ylim=c(0,25), pch=20, cex=0.3)

# Now remove all of those images from RIL mean that are no good
ril.mean<-ril.mean[ril.mean$plant_id %in% filter_2$plant_id, ]
plot(ril.mean$height_above_bound~ril.mean$dap_i, pch=20, cex=0.3, ylab="height", xlab="day")

# Get only measures from the final reliable timepoint (DAP 30 & 31)
bp14<-ril.mean[ril.mean$dap_i == 30 | ril.mean$dap_i == 31,]
bp14$experiment<-rep('BP14', nrow(bp14))
bp14$year<-rep('2014', nrow(bp14))
bp14$plot<-rep('bellweather', nrow(bp14))
bp14$subplot_id<-rep('bellweather', nrow(bp14))
bp14<-bp14[,c(11,12,4,13,14,3,6)]
colnames(bp14)[c(6,7)]<-c('id','height')

#############################################################################
# Begin ANOVA
#############################################################################
# Perform ANOVA using Anova from cars so we can specify type of ANOVA
# Performed on raw values not on values fitted using BLUP
##################################################################################
# Lets see what major components of variance are between trials
##################################################################################
# Lets to an ANOVA to test for differences in location, treatment 
dr13.lastday<-DR13[,c(2:7,15)]
colnames(dr13.lastday)[7]<-c('height')
dn13.lastday<-DN13[,c(2:7,16)]
colnames(dn13.lastday)[7]<-c('height')
dr14.lastday<-DR14[,c(2:7,13)]
colnames(dr14.lastday)[7]<-c('height')
dn14.lastday<-DN14[,c(2:7,11)]
colnames(dn14.lastday)[7]<-c('height')

# Get the measurements for the last day for all drought experiments, all field experiments, drought and density
dr.all_last_day<-rbind(dr13.lastday,dr14.lastday,bp14, dl.height)
# All measurements of last day in field experiments
all.field_last_day<-rbind(dr13.lastday,dn13.lastday,dr14.lastday,dn14.lastday)
# Pool field measurments in field Drought study
dr.height_last_day<-rbind(dr13.lastday, dr14.lastday)
# Pool field measurments in field Density study
dn.height_last_day<-rbind(dn13.lastday, dn14.lastday)


# All drought
# Lets make a column for drought experiment location
dr.all_last_day$location<-rep('NA', nrow(dr.all_last_day))
# Add location as a field
dr.all_last_day[dr.all_last_day$experiment == 'DR13', 'location']<-c("illinois")
dr.all_last_day[dr.all_last_day$experiment == 'DR14', 'location']<-c("illinois")
dr.all_last_day[dr.all_last_day$experiment == 'DL13', 'location']<-c("stanford")
dr.all_last_day[dr.all_last_day$experiment == 'BP14', 'location']<-c("ddpsc")


# Do an ANOVA to find the most significant factor in all drought experiemnts

# All drought

dr.all_last_day$treatment<-as.factor(dr.all_last_day$treatment)
# Notice intercept term is removed form the ANOVA style model (y ~ -1 + ...)
model.all.dr<-lm(height~0+location+treatment%in%location+id, data=dr.all_last_day)
summary(model.all.dr)
all.dr.aov.table<-Anova(model.all.dr, type=c(3), singular.ok = T)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/table_data/anova_tables/")
write.csv(all.dr.aov.table, file="all_drought_treatment_aov.csv", quote=F, row.names=T)

# Now do the ANOVA in drought for the field
############################## 
# Drought field
############################## 
dr.height_last_day$year<-as.factor(dr.height_last_day$year)
dr.height_last_day$plot<-as.factor(dr.height_last_day$plot)
model.dr<-lm(height~0+year+plot%in%treatment+treatment+id, data=dr.height_last_day)
dr.field.aov.table<-Anova(model.dr, type=c(3), singular.ok = T)

summary(model.dr)
write.csv(dr.field.aov.table, "field_drought_treatment_aov.csv", quote=F, row.names=T)

# Now do the ANOVA in drought for the field drought 2013
# DR13
dr13.lastday$plot<-as.factor(dr13.lastday$plot)
dr13.lastday$treatment<-as.factor(dr13.lastday$treatment)
model.dr13<-lm(height~0+treatment+plot%in%treatment+id, data=dr13.lastday)
dr13.aov.table<-Anova(model.dr13, type=c(3), singular.ok = T)
dr13.model.summary<-summary(model.dr13)
write.csv(dr13.aov.table, "dr13_aov.table.csv", quote=F, row.names=T)

# Now do the ANOVA in drought for the field drought 2014
# DR14
dr14.lastday$plot<-as.factor(dr14.lastday$plot)
dr14.lastday$treatment<-as.factor(dr14.lastday$treatment)
model.dr14<-lm(height~0+plot%in%treatment+treatment+id, data=dr14.lastday)
Anova(model.dr14, type=c(3), singular.ok = T)
summary(model.dr14)
dr14.aov.table<-Anova(model.dr14)
dr14.model.summary<-summary(model.dr14)
write.csv(dr14.aov.table, "dr14_aov.table.csv", quote=F, row.names=T)


# Do ANOVA for Dinneny drought experiment
############################## 
# DL
############################## 
model.dl<-lm(height~0+treatment+id, data=dl.height)
dl.aov.table<-Anova(model.dl, type=c(3), singular.ok = T)
summary(model.dl)
write.csv(dl.aov.table, 'carnegie_drought_aov.csv', quote=F, row.names=T)

# Do ANOVA for Bellweather data
############################## 
# BP
############################## 
model.bp<-lm(height~0+treatment+id, data=bp14)
bp.aov.table<-Anova(model.bp, type=c(3), singular.ok = T)
summary(model.bp)
write.csv(bp.aov.table, file="bellweather_ril_aov.table.csv", quote=F, row.names=T)



###### What day does treatment become a signficiant factor?
# DR13
t.points<-c(9:ncol(DR13))
DR13_treatment_sig<-get.field_anova.qtl_format(DR13, t.points)

# DR14
t.points<-c(9:ncol(DR14))
DR14_treatment_sig<-get.field_anova.qtl_format(DR14, t.points)

# BP
ril_height_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril/ril_height_qtl_raw.csv")
t.points<-c(9:ncol(ril_height_qtl))
BP14_treatment_sig<-get.bp_anova.qtl_format(ril_height_qtl, t.points)



############################## 
# Now Density
############################## 

# All DN
dn.height_last_day<-all.field_last_day[all.field_last_day$experiment == 'DN13' | all.field_last_day$experiment == 'DN14',]
dn.height_last_day$year<-as.factor(dn.height_last_day$year)
dn.height_last_day$plot<-as.factor(dn.height_last_day$plot)
dn.height_last_day$treatment<-as.factor(dn.height_last_day$treatment)
model.dn.height<-lm(height~0+treatment+plot%in%treatment+year+id, data=dn.height_last_day)
all_density_treatment_aov<-Anova(model.dn.height, type=c(3), singular.ok = T)
write.csv(all_density_treatment_aov, file="all_density_treatment_aov.csv", quote=F, row.names=T)


# DN13
dn13.lastday$plot<-as.factor(dn13.lastday$plot)
model.dn13<-lm(height~0+treatment+plot%in%treatment+id, data=dn13.lastday)
dn13_aov.table<-Anova(model.dn13, type=c(3), singular.ok = T)
write.csv(dn13_aov.table, file="dn13_aov.table.csv", quote=F, row.names=T)


# DN14
dn14.lastday$plot<-as.factor(dn14.lastday$plot)
model.dn14<-lm(height~0+treatment+plot%in%treatment+id, data=dn14.lastday)
dn14_aov.table<-Anova(model.dn14, type=c(3), singular.ok = T)
write.csv(dn14_aov.table, file="dn14_aov.table.csv", quote=F, row.names=T)


# All field 
all.field_last_day$year<-as.factor(all.field_last_day$year)
all.field_last_day$plot<-as.factor(all.field_last_day$plot)
model.all.field<-lm(height~0+year+plot%in%treatment+treatment+id, data=all.field_last_day)
all.field.aov<-Anova(model.all.field, type=c(3), singular.ok = T)
write.csv(all.field.aov, file="all_field_aov.csv", quote=F, row.names=T)



###### What day does treatment become a signficiant factor?
# DN13
t.points<-c(9:ncol(DN13))
DN13_treatment_sig<-get.field_anova.qtl_format(DN13, t.points)

# DN14
t.points<-c(9:ncol(DN14))
DN14_treatment_sig<-get.field_anova.qtl_format(DN14, t.points)

save.image('anova_analysis_ril.Rdata')

############################## 
# Define fxns
############################## 

get.field_anova.qtl_format<-function(dataset, times) {
  # Get vector of individual treatments and height types
  
  pheno<-dataset
  t.points<-times
  pval<-c()
  
  # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with all factors.
  for(d in t.points) {
    # Use only RILs with all measurements for each treatment.phenotype
    cc.pheno<-pheno[complete.cases(pheno[,d]),c(4,5,7,d)]
    # Build linear model each cofactor is a random effect
    model<-lm(cc.pheno[,4]~0+id+treatment+plot%in%treatment, data=cc.pheno)
    # Extract p-value from model object, save individual 
    df<-Anova(model, type=c(3), singular.ok = T)
    p<-df[1:3,4]
    p<-c(colnames(pheno)[d], p)
    names(p)<-c('date',rownames(df)[1:3])
    pval<-rbind(pval, p)
  }
  
  return(pval)
}


get.bp_anova.qtl_format<-function(dataset, times) {
  # Get vector of individual treatments and height types
  
  pheno<-dataset
  t.points<-times
  pval<-c()
  
  # For each treatment phenotype(date), and height_type calculate the proportion of variance associated with all factors.
  for(d in t.points) {
    # Use only RILs with all measurements for each treatment.phenotype
    cc.pheno<-pheno[complete.cases(pheno[,d]),c(4,7,d)]
    # Build linear model each cofactor is a random effect
    model<-lm(cc.pheno[,3]~0+id+treatment, data=cc.pheno)
    # Extract p-value from model object, save individual 
    df<-Anova(model, type=c(3), singular.ok = T)
    p<-df[1:3,4]
    p<-c(colnames(pheno)[d], p)
    names(p)<-c('date',rownames(df)[1:3])
    pval<-rbind(pval, p)
  }
  
  return(pval)
}

