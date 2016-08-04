# Calculate broad sense heritability and partition variance of height 
# for each of the field experiments and the Dinneny Lab growout

# Load packages
library(lme4)
library(ggplot2)
library(lattice)

# Load in the datasets
# Illinois
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")
DN13<-read.csv('DN13_qtl.csv')
DN14<-read.csv('DN14_qtl.csv')
DR13<-read.csv('DR13_qtl.csv')
DR14<-read.csv('DR14_qtl.csv')

# Dinneny Lab
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/dinneny")
DL14<-read.csv("dinneny_qtl.csv")


###### Lets first do DN13
################################# DN13 ################################# 
# Total variance first
i.treat<-unique(DN13$treatment)

pheno<-DN13
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")

# Full variance partitioning
if (((length(unique(pheno[,2])) > 1 )) ) {
  # Create variables to store values
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  # For each phenotype calculate variance
  for(i in 5:length(colnames(pheno))){
    # Use only RILs with all measurements for each phenotype
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  colnames(variance)<-colnames(pheno)[5:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')

}
  
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")
traits<-colnames(variance)
types<-rownames(variance)
DN13_height_total_variance<-variance
write.csv(DN13_height_total_variance, file="DN13_height_total_variance.csv", quote=F)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}
  
colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))
  H2
DN13_height_total_variance.l<-variance.l  
write.csv(DN13_height_total_variance.l, file="DN13_height_total_variance.l.csv", quote=F)
  
pdf("DN13_height_total_variance.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

#### Now partition just genetic variance in each treatment
# DN13
i.treat<-unique(DN13$treatment)

pheno<-DN13
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")


for (t in 1:length(i.treat)) {
  # Create variables to store values
  treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
  H2<-c()
  e2<-c()
  
  # For each treatment.phenotype calculate variance
  for(i in 5:length(colnames(treatment.pheno))){
    # Use only RILs with all measurements for each treatment.phenotype
    cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  
  variance<-rbind(H2, e2)
  colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
  colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
  rownames(variance)<-c('Genotype', 'Error')
  assign(paste('variance', i.treat[t], sep="_"), variance)
}

traits<-colnames(variance)
types<-rownames(variance)
DN13_height_H2<-variance

  variance.l<-c()
  t<-c()
  for(i in 1:length(traits)) {
    t<-rep(traits[i], length(types))
    rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
    variance.l<-rbind(variance.l, rf)
  }
  colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
  rownames(variance.l)<-c(1:nrow(variance.l))
  DN13_height_H2.l<-variance.l

  setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")

  write.csv(DN13_height_H2.l, "DN13_height_heritability.csv", quote=F)

pdf("DN13_height_heritability.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())


###### Now DN14
################################# DN14 ################################# 
# Total variance first
i.treat<-unique(DN14$treatment)

pheno<-DN14
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")

# Full variance partitioning
if (((length(unique(pheno[,2])) > 1 )) & (length(unique(pheno[,3])) > 1)) {
  # Create variables to store values
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  # For each phenotype calculate variance
  for(i in 5:length(colnames(pheno))){
    # Use only RILs with all measurements for each phenotype
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  colnames(variance)<-colnames(pheno)[5:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  #her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
  #write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA) 
}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")

traits<-colnames(variance)
types<-rownames(variance)
DN14_height_total_variance<-variance
write.csv(DN14_height_total_variance, "DN14_height_total_variance.csv", quote=F)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

DN14_height_total_variance.l<-variance.l  
write.csv(DN14_height_total_variance.l, "DN14_height_total_variance.l.csv", quote=F)

pdf("DN14_height_total_variance.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

#### Now partition just genetic variance in each treatment
# DN14 
i.treat<-unique(DN14$treatment)

pheno<-DN14
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")


for (t in 1:length(i.treat)) {
  # Create variables to store values
  treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
  H2<-c()
  e2<-c()
  
  # For each treatment.phenotype calculate variance
  for(i in 5:length(colnames(treatment.pheno))){
    # Use only RILs with all measurements for each treatment.phenotype
    cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  
  variance<-rbind(H2, e2)
  colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
  colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
  rownames(variance)<-c('Genotype', 'Error')
  DN14_height_H2<-variance
}

traits<-colnames(variance)
types<-rownames(variance)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

DN14_height_heritability.l<-variance.l  

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")

write.csv(variance.l, "DN14_height_heritability.csv", quote=F)

pdf("DN14_height_heritability.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())


###### Now DR13
################################# DR13 ################################# 
# Total variance first
# DR13

i.treat<-unique(DR13$treatment)

pheno<-DR13
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")

# Full variance partitioning
if (((length(unique(pheno[,2])) > 1 )) & (length(unique(pheno[,3])) > 1)) {
  # Create variables to store values
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  # For each phenotype calculate variance
  for(i in 5:length(colnames(pheno))){
    # Use only RILs with all measurements for each phenotype
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  colnames(variance)<-colnames(pheno)[5:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
 
}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")
traits<-colnames(variance)
types<-rownames(variance)
DR13_height_total_variance<-variance
write.csv(DR13_height_total_variance,"DR13_height_total_variance.csv", quote=F)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

DR13_height_total_variance.l<-variance.l  
write.csv(DR13_height_total_variance.l,"DR13_height_total_variance.l.csv", quote=F)

pdf("DR13_height_total_variance.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

#### Now partition just genetic variance in each treatment
# DR13
i.treat<-unique(DR13$treatment)

pheno<-DR13
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")


for (t in 1:length(i.treat)) {
  # Create variables to store values
  treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
  H2<-c()
  e2<-c()
  
  # For each treatment.phenotype calculate variance
  for(i in 5:length(colnames(treatment.pheno))){
    # Use only RILs with all measurements for each treatment.phenotype
    cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  
  variance<-rbind(H2, e2)
  colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
  colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
  rownames(variance)<-c('Genotype', 'Error')
  DR13_height_H2<-variance
}

traits<-colnames(variance)
types<-rownames(variance)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")


DR13_height_heritability.l<-variance.l  
write.csv(DR13_height_heritability.l,"DR13_height_heritability.csv",quote=F)

pdf("DR13_height_heritability.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

# DR14
################################# DR14 ################################# 
###### Now DR14
i.treat<-unique(DR14$treatment)

pheno<-DR14
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")

# Full variance partitioning
if (((length(unique(pheno[,2])) > 1 )) & (length(unique(pheno[,3])) > 1)) {
  # Create variables to store values
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  # For each phenotype calculate variance
  for(i in 5:length(colnames(pheno))){
    # Use only RILs with all measurements for each phenotype
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  colnames(variance)<-colnames(pheno)[5:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')

}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")
traits<-colnames(variance)
types<-rownames(variance)
DR14_height_total_variance<-variance
write.csv(DR14_height_total_variance, "DR14_height_total_variance.csv", quote=F)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

DR14_height_total_variance.l<-variance.l  
write.csv(DR14_height_total_variance.l, "DR14_height_total_variance.l.csv", quote=F)


pdf("DR14_height_total_variance.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

#### Now partition just genetic variance in each treatment
# DR14
i.treat<-unique(DR14$treatment)

pheno<-DR14
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")


for (t in 1:length(i.treat)) {
  # Create variables to store values
  treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
  H2<-c()
  e2<-c()
  
  # For each treatment.phenotype calculate variance
  for(i in 5:length(colnames(treatment.pheno))){
    # Use only RILs with all measurements for each treatment.phenotype
    cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  
  variance<-rbind(H2, e2)
  colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
  colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
  rownames(variance)<-c('Genotype', 'Error')
  DR14_height_H2<-variance
}

traits<-colnames(variance)
types<-rownames(variance)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
 
variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))


DR14_height_heritability.l<-variance.l

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")

write.csv(variance.l, "DR14_height_heritability.csv", quote=F)

pdf("DR14_height_heritability.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())



# DL14
################################# DL14 ################################# 
###### Now DL14
i.treat<-unique(DL14$treatment)
pheno<-DL14
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")
pheno<-pheno[,c(1,2,3,4,6)]

# Full variance partitioning
if (((length(unique(pheno[,2])) > 1 ))) {
  # Create variables to store values
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  # For each phenotype calculate variance
  for(i in 5:length(colnames(pheno))){
    # Use only RILs with all measurements for each phenotype
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  colnames(variance)<-colnames(pheno)[5:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  
}

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")
traits<-colnames(variance)
types<-rownames(variance)
DL14_height_total_variance<-variance
write.csv(DL14_height_total_variance, "DL14_height_total_variance.csv", quote=F)

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))

DL14_height_total_variance.l<-variance.l  
write.csv(DL14_height_total_variance.l, "DL14_height_total_variance.l.csv", quote=F)


#pdf("DL14_height_total_variance.pdf")
#p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p
#invisible(dev.off())

p = ggplot(data=variance.l, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
p=p + geom_bar(width = 1, stat="identity")
p=p+facet_grid(facets=. ~ Trait)
p=p+xlab("Trait")
pdf("DL14_height_total_variance.pdf")
p
dev.off()

#### Now partition just genetic variance in each treatment
# DL14
i.treat<-unique(DL14$treatment)

pheno<-DL14
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")


for (t in 1:length(i.treat)) {
  # Create variables to store values
  treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
  H2<-c()
  e2<-c()
  
  # For each treatment.phenotype calculate variance
  for(i in 5:length(colnames(treatment.pheno))){
    # Use only RILs with all measurements for each treatment.phenotype
    cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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
  
  variance<-rbind(H2, e2)
  colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
  colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
  rownames(variance)<-c('Genotype', 'Error')
  DL14_height_H2<-variance
}

traits<-colnames(variance)
types<-rownames(variance)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")

variance.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
  variance.l<-rbind(variance.l, rf)
}

colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(variance.l)<-c(1:nrow(variance.l))
}

DL14_height_heritability.l<-variance.l

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")

write.csv(variance.l, "DL14_height_heritability.csv", quote=F)

pdf("DL14_height_heritability.pdf")
p<-ggplot(variance.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())



################################# LOAD IN BELLWEATHER RIL ################################# 
BP14<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance/ril_height_above_bound_variance.l.csv")
BP14$Trait<-as.character(BP14$Trait)

################################# ALL VARIANCE ################################# 
all_variance<-rbind(DN13_height_total_variance.l[29:32,], DN14_height_total_variance.l[9:12,], DL14_height_total_variance.l, DR13_height_total_variance.l[29:32,], DR14_height_total_variance.l[17:20,], BP14[89:92,]) 
all_variance$Trait<-c(rep("DN13", 4), rep("DN14", 4), rep("DL14", 4), rep("DR13", 4), rep("DR14", 4), rep("BP14", 4))
#p<-ggplot(all_variance, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'G x Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
all_variance$Trait<-factor(all_variance$Trait, levels = c('DN13', 'DN14', 'DL14', 'DR13', 'DR14', 'BP14'), ordered=TRUE)

p = ggplot(data=all_variance, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
p=p + geom_bar(width = 1, stat="identity")
p=p+facet_grid(facets=. ~ Trait)
p=p+xlab("Trait") + theme(legend.title=element_blank())
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning")
pdf("all_ril_variance_at_final_timepoint.pdf")
p
dev.off()

save.image(file="/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/ril_correlation_H2_and_total_variance.Rdata")




