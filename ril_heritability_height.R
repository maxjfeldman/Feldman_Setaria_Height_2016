# Calculate broad sense heritability and partition variance of height 
# for each of the experiments

# Load packages
library(lme4)
library(ggplot2)
library(lattice)
library(lmomco)

# Load in the datasets
# Illinois
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")
DN13<-read.csv('DN13_qtl.csv')
DN14<-read.csv('DN14_qtl.csv')
DR13<-read.csv('DR13_qtl.csv')
DR14<-read.csv('DR14_qtl.csv')
 
# Dinneny Lab
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/dinneny")
DL13<-read.csv("dinneny_qtl.csv")


setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril")
BP14<-read.csv('ril_height_qtl_raw.csv')


################################################################################
# Define functions
################################################################################

# Broad sense heritability
get_h2<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")
  H2<-c()
  for (i in 5:length(colnames(pheno))){
    # Get complete cases
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
  
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno)
  
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.numeric(VarCorr(model))
    res<-attr(VarCorr(model), "sc")^2
  
    # Extract individual components (order will remain the same)
    gxt.var<-re[1]
    geno.var<-re[2]
    treat.var<-re[3]
    # Total variance is sum of all variances
    tot.var<-sum(re, res)
    
    reps.t1<-table(pheno[pheno$treatment == i.treat[1], 'id'])
    reps.t2<-table(pheno[pheno$treatment == i.treat[2], 'id'])
    reps.treatment<-c(reps.t1, reps.t2)
    
    reps.t1<-as.character(unique(pheno[pheno$treatment == i.treat[1], 'id']))
    reps.t2<-as.character(unique(pheno[pheno$treatment == i.treat[2], 'id']))
    unique.combined <- c(as.character(reps.t1), as.character(reps.t2))
    
    freq.unique.combined <- table(unique.combined)
    
    # Calculate the harmonic mean replication within treatment blocks
    hm_treatment<-harmonic.mean(freq.unique.combined)$harmean
    
    # Now get a count of total genotypic replication
    reps.total<-table(pheno[,'id'])
    # Get the harmonic mean of this quantity
    hm_total<-harmonic.mean(reps.total)$harmean
    
    # Calculate heritability as described by AEL 
    # H2 = geno.var/(geno.var + (gxt.var/harmonic mean of treatment block replication) + (residual.var/harmonic mean of total genotype replication) )
    h2<-((geno.var)/(geno.var + (gxt.var/hm_treatment) + (res/hm_total)))
    
    # This is the heritability
    H2<-c(H2,h2)
     
  }
  names(H2)<-colnames(pheno)[5:ncol(pheno)]
  return(H2)
}

# Heritability within treatment
get_h2_in_treatment<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")
  variance.out<-c()
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
    variance.out<-cbind(variance.out, variance)
  }
  return(variance.out)
}  

# Total variance partition
get_total_var<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")
  
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  for (i in 5:length(colnames(pheno))){
    # Get complete cases
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno)
    
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
  return(variance)
}


# Total variance partition field (includes plot)
get_total_var_field<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")
  
  H2<-c()
  t2<-c()
  e2<-c()
  p2<-c()
  gxt2<-c()
  
  for (i in 5:length(colnames(pheno))){
    # Get complete cases
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
    
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|plot)+(1|id:treatment), data=cc.pheno)
    
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.numeric(VarCorr(model))
    res<-attr(VarCorr(model), "sc")^2
    
    # Extract individual components (order will remain the same)
    gxt.var<-re[1]
    geno.var<-re[2]
    plot.var<-re[3]
    treat.var<-re[4]
    # Total variance is sum of all variances
    tot.var<-sum(re, res)
    
    # Get proportion of variance for all factors
    h<-geno.var/tot.var
    t<-treat.var/tot.var
    e<-res/tot.var
    gxt<-gxt.var/tot.var
    p<-plot.var/tot.var
    # Append variables to a vector of variables
    H2<-c(H2,h)
    t2<-c(t2,t)
    e2<-c(e2,e)
    gxt2<-c(gxt2, gxt)
    p2<-c(p2, p)
    
  }
  variance<-rbind(H2, t2, p2, gxt2, e2)
  colnames(variance)<-colnames(pheno)[5:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment','Plot', 'G x Treatment', 'Error')
  return(variance)
}


################################################################################
# END define Fxn
################################################################################

# Do analysis...

################
# DN13
################

# First lets do broad sense heritability
DN13.h2<-get_h2(DN13)
DN13.h2_in_treatment<-get_h2_in_treatment(DN13)

# Lets write these files out 
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
write.csv(DN13.h2, file="DN13_height_h2.csv", quote=F)
write.csv(DN13.h2_in_treatment, file="DN13_height_h2_in_treatment.csv", quote=F)

# Lets also make a combined file in long form for plotting
treatments<-c()
entries<-c()
gs<-c()
es<-c()
h2<-c()
days<-c()
for(i in 1:ncol(DN13.h2_in_treatment)){
  c.name<-colnames(DN13.h2_in_treatment)[i]
  treatment<-strsplit(c.name, '_')[[1]][1]
  trait<-strsplit(c.name, '_')[[1]][2]
  day<-strsplit(c.name, '_')[[1]][3]
  experiment<-strsplit(c.name, '_')[[1]][4]
  entry<-paste(trait, day, experiment, sep='_')
  # Get individual values
  g<-DN13.h2_in_treatment[1,i]
  e<-DN13.h2_in_treatment[2,i]
  treatments<-c(treatments, treatment)
  entries<-c(entries, entry)
  gs<-c(gs,g)
  es<-c(es,e)
  days<-c(days, day)
}

DN13.h2.all<-cbind(entries, treatments, days, gs, es)
colnames(DN13.h2.all)<-c('trait', 'treatment', 'dap', 'h2', 'residual')

# Make small data.frame of only H2
trait<-names(DN13.h2)
treatment<-rep('all', length(DN13.h2))
dap<-sort(unique(days))
h2<-DN13.h2
residual<-1-h2
all<-cbind(trait, treatment, dap, h2, residual)

# Combine h2 and partition of genetic variance within treatment
DN13.h2.all<-rbind(DN13.h2.all, all)
rownames(DN13.h2.all)<-c(1:nrow(DN13.h2.all))
DN13.h2.all<-as.data.frame(DN13.h2.all)

# Change to numeric
DN13.h2.all$dap<-as.numeric(as.character(DN13.h2.all$dap))
DN13.h2.all$h2<-as.numeric(as.character(DN13.h2.all$h2))
DN13.h2.all$residual<-as.numeric(as.character(DN13.h2.all$residual))

# Lets make a plot of these heritabilities...
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
p<-ggplot(DN13.h2.all, aes(x=dap, y=h2, group=treatment)) + geom_line(aes(colour = treatment)) + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("DN13_height_heritability_all.pdf")
p
dev.off()
# Write to file
write.csv(DN13.h2.all, file="DN13_height_h2_all.csv", quote=F)

#####
# Now lets write out total variance
DN13.total_var<-get_total_var(DN13)
DN13.total_var_field<-get_total_var_field(DN13)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")

# Write output for all variance components
traits<-colnames(DN13.total_var)
types<-rownames(DN13.total_var)
DN13_height_total_variance<-DN13.total_var
write.csv(DN13_height_total_variance, file="DN13_height_total_variance_np.csv", quote=F)

# Convert into long form...
DN13.total_var.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DN13.total_var[,i], stringsAsFactors=F))
  DN13.total_var.l<-rbind(DN13.total_var.l, rf)
}

colnames(DN13.total_var.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DN13.total_var.l)<-c(1:nrow(DN13.total_var.l))
DN13_height_total_variance.l<-DN13.total_var.l 
write.csv(DN13_height_total_variance.l, file="DN13_height_total_variance_np.l.csv", quote=F)

pdf("DN13_height_total_variance_np.pdf")
p<-ggplot(DN13.total_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

# Write output for all variance components including field plot as a factor
traits<-colnames(DN13.total_var_field)
types<-rownames(DN13.total_var_field)
DN13_height_total_variance_field<-DN13.total_var_field
write.csv(DN13_height_total_variance_field, file="DN13_height_total_variance_field.csv", quote=F)

# Convert into long form...
DN13.total_var_field.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DN13.total_var_field[,i], stringsAsFactors=F))
  DN13.total_var_field.l<-rbind(DN13.total_var_field.l, rf)
}

colnames(DN13.total_var_field.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DN13.total_var_field.l)<-c(1:nrow(DN13.total_var_field.l))
DN13_height_total_variance_field.l<-DN13.total_var_field.l 
write.csv(DN13_height_total_variance_field.l, file="DN13_height_total_variance_field.l.csv", quote=F)

pdf("DN13_height_total_variance_field.pdf")
p<-ggplot(DN13.total_var_field.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())



################
# DN14
################

# First lets do broad sense heritability
DN14.h2<-get_h2(DN14)
DN14.h2_in_treatment<-get_h2_in_treatment(DN14)

# Lets write these files out 
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
write.csv(DN14.h2, file="DN14_height_h2.csv", quote=F)
write.csv(DN14.h2_in_treatment, file="DN14_height_h2_in_treatment.csv", quote=F)

# Lets also make a combined file in long form for plotting
treatments<-c()
entries<-c()
gs<-c()
es<-c()
h2<-c()
days<-c()
for(i in 1:ncol(DN14.h2_in_treatment)){
  c.name<-colnames(DN14.h2_in_treatment)[i]
  treatment<-strsplit(c.name, '_')[[1]][1]
  trait<-strsplit(c.name, '_')[[1]][2]
  day<-strsplit(c.name, '_')[[1]][3]
  experiment<-strsplit(c.name, '_')[[1]][4]
  entry<-paste(trait, day, experiment, sep='_')
  # Get individual values
  g<-DN14.h2_in_treatment[1,i]
  e<-DN14.h2_in_treatment[2,i]
  treatments<-c(treatments, treatment)
  entries<-c(entries, entry)
  gs<-c(gs,g)
  es<-c(es,e)
  days<-c(days, day)
}

DN14.h2.all<-cbind(entries, treatments, days, gs, es)
colnames(DN14.h2.all)<-c('trait', 'treatment', 'dap', 'h2', 'residual')

# Make small data.frame of only H2
trait<-names(DN14.h2)
treatment<-rep('all', length(DN14.h2))
dap<-sort(unique(days))
h2<-DN14.h2
residual<-1-h2
all<-cbind(trait, treatment, dap, h2, residual)

# Combine h2 and partition of genetic variance within treatment
DN14.h2.all<-rbind(DN14.h2.all, all)
rownames(DN14.h2.all)<-c(1:nrow(DN14.h2.all))
DN14.h2.all<-as.data.frame(DN14.h2.all)

# Change to numeric
DN14.h2.all$dap<-as.numeric(as.character(DN14.h2.all$dap))
DN14.h2.all$h2<-as.numeric(as.character(DN14.h2.all$h2))
DN14.h2.all$residual<-as.numeric(as.character(DN14.h2.all$residual))

# Lets make a plot of these heritabilities...
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
p<-ggplot(DN14.h2.all, aes(x=dap, y=h2, group=treatment)) + geom_line(aes(colour = treatment)) + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("DN14_height_heritability_all.pdf")
p
dev.off()
# Write to file
write.csv(DN14.h2.all, file="DN14_height_h2_all.csv", quote=F)

#####
# Now lets write out total variance
DN14.total_var<-get_total_var(DN14)
DN14.total_var_field<-get_total_var_field(DN14)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")

# Write output for all variance components
traits<-colnames(DN14.total_var)
types<-rownames(DN14.total_var)
DN14_height_total_variance<-DN14.total_var
write.csv(DN14_height_total_variance, file="DN14_height_total_variance_np.csv", quote=F)

# Convert into long form...
DN14.total_var.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DN14.total_var[,i], stringsAsFactors=F))
  DN14.total_var.l<-rbind(DN14.total_var.l, rf)
}

colnames(DN14.total_var.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DN14.total_var.l)<-c(1:nrow(DN14.total_var.l))
DN14_height_total_variance.l<-DN14.total_var.l 
write.csv(DN14_height_total_variance.l, file="DN14_height_total_variance_np.l.csv", quote=F)

pdf("DN14_height_total_variance_np.pdf")
p<-ggplot(DN14.total_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

# Write output for all variance components including field plot as a factor
traits<-colnames(DN14.total_var_field)
types<-rownames(DN14.total_var_field)
DN14_height_total_variance_field<-DN14.total_var_field
write.csv(DN14_height_total_variance_field, file="DN14_height_total_variance_field.csv", quote=F)

# Convert into long form...
DN14.total_var_field.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DN14.total_var_field[,i], stringsAsFactors=F))
  DN14.total_var_field.l<-rbind(DN14.total_var_field.l, rf)
}

colnames(DN14.total_var_field.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DN14.total_var_field.l)<-c(1:nrow(DN14.total_var_field.l))
DN14_height_total_variance_field.l<-DN14.total_var_field.l 
write.csv(DN14_height_total_variance_field.l, file="DN14_height_total_variance_field.l.csv", quote=F)

pdf("DN14_height_total_variance_field.pdf")
p<-ggplot(DN14.total_var_field.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())



################
# DR13
################

# First lets do broad sense heritability
DR13.h2<-get_h2(DR13)
DR13.h2_in_treatment<-get_h2_in_treatment(DR13)

# Lets write these files out 
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
write.csv(DR13.h2, file="DR13_height_h2.csv", quote=F)
write.csv(DR13.h2_in_treatment, file="DR13_height_h2_in_treatment.csv", quote=F)

# Lets also make a combined file in long form for plotting
treatments<-c()
entries<-c()
gs<-c()
es<-c()
h2<-c()
days<-c()
for(i in 1:ncol(DR13.h2_in_treatment)){
  c.name<-colnames(DR13.h2_in_treatment)[i]
  treatment<-strsplit(c.name, '_')[[1]][1]
  trait<-strsplit(c.name, '_')[[1]][2]
  day<-strsplit(c.name, '_')[[1]][3]
  experiment<-strsplit(c.name, '_')[[1]][4]
  entry<-paste(trait, day, experiment, sep='_')
  # Get individual values
  g<-DR13.h2_in_treatment[1,i]
  e<-DR13.h2_in_treatment[2,i]
  treatments<-c(treatments, treatment)
  entries<-c(entries, entry)
  gs<-c(gs,g)
  es<-c(es,e)
  days<-c(days, day)
}

DR13.h2.all<-cbind(entries, treatments, days, gs, es)
colnames(DR13.h2.all)<-c('trait', 'treatment', 'dap', 'h2', 'residual')

# Make small data.frame of only H2
trait<-names(DR13.h2)
treatment<-rep('all', length(DR13.h2))
dap<-sort(unique(days))
h2<-DR13.h2
residual<-1-h2
all<-cbind(trait, treatment, dap, h2, residual)

# Combine h2 and partition of genetic variance within treatment
DR13.h2.all<-rbind(DR13.h2.all, all)
rownames(DR13.h2.all)<-c(1:nrow(DR13.h2.all))
DR13.h2.all<-as.data.frame(DR13.h2.all)

# Change to numeric
DR13.h2.all$dap<-as.numeric(as.character(DR13.h2.all$dap))
DR13.h2.all$h2<-as.numeric(as.character(DR13.h2.all$h2))
DR13.h2.all$residual<-as.numeric(as.character(DR13.h2.all$residual))

# Lets make a plot of these heritabilities...
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
p<-ggplot(DR13.h2.all, aes(x=dap, y=h2, group=treatment)) + geom_line(aes(colour = treatment)) + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("DR13_height_heritability_all.pdf")
p
dev.off()
# Write to file
write.csv(DR13.h2.all, file="DR13_height_h2_all.csv", quote=F)

#####
# Now lets write out total variance
DR13.total_var<-get_total_var(DR13)
DR13.total_var_field<-get_total_var_field(DR13)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")

# Write output for all variance components
traits<-colnames(DR13.total_var)
types<-rownames(DR13.total_var)
DR13_height_total_variance<-DR13.total_var
write.csv(DR13_height_total_variance, file="DR13_height_total_variance_np.csv", quote=F)

# Convert into long form...
DR13.total_var.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DR13.total_var[,i], stringsAsFactors=F))
  DR13.total_var.l<-rbind(DR13.total_var.l, rf)
}

colnames(DR13.total_var.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DR13.total_var.l)<-c(1:nrow(DR13.total_var.l))
DR13_height_total_variance.l<-DR13.total_var.l 
write.csv(DR13_height_total_variance.l, file="DR13_height_total_variance_np.l.csv", quote=F)

pdf("DR13_height_total_variance_np.pdf")
p<-ggplot(DR13.total_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

# Write output for all variance components including field plot as a factor
traits<-colnames(DR13.total_var_field)
types<-rownames(DR13.total_var_field)
DR13_height_total_variance_field<-DR13.total_var_field
write.csv(DR13_height_total_variance_field, file="DR13_height_total_variance_field.csv", quote=F)

# Convert into long form...
DR13.total_var_field.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DR13.total_var_field[,i], stringsAsFactors=F))
  DR13.total_var_field.l<-rbind(DR13.total_var_field.l, rf)
}

colnames(DR13.total_var_field.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DR13.total_var_field.l)<-c(1:nrow(DR13.total_var_field.l))
DR13_height_total_variance_field.l<-DR13.total_var_field.l 
write.csv(DR13_height_total_variance_field.l, file="DR13_height_total_variance_field.l.csv", quote=F)

pdf("DR13_height_total_variance_field.pdf")
p<-ggplot(DR13.total_var_field.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())


################
# DR14
################

# First lets do broad sense heritability
DR14.h2<-get_h2(DR14)
DR14.h2_in_treatment<-get_h2_in_treatment(DR14)

# Lets write these files out 
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
write.csv(DR14.h2, file="DR14_height_h2.csv", quote=F)
write.csv(DR14.h2_in_treatment, file="DR14_height_h2_in_treatment.csv", quote=F)

# Lets also make a combined file in long form for plotting
treatments<-c()
entries<-c()
gs<-c()
es<-c()
h2<-c()
days<-c()
for(i in 1:ncol(DR14.h2_in_treatment)){
  c.name<-colnames(DR14.h2_in_treatment)[i]
  treatment<-strsplit(c.name, '_')[[1]][1]
  trait<-strsplit(c.name, '_')[[1]][2]
  day<-strsplit(c.name, '_')[[1]][3]
  experiment<-strsplit(c.name, '_')[[1]][4]
  entry<-paste(trait, day, experiment, sep='_')
  # Get individual values
  g<-DR14.h2_in_treatment[1,i]
  e<-DR14.h2_in_treatment[2,i]
  treatments<-c(treatments, treatment)
  entries<-c(entries, entry)
  gs<-c(gs,g)
  es<-c(es,e)
  days<-c(days, day)
}

DR14.h2.all<-cbind(entries, treatments, days, gs, es)
colnames(DR14.h2.all)<-c('trait', 'treatment', 'dap', 'h2', 'residual')

# Make small data.frame of only H2
trait<-names(DR14.h2)
treatment<-rep('all', length(DR14.h2))
dap<-sort(unique(days))
h2<-DR14.h2
residual<-1-h2
all<-cbind(trait, treatment, dap, h2, residual)

# Combine h2 and partition of genetic variance within treatment
DR14.h2.all<-rbind(DR14.h2.all, all)
rownames(DR14.h2.all)<-c(1:nrow(DR14.h2.all))
DR14.h2.all<-as.data.frame(DR14.h2.all)

# Change to numeric
DR14.h2.all$dap<-as.numeric(as.character(DR14.h2.all$dap))
DR14.h2.all$h2<-as.numeric(as.character(DR14.h2.all$h2))
DR14.h2.all$residual<-as.numeric(as.character(DR14.h2.all$residual))

# Lets make a plot of these heritabilities...
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
p<-ggplot(DR14.h2.all, aes(x=dap, y=h2, group=treatment)) + geom_line(aes(colour = treatment)) + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("DR14_height_heritability_all.pdf")
p
dev.off()
# Write to file
write.csv(DR14.h2.all, file="DR14_height_h2_all.csv", quote=F)

#####
# Now lets write out total variance
DR14.total_var<-get_total_var(DR14)
DR14.total_var_field<-get_total_var_field(DR14)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_field")

# Write output for all variance components
traits<-colnames(DR14.total_var)
types<-rownames(DR14.total_var)
DR14_height_total_variance<-DR14.total_var
write.csv(DR14_height_total_variance, file="DR14_height_total_variance_np.csv", quote=F)

# Convert into long form...
DR14.total_var.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DR14.total_var[,i], stringsAsFactors=F))
  DR14.total_var.l<-rbind(DR14.total_var.l, rf)
}

colnames(DR14.total_var.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DR14.total_var.l)<-c(1:nrow(DR14.total_var.l))
DR14_height_total_variance.l<-DR14.total_var.l 
write.csv(DR14_height_total_variance.l, file="DR14_height_total_variance_np.l.csv", quote=F)

pdf("DR14_height_total_variance_np.pdf")
p<-ggplot(DR14.total_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())

# Write output for all variance components including field plot as a factor
traits<-colnames(DR14.total_var_field)
types<-rownames(DR14.total_var_field)
DR14_height_total_variance_field<-DR14.total_var_field
write.csv(DR14_height_total_variance_field, file="DR14_height_total_variance_field.csv", quote=F)

# Convert into long form...
DR14.total_var_field.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DR14.total_var_field[,i], stringsAsFactors=F))
  DR14.total_var_field.l<-rbind(DR14.total_var_field.l, rf)
}

colnames(DR14.total_var_field.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DR14.total_var_field.l)<-c(1:nrow(DR14.total_var_field.l))
DR14_height_total_variance_field.l<-DR14.total_var_field.l 
write.csv(DR14_height_total_variance_field.l, file="DR14_height_total_variance_field.l.csv", quote=F)

pdf("DR14_height_total_variance_field.pdf")
p<-ggplot(DR14.total_var_field.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())







################
# DL13
################

# First lets do broad sense heritability
DL13.h2<-get_h2(DL13)
DL13.h2_in_treatment<-get_h2_in_treatment(DL13)

# Lets write these files out 
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_bellweather")
write.csv(DL13.h2, file="DL13_height_h2.csv", quote=F)
write.csv(DL13.h2_in_treatment, file="DL13_height_h2_in_treatment.csv", quote=F)

# Lets also make a combined file in long form for plotting
treatments<-c()
entries<-c()
gs<-c()
es<-c()
h2<-c()
days<-c()
for(i in 1:ncol(DL13.h2_in_treatment)){
  c.name<-colnames(DL13.h2_in_treatment)[i]
  treatment<-strsplit(c.name, '_')[[1]][1]
  # Lets do something a bit fancy because trait names contain underscore '_'
  trait<-strsplit(c.name, '_')[[1]]
  trait<-trait[-c(1,length(trait))]
  trait<-paste0(trait, collapse = "_")
  experiment<-strsplit(c.name, '_')[[1]][length(strsplit(c.name, '_')[[1]])]
  entry<-paste(trait, experiment, sep='_')
  # Get individual values
  g<-DL13.h2_in_treatment[1,i]
  e<-DL13.h2_in_treatment[2,i]
  treatments<-c(treatments, treatment)
  entries<-c(entries, entry)
  gs<-c(gs,g)
  es<-c(es,e)
}

DL13.h2.all<-cbind(entries, treatments, gs, es)
colnames(DL13.h2.all)<-c('trait', 'treatment', 'h2', 'residual')

# Make small data.frame of only H2
trait<-names(DL13.h2)
treatment<-rep('all', length(DL13.h2))
h2<-DL13.h2
residual<-1-h2
all<-cbind(trait, treatment, h2, residual)

# Combine h2 and partition of genetic variance within treatment
DL13.h2.all<-rbind(DL13.h2.all, all)
rownames(DL13.h2.all)<-c(1:nrow(DL13.h2.all))
DL13.h2.all<-as.data.frame(DL13.h2.all)

# Change to numeric
DL13.h2.all$h2<-as.numeric(as.character(DL13.h2.all$h2))
DL13.h2.all$residual<-as.numeric(as.character(DL13.h2.all$residual))

# Lets make a plot of these heritabilities...
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_bellweather")
p<-ggplot(DL13.h2.all, aes(x=trait, y=h2, group=treatment)) + geom_line(aes(colour = treatment)) + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("DL13_height_heritability_all.pdf")
p
dev.off()
# Write to file
write.csv(DL13.h2.all, file="DL13_height_h2_all.csv", quote=F)

#####
# Now lets write out total variance
DL13.total_var<-get_total_var(DL13)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_bellweather")

# Write output for all variance components
traits<-colnames(DL13.total_var)
types<-rownames(DL13.total_var)
DL13_height_total_variance<-DL13.total_var
write.csv(DL13_height_total_variance, file="DL13_height_total_variance.csv", quote=F)

# Convert into long form...
DL13.total_var.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(DL13.total_var[,i], stringsAsFactors=F))
  DL13.total_var.l<-rbind(DL13.total_var.l, rf)
}

colnames(DL13.total_var.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(DL13.total_var.l)<-c(1:nrow(DL13.total_var.l))
DL13_height_total_variance.l<-DL13.total_var.l 
write.csv(DL13_height_total_variance.l, file="DL13_height_total_variance.l.csv", quote=F)

pdf("DL13_height_total_variance.pdf")
p<-ggplot(DL13.total_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())



################
# BP14
################

# First lets do broad sense heritability
BP14.h2<-get_h2(BP14)
BP14.h2_in_treatment<-get_h2_in_treatment(BP14)

# Lets write these files out 
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_bellweather")
write.csv(BP14.h2, file="BP14_height_h2.csv", quote=F)
write.csv(BP14.h2_in_treatment, file="BP14_height_h2_in_treatment.csv", quote=F)

# Lets also make a combined file in long form for plotting
treatments<-c()
entries<-c()
gs<-c()
es<-c()
h2<-c()
days<-c()
for(i in 1:ncol(BP14.h2_in_treatment)){
  c.name<-colnames(BP14.h2_in_treatment)[i]
  treatment<-strsplit(c.name, '_')[[1]][1]
  trait<-strsplit(c.name, '_')[[1]][2]
  day<-strsplit(c.name, '_')[[1]][3]
  experiment<-strsplit(c.name, '_')[[1]][4]
  entry<-paste(trait, day, experiment, sep='_')
  # Get individual values
  g<-BP14.h2_in_treatment[1,i]
  e<-BP14.h2_in_treatment[2,i]
  treatments<-c(treatments, treatment)
  entries<-c(entries, entry)
  gs<-c(gs,g)
  es<-c(es,e)
  days<-c(days, day)
}

BP14.h2.all<-cbind(entries, treatments, days, gs, es)
colnames(BP14.h2.all)<-c('trait', 'treatment', 'dap', 'h2', 'residual')

# Make small data.frame of only H2
trait<-names(BP14.h2)
treatment<-rep('all', length(BP14.h2))
dap<-sort(unique(days))
h2<-BP14.h2
residual<-1-h2
all<-cbind(trait, treatment, dap, h2, residual)

# Combine h2 and partition of genetic variance within treatment
BP14.h2.all<-rbind(BP14.h2.all, all)
rownames(BP14.h2.all)<-c(1:nrow(BP14.h2.all))
BP14.h2.all<-as.data.frame(BP14.h2.all)

# Change to numeric
BP14.h2.all$dap<-as.numeric(as.character(BP14.h2.all$dap))
BP14.h2.all$h2<-as.numeric(as.character(BP14.h2.all$h2))
BP14.h2.all$residual<-as.numeric(as.character(BP14.h2.all$residual))

# Lets make a plot of these heritabilities...
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_bellweather")
p<-ggplot(BP14.h2.all, aes(x=dap, y=h2, group=treatment)) + geom_line(aes(colour = treatment)) + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("BP14_height_heritability_all.pdf")
p
dev.off()
# Write to file
write.csv(BP14.h2.all, file="BP14_height_h2_all.csv", quote=F)

#####
# Now lets write out total variance
BP14.total_var<-get_total_var(BP14)

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/total_variance_bellweather")

# Write output for all variance components
traits<-colnames(BP14.total_var)
types<-rownames(BP14.total_var)
BP14_height_total_variance<-BP14.total_var
write.csv(BP14_height_total_variance, file="BP14_height_total_variance.csv", quote=F)

# Convert into long form...
BP14.total_var.l<-c()
t<-c()
for(i in 1:length(traits)) {
  t<-rep(traits[i], length(types))
  rf<-cbind(t, types, data.frame(BP14.total_var[,i], stringsAsFactors=F))
  BP14.total_var.l<-rbind(BP14.total_var.l, rf)
}

colnames(BP14.total_var.l)<-c('Trait', 'Type', 'Proportion_explained')
rownames(BP14.total_var.l)<-c(1:nrow(BP14.total_var.l))
BP14_height_total_variance.l<-BP14.total_var.l 
write.csv(BP14_height_total_variance.l, file="BP14_height_total_variance.l.csv", quote=F)

# This is Figure 1c.
pdf("BP14_height_total_variance.pdf")
p<-ggplot(BP14.total_var.l, aes(x=Trait, y=Proportion_explained, group=Type)) + geom_line(aes(colour = Type)) + theme_bw() + scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
invisible(dev.off())



# Lets test if method used to measure PHT in the field in 2014 is more heritable than the method used in 2013

# Calculate the mean variance between years
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/variance_partitioning/ril/H2_field")
DN13_H2<-read.csv("DN13_height_h2.csv")
DN14_H2<-read.csv("DN14_height_h2.csv")
DR13_H2<-read.csv("DR13_height_h2.csv")
DR14_H2<-read.csv("DR14_height_h2.csv")

# Lets take the mean of H2 from each experiment and each year
# 2013
DR13_H2.ave<-mean(DR13_H2[,2])
DN13_H2.ave<-mean(DN13_H2[,2])
# 2014
DR14_H2.ave<-mean(DR14_H2[,2])
DN14_H2.ave<-mean(DN14_H2[,2])

# Get average heritability for each year
year2013<-mean(DR13_H2.ave, DN13_H2.ave)
year2014<-mean(DR14_H2.ave, DN14_H2.ave)

# Is one year more heritable than the other?
year2014 - year2013
# Yes 2014 is ~11.3% more heritable than 2013
















