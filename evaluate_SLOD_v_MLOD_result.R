###################
# In this script we examine the properties of SLOD and MLOD using their summary tables
# The most powerful part of this analysis was the examination of pLOD model values / model permutation threshold
# This may go into the methods section


################### 
library(ggplot2)
library(stringr)

################### Lets make some fxns
merged_table<-c()
unify_chr<-function(temp){
  all_qtl<-sort(table(temp$marker), decreasing=T)
  if(length(all_qtl) > 1){
    m.name<-names(all_qtl)[1]
    ave.pos<-mean(temp[temp$marker == m.name, 'pos'])
    cr<-unique(temp[temp$marker == names(all_qtl)[1],'chr'])
    po<-unique(temp[temp$marker == names(all_qtl)[1],'pos'])
    max.pos<-ave.pos+10
    min.pos<-ave.pos-10
    subset<-temp[temp$pos > min.pos & temp$pos < max.pos,]
    subset$marker<-rep(m.name, length(subset$marker))
    subset$chr<-rep(cr, length(subset$chr))
    subset$pos<-rep(po, length(subset$pos))
    temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
    merged_table<<-rbind(merged_table, subset)
    unify_chr(temp)
  } 
  if(length(all_qtl) == 1){
    merged_table<<-rbind(merged_table, temp)
    unify_chr(temp)
  }
}


unify_marker<-function(input){
  
  chrs<-sort(unique(input$chr))
  merged_table<<-c()
  for(ch in chrs) {
    temp<-input[input$chr == ch,]
    temp$marker<-as.character(temp$marker)
    unify_chr(temp)
  }
  return(merged_table)
}


#####################################################################
# Start with SLOD
#####################################################################

# Now lets read in data
setwd("/Users/mfeldman/Desktop/current_project/temp/summary.table")

# SLOD
S.DN13.d<-read.csv("summary.table.slod.dense_DN13_height_blup.csv")
S.DN13.s<-read.csv("summary.table.slod.sparse_DN13_height_blup.csv")
S.DN14.d<-read.csv("summary.table.slod.dense_DN14_height_blup.csv")
S.DN14.s<-read.csv("summary.table.slod.sparse_DN14_height_blup.csv")
S.DR13.d<-read.csv("summary.table.slod.dry_DR13_height_blup.csv")
S.DR13.w<-read.csv("summary.table.slod.wet_DR13_height_blup.csv")
S.DR14.d<-read.csv("summary.table.slod.dry_DR14_height_blup.csv")
S.DR14.w<-read.csv("summary.table.slod.wet_DR14_height_blup.csv")
S.BP14.d<-read.csv("summary.table.slod.dry_BP14_ril_height_above_bound.csv")
S.BP14.w<-read.csv("summary.table.slod.wet_BP14_ril_height_above_bound.csv")

slod.st<-rbind(S.DN13.d,S.DN13.s, S.DN14.d, S.DN14.s, S.DR13.d, S.DR13.w, S.DR14.d, S.DR14.w, S.BP14.d, S.BP14.w)
slod.st$trait<-as.character(slod.st$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
for(r in 1:nrow(slod.st)){
  treat<-strsplit(slod.st$trait[r], "\\.")[[1]][1]
  trait<-strsplit(slod.st$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(slod.st))
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  type<-c(type, t)
}

exp<-experiment
slod.st<-cbind(slod.st, treatment, exp, year, type)
colnames(slod.st)[1]<-c('marker')


slod.st_uni<-unify_marker(slod.st)
slod.st_uni<-unique(slod.st_uni)

# Lets add growth day
day<-c()
day<-str_extract(slod.st_uni[,'trait'], '_[0-9]+_')
day<-str_extract(day, '[0-9]+')

slod.st_uni$day<-rep("NA", nrow(slod.st_uni))
slod.st_uni$day<-day
slod.st_uni$day<-as.numeric(as.character(slod.st_uni$day))



#####################################################################
# MLOD
#####################################################################
# Now lets read in data
setwd("/Users/mfeldman/Desktop/current_project/temp/summary.table")
# Now lets load MLOD
# mlod
M.DN13.d<-read.csv("summary.table.mlod.dense_DN13_height_blup.csv")
M.DN13.s<-read.csv("summary.table.mlod.sparse_DN13_height_blup.csv")
M.DN14.d<-read.csv("summary.table.mlod.dense_DN14_height_blup.csv")
M.DN14.s<-read.csv("summary.table.mlod.sparse_DN14_height_blup.csv")
M.DR13.d<-read.csv("summary.table.mlod.dry_DR13_height_blup.csv")
M.DR13.w<-read.csv("summary.table.mlod.wet_DR13_height_blup.csv")
M.DR14.d<-read.csv("summary.table.mlod.dry_DR14_height_blup.csv")
M.DR14.w<-read.csv("summary.table.mlod.wet_DR14_height_blup.csv")
M.BP14.d<-read.csv("summary.table.mlod.dry_BP14_ril_height_above_bound.csv")
M.BP14.w<-read.csv("summary.table.mlod.wet_BP14_ril_height_above_bound.csv")

# insert S.DR13.w and bellweather
mlod.st<-rbind(M.DN13.d,M.DN13.s, M.DN14.d, M.DN14.s, M.DR13.d, M.DR13.w, M.DR14.d, M.DR14.w, M.BP14.d, M.BP14.w)

mlod.st$trait<-as.character(mlod.st$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
for(r in 1:nrow(mlod.st)){
  treat<-strsplit(mlod.st$trait[r], "\\.")[[1]][1]
  trait<-strsplit(mlod.st$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(mlod.st))
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  type<-c(type, t)
}

exp<-experiment
mlod.st<-cbind(mlod.st, treatment, exp, year, type)
colnames(mlod.st)[1]<-c('marker')


mlod.st_uni<-unify_marker(mlod.st)
mlod.st_uni<-unique(mlod.st_uni)

# Lets add growth day
day<-c()
day<-str_extract(mlod.st_uni[,'trait'], '_[0-9]+_')
day<-str_extract(day, '[0-9]+')


mlod.st_uni$day<-rep("NA", nrow(mlod.st_uni))
mlod.st_uni$day<-day
mlod.st_uni$day<-as.numeric(as.character(mlod.st_uni$day))
mlod.st_uni$day<-rep("NA", nrow(mlod.st_uni))
mlod.st_uni$day<-day
mlod.st_uni$day<-as.numeric(as.character(mlod.st_uni$day))




###### BP14
BP14.uni.s<-slod.st_uni[slod.st_uni$exp == 'BP14',]
BP14.uni.m<-mlod.st_uni[mlod.st_uni$exp == 'BP14',]
###### Lets break them out by treatment
# dry
BP14.s.d<-unique(BP14.uni.s[BP14.uni.s$treatment == 'dry', 'marker'])
BP14.m.d<-unique(BP14.uni.m[BP14.uni.m$treatment == 'dry', 'marker'])

# wet
BP14.s.w<-unique(BP14.uni.s[BP14.uni.s$treatment == 'wet', 'marker'])
BP14.m.w<-unique(BP14.uni.m[BP14.uni.m$treatment == 'wet', 'marker'])

BP14.d.count<-c(length(BP14.s.d), length(BP14.m.d))
BP14.w.count<-c(length(BP14.s.w), length(BP14.m.w))

BP14.count<-rbind(BP14.d.count, BP14.w.count)
rownames(BP14.count)<-c("dry", "wet")
colnames(BP14.count)<-c("slod", "mlod")

###### DR14
DR14.uni.s<-slod.st_uni[slod.st_uni$exp == 'DR14',]
DR14.uni.m<-mlod.st_uni[mlod.st_uni$exp == 'DR14',]

###### Lets break them out by treatment
# dry
DR14.s.d<-unique(DR14.uni.s[DR14.uni.s$treatment == 'dry', 'marker'])
DR14.m.d<-unique(DR14.uni.m[DR14.uni.m$treatment == 'dry', 'marker'])

# wet
DR14.s.w<-unique(DR14.uni.s[DR14.uni.s$treatment == 'wet', 'marker'])
DR14.m.w<-unique(DR14.uni.m[DR14.uni.m$treatment == 'wet', 'marker'])

DR14.d.count<-c(length(DR14.s.d), length(DR14.m.d))
DR14.w.count<-c(length(DR14.s.w), length(DR14.m.w))

DR14.count<-rbind(DR14.d.count, DR14.w.count)
rownames(DR14.count)<-c("dry", "wet")
colnames(DR14.count)<-c("slod", "mlod")


###### DR13
DR13.uni.s<-slod.st_uni[slod.st_uni$exp == 'DR13',]
DR13.uni.m<-mlod.st_uni[mlod.st_uni$exp == 'DR13',]

###### Lets break them out by treatment
# dry
DR13.s.d<-unique(DR13.uni.s[DR13.uni.s$treatment == 'dry', 'marker'])
DR13.m.d<-unique(DR13.uni.m[DR13.uni.m$treatment == 'dry', 'marker'])

# wet
DR13.s.w<-unique(DR13.uni.s[DR13.uni.s$treatment == 'wet', 'marker'])
DR13.m.w<-unique(DR13.uni.m[DR13.uni.m$treatment == 'wet', 'marker'])

DR13.d.count<-c(length(DR13.s.d), length(DR13.m.d))
DR13.w.count<-c(length(DR13.s.w), length(DR13.m.w))

DR13.count<-rbind(DR13.d.count, DR13.w.count)
rownames(DR13.count)<-c("dry", "wet")
colnames(DR13.count)<-c("slod", "mlod")

###### DN14
DN14.uni.s<-slod.st_uni[slod.st_uni$exp == 'DN14',]
DN14.uni.m<-mlod.st_uni[mlod.st_uni$exp == 'DN14',]

###### Lets break them out by treatment
# dry
DN14.s.d<-unique(DN14.uni.s[DN14.uni.s$treatment == 'dense', 'marker'])
DN14.m.d<-unique(DN14.uni.m[DN14.uni.m$treatment == 'dense', 'marker'])

# wet
DN14.s.s<-unique(DN14.uni.s[DN14.uni.s$treatment == 'sparse', 'marker'])
DN14.m.s<-unique(DN14.uni.m[DN14.uni.m$treatment == 'sparse', 'marker'])

DN14.d.count<-c(length(DN14.s.d), length(DN14.m.d))
DN14.s.count<-c(length(DN14.s.s), length(DN14.m.s))

DN14.count<-rbind(DN14.d.count, DN14.s.count)
rownames(DN14.count)<-c("dense", "sparse")
colnames(DN14.count)<-c("slod", "mlod")


###### DN13
DN13.uni.s<-slod.st_uni[slod.st_uni$exp == 'DN13',]
DN13.uni.m<-mlod.st_uni[mlod.st_uni$exp == 'DN13',]

###### Lets break them out by treatment
# dry
DN13.s.d<-unique(DN13.uni.s[DN13.uni.s$treatment == 'dense', 'marker'])
DN13.m.d<-unique(DN13.uni.m[DN13.uni.m$treatment == 'dense', 'marker'])

# wet
DN13.s.s<-unique(DN13.uni.s[DN13.uni.s$treatment == 'sparse', 'marker'])
DN13.m.s<-unique(DN13.uni.m[DN13.uni.m$treatment == 'sparse', 'marker'])

DN13.d.count<-c(length(DN13.s.d), length(DN13.m.d))
DN13.s.count<-c(length(DN13.s.s), length(DN13.m.s))

DN13.count<-rbind(DN13.d.count, DN13.s.count)
rownames(DN13.count)<-c("dense", "sparse")
colnames(DN13.count)<-c("slod", "mlod")

rownames(BP14.count)<-paste('BP14', rownames(BP14.count), sep="_")
rownames(DR14.count)<-paste('DR14', rownames(DR14.count), sep="_")
rownames(DR13.count)<-paste('DR13', rownames(DR13.count), sep="_")
rownames(DN14.count)<-paste('DN14', rownames(DN14.count), sep="_")
rownames(DN13.count)<-paste('DN13', rownames(DN13.count), sep="_")

# Make table of this data and see if one method is more variable than the other
sm.qtl.count.table<-rbind(BP14.count, DR14.count, DR13.count, DN14.count, DN13.count)
slod.var<-var(sm.qtl.count.table[,1])
mlod.var<-var(sm.qtl.count.table[,2])

slod.cv<-sqrt(slod.var)/mean(sm.qtl.count.table[,1])
mlod.cv<-sqrt(mlod.var)/mean(sm.qtl.count.table[,2])

#### BP14

# Start with mean
BP14.uni.s.d.ave<-mean(BP14.uni.s[BP14.uni.s$treatment == 'dry','lod'])
BP14.uni.s.w.ave<-mean(BP14.uni.s[BP14.uni.s$treatment == 'wet','lod'])

BP14.uni.m.d.ave<-mean(BP14.uni.m[BP14.uni.m$treatment == 'dry','lod'])
BP14.uni.m.w.ave<-mean(BP14.uni.m[BP14.uni.m$treatment == 'wet','lod'])

# Now do median
BP14.uni.s.d.med<-median(BP14.uni.s[BP14.uni.s$treatment == 'dry','lod'])
BP14.uni.s.w.med<-median(BP14.uni.s[BP14.uni.s$treatment == 'wet','lod'])

BP14.uni.m.d.med<-median(BP14.uni.m[BP14.uni.m$treatment == 'dry','lod'])
BP14.uni.m.w.med<-median(BP14.uni.m[BP14.uni.m$treatment == 'wet','lod'])

#### DR14

# Start with mean
DR14.uni.s.d.ave<-mean(DR14.uni.s[DR14.uni.s$treatment == 'dry','lod'])
DR14.uni.s.w.ave<-mean(DR14.uni.s[DR14.uni.s$treatment == 'wet','lod'])

DR14.uni.m.d.ave<-mean(DR14.uni.m[DR14.uni.m$treatment == 'dry','lod'])
DR14.uni.m.w.ave<-mean(DR14.uni.m[DR14.uni.m$treatment == 'wet','lod'])

# Now do median
DR14.uni.s.d.med<-median(DR14.uni.s[DR14.uni.s$treatment == 'dry','lod'])
DR14.uni.s.w.med<-median(DR14.uni.s[DR14.uni.s$treatment == 'wet','lod'])

DR14.uni.m.d.med<-median(DR14.uni.m[DR14.uni.m$treatment == 'dry','lod'])
DR14.uni.m.w.med<-median(DR14.uni.m[DR14.uni.m$treatment == 'wet','lod'])

#### DR13

# Start with mean
DR13.uni.s.d.ave<-mean(DR13.uni.s[DR13.uni.s$treatment == 'dry','lod'])
DR13.uni.s.w.ave<-mean(DR13.uni.s[DR13.uni.s$treatment == 'wet','lod'])

DR13.uni.m.d.ave<-mean(DR13.uni.m[DR13.uni.m$treatment == 'dry','lod'])
DR13.uni.m.w.ave<-mean(DR13.uni.m[DR13.uni.m$treatment == 'wet','lod'])

# Now do median
DR13.uni.s.d.med<-median(DR13.uni.s[DR13.uni.s$treatment == 'dry','lod'])
DR13.uni.s.w.med<-median(DR13.uni.s[DR13.uni.s$treatment == 'wet','lod'])

DR13.uni.m.d.med<-median(DR13.uni.m[DR13.uni.m$treatment == 'dry','lod'])
DR13.uni.m.w.med<-median(DR13.uni.m[DR13.uni.m$treatment == 'wet','lod'])

#### DN14

# Start with mean
DN14.uni.s.d.ave<-mean(DN14.uni.s[DN14.uni.s$treatment == 'dense','lod'])
DN14.uni.s.s.ave<-mean(DN14.uni.s[DN14.uni.s$treatment == 'sparse','lod'])

DN14.uni.m.d.ave<-mean(DN14.uni.m[DN14.uni.m$treatment == 'dense','lod'])
DN14.uni.m.s.ave<-mean(DN14.uni.m[DN14.uni.m$treatment == 'sparse','lod'])

# Now do median
DN14.uni.s.d.med<-median(DN14.uni.s[DN14.uni.s$treatment == 'dense','lod'])
DN14.uni.s.s.med<-median(DN14.uni.s[DN14.uni.s$treatment == 'sparse','lod'])

DN14.uni.m.d.med<-median(DN14.uni.m[DN14.uni.m$treatment == 'dense','lod'])
DN14.uni.m.s.med<-median(DN14.uni.m[DN14.uni.m$treatment == 'sparse','lod'])

#### DN13

# Start with mean
DN13.uni.s.d.ave<-mean(DN13.uni.s[DN13.uni.s$treatment == 'dense','lod'])
DN13.uni.s.s.ave<-mean(DN13.uni.s[DN13.uni.s$treatment == 'sparse','lod'])

DN13.uni.m.d.ave<-mean(DN13.uni.m[DN13.uni.m$treatment == 'dense','lod'])
DN13.uni.m.s.ave<-mean(DN13.uni.m[DN13.uni.m$treatment == 'sparse','lod'])

# Now do median
DN13.uni.s.d.med<-median(DN13.uni.s[DN13.uni.s$treatment == 'dense','lod'])
DN13.uni.s.s.med<-median(DN13.uni.s[DN13.uni.s$treatment == 'sparse','lod'])

DN13.uni.m.d.med<-median(DN13.uni.m[DN13.uni.m$treatment == 'dense','lod'])
DN13.uni.m.s.med<-median(DN13.uni.m[DN13.uni.m$treatment == 'sparse','lod'])

######### Add mean and median values to the table
mean.slod<-c(BP14.uni.s.d.ave, BP14.uni.s.w.ave, DR14.uni.s.d.ave, DR14.uni.s.w.ave, DR13.uni.s.d.ave, DR13.uni.s.w.ave, DN14.uni.s.d.ave, DN14.uni.s.s.ave, DN13.uni.s.d.ave, DN13.uni.s.s.ave)
mean.mlod<-c(BP14.uni.m.d.ave, BP14.uni.m.w.ave, DR14.uni.m.d.ave, DR14.uni.m.w.ave, DR13.uni.m.d.ave, DR13.uni.m.w.ave, DN14.uni.m.d.ave, DN14.uni.m.s.ave, DN13.uni.m.d.ave, DN13.uni.m.s.ave)

median.slod<-c(BP14.uni.s.d.med, BP14.uni.s.w.med, DR14.uni.s.d.med, DR14.uni.s.w.med, DR13.uni.s.d.med, DR13.uni.s.w.med, DN14.uni.s.d.med, DN14.uni.s.s.med, DN13.uni.s.d.med, DN13.uni.s.s.med)
median.mlod<-c(BP14.uni.m.d.med, BP14.uni.m.w.med, DR14.uni.m.d.med, DR14.uni.m.w.med, DR13.uni.m.d.med, DR13.uni.m.w.med, DN14.uni.m.d.med, DN14.uni.m.s.med, DN13.uni.m.d.med, DN13.uni.m.s.med)




temp<-cbind(mean.slod, mean.mlod, median.slod, median.mlod)
temp2<-cbind(sm.qtl.count.table, temp)

# Load in a manually entered data.frame that contains data on the penalized LOD score of
# the resulting QTL model and on the penalty used to generate the model
setwd("/Users/mfeldman/Desktop/current_project/temp/")
model.scores<-read.csv("model.scores.csv")

slod.ratio<-model.scores$plod.slod/model.scores$m.perm.slod
mlod.ratio<-model.scores$plod.mlod/model.scores$m.perm.mlod

model.scores<-cbind(model.scores, slod.ratio, mlod.ratio)
model.scores<-as.data.frame(model.scores)

# Correlation between plod.slod and slod.ratio
cor(model.scores[,'plod.slod'], model.scores[,'slod.ratio'])
cor(model.scores[,'plod.mlod'], model.scores[,'mlod.ratio'])
cor(model.scores[,'plod.slod'], model.scores[,'plod.mlod'])

temp2<-cbind(model.scores, temp2)
# slod number v. slod ratio
cor(temp2[,'slod'], temp2[,'slod.ratio'])
# slod number v. slod plod
cor(temp2[,'slod'], temp2[,'plod.slod'])

# mlod number v. mlod ratio
cor(temp2[,'mlod'], temp2[,'mlod.ratio'])
# mlod number v. mlod plod
cor(temp2[,'mlod'], temp2[,'plod.mlod'])

# Time points
BP14.tps<-length(unique(BP14.uni.s$trait))/2
DR14.tps<-length(unique(DR14.uni.s$trait))/2
DR13.tps<-length(unique(DR13.uni.s$trait))/2
DN14.tps<-length(unique(DN14.uni.s$trait))/2
DN13.tps<-length(unique(DN13.uni.s$trait))/2

timepoint<-c(BP14.tps, BP14.tps, DR14.tps, DR14.tps, DR13.tps, DR13.tps, DN14.tps, DN14.tps, DN13.tps, DN13.tps)

temp2<-cbind(temp2, timepoint)
# Check ratio v. timepoint
cor(temp2[,'slod.ratio'], temp2[,'timepoint'])
cor(temp2[,'mlod.ratio'], temp2[,'timepoint'])

# Check pLOD vs # of QTL
cor(temp2[,'plod.slod'], temp2[,'slod'])
cor(temp2[,'plod.mlod'], temp2[,'mlod'])

cor(temp2[,'plod.slod'], temp2[,'timepoint'])
cor(temp2[,'plod.mlod'], temp2[,'timepoint'])

cor(temp2[,'slod'], temp2[,'timepoint'])
cor(temp2[,'mlod'], temp2[,'timepoint'])

mean(temp2$slod.ratio)
mean(temp2$mlod.ratio)














