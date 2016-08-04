
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
    #unify_chr(temp)
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


##########################
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
S.BP14_rate.d<-read.csv("summary.table.slod.dry_BP14_AGR.csv")
S.BP14_rate.w<-read.csv("summary.table.slod.wet_BP14_AGR.csv")

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

exp<-unique(slod.st_uni$exp)
exp<-as.character(exp)

########################################################################################################
# Make a plot of the % variance explained and allelic effects through time of the SLOD results
########################################################################################################

# Lets make a plot of proportional contribution of each QTL through time in each experiment
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
pdf('Figure_S9.pdf')
for(e in exp){
  temp<-slod.st_uni[slod.st_uni$exp == e,]
  p<-ggplot(temp,aes(x=day, y=prop.var, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ggtitle(e)
  print(p)
}  
dev.off()

# Lets make a plot of alleleic effect size of each QTL through time in each experiment
exp<-unique(slod.st_uni$exp)
exp<-as.character(exp)
pdf('Figure_S10.pdf')
for(e in exp){
  temp<-slod.st_uni[slod.st_uni$exp == e,]
  p<-ggplot(temp,aes(x=day, y=additive.fx, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ggtitle(e)
  print(p + geom_pointrange(aes(ymin = (additive.fx-additive.fx_se), ymax=(additive.fx+additive.fx_se)))) 
}  
dev.off()


# Lets make a figure for the Bellweather growouts (Figure_7a and Figure_7b). This will represent the supplemental Figures.

temp<-slod.st_uni[slod.st_uni$exp == 'BP14',]
# Lets add a color column
temp$t.col<-rep('NA', nrow(temp))
temp[temp$treatment == 'dry', 't.col']<-c("orange")
temp[temp$treatment == 'wet', 't.col']<-c("blue")

# Plot % variance over time on phenotyper
pdf("Figure_7a.pdf")
p<-ggplot(temp,aes(x=day, y=prop.var, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ylab("% variance explained") + xlab("Days after planting") + scale_colour_manual(values=c("orange", "blue"))
print(p)
dev.off()

# Plot fx size
pdf("Figure_7b.pdf")
p<-ggplot(temp,aes(x=day, y=additive.fx, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ylab("Alleic effect") + xlab("Days after planting") + scale_colour_manual(values=c("orange", "blue"))
print(p + geom_pointrange(aes(ymin = (additive.fx-additive.fx_se), ymax=(additive.fx+additive.fx_se)))) 
dev.off()



#####################################################################
# Lets look at the rates
#####################################################################


##########################
# Now lets read in data
setwd("/Users/mfeldman/Desktop/temp/summary.table")

# SLOD
S.BP14_rate.d<-read.csv("summary.table.slod.dry_BP14_AGR.csv")
S.BP14_rate.w<-read.csv("summary.table.slod.wet_BP14_AGR.csv")

S.BP14_rate<-rbind(S.BP14_rate.d, S.BP14_rate.w)

S.BP14_rate$trait<-as.character(S.BP14_rate$trait)

treatment<-c()
experiment<-c()
year<-c()
type<-c()
for(r in 1:nrow(S.BP14_rate)){
  treat<-strsplit(S.BP14_rate$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.BP14_rate$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.BP14_rate))
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  type<-c(type, t)
}

exp<-experiment
S.BP14_rate<-cbind(S.BP14_rate, treatment, exp, year, type)
colnames(S.BP14_rate)[1]<-c('marker')


S.BP14_rate_uni<-unify_marker(S.BP14_rate)
S.BP14_rate_uni<-unique(S.BP14_rate_uni)

# Lets add growth day
day<-c()
day<-str_extract(S.BP14_rate_uni[,'trait'], '_[0-9]+_')
day<-str_extract(day, '[0-9]+')

S.BP14_rate_uni$day<-rep("NA", nrow(S.BP14_rate_uni))
S.BP14_rate_uni$day<-day
S.BP14_rate_uni$day<-as.numeric(as.character(S.BP14_rate_uni$day))

# Lets plot proportion of variance in each experiment into a .pdf file 
exp<-unique(S.BP14_rate_uni$exp)
exp<-as.character(exp)
pdf('BP14_AGR_prop.var.pdf')
for(e in exp){
  temp<-S.BP14_rate_uni[S.BP14_rate_uni$exp == e,]
  p<-ggplot(temp,aes(x=day, y=prop.var, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ggtitle(e)
  print(p)
}  
dev.off()


# Lets plot fx size in each experiment into a .pdf file 
exp<-unique(S.BP14_rate_uni$exp)
exp<-as.character(exp)
pdf('BP14_AGR_fx.size.pdf')
for(e in exp){
  temp<-S.BP14_rate_uni[S.BP14_rate_uni$exp == e,]
  p<-ggplot(temp,aes(x=day, y=additive.fx, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ggtitle(e)
  print(p + geom_pointrange(aes(ymin = (additive.fx-additive.fx_se), ymax=(additive.fx+additive.fx_se)))) 
}  
dev.off()








#####################################################################
# Lets look at the rates
#####################################################################


##########################
# Now lets read in data
setwd("/Users/mfeldman/Desktop/current_project/temp/summary.table")

# SLOD
S.BP14_rate.d<-read.csv("summary.table.slod.dry_BP14_AGR.csv")
S.BP14_rate.w<-read.csv("summary.table.slod.wet_BP14_AGR.csv")

S.BP14_rate<-rbind(S.BP14_rate.d, S.BP14_rate.w)

S.BP14_rate$trait<-as.character(S.BP14_rate$trait)

treatment<-c()
experiment<-c()
year<-c()
type<-c()
for(r in 1:nrow(S.BP14_rate)){
  treat<-strsplit(S.BP14_rate$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.BP14_rate$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.BP14_rate))
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  type<-c(type, t)
}

exp<-experiment
S.BP14_rate<-cbind(S.BP14_rate, treatment, exp, year, type)
colnames(S.BP14_rate)[1]<-c('marker')


S.BP14_rate_uni<-unify_marker(S.BP14_rate)
S.BP14_rate_uni<-unique(S.BP14_rate_uni)

# Lets add growth day
day<-c()
day<-str_extract(S.BP14_rate_uni[,'trait'], '_[0-9]+_')
day<-str_extract(day, '[0-9]+')

S.BP14_rate_uni$day<-rep("NA", nrow(S.BP14_rate_uni))
S.BP14_rate_uni$day<-day
S.BP14_rate_uni$day<-as.numeric(as.character(S.BP14_rate_uni$day))

# Lets plot proportion of variance in each experiment into a .pdf file 
exp<-unique(S.BP14_rate_uni$exp)
exp<-as.character(exp)
pdf('BP14_AGR_prop.var.pdf')
for(e in exp){
  temp<-S.BP14_rate_uni[S.BP14_rate_uni$exp == e,]
  p<-ggplot(temp,aes(x=day, y=prop.var, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ggtitle(e)
  print(p)
}  
dev.off()


# Lets plot fx size in each experiment into a .pdf file 
exp<-unique(S.BP14_rate_uni$exp)
exp<-as.character(exp)
pdf('BP14_AGR_fx.size.pdf')
for(e in exp){
  temp<-S.BP14_rate_uni[S.BP14_rate_uni$exp == e,]
  p<-ggplot(temp,aes(x=day, y=additive.fx, color=treatment)) + geom_line() + facet_wrap(~marker, scales = "free_y") + theme_bw() + ggtitle(e)
  print(p + geom_pointrange(aes(ymin = (additive.fx-additive.fx_se), ymax=(additive.fx+additive.fx_se)))) 
}  
dev.off()














