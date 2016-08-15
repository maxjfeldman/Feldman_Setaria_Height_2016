# This is a script to reformat all height data from Illinois for QTL analysis

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/illinois_data")

# Lets start with DN13
DN13<-read.csv("13DN_height_step_2.csv", colClasses=c('treatment'='character', 'DAS'='character'))
# First lets remove all data that is known to be bad (flag == 2)
DN13<-DN13[DN13$flag != 2, ]
# Remove all 'artificial' entries
DN13<-DN13[DN13$trait != 'artificial',]
DN13$data<-as.numeric(as.character(DN13$data))
# Change density measure to 'dense' and 'sparse'
DN13[DN13$treatment == 'high_density', 'treatment']<-c('dense')
DN13[DN13$treatment == 'low_density', 'treatment']<-c('sparse')
# Remove all unneccessary columns
DN13<-DN13[,c(1,2,8,9,10,11,13,16)]
# Now get mean of the data by subplot_ID and DAS
DN13.ag<-aggregate(as.numeric(DN13[,8]), by=list(DN13$year, DN13$experiment, DN13$genotype, DN13$treatment, DN13$plot, DN13$subplot_ID, DN13$DAS), mean, na.action=na.pass, na.rm=TRUE)
# Need to rename columns
colnames(DN13.ag)<-c('year','experiment','id','treatment','plot','subplot_id','DAS','height')

# Lets make a height trait for each day
days<-sort(unique(DN13$DAS))
DN13_qtl<-c()
for(d in 1:length(days)) {
  temp<-DN13.ag[DN13.ag$DAS == days[d],]
  colnames(temp)[8]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DN13_qtl<-temp[,-c(7)]
  }
  if (d > 1 ) {
    DN13_qtl<-merge(DN13_qtl, temp[,-c(7)], by=c('year', 'experiment', 'id', 'treatment', 'plot', 'subplot_id'), all=T)
  }
}

# Okay now lets finish this other crap...
DN13_qtl$Obs<-c(1:nrow(DN13_qtl))
DN13_qtl$experiment<-rep("DN13", nrow(DN13_qtl))
DN13_qtl$measurement<-rep("unknown", nrow(DN13_qtl))
DN13_qtl<-DN13_qtl[,c(15,2,1,4,5,6,3,16,7:14)]
DN13_qtl$plot<-substr(DN13_qtl$subplot_id, 1,2)

# write out 
write.csv(DN13_qtl, file="DN13_qtl.csv", row.names=F, quote=F)

# DR13 next
DR13<-read.csv("13DR_height_step_2.csv", colClasses=c('treatment'='character', 'DAS'='character'))
# First lets remove all data that is known to be bad (flag == 2)
DR13$data<-as.numeric(as.character(DR13$data))
DR13<-DR13[DR13$flag != 2, ]
# Remove all 'artificial' entries
DR13<-DR13[DR13$trait != 'artificial',]

# Divide into plant height through time and final measure of culm height
DR13_culm<-DR13[DR13$trait == "culm_height",]

# Remove all unneccessary columns
DR13<-DR13[,c(1,2,9,10,11,12,15,18)]
DR13_culm<-DR13_culm[,c(1,2,9,10,11,12,15,18)]

# Lets make a height trait for each day
days<-sort(unique(DR13$DAS))
DR13_qtl<-c()
for(d in 1:length(days)) {
  temp<-DR13[DR13$DAS == days[d],]
  colnames(temp)[8]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DR13_qtl<-temp[,-c(7)]
  }
  if (d > 1 ) {
    DR13_qtl<-merge(DR13_qtl, temp[,-c(7)], by=c('year', 'experiment', 'genotype', 'treatment', 'plot', 'subplot_ID'), all=T)
  }
}


# Add the rest of the shit...
DR13_qtl$Obs<-c(1:nrow(DR13_qtl))
DR13_qtl$experiment<-rep("DR13", nrow(DR13_qtl))
DR13_qtl$measurement<-rep("unknown", nrow(DR13_qtl))
DR13_qtl<-DR13_qtl[,c(15,2,1,4,5,6,3,16,7:14)]
colnames(DR13_qtl)[c(6,7)]<-c("subplot_id", "id")

# Write to file
write.csv(DR13_qtl, file="DR13_qtl.csv", row.names=F, quote=F)

# DN14 next
DN14<-read.csv("14DN_height_step_2.csv", colClasses=c('treatment'='character'))

# First lets remove all data that is known to be bad (flag == 2)
DN14<-DN14[DN14$flag != 2, ]
DN14$data<-as.numeric(as.character(DN14$data))
# Remove all 'artificial' entries
DN14<-DN14[DN14$trait != 'artificial',]
DN14[DN14$treatment == 'high_density', 'treatment']<-c('dense')
DN14[DN14$treatment == 'low_density', 'treatment']<-c('sparse')

# Add a DAS column
DN14$DAS<-DN14$DOY - DN14$sowing_doy

# DN14 has a lot of different DAS, need to group into a single day
# Used the first day in each group as a reference point
DN14[DN14$measurement_number == '1', 'DAS']<-c("25")
DN14[DN14$measurement_number == '2', 'DAS']<-c("46")
DN14[DN14$measurement_number == '3', 'DAS']<-c("67")

# Divide into plant height through time and final measure of culm height
DN14_culm<-DN14[DN14$trait == "culm_height",]

# Remove all unneccessary columns
DN14<-DN14[,c(1,2,8,9,10,11,17,19)]
#DN14_culm<-DN14_culm[,c(1,2,8,9,10,11,17,19)]
DN14$data<-as.numeric(as.character(DN14$data))

# Need to aggregate height by replicate
DN14.ag<-aggregate(as.numeric(DN14[,7]), by=list(DN14$year, DN14$experiment, DN14$genotype, DN14$treatment, DN14$plot, DN14$subplot_ID, DN14$DAS), mean, na.action=na.pass, na.rm=TRUE)
colnames(DN14.ag)<-c('year', 'experiment', 'id', 'treatment', 'plot', 'subplot_id', 'DAS', 'height')

# Lets make a height trait for each day
days<-sort(unique(DN14$DAS))
DN14_qtl<-c()
for(d in 1:length(days)) {
  temp<-DN14.ag[DN14.ag$DAS == days[d],]
  colnames(temp)[8]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DN14_qtl<-temp[,-c(7)]
  }
  if (d > 1 ) {
    DN14_qtl<-merge(DN14_qtl, temp[,-c(7)], by=c('year', 'experiment', 'id', 'treatment', 'plot', 'subplot_id'), all=T)
  }
}


# Add the rest of the shit...
DN14_qtl$Obs<-c(1:nrow(DN14_qtl))
DN14_qtl$experiment<-rep("DN14", nrow(DN14_qtl))
DN14_qtl$measurement<-rep("unknown", nrow(DN14_qtl))
DN14_qtl<-DN14_qtl[,c(10,2,1,4,5,6,3,11,7:9)]
#colnames(DN14_qtl)[c(11)]<-c("culm_height")

# Write to file
write.csv(DN14_qtl, file="DN14_qtl.csv", row.names=F, quote=F)

# DR14
DR14<-read.csv("14DR_height_step_2.csv", colClasses=c('treatment'='character'))
# First lets remove all data that is known to be bad (flag == 2)
DR14<-DR14[DR14$flag != 2, ]
# Remove all 'artificial' entries
DR14<-DR14[DR14$DAS != '0',]
DR14$data<-as.numeric(as.character(DR14$data))

# Divide into plant height through time and final measure of culm height
DR14_culm<-DR14[DR14$trait == "culm_height",]

# Remove all unneccessary columns
DR14<-DR14[,c(1,2,8,9,10,11,15,19)]
DR14_culm<-DR14_culm[,c(1,2,9,10,11,12,15,19)]

# Aggregate replicates
DR14.ag<-aggregate(as.numeric(DR14[,8]), by=list(DR14$year, DR14$experiment, DR14$genotype, DR14$treatment, DR14$plot, DR14$subplot_ID, DR14$DAS), mean, na.action=na.pass, na.rm=TRUE)
colnames(DR14.ag)<-c('year', 'experiment', 'id', 'treatment', 'plot', 'subplot_id', 'DAS', 'height')

# Lets make a height trait for each day
days<-sort(unique(DR14.ag$DAS))
DR14_qtl<-c()
for(d in 1:length(days)) {
  temp<-DR14.ag[DR14.ag$DAS == days[d],]
  colnames(temp)[8]<-paste('height', days[d], sep="_")
  if (d == 1) {
    DR14_qtl<-temp[,-c(7)]
  }
  if (d > 1 ) {
    DR14_qtl<-merge(DR14_qtl, temp[,-c(7)], by=c('year', 'experiment', 'id', 'treatment', 'plot', 'subplot_id'), all=T)
  }
}


# Add the rest of the shit...
DR14_qtl$Obs<-c(1:nrow(DR14_qtl))
DR14_qtl$experiment<-rep("DR14", nrow(DR14_qtl))
DR14_qtl$measurement<-rep("unknown", nrow(DR14_qtl))
DR14_qtl<-DR14_qtl[,c(12,2,1,4,5,6,3,13,7:11)]

#colnames(DR14_qtl)[c(13)]<-c("culm_height")
#DR14_qtl[,8:ncol(DR14_qtl)]<-DR14_qtl[,9:ncol(DR14_qtl)]

# Write to file
write.csv(DR14_qtl, file="DR14_qtl.csv", row.names=F, quote=F)


