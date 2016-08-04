setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/qtl/summary_tables_known_qtl")

# Read in .csv
DN13.k.qtl<-read.csv("DN13_height_blup_known_qtl_concatenated_summary_table.csv", header=F)
DN14.k.qtl<-read.csv("DN14_height_blup_known_qtl_concatenated_summary_table.csv", header=F)
DR13.k.qtl<-read.csv("DR13_height_blup_known_qtl_concatenated_summary_table.csv", header=F)
DR14.k.qtl<-read.csv("DR14_height_blup_known_qtl_concatenated_summary_table.csv", header=F)
DL13.k.qtl<-read.csv("DL13_traits_qtl_concatenated_summary_table.csv", header=F)
DL13.k.qtl<-DL13.k.qtl[DL13.k.qtl$V9 == 'plant_height_DL13',]
BP14.k.qtl<-read.csv("BP14_height_above_bound_qtl_concatenated_summary_table.csv", header=F)

# Put all the QTL together
all.k.qtl<-rbind(DN13.k.qtl, DN14.k.qtl, DR13.k.qtl, DR14.k.qtl, BP14.k.qtl)
all.k.qtl<-all.k.qtl[,-c(1)]
colnames(all.k.qtl)<-c("marker", "chr", "pos", "lod", "prop.var", "fx.size", "fx.se", "trait", "treatment", "exp", "year", "type")

# Lets prepare one just for a static final timepoint analysis
all.k.qtl_static<-rbind(DN13.k.qtl, DN14.k.qtl, DR13.k.qtl, DR14.k.qtl, BP14.k.qtl, DL13.k.qtl)
all.k.qtl_static<-all.k.qtl_static[,-c(1)]
colnames(all.k.qtl_static)<-c("marker", "chr", "pos", "lod", "prop.var", "fx.size", "fx.se", "trait", "treatment", "exp", "year", "type")


# Add column for day
day<-c()
for(i in 1:nrow(all.k.qtl)) {
  r<-all.k.qtl[i,]
  d<-strsplit(as.character(r$trait), "_")[[1]][2]
  day<-c(day, d)
}

all.k.qtl<-cbind(all.k.qtl, day)

#################################################################################
# Lets examine the architecture of the 9 QTL on the 'last day' 
# Same as listed in the % variance section
#################################################################################
last_days<-c("height_67_DN13", "height_67_DN14", "height_67_DR13", "height_47_DR14", "height_30_BP14", "plant_height_DL13")

all.k.qtl_static_final<-all.k.qtl_static[all.k.qtl_static$trait %in% last_days, ]
known.qtl_prop.table_final<-c()
for(tr in last_days){
  temp<-all.k.qtl_static_final[all.k.qtl_static_final$trait == tr, ]
  treatments<-unique(temp$treatment)
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    #markers<-temp2$marker
    prop.var<-temp2$prop.var
    logical<-temp2$fx.size < 0
    for(l in 1:length(logical)){
      if (logical[l] == "TRUE") {
        prop.var[l]<-prop.var[l] * -1
      }
    }
    line_entry<-c(paste(unique(temp2$exp), unique(temp2$year), sep=""), t, prop.var)
    known.qtl_prop.table_final<-rbind(known.qtl_prop.table_final, line_entry)
  }
}

markers<-unique(as.character(all.k.qtl_static_final$marker))
colnames(known.qtl_prop.table_final)<-c("Experiment", "Treatment", markers)
row.names(known.qtl_prop.table_final)<-c(1:nrow(known.qtl_prop.table_final))

known.qtl_prop.table_final<-as.data.frame(known.qtl_prop.table_final)
for(i in 3:ncol(known.qtl_prop.table_final)){
  known.qtl_prop.table_final[,i]<-as.numeric(as.character(known.qtl_prop.table_final[,i]))
}

total_exp<-c()
for(i in 1:nrow(known.qtl_prop.table_final)){
  pve<-sum(abs(known.qtl_prop.table_final[i,3:ncol(known.qtl_prop.table_final)]))
  total_exp<-c(total_exp, pve)
}

known.qtl_prop.table_final<-cbind(known.qtl_prop.table_final, total_exp)

treatments<-as.character(unique(known.qtl_prop.table_final$Treatment))
ave_pve_marker<-c()
for(t in treatments) {
  temp<-known.qtl_prop.table_final[known.qtl_prop.table_final$Treatment == t,]
  tr_ave_pve_marker<-c()
  for(i in 3:ncol(temp)){
    pve<-mean(abs(as.numeric(as.character(temp[,i]))))
    tr_ave_pve_marker<-c(tr_ave_pve_marker, pve)
  }
  tr_ave_pve_marker<-c("Average", t, tr_ave_pve_marker)
  ave_pve_marker<-rbind(ave_pve_marker, tr_ave_pve_marker)
}

grand_mean<-c()
for(i in 3:ncol(known.qtl_prop.table_final)){
  pve<-mean(abs(known.qtl_prop.table_final[,i]))
  grand_mean<-c(grand_mean, pve)
}

grand_mean<-c("Average", "Total", grand_mean)
grand_mean<-data.frame(matrix(grand_mean, nrow=1))
colnames(grand_mean)<-colnames(known.qtl_prop.table_final)
ave_pve_marker<-data.frame(matrix(ave_pve_marker, nrow=4))
colnames(ave_pve_marker)<-colnames(known.qtl_prop.table_final)
known.qtl_prop.table_final<-rbind(known.qtl_prop.table_final, ave_pve_marker, grand_mean)

library(gplots)
qtl_prop.table.4.plot<-known.qtl_prop.table_final[,3:ncol(known.qtl_prop.table_final)]
row.names(qtl_prop.table.4.plot)<-paste(known.qtl_prop.table_final$Experiment, known.qtl_prop.table_final$Treatment, sep=".")
for(i in 1:ncol(qtl_prop.table.4.plot)){
  qtl_prop.table.4.plot[,i]<-as.numeric(qtl_prop.table.4.plot[,i])
}

qtl_prop.table.4.plot[13:17,c(5:7)]<-qtl_prop.table.4.plot[13:17,c(5:7)]*-1
qtl_prop.table.4.plot<-as.matrix(qtl_prop.table.4.plot)
qtl_prop.table.4.plot<-round(qtl_prop.table.4.plot, 1)
colnames(qtl_prop.table.4.plot)[10]<-c("All Markers")

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")
pdf("Figure_5.pdf")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 200)
heatmap.2(qtl_prop.table.4.plot, cellnote=qtl_prop.table.4.plot,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, margins=c(10,10))
dev.off()

#################################################################################
# Two different objective: 
# 1) Get total varaince explained by 9 known QTL 
# 2) Plot % variance explained for each QTL through each growout
#################################################################################
# Convert character columns from factor to character
for(i in c(1:4,7:12)){
  all.k.qtl[,i]<-as.character(all.k.qtl[,i])
}

# Convert numeric columns from factor to numeric
for(i in c(5:6,13)){
  all.k.qtl[,i]<-as.numeric(as.character(all.k.qtl[,i]))
}

#################################################################################
# Objective 1
#################################################################################
total.prop.var<-c()
traits<-unique(all.k.qtl$trait)
for(t in traits) {
  temp1<-all.k.qtl[all.k.qtl$trait == t,]
  treats<-unique(temp1$treatment)
  for (tr in treats){
    temp2<-temp1[temp1$treatment == tr,]
    newline<-c('total', 'NA', 'NA', 'NA', sum(temp2$prop.var, na.rm=T), sum(temp2$fx.size, na.rm=T), 'NA', as.character(t), as.character(tr), as.character(unique(temp2$exp)), unique(temp2$year), 'raw', unique(temp2$day))
    temp2<-rbind(temp2, newline)
    total.prop.var<-rbind(total.prop.var, temp2)
  }
}


total.only<-total.prop.var[total.prop.var$marker == 'total',]
total.only$day<-as.numeric(total.only$day)
total.only$exp<-paste(total.only$exp, total.only$year, sep="")
total.only<-total.only[order(total.only$treatment, total.only$exp, total.only$day),]
total.only<-total.only[,c(5,6,8:13)]

write.csv(total.only, file="total_variance_explained_by_major_QTL.csv", quote=F, row.names=F)


#################################################################################
# Objective 2
#################################################################################
markers<-unique(all.k.qtl$marker)

m.var<-c()
m.fx<-c()
for(m in 1:length(markers)){
  mark<-markers[m]
  temp<-all.k.qtl[all.k.qtl$marker == mark,]
  temp$exp<-paste(temp$exp, temp$year, sep="")
  temp<-temp[order(temp$treatment, temp$exp, temp$day),]
  colnames(temp)[c(5,6)]<-c(paste(as.character(mark), 'prop.var', sep="_"), paste(as.character(mark), 'fx.size', sep="_"))
  temp2<-temp[,c(5,8:13)]
  temp3<-temp[,c(6,8:13)]
  if(m == 1) {
    m.var<-temp2
    m.fx<-temp3
  }
  if(m > 1) {
    m.var<-merge(m.var, temp2, by=c('trait', 'treatment', 'exp', 'year','type','day'))
    m.fx<-merge(m.fx, temp3, by=c('trait', 'treatment', 'exp', 'year','type','day'))
  }
}

colnames(m.var)[7:15]<-as.character(markers)
rownames(m.var)<-paste(m.var$trait, m.var$treatment, sep="_")
m.var<-m.var[order(m.var$treatment, m.var$exp, m.var$day),]
write.csv(m.var, file="variance_attributed_to_major_QTL_by_trait.csv", quote=F, row.names=F)
colnames(m.fx)[7:15]<-as.character(markers)
rownames(m.fx)<-paste(m.fx$trait, m.fx$treatment, sep="_")
m.fx<-m.fx[order(m.fx$treatment, m.fx$exp, m.fx$day),]

pdf("variance_attributed_to_major_QTL_heatmap.pdf")
heatmap(as.matrix(m.var[,7:15]), Rowv = NA, Colv = NA)
dev.off()
heatmap(as.matrix(m.fx[,7:15]), Rowv = NA, Colv = NA)



