# Read in libraries
library(ggplot2)
library(stringr)

############################################################################
# Build plotting functions
############################################################################

plot_genotype_3_loci<-function(input, AAA, BBB, big, small, name1, name2, treatment){
  t<-treatment
  smallies<-input[input$id %in% small & input$treatment == t,]
  smallies$cat<-rep(name2, nrow(smallies))
  biggie<-input[input$id %in% big & input$treatment == t,]
  biggie$cat<-rep(name1, nrow(biggie)) 
  a10<-input[input$id %in% AAA & input$treatment == t,]
  a10$cat<-rep("AAA", nrow(a10)) 
  b100<-input[input$id %in% BBB & input$treatment == t,]
  b100$cat<-rep("BBB", nrow(b100)) 
  plotter<-rbind(smallies,biggie,a10, b100)
  plotter<-plotter[,c(ncol(plotter),1:(ncol(plotter)-1))]
  out<-c()
  for(i in (10:ncol(plotter))){
    col<-colnames(plotter)[i]
    day<-strsplit(col, '_')[[1]][2]
    height<-plotter[,i]
    temp<-plotter[,c(1:9)]
    #temp<-plotter[,c(1:8,14)]
    day<-rep(day, nrow(temp))
    temp2<-cbind(temp, day, height)
    out<-rbind(out, temp2)
  }
  out$day<-as.numeric(as.character(out$day))
  out$height<-as.numeric(as.character(out$height))
  return(out)
}

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



############################################################################
# Load in data
############################################################################

BP14_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound/ril_height_above_bound_qtl.csv")
DR14_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field/DR14.height.blup_qtl.csv")
DR13_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field/DR13.height.blup_qtl.csv")
DN14_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field/DN14.height.blup_qtl.csv")
DN13_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/illinois_field/DN13.height.blup_qtl.csv")

BP14_rate_qtl<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound/BP14_AGR_qtl.csv")


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


############################################################################
# Lets plot how the height of different genotypes change through time in different treatment blocks
############################################################################

########################### Start with Bellweather experiment

map<-read.csv("/Users/mfeldman/Dropbox/setaria_height_paper/results/genetic_map/GBS_map_A10xB100_v0.96.csv")

marker1<-c("S5_41999990")
marker2<-c("S2_45043544")
marker3<-c("S9_6724364")

B100.m1<-as.character(map[map[,marker1] == 'BB','id'])
A10.m1<-as.character(map[map[,marker1] == 'AA','id'])

B100.m2<-as.character(map[map[,marker2] == 'BB','id'])
A10.m2<-as.character(map[map[,marker2] == 'AA','id'])

B100.m3<-as.character(map[map[,marker3] == 'BB','id'])
A10.m3<-as.character(map[map[,marker3] == 'AA','id'])

# 3 QTL
AAA<-A10.m3[A10.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]
BBB<-B100.m3[B100.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]
# Largest (A10 on 5@100, B100 on 2@91, B100 on 9@35)
BAB<-B100.m3[B100.m3 %in% (A10.m1[A10.m1 %in% B100.m2])]
# Smallest (B100 on 5@100, A10 on 2@91, A10 on 9@35)
ABA<-A10.m3[A10.m3 %in% (B100.m1[B100.m1 %in% A10.m2])]

# Assign directory to write to
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")

# First wet
t<-c('wet')
output<-plot_genotype_3_loci(BP14_qtl, AAA, BBB, BAB, ABA, 'ABB', 'BAA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "ABA", "BAB", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
pdf("Figure_9a.pdf")
print(p)
dev.off()

# Second dry
t<-c('dry')
output<-plot_genotype_3_loci(BP14_qtl, AAA, BBB, BAB, ABA, 'ABB', 'BAA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "ABA", "BAB", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
pdf("Figure_9b.pdf")
print(p)
dev.off()

# Now rate in wet condition
t<-c('wet')
output<-plot_genotype_3_loci(BP14_rate_qtl, AAA, BBB, BAB, ABA, 'ABB', 'BAA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "ABA", "BAB", "BBB")) + ylab("Vertical growth rate (cm/day)") + xlab("Days after planting") + theme_bw() 
pdf("Figure_9c.pdf")
print(p)
dev.off()

# Now rate in dry condition
t<-c('dry')
output<-plot_genotype_3_loci(BP14_rate_qtl, AAA, BBB, BAB, ABA, 'ABB', 'BAA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "ABA", "BAB", "BBB")) + ylab("Vertical growth rate (cm/day)") + xlab("Days after planting") + theme_bw() 
pdf("Figure_9d.pdf")
print(p)
dev.off()



########################### Drought 2014
# DR14
colnames(S.DR14.w)[1]<-c('marker')
colnames(S.DR14.d)[1]<-c('marker')
S.DR14<-rbind(S.DR14.d, S.DR14.w)

S.DR14.st<-unify_marker(S.DR14)
S.DR14.mark<-unique(S.DR14.st$marker)
lods<-c()
for(m in S.DR14.mark){
  lod<-mean(S.DR14.st[S.DR14.st$marker == m, 'lod'])
  lods<-c(lods, lod)
}
names(lods)<-S.DR14.mark
lods<-sort(lods, decreasing=T)

# Examine average significance of QTL 
lods

# Assign markers
marker1<-names(lods)[1] # B100 -
marker2<-names(lods)[2] # B100 + 
marker3<-names(lods)[3] # B100 -

B100.m1<-as.character(map[map[,marker1] == 'BB','id'])
A10.m1<-as.character(map[map[,marker1] == 'AA','id'])

B100.m2<-as.character(map[map[,marker2] == 'BB','id'])
A10.m2<-as.character(map[map[,marker2] == 'AA','id'])

B100.m3<-as.character(map[map[,marker3] == 'BB','id'])
A10.m3<-as.character(map[map[,marker3] == 'AA','id'])

# 3 QTL
AAA<-A10.m3[A10.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]
BBB<-B100.m3[B100.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]
# Largest (A10 on 5@100, B100 on 9@37, A10 on 7@99)
big<-A10.m3[A10.m3 %in% (A10.m1[A10.m1 %in% B100.m2])]
# Smallest (B100 on 5@100, A10 on 9@37, B100 on 7@99)
small<-B100.m3[B100.m3 %in% (B100.m1[B100.m1 %in% A10.m2])]

t<-c('wet')
output<-plot_genotype_3_loci(DR14_qtl, AAA, BBB, big, small, 'ABA', 'BAB', t)
p<-ggplot(output, aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"))  + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)

t<-c('dry')
output<-plot_genotype_3_loci(DR14_qtl, AAA, BBB, big, small, 'ABA', 'BAB', t)
p<-ggplot(output, aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)


########################### Drought 2013
# DR13
colnames(S.DR13.w)[1]<-c('marker')
colnames(S.DR13.d)[1]<-c('marker')
S.DR13<-rbind(S.DR13.d, S.DR13.w)

S.DR13.st<-unify_marker(S.DR13)
S.DR13.mark<-unique(S.DR13.st$marker)
lods<-c()
for(m in S.DR13.mark){
  lod<-mean(S.DR13.st[S.DR13.st$marker == m, 'lod'])
  lods<-c(lods, lod)
}

names(lods)<-S.DR13.mark
lods<-sort(lods, decreasing=T)

# Examine average significance of QTL 
lods

# Assign markers
marker1<-names(lods)[1]
marker2<-names(lods)[2]
marker3<-names(lods)[3]


B100.m1<-as.character(map[map[,marker1] == 'BB','id'])
A10.m1<-as.character(map[map[,marker1] == 'AA','id'])

B100.m2<-as.character(map[map[,marker2] == 'BB','id'])
A10.m2<-as.character(map[map[,marker2] == 'AA','id'])

B100.m3<-as.character(map[map[,marker3] == 'BB','id'])
A10.m3<-as.character(map[map[,marker3] == 'AA','id'])


# 3 QTL
AAA<-A10.m3[A10.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]
BBB<-B100.m3[B100.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]
# Largest (A10 on 5@100, B100 on 2@78, B100 on 2@1.2)
big<-A10.m3[A10.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]
# Smallest
small<-B100.m3[B100.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]


t<-c('wet')
output<-plot_genotype_3_loci(DR13_qtl, AAA, BBB, big, small, 'ABB', 'BAA', t)
p<-ggplot(output[output$day != 92,], aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "AAB", "BBA", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)

t<-c('dry')
output<-plot_genotype_3_loci(DR13_qtl, AAA, BBB, big, small, 'BBA', 'AAB', t)
p<-ggplot(output[output$day != 92,], aes(day, height,colour=cat)) + geom_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "AAB", "BBA", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)



########################### Density 2014
# DN14
colnames(S.DN14.s)[1]<-c('marker')
colnames(S.DN14.d)[1]<-c('marker')
S.DN14<-rbind(S.DN14.d, S.DN14.s)

S.DN14.st<-unify_marker(S.DN14)
S.DN14.mark<-unique(S.DN14.st$marker)

S.DN14.st<-unify_marker(S.DN14)
S.DN14.mark<-unique(S.DN14.st$marker)
lods<-c()
for(m in S.DN14.mark){
  lod<-mean(S.DN14.st[S.DN14.st$marker == m, 'lod'])
  lods<-c(lods, lod)
}
names(lods)<-S.DN14.mark
lods<-sort(lods, decreasing=T)

# Examine average significance of QTL 
lods

# Assign markers
marker1<-names(lods)[1]
marker2<-names(lods)[2]
marker3<-names(lods)[3]

B100.m1<-as.character(map[map[,marker1] == 'BB','id'])
A10.m1<-as.character(map[map[,marker1] == 'AA','id'])

B100.m2<-as.character(map[map[,marker2] == 'BB','id'])
A10.m2<-as.character(map[map[,marker2] == 'AA','id'])

B100.m3<-as.character(map[map[,marker3] == 'BB','id'])
A10.m3<-as.character(map[map[,marker3] == 'AA','id'])


# 3 QTL
AAA<-A10.m3[A10.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]
BBB<-B100.m3[B100.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]
# Largest (A10 on 5@100, A10 on 7@99, B100 on 9@37)
big<-B100.m3[B100.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]
# Smallest
small<-A10.m3[A10.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]


t<-c('sparse')
output<-plot_genotype_3_loci(DN14_qtl, AAA, BBB, big, small, 'AAB', 'BBA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + stat_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "AAB", "BBA", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)

t<-c('dense')
output<-plot_genotype_3_loci(DN14_qtl, AAA, BBB, big, small, 'AAB', 'BBA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + stat_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "AAB", "BBA", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)

########################### Density 2013
# DN13
colnames(S.DN13.s)[1]<-c('marker')
colnames(S.DN13.d)[1]<-c('marker')
S.DN13<-rbind(S.DN13.d, S.DN13.s)

S.DN13.st<-unify_marker(S.DN13)
S.DN13.mark<-unique(S.DN13.st$marker)

lods<-c()
for(m in S.DN13.mark){
  lod<-mean(S.DN13.st[S.DN13.st$marker == m, 'lod'])
  lods<-c(lods, lod)
}
names(lods)<-S.DN13.mark
lods<-sort(lods, decreasing=T)

# Examine average significance of QTL 
lods

# Assign markers
marker1<-names(lods)[1]
marker2<-names(lods)[2]
marker3<-names(lods)[3]

B100.m1<-as.character(map[map[,marker1] == 'BB','id'])
A10.m1<-as.character(map[map[,marker1] == 'AA','id'])

B100.m2<-as.character(map[map[,marker2] == 'BB','id'])
A10.m2<-as.character(map[map[,marker2] == 'AA','id'])

B100.m3<-as.character(map[map[,marker3] == 'BB','id'])
A10.m3<-as.character(map[map[,marker3] == 'AA','id'])


# 3 QTL
AAA<-A10.m3[A10.m3 %in% (A10.m1[A10.m1 %in% A10.m2])]
BBB<-B100.m3[B100.m3 %in% (B100.m1[B100.m1 %in% B100.m2])]
# Largest (A10 on 5@100, B100 on 2@77, B100 on 8@36)
big<-B100.m3[B100.m3 %in% (A10.m1[A10.m1 %in% B100.m2])]
# Smallest
small<-A10.m3[A10.m3 %in% (B100.m1[B100.m1 %in% A10.m2])]



t<-c('sparse')
output<-plot_genotype_3_loci(DN13_qtl, AAA, BBB, big, small, 'ABB', 'BAA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + stat_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "AAB", "BBA", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)

t<-c('dense')
output<-plot_genotype_3_loci(DN13_qtl, AAA, BBB, big, small, 'ABB', 'BAA', t)
p<-ggplot(output, aes(day, height,colour=cat)) + stat_smooth(size = 1.5) + scale_color_manual(name=c("Genotype"), values=c("violet", "purple4", "dark green", "green"), labels=c("AAA", "AAB", "BBA", "BBB")) + ylab("Height (cm)") + xlab("Days after planting") + theme_bw() 
print(p)



############################################################################
# Lets plot % positive or negative effect size at each time point (From funqtl package results)
############################################################################

# Lets start with Bellweather
##################

S.BP14<-rbind(S.BP14.d, S.BP14.w)
colnames(S.BP14)[1]<-c('marker')


S.BP14$trait<-as.character(S.BP14$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
dap<-c()
for(r in 1:nrow(S.BP14)){
  treat<-strsplit(S.BP14$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.BP14$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.BP14))
  d<-strsplit(trait, "_")[[1]][2]
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  #type<-c(type, t)
  dap<-c(dap, d)
}

exp<-experiment
S.BP14<-cbind(S.BP14, treatment, exp, year, dap)

S.BP14.um<-unify_marker(S.BP14)

treatments<-unique(S.BP14.um$treatment)
days<-unique(S.BP14.um$dap)

output<-c()
for(d in days){
  temp1<-S.BP14.um[S.BP14.um$dap == d,]
  for(t in treatments){
    temp2<-temp1[temp1$treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

output.dry<-output[output$treatment == 'dry',]
output.wet<-output[output$treatment == 'wet',]

pdf("BP14_parental_QTL_contribution.pdf")
p<-ggplot(output.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry")
print(p)
p<-ggplot(output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()

days<-unique(output$dap)
treatments<-unique(output$treatment)
genotypes<-unique(output$genotype)
prop_of_genetic<-c()
for (d in days){
  temp<-output[output$dap == d,]
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    total.var<-sum(temp2$prop.var)
    temp2$prop.var<-temp2$prop.var/total.var
    prop_of_genetic<-rbind(prop_of_genetic, temp2)
  }
}

prop_of_genetic.wet<-prop_of_genetic[prop_of_genetic$treatment == 'wet',]
prop_of_genetic.dry<-prop_of_genetic[prop_of_genetic$treatment == 'dry',]

pdf("Figure_10a.pdf")
p<-ggplot(prop_of_genetic.wet, aes(dap, prop.var, colour=genotype)) + geom_line(size=1.5) + ggtitle("Wet") + ylab("% Variance explained") + xlab("Days after planting") + scale_colour_manual(name=c("Genotype"), values=c("red", "blue"), labels=c("A10", "B100"))
print(p)
dev.off()

pdf("Figure_10b.pdf")
p<-ggplot(prop_of_genetic.dry, aes(dap, prop.var, colour=genotype)) + geom_line(size=1.5) + ggtitle("Dry") + ylab("% Variance explained") + xlab("Days after planting") + scale_colour_manual(name=c("Genotype"), values=c("red", "blue"), labels=c("A10", "B100"))
print(p)
dev.off()

# Now BP14_rate
##################

S.BP14_rate<-rbind(S.BP14_rate.d, S.BP14_rate.w)
colnames(S.BP14_rate)[1]<-c('marker')


S.BP14_rate$trait<-as.character(S.BP14_rate$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
dap<-c()
for(r in 1:nrow(S.BP14_rate)){
  treat<-strsplit(S.BP14_rate$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.BP14_rate$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.BP14_rate))
  d<-strsplit(trait, "_")[[1]][2]
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  #type<-c(type, t)
  dap<-c(dap, d)
}

exp<-experiment
S.BP14_rate<-cbind(S.BP14_rate, treatment, exp, year, dap)

S.BP14_rate.um<-unify_marker(S.BP14_rate)

treatments<-unique(S.BP14_rate.um$treatment)
days<-unique(S.BP14_rate.um$dap)

output<-c()
for(d in days){
  temp1<-S.BP14_rate.um[S.BP14_rate.um$dap == d,]
  for(t in treatments){
    temp2<-temp1[temp1$treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

output.dry<-output[output$treatment == 'dry',]
output.wet<-output[output$treatment == 'wet',]

pdf("BP14_rate_parental_QTL_contribution.pdf")
p<-ggplot(output.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("Dry") + ylab("% Variance explained") + xlab("Days after planting") + scale_colour_manual(name=c("Genotype"), values=c("red", "blue"), labels=c("A10", "B100"))
print(p)
p<-ggplot(output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("Wet") + ylab("% Variance explained") + xlab("Days after planting") + scale_colour_manual(name=c("Genotype"), values=c("red", "blue"), labels=c("A10", "B100"))
print(p)
dev.off()

days<-unique(output$dap)
treatments<-unique(output$treatment)
genotypes<-unique(output$genotype)
prop_of_genetic<-c()
for (d in days){
  temp<-output[output$dap == d,]
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    total.var<-sum(temp2$prop.var)
    temp2$prop.var<-temp2$prop.var/total.var
    prop_of_genetic<-rbind(prop_of_genetic, temp2)
  }
}

prop_of_genetic.wet<-prop_of_genetic[prop_of_genetic$treatment == 'wet',]
prop_of_genetic.dry<-prop_of_genetic[prop_of_genetic$treatment == 'dry',]

pdf("BP14_rate_parental_prop_of_genetic_dry.pdf")
p<-ggplot(prop_of_genetic.dry, aes(dap, prop.var, colour=genotype)) + geom_line(size = 1.5) + ggtitle("Dry") + ylab("% Variance explained") + xlab("Days after planting") + scale_colour_manual(name=c("Genotype"), values=c("red", "blue"), labels=c("A10", "B100"))
print(p)
dev.off()
pdf("BP14_rate_parental_prop_of_genetic_wet.pdf")
p<-ggplot(prop_of_genetic.wet, aes(dap, prop.var, colour=genotype)) + geom_line(size = 1.5) + ggtitle("Wet") + ylab("% Variance explained") + xlab("Days after planting") + scale_colour_manual(name=c("Genotype"), values=c("red", "blue"), labels=c("A10", "B100"))
print(p)
dev.off()


# Now DR14
##################

S.DR14<-rbind(S.DR14.d, S.DR14.w)
colnames(S.DR14)[1]<-c('marker')


S.DR14$trait<-as.character(S.DR14$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
dap<-c()
for(r in 1:nrow(S.DR14)){
  treat<-strsplit(S.DR14$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.DR14$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.DR14))
  d<-strsplit(trait, "_")[[1]][2]
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  #type<-c(type, t)
  dap<-c(dap, d)
}

exp<-experiment
S.DR14<-cbind(S.DR14, treatment, exp, year, dap)

S.DR14.um<-unify_marker(S.DR14)

treatments<-unique(S.DR14.um$treatment)
days<-unique(S.DR14.um$dap)

output<-c()
for(d in days){
  temp1<-S.DR14.um[S.DR14.um$dap == d,]
  for(t in treatments){
    temp2<-temp1[temp1$treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

output.dry<-output[output$treatment == 'dry',]
output.wet<-output[output$treatment == 'wet',]

pdf("DR14_parental_QTL_contribution.pdf")
p<-ggplot(output.dry, aes(dap, prop.var, colour=genotype)) + geom_line(size = 1.5) + ggtitle("Dry") + scale_colour_manual(name=c("Genotype"), values = c("red", "blue"), labels=c("A10", "B100")) + xlab("Days after planting") + ylab("% of total variance explained") + scale_y_continuous(limits=c(0,1))

print(p)
p<-ggplot(output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()

days<-unique(output$dap)
treatments<-unique(output$treatment)
genotypes<-unique(output$genotype)
prop_of_genetic<-c()
for (d in days){
  temp<-output[output$dap == d,]
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    total.var<-sum(temp2$prop.var)
    temp2$prop.var<-temp2$prop.var/total.var
    prop_of_genetic<-rbind(prop_of_genetic, temp2)
  }
}

prop_of_genetic.wet<-prop_of_genetic[prop_of_genetic$treatment == 'wet',]
prop_of_genetic.dry<-prop_of_genetic[prop_of_genetic$treatment == 'dry',]

pdf("DR14_parental_prop_of_genetic.pdf")
p<-ggplot(prop_of_genetic.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry")
print(p)
p<-ggplot(prop_of_genetic.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()


# Now DR13
##################

S.DR13<-rbind(S.DR13.d, S.DR13.w)
colnames(S.DR13)[1]<-c('marker')


S.DR13$trait<-as.character(S.DR13$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
dap<-c()
for(r in 1:nrow(S.DR13)){
  treat<-strsplit(S.DR13$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.DR13$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.DR13))
  d<-strsplit(trait, "_")[[1]][2]
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  #type<-c(type, t)
  dap<-c(dap, d)
}

exp<-experiment
S.DR13<-cbind(S.DR13, treatment, exp, year, dap)

S.DR13.um<-unify_marker(S.DR13)

treatments<-unique(S.DR13.um$treatment)
days<-unique(S.DR13.um$dap)

output<-c()
for(d in days){
  temp1<-S.DR13.um[S.DR13.um$dap == d,]
  for(t in treatments){
    temp2<-temp1[temp1$treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

output.dry<-output[output$treatment == 'dry',]
output.wet<-output[output$treatment == 'wet',]

pdf("DR13_parental_QTL_contribution.pdf")
p<-ggplot(output.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry")
print(p)
p<-ggplot(output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()

days<-unique(output$dap)
treatments<-unique(output$treatment)
genotypes<-unique(output$genotype)
prop_of_genetic<-c()
for (d in days){
  temp<-output[output$dap == d,]
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    total.var<-sum(temp2$prop.var)
    temp2$prop.var<-temp2$prop.var/total.var
    prop_of_genetic<-rbind(prop_of_genetic, temp2)
  }
}

prop_of_genetic.wet<-prop_of_genetic[prop_of_genetic$treatment == 'wet',]
prop_of_genetic.dry<-prop_of_genetic[prop_of_genetic$treatment == 'dry',]

pdf("DR13_parental_prop_of_genetic.pdf")
p<-ggplot(prop_of_genetic.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry")
print(p)
p<-ggplot(prop_of_genetic.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()


# Now DN14
##################

S.DN14<-rbind(S.DN14.d, S.DN14.s)
colnames(S.DN14)[1]<-c('marker')


S.DN14$trait<-as.character(S.DN14$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
dap<-c()
for(r in 1:nrow(S.DN14)){
  treat<-strsplit(S.DN14$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.DN14$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.DN14))
  d<-strsplit(trait, "_")[[1]][2]
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  #type<-c(type, t)
  dap<-c(dap, d)
}

exp<-experiment
S.DN14<-cbind(S.DN14, treatment, exp, year, dap)

S.DN14.um<-unify_marker(S.DN14)

treatments<-unique(S.DN14.um$treatment)
days<-unique(S.DN14.um$dap)

output<-c()
for(d in days){
  temp1<-S.DN14.um[S.DN14.um$dap == d,]
  for(t in treatments){
    temp2<-temp1[temp1$treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

output.dense<-output[output$treatment == 'dense',]
output.sparse<-output[output$treatment == 'sparse',]

pdf("DN14_parental_QTL_contribution.pdf")
p<-ggplot(output.dense, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dense")
print(p)
p<-ggplot(output.sparse, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("sparse")
print(p)
dev.off()

days<-unique(output$dap)
treatments<-unique(output$treatment)
genotypes<-unique(output$genotype)
prop_of_genetic<-c()
for (d in days){
  temp<-output[output$dap == d,]
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    total.var<-sum(temp2$prop.var)
    temp2$prop.var<-temp2$prop.var/total.var
    prop_of_genetic<-rbind(prop_of_genetic, temp2)
  }
}

prop_of_genetic.dense<-prop_of_genetic[prop_of_genetic$treatment == 'dense',]
prop_of_genetic.sparse<-prop_of_genetic[prop_of_genetic$treatment == 'sparse',]

pdf("DN14_parental_prop_of_genetic.pdf")
p<-ggplot(prop_of_genetic.dense, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry")
print(p)
p<-ggplot(prop_of_genetic.sparse, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()


# Now DN13
##################

S.DN13<-rbind(S.DN13.d, S.DN13.s)
colnames(S.DN13)[1]<-c('marker')


S.DN13$trait<-as.character(S.DN13$trait)
treatment<-c()
experiment<-c()
year<-c()
type<-c()
dap<-c()
for(r in 1:nrow(S.DN13)){
  treat<-strsplit(S.DN13$trait[r], "\\.")[[1]][1]
  trait<-strsplit(S.DN13$trait[r], "\\.")[[1]][2]
  exp<-substr(trait, nchar(trait)-3, nchar(trait))
  yr<-substr(trait, nchar(trait)-1, nchar(trait))
  t<-rep('raw', nrow(S.DN13))
  d<-strsplit(trait, "_")[[1]][2]
  
  treatment<-c(treatment, treat)
  experiment<-c(experiment, exp)
  year<-c(year, yr)
  #type<-c(type, t)
  dap<-c(dap, d)
}

exp<-experiment
S.DN13<-cbind(S.DN13, treatment, exp, year, dap)

S.DN13.um<-unify_marker(S.DN13)

treatments<-unique(S.DN13.um$treatment)
days<-unique(S.DN13.um$dap)

output<-c()
for(d in days){
  temp1<-S.DN13.um[S.DN13.um$dap == d,]
  for(t in treatments){
    temp2<-temp1[temp1$treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

output.dense<-output[output$treatment == 'dense',]
output.sparse<-output[output$treatment == 'sparse',]

pdf("DN13_parental_QTL_contribution.pdf")
p<-ggplot(output.dense, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dense")
print(p)
p<-ggplot(output.sparse, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("sparse")
print(p)
dev.off()

days<-unique(output$dap)
treatments<-unique(output$treatment)
genotypes<-unique(output$genotype)
prop_of_genetic<-c()
for (d in days){
  temp<-output[output$dap == d,]
  for(t in treatments){
    temp2<-temp[temp$treatment == t,]
    total.var<-sum(temp2$prop.var)
    temp2$prop.var<-temp2$prop.var/total.var
    prop_of_genetic<-rbind(prop_of_genetic, temp2)
  }
}

prop_of_genetic.dense<-prop_of_genetic[prop_of_genetic$treatment == 'dense',]
prop_of_genetic.sparse<-prop_of_genetic[prop_of_genetic$treatment == 'sparse',]

pdf("DN13_parental_prop_of_genetic.pdf")
p<-ggplot(prop_of_genetic.dense, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry")
print(p)
p<-ggplot(prop_of_genetic.sparse, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet")
print(p)
dev.off()