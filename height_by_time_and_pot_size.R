library(ggplot2)

setwd("/Users/mfeldman/Desktop/jurkowski")

heights<-read.csv('mj_plant_heights.csv')
# Lets look at the structure of the data...
heights[1:5,]

# First lets remove the last three columns
dim(heights)
heights<-heights[,c(1:19)]

# Lets rename the columns...
c.names<-c("plant_id", "pot_size", "7", "10", "12", "14", "17", "19","21","24","26","28","31","33","35","39","46","48","50")
colnames(heights)<-c.names

# Remove first row
heights<-heights[-c(1),]

# Get genotype and plant id from the 'plant_id' field
genotype<-c()
plant<-c()
for (r in 1:nrow(heights)) {
  id<-as.character(heights[r,'plant_id'])
  s.id<-strsplit(id, " ")
  g<-s.id[[1]][1]
  p<-s.id[[1]][3]
  
  genotype<-c(genotype, g)
  plant<-c(plant, p)
}

heights$plant_id<-plant
heights$genotype<-genotype

heights$pot_size<-as.character(heights$pot_size)
heights[heights$pot_size == "Mini Pot",'pot_size']<-c("mini_pot")
heights[heights$pot_size == "Phenotyper Pot",'pot_size']<-c("bellweather_pot")
heights[heights$pot_size == "Tall Pot",'pot_size']<-c("tree_pot")
heights[heights$pot_size == "2 Gallon Pot",'pot_size']<-c("two_gallon_pot")

heights[heights$genotype == "A10",'genotype']<-c("S. viridis")
heights[heights$genotype == "B100",'genotype']<-c("S. italica")


# Lets re-order columns and then convert to long form
heights<-heights[,c(20,1:19)]

heights.l<-c()
for(i in 4:ncol(heights)){
  h<-heights[,c(1:3,i)]
  d<-colnames(heights)[i]
  dap<-rep(d, nrow(h))
  colnames(h)[4]<-c('height')
  height_day<-cbind(h,dap)
  heights.l<-rbind(heights.l, height_day)
}

heights.l$height<-as.numeric(as.character(heights.l$height))
heights.l$dap<-as.numeric(as.character(heights.l$dap))
heights.l$pot_size<-as.factor(heights.l$pot_size)
heights.l$pot_size<-factor(heights.l$pot_size,levels(heights.l$pot_size)[c(2,3,4,1)])

p <- ggplot(heights.l, aes(dap, height, color=factor(genotype))) + geom_point() + facet_wrap(~pot_size)
p <- ggplot(heights.l, aes(dap, height, color=factor(pot_size))) + geom_point() + facet_wrap(~genotype)


p <- ggplot(heights.l, aes(dap, height, color=factor(genotype))) + geom_smooth() + facet_wrap(~pot_size)  + ylab("Height (mm)") + xlab("Days after planting") + theme_bw()
pdf("Figure_10.pdf")
p + labs(color = "Genotype")
dev.off()



p <- ggplot(heights.l, aes(dap, height, color=factor(pot_size))) + geom_smooth() + facet_wrap(~genotype) + theme_bw() + ylab("Height (mm)") + xlab("Days after planting") 
pdf("Figure_S12.pdf")
p + labs(color = "Pot size")
dev.off()