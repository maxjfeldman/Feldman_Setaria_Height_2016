# This is a script to subset a genetic map file by a QTL confidence interval
# and sort the map by phenotypic value

# Load dependencies
library(qtl)

# Read in data for trait of interest from summary table from the wet/well-watered condition (we'll call it cr_table for crown root table)
cr_table<-read.csv("/Users/mfeldman/Desktop/summary.table.mqm.wet.height_10_BP14.csv")

# Here are two ways to do this... 
# I recommend the second method.

# Lets look at both...

#############################################################################################
# Method 1: Use external map and trait file
#############################################################################################

# Read in genetic map and trait file
# The path will be different on your computer

map<-read.csv("/Users/mfeldman/Desktop/example_for_MC/GBS_map_A10xB100_v0.96.csv")
trait<-read.csv("/Users/mfeldman/Desktop/ril_height_above_bound_qtl.csv")

# It looks like the pipeline returned both marker names on each side of the confidence interval 
# The markers that denote the confidence interval can be found in fields: L.CI_marker and R.CI_marker
# This is ideal... There is a bug in R/qtl in which sometimes you only get a cM position and it gives a strange value marker name
# We'll cover how to deal with this in Method 2 below, but this time it is not a problem

# Lets make the chart for the QTL on Chromosome 6
# Get the marker on the left border, it is represented as a factor but you want it represented as a character string.. as.character()
lb_marker<-as.character(cr_table[cr_table$chr == 5,'L.CI_marker'])
# Now do the same thing for the right border
rb_marker<-as.character(cr_table[cr_table$chr == 5,'R.CI_marker'])

# Now find the index that is represented by the left border marker and right border marker
lb_index<-which(colnames(map)==lb_marker)
rb_index<-which(colnames(map)==rb_marker)

# Lets subset the map to get only markers spanning these indices... Need the first column as this contains the genotype ids
qtl_interval_on_5<-map[,c(1, lb_index:rb_index)]

# This new map subset is 16 markers long
dim(qtl_interval_on_5)

# Okay, lets add in the phenotypes now, remember we're only dealing with the wet treatment results so lets remove the rest
# We'll need the id/genotype and the trait of interest id column is column number 7 and height on day 10 is number 12
colnames(trait)
cr_number.wet<-trait[trait$treatment == 'wet', c(7,11)]
# Lets take an average based upon genotype
cr_number.wet<-aggregate(cr_number.wet[,2], by=list(cr_number.wet$id), FUN=mean, na.action=na.omit)
# Re-add the column names
colnames(cr_number.wet)<-c("id", "height")

# merge the markers and phenotypes
qtl_interval_on_5<-merge(qtl_interval_on_5, cr_number.wet, by="id")

# Finally lets sort the table by phenotype...
qtl_interval_on_5<-qtl_interval_on_5[order(qtl_interval_on_5$height), ]

# Write to file
setwd("/Users/mfeldman/Desktop")
write.csv(qtl_interval_on_5, "Method_1.qtl_interval_on_5.csv", quote=F, row.names=F)

qtl_interval_on_5<-qtl_interval_on_5[-is.na(qtl_interval_on_5$height),]

# Clear workspace
rm(list=ls())
