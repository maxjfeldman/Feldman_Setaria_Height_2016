library(ggplot2, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(MASS, warn.conflicts = FALSE)
library(car, warn.conflicts = FALSE)
library(nlme, warn.conflicts = FALSE)
library(mvtnorm, warn.conflicts = FALSE)
library(grid, warn.conflicts = FALSE)

# Set active directory
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril")

# Add planting data
planting_date = as.POSIXct("2014-1-13")
data<-read.csv("sv.vis_data_ril.csv")
data$date<-as.POSIXct(data$timestamp, origin = "1970-01-01")
data$dap <- as.numeric(data$date - planting_date)
data$dap_i<-as.integer(data$dap)


# Read a key that associates phenotyper ID (aka 'plant_id') to RIL id
key<-read.csv('/Users/mfeldman/Desktop/setaria_ril/phe2ril2.csv', header=F)
colnames(key)<-c('phe', 'ril')

# Make the RIL name same as in QTL map file
key$ril<-sprintf("%03d", key$ril)
key$ril<-paste('RIL_', key$ril, sep="")

# Add a new (empty column called genotype)
data$genotype<-NA
# Add genotypes 
for(i in 1:nrow(key)) {
  phe_name<-key[i,1]
  ril_name<-key[i,2]
  data$genotype[grep(paste('Dr', phe_name, 'A', sep=""), data$plantbarcode)]<-ril_name
}

# Add parental genotypes
data$genotype[grep('Dp1', data$plantbarcode)]<-c("A10")
data$genotype[grep('Dp2', data$plantbarcode)]<-c("B100")

# Add treatment
data$treatment<-c('none')
data$treatment[grep("AA", data$plantbarcode)]<-c('wet')
data$treatment[grep("AB", data$plantbarcode)]<-c('dry')

data<-subset(data, select=-c(image_id,run_id,timestamp,lifter,gain,measurementlabel,id,exposure,other,image_id.1))


# Do calibration for area and height 
# Directly from Fahlgren et al., 2015

zoom.lm = lm(zoom.camera ~ zoom, data=data.frame(zoom=c(1,6000),
                                                 zoom.camera=c(1,6)))
summary(zoom.lm)

# Get calibration data
if (!file.exists('/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril/zoom_calibration_data.txt')) {
  download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                'zoom_calibration_data.txt')
}

z.data = read.table(file="/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril/zoom_calibration_data.txt", sep="\t", header=TRUE)
z.data$px_cm = z.data$length_px / z.data$length_cm
z.data$zoom.camera = predict(object = zoom.lm, newdata=z.data)
data$zoom = sub('z', '', data$zoom)
data$zoom<-as.numeric(data$zoom)
data$sv.zoom.camera = predict(object = zoom.lm, newdata=data)


# Look for model that works 
# Exponential first
area.coef = coef(nls(log(rel_area) ~ log(a * exp(b * zoom.camera)),
                     z.data, start = c(a = 1, b = 0.01)))
area.coef = data.frame(a=area.coef[1], b=area.coef[2])
area.nls = nls(rel_area ~ a * exp(b * zoom.camera),
               data = z.data, start=c(a=area.coef$a, b=area.coef$b))
summary(area.nls)

# Polynomial
area.pol = lm(rel_area ~ zoom.camera + I(zoom.camera^2), z.data)
summary(area.pol)

# Exponential
len.coef = coef(nls(log(px_cm) ~ log(a * exp(b * zoom.camera)),
                    z.data[z.data$camera == 'VIS SV',], start = c(a = 1, b = 0.01)))
len.coef = data.frame(a=len.coef[1], b=len.coef[2])
len.nls = nls(px_cm ~ a * exp(b * zoom.camera),
              data = z.data[z.data$camera == 'VIS SV',],
              start=c(a=len.coef$a, b=len.coef$b))
summary(len.nls)

# Polynomial
len.poly = lm(px_cm ~ zoom.camera + I(zoom.camera^2),
              data=z.data[z.data$camera == 'VIS SV',])
summary(len.poly)


########## CONVERT THE VIS DATA ##########
data.zoom = data
data$zoom.camera = data$sv.zoom.camera

data$zoom.camera = data$sv.zoom.camera
data$px_cm = predict(object = len.poly, newdata=data)
data.zoom$extent_x = data$width / data$px_cm
data.zoom$extent_y = data$height / data$px_cm
data.zoom$height_above_bound = data$height_above_bound / data$px_cm
data.zoom$center.of.mass.y = data$center.of.mass.y / data$px_cm
data.zoom$ellipse_center_y = data$ellipse_center_y / data$px_cm
data.zoom$canopy_height = data$canopy_height / data$px_cm
data.zoom$canopy_width = data$canopy_width / data$px_cm



# Lets save this file as a .csv
write.csv(data.zoom, file="sv.vis_data_ril.csv", row.names=F, quote=F)
