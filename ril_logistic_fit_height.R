# Script to model hight of RILs
library(ggplot2)
library(lattice)
library(lme4)
library(nlme)
library(mvtnorm)

# Read in .csv
setwd("/Users/mfeldman/Desktop/setaria_height_paper/data/bellweather_ril")
ril<-read.csv("sv.vis_data_ril.csv", colClasses=c('treatment'='character'))

# Redo first few steps from variance scripts

# Need to remove plants out of bound
ril<-ril[ril$in_bound == 'True',]

# Lets remove the entry named '0370' as it has no treatment or information about genotype
ril<-ril[ril$plantbarcode != '0370',]

# Lets remove all columns that arn't to do with height
ril<-subset(ril, select=-c(camera, imgtype, area, hull.area, solidity, perimeter, width, height, longest_axis, center.of.mass.x, hull_vertices, ellipse_center_x, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity, y.position, height_below_bound,above_bound_area,percent_above_bound_area, below_bound_area,percent_below_bound_area,date,dap,extent_x))

# This order makes more sense
ril<-ril[,c(4,2,13,5,12,3,12,1,11,15,6,8,9,10,7)]

# From the previous section we believe that the mean of 'height_above_bound' is the best measure
# Perform this operation
ril.mean<-aggregate(ril[,9:13], by=list(ril$plantbarcode, ril$cartag, ril$genotype, ril$treatment, ril$dap_i), mean)
colnames(ril.mean)[1:5]<-c('plant_id', 'cartag','genotype','treatment','dap_i')

# Lets do our best to remove empty cars...
plot(ril.mean$height_above_bound~ril.mean$dap_i, pch=20, cex=0.3, ylab="height", xlab="day")

# Need to zoom in 17 - 32 dap_i
close_up<-ril.mean[ril.mean$dap_i > 16,]
plot(close_up$height_above_bound~close_up$dap_i, ylim=c(0,25), pch=20, cex=0.3)
filter_1<-close_up[close_up$height_above_bound > 5,]
plot(filter_1$height_above_bound~filter_1$dap_i, ylim=c(0,25), pch=20, cex=0.3)
filter_2<-filter_1[(filter_1$height_above_bound > 12 & (filter_1$dap_i == 28 | filter_1$dap_i == 29)),'plant_id']
filter_2<-filter_1[filter_1$plant_id %in% filter_2, ]
plot(filter_2$height_above_bound~filter_2$dap_i, ylim=c(0,25), pch=20, cex=0.3)

# Now remove all of those images from RIL mean that are no good
ril.mean<-ril.mean[ril.mean$plant_id %in% filter_2$plant_id, ]
plot(ril.mean$height_above_bound~ril.mean$dap_i, pch=20, cex=0.3, ylab="height", xlab="day")

write.csv(ril.mean, file='ril_bellweather.mean.csv',quote=F,row.names=F)
ril.mean$treatment<-as.character(ril.mean$treatment)

# Lets try modeling this
rils<-unique(ril.mean$genotype)
treatments<-unique(ril.mean$treatment)

#####################################################################################
# Use Gompertz
#####################################################################################
ril.gomp<-ril.mean

rils<-sort(unique(ril.gomp$genotype))
treatments<-unique(ril.gomp$treatment)
dap_i<-unique(ril.gomp$dap_i)

ril_final_parameters.gomp<-c()
ril_final_model_fit.gomp<-c()
ril_final_rates.gomp<-c()
problem_lines_gomp<-c()

for (i in 1:length(rils)) {
  r<-rils[i]
  temp1<-ril.gomp[ril.gomp$genotype == as.character(r),]
  per.ril<-c()
  for (j in 1:length(treatments)) {
    t<-treatments[j]
    per.t<-c()
    temp2<-temp1[temp1$treatment == as.character(t),]
    if (nrow(temp2) < 1) {next;}
    #temp2<-temp2[complete.cases(temp2),]
    # Log of 0 is INF need to replace with another small #
    temp2$height_above_bound<-replace(temp2$height_above_bound, temp2$height==0, 1)
    #tmp.logis  <- getInitial(biomass ~ SSfpl(dap_i, Asym, xmid, scal), data = temp2, control= list(warnOnly = T))
    #rough.logis <- nls(biomass ~ SSlogis(dap_i, Asym, xmid, scal), trace = F, start = tmp.logis, na.action=na.omit, control = list(maxiter=500, warnOnly = T), data = temp2)
    #rough.logis <- nls(height ~ SSlogis(dap_i, Asym, xmid, scal), trace = F, na.action=na.omit, control = list(maxiter=500, warnOnly = T), data = temp2)
    #fit.logis <- gnls(height ~ SSlogis(dap_i, Asym, xmid, scal), data = temp2, na.action=na.omit, start = coef(rough.logis), weights= varExp(form = ~ fitted(.)), control=list(tolerance = 100))
    #out.gomp <- output.logis.gnls(fit.logis, sort(temp2$dap_i), T)
    #tmp.gomp <- getInitial(height ~ SSgompertz(dap_i, Asym, b2, b3), data = temp2)
    #fit.gomp <- gnls(height ~ SSgompertz(dap_i, Asym, b2, b3), data = temp2, weights= varExp(form = ~ fitted(.)))
    error_check<-tryCatch(fit.gomp <- gnls(height_above_bound ~ SSgompertz(dap_i, Asym, b2, b3), data = temp2, weights= varExp(form = ~ fitted(.))) ,
                           warning = function(w) {print('Warning!'); return('warning')},
                           error = function(e) {print(paste("Non-convergence of line:", r, "treatment ", t, sep=" ")); problem_lines_gomp<-c(problem_lines_gomp, paste(r,t)); return('error')},
                           paste("okay"))
    if (error_check == 'error') {next;}
    #out.gomp <- output.gomp.gnls(fit.gomp, sort(temp2$dap_i), CI = T)
    out.gomp <- output.gomp.gnls(fit.gomp, sort(dap_i), CI = T)
    #per.t$ril<-rep(as.character(temp2$genotype[1]), length(unique(temp2$dap_i)))
    per.t$ril<-rep(as.character(temp2$genotype[1]), length(unique(dap_i)))
    #per.t$treatment<-rep(t, length(unique(temp2$dap_i)))
    per.t$treatment<-rep(t, length(unique(dap_i)))
    #per.t$dap_i<-sort((unique(temp2$dap_i)))
    per.t$dap_i<-sort((unique(dap_i)))
    if (length(unique(out.gomp$rates$M))  == length(per.t$ril)) {
      per.t$M<-unique(out.gomp$rates$M)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$M<-length(out.gomp$rates$M[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$M.lo))  == length(per.t$ril)) {
      per.t$M.lo<-unique(out.gomp$rates$M.lo)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$M.lo<-length(out.gomp$rates$M.lo[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$M.hi))  == length(per.t$ril)) {
      per.t$M.hi<-unique(out.gomp$rates$M.hi)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$M.hi<-length(out.gomp$rates$M.hi[seq(1, nrow(out.gomp$rates), count)])
    }
    
    if (length(unique(out.gomp$rates$AGR))  == length(per.t$ril)) {
      per.t$AGR<-unique(out.gomp$rates$AGR)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$AGR<-length(out.gomp$rates$AGR[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$AGR.lo))  == length(per.t$ril)) {
      per.t$AGR.lo<-unique(out.gomp$rates$AGR.lo)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$AGR.lo<-length(out.gomp$rates$AGR.lo[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$AGR.hi))  == length(per.t$ril)) {
      per.t$AGR.hi<-unique(out.gomp$rates$AGR.hi)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$AGR.hi<-length(out.gomp$rates$AGR.hi[seq(1, nrow(out.gomp$rates), count)])
    }
    
    if (length(unique(out.gomp$rates$RGRt))  == length(per.t$ril)) {
      per.t$RGRt<-unique(out.gomp$rates$RGRt)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$RGRt<-length(out.gomp$rates$RGRt[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$RGRt.lo))  == length(per.t$ril)) {
      per.t$RGRt.lo<-unique(out.gomp$rates$RGRt.lo)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$RGRt.lo<-length(out.gomp$rates$RGRt.lo[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$RGRt.hi))  == length(per.t$ril)) {
      per.t$RGRt.hi<-unique(out.gomp$rates$RGRt.hi)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$RGRt.hi<-length(out.gomp$rates$RGRt.hi[seq(1, nrow(out.gomp$rates), count)])
    }
    
    if (length(unique(out.gomp$rates$RGRm))  == length(per.t$ril)) {
      per.t$RGRm<-unique(out.gomp$rates$RGRm)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$RGRm<-length(out.gomp$rates$RGRm[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$RGRm.lo))  == length(per.t$ril)) {
      per.t$RGRm.lo<-unique(out.gomp$rates$RGRm.lo)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$RGRm.lo<-length(out.gomp$rates$RGRm.lo[seq(1, nrow(out.gomp$rates), count)])
    }
    if (length(unique(out.gomp$rates$RGRm.hi))  == length(per.t$ril)) {
      per.t$RGRm.hi<-unique(out.gomp$rates$RGRm.hi)
    } else {
      count<-table(out.gomp$rates$times)[1]
      per.t$RGRm.hi<-length(out.gomp$rates$RGRm.hi[seq(1, nrow(out.gomp$rates), count)])
    } 
    per.t<-as.data.frame(per.t)
    ril_final_rates.gomp<-rbind(ril_final_rates.gomp, per.t)
    # Get parameters
    parmas<-out.gomp$params
    parmas<-c(as.character(r),as.character(t),parmas)
    names(parmas)[1:2]<-c("genotype", "treatment")
    ril_final_parameters.gomp<-rbind(ril_final_parameters.gomp, parmas)
    # Get model fit summary
    fit_sum<-out.gomp$summary
    fit_sum<-c(as.character(r),as.character(t),fit_sum)
    names(fit_sum)[1:2]<-c("genotype", "treatment")
    ril_final_model_fit.gomp<-rbind(ril_final_model_fit.gomp, fit_sum)
  }
}

setwd("/Users/mfeldman/Desktop/setaria_height_paper/results/predicted_height/bellweather_ril")
write.csv(ril_final_model_fit.gomp, 'ril_final_model_fit.gomp.csv', row.names=F, quote=F)
write.csv(ril_final_parameters.gomp, 'ril_final_parameters.gomp.csv', row.names=F, quote=F)
write.csv(ril_final_rates.gomp, 'ril_final_rates.gomp.csv', row.names=F, quote=F)


#####################################################################################
# Use 3-parameter logistic
#####################################################################################

ril.logis<-ril.mean

# Lets try modeling this
rils<-sort(unique(ril.logis$genotype))
treatments<-unique(ril.logis$treatment)
problem_lines_logis<-c()
ril_final_rates.logis<-c()
ril_final_parameters.logis<-c()
ril_final_model_fit.logis<-c()
dap_i<-ril.mean$dap_i

for (i in 1:length(rils)) {
  r<-rils[i]
  temp1<-ril.logis[ril.logis$genotype == as.character(r),]
  per.ril<-c()
  for (j in 1:length(treatments)) {
    t<-treatments[j]
    per.t<-c()
    temp2<-temp1[temp1$treatment == as.character(t),]
    if (nrow(temp2) < 1) {next;}
    #temp2<-temp2[complete.cases(temp2),]
    # Log of 0 is INF need to replace with another small #
    temp2$height_above_bound<-replace(temp2$height_above_bound, temp2$height_above_bound==0, 1)
    error_check<-c('empty')
    #tmp.logis  <- getInitial(biomass ~ SSfpl(dap_i, Asym, xmid, scal), data = temp2, control= list(warnOnly = T))
    #rough.logis <- nls(biomass ~ SSlogis(dap_i, Asym, xmid, scal), trace = F, start = tmp.logis, na.action=na.omit, control = list(maxiter=500, warnOnly = T), data = temp2)
    #rough.logis <- nls(height ~ SSfpl(dap_i, A, B, xmid, scal), trace = F, na.action=na.omit, control = list(maxiter=500, warnOnly = T), data = temp2)
    #fit.logis <- gnls(height_above_bound ~ SSlogis(dap_i, Asym, xmid, scal), data = temp2, na.action=na.omit, start = coef(rough.logis.a10.d), weights= varExp(form = ~ fitted(.)), control=list(tolerance = 100))
    #fit.logis <- gnls(height ~ SSlogis(dap_i, Asym, xmid, scal), data = temp2, na.action=na.omit, weights= varExp(form = ~ fitted(.)), control=list(tolerance = 100))
    print(c(as.character(r),i))
    print(c(t,j))
    error_check<-tryCatch( fit.logis <- gnls(height_above_bound ~ SSlogis(dap_i, Asym, xmid, scal), data = temp2, weights= varExp(form = ~ fitted(.))) ,
                           warning = function(w) {print('Warning!'); return('warning')},
                           error = function(e) {print(paste("Non-convergence of line:", r, "treatment ", t, sep=" ")); problem_lines_logis<-c(problem_lines_logis, paste(r,t)); return('error')},
                           #finally = print("okay"))
                           paste("okay"))
    if (error_check == 'error') {next;}
    
    out.logis <- output.logis.gnls(fit.logis, sort(temp2$dap_i), T)
    out.logis <- output.logis.gnls(fit.logis, sort(dap_i), T)
    #per.t$ril<-rep(as.character(temp2$genotype[1]), length(unique(temp2$dap_i)))
    per.t$ril<-rep(as.character(temp2$genotype[1]), length(unique(dap_i)))
    #per.t$treatment<-rep(t, length(unique(temp2$dap_i)))
    per.t$treatment<-rep(t, length(unique(dap_i)))
    #per.t$dap_i<-sort((unique(temp2$dap_i)))
    per.t$dap_i<-sort((unique(dap_i)))
    if (length(unique(out.logis$rates$M))  == length(per.t$ril)) {
      per.t$M<-unique(out.logis$rates$M)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$M<-length(out.logis$rates$M[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$M.lo))  == length(per.t$ril)) {
      per.t$M.lo<-unique(out.logis$rates$M.lo)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$M.lo<-length(out.logis$rates$M.lo[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$M.hi))  == length(per.t$ril)) {
      per.t$M.hi<-unique(out.logis$rates$M.hi)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$M.hi<-length(out.logis$rates$M.hi[seq(1, nrow(out.logis$rates), count)])
    }
    
    if (length(unique(out.logis$rates$AGR))  == length(per.t$ril)) {
      per.t$AGR<-unique(out.logis$rates$AGR)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$AGR<-length(out.logis$rates$AGR[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$AGR.lo))  == length(per.t$ril)) {
      per.t$AGR.lo<-unique(out.logis$rates$AGR.lo)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$AGR.lo<-length(out.logis$rates$AGR.lo[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$AGR.hi))  == length(per.t$ril)) {
      per.t$AGR.hi<-unique(out.logis$rates$AGR.hi)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$AGR.hi<-length(out.logis$rates$AGR.hi[seq(1, nrow(out.logis$rates), count)])
    }
    
    if (length(unique(out.logis$rates$RGRt))  == length(per.t$ril)) {
      per.t$RGRt<-unique(out.logis$rates$RGRt)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$RGRt<-length(out.logis$rates$RGRt[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$RGRt.lo))  == length(per.t$ril)) {
      per.t$RGRt.lo<-unique(out.logis$rates$RGRt.lo)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$RGRt.lo<-length(out.logis$rates$RGRt.lo[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$RGRt.hi))  == length(per.t$ril)) {
      per.t$RGRt.hi<-unique(out.logis$rates$RGRt.hi)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$RGRt.hi<-length(out.logis$rates$RGRt.hi[seq(1, nrow(out.logis$rates), count)])
    }
    
    if (length(unique(out.logis$rates$RGRm))  == length(per.t$ril)) {
      per.t$RGRm<-unique(out.logis$rates$RGRm)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$RGRm<-length(out.logis$rates$RGRm[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$RGRm.lo))  == length(per.t$ril)) {
      per.t$RGRm.lo<-unique(out.logis$rates$RGRm.lo)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$RGRm.lo<-length(out.logis$rates$RGRm.lo[seq(1, nrow(out.logis$rates), count)])
    }
    if (length(unique(out.logis$rates$RGRm.hi))  == length(per.t$ril)) {
      per.t$RGRm.hi<-unique(out.logis$rates$RGRm.hi)
    } else {
      count<-table(out.logis$rates$times)[1]
      per.t$RGRm.hi<-length(out.logis$rates$RGRm.hi[seq(1, nrow(out.logis$rates), count)])
    }
    per.t<-as.data.frame(per.t)
    ril_final_rates.logis<-rbind(ril_final_rates.logis, per.t)
    # Get parameters
    parmas<-out.logis$params
    parmas<-c(as.character(r),as.character(t),parmas)
    names(parmas)[1:2]<-c("genotype", "treatment")
    ril_final_parameters.logis<-rbind(ril_final_parameters.logis, parmas)
    # Get model fit summary
    fit_sum<-out.logis$summary
    fit_sum<-c(as.character(r),as.character(t),fit_sum)
    names(fit_sum)[1:2]<-c("genotype", "treatment")
    ril_final_model_fit.logis<-rbind(ril_final_model_fit.logis, fit_sum)
  }
}

setwd("/Users/mfeldman/Desktop/setaria_height_paper/results/predicted_height/bellweather_ril")
write.csv(ril_final_model_fit.logis, 'ril_final_model_fit.logis.csv', row.names=F, quote=F)
write.csv(ril_final_parameters.logis, 'ril_final_parameters.logis.csv', row.names=F, quote=F)
write.csv(ril_final_rates.logis, 'ril_final_rates.logis.csv', row.names=F, quote=F)



#####################################################################################
# Use 4-parameter logistic
#####################################################################################


ril.fpl<-ril.mean
#ril.fpl<-ril.gomp

# Lets try modeling this
rils<-sort(unique(ril.fpl$genotype))
treatments<-unique(ril.fpl$treatment)
problem_lines_fpl<-c()
ril_final_rates.fpl<-c()
ril_final_parameters.fpl<-c()
ril_final_model_fit.fpl<-c()

for (i in 1:length(rils)) {
  r<-rils[i]
  temp1<-ril.fpl[ril.fpl$genotype == r,]
  per.ril<-c()
  for (j in 1:length(treatments)) {
    t<-treatments[j]
    per.t<-c()
    temp2<-temp1[temp1$treatment == t,]
    if (nrow(temp2) < 1) {next;}
    temp2<-temp2[complete.cases(temp2),]
    # Log of 0 is INF need to replace with another small #
    temp2$height_above_bound<-replace(temp2$height_above_bound, temp2$height_above_bound==0, 1)
    #tmp.fpl  <- getInitial(biomass ~ SSfpl(dap_i, Asym, xmid, scal), data = temp2, control= list(warnOnly = T))
    #rough.fpl <- nls(biomass ~ SSlogis(dap_i, Asym, xmid, scal), trace = F, start = tmp.fpl, na.action=na.omit, control = list(maxiter=500, warnOnly = T), data = temp2)
    #rough.fpl <- nls(height ~ SSfpl(dap_i, A, B, xmid, scal), trace = F, na.action=na.omit, control = list(maxiter=500, warnOnly = T), data = temp2)
    error_check<-c('empty')
    error_check<-tryCatch( fit.fpl <- gnls(height_above_bound ~ SSfpl(dap_i, A, B, xmid, scal), data = temp2, na.action=na.omit, weights= varExp(form = ~ fitted(.)), control=list(tolerance = 100)) ,
                           warning = function(w) {print('Warning!'); return('warning')},
                           error = function(e) {print(paste("Non-convergence of line:", r, "treatment ", t, sep=" ")); problem_lines_fpl<-c(problem_lines_fpl, paste(r,t)); return('error')},
                           #finally = print("okay"))
                           paste("okay"))
    if (error_check == 'error') {next;}
    #fit.fpl <- gnls(height ~ SSfpl(dap_i, A, B, xmid, scal), data = temp2, na.action=na.omit, weights= varExp(form = ~ fitted(.)), control=list(tolerance = 100))
    out.fpl <- output.fpl.gnls(fit.fpl, sort(temp2$dap_i), T)
    per.t$ril<-rep(as.character(temp2$genotype[1]), length(unique(temp2$dap_i)))
    per.t$treatment<-rep(t, length(unique(temp2$dap_i)))
    per.t$dap_i<-sort((unique(temp2$dap_i)))
    if (length(unique(out.fpl$rates$M))  == length(per.t$ril)) {
      per.t$M<-unique(out.fpl$rates$M)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$M<-length(out.fpl$rates$M[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$M.lo))  == length(per.t$ril)) {
      per.t$M.lo<-unique(out.fpl$rates$M.lo)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$M.lo<-length(out.fpl$rates$M.lo[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$M.hi))  == length(per.t$ril)) {
      per.t$M.hi<-unique(out.fpl$rates$M.hi)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$M.hi<-length(out.fpl$rates$M.hi[seq(1, nrow(out.fpl$rates), count)])
    }
    
    if (length(unique(out.fpl$rates$AGR))  == length(per.t$ril)) {
      per.t$AGR<-unique(out.fpl$rates$AGR)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$AGR<-length(out.fpl$rates$AGR[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$AGR.lo))  == length(per.t$ril)) {
      per.t$AGR.lo<-unique(out.fpl$rates$AGR.lo)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$AGR.lo<-length(out.fpl$rates$AGR.lo[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$AGR.hi))  == length(per.t$ril)) {
      per.t$AGR.hi<-unique(out.fpl$rates$AGR.hi)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$AGR.hi<-length(out.fpl$rates$AGR.hi[seq(1, nrow(out.fpl$rates), count)])
    }
    
    if (length(unique(out.fpl$rates$RGRt))  == length(per.t$ril)) {
      per.t$RGRt<-unique(out.fpl$rates$RGRt)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$RGRt<-length(out.fpl$rates$RGRt[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$RGRt.lo))  == length(per.t$ril)) {
      per.t$RGRt.lo<-unique(out.fpl$rates$RGRt.lo)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$RGRt.lo<-length(out.fpl$rates$RGRt.lo[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$RGRt.hi))  == length(per.t$ril)) {
      per.t$RGRt.hi<-unique(out.fpl$rates$RGRt.hi)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$RGRt.hi<-length(out.fpl$rates$RGRt.hi[seq(1, nrow(out.fpl$rates), count)])
    }
    
    if (length(unique(out.fpl$rates$RGRm))  == length(per.t$ril)) {
      per.t$RGRm<-unique(out.fpl$rates$RGRm)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$RGRm<-length(out.fpl$rates$RGRm[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$RGRm.lo))  == length(per.t$ril)) {
      per.t$RGRm.lo<-unique(out.fpl$rates$RGRm.lo)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$RGRm.lo<-length(out.fpl$rates$RGRm.lo[seq(1, nrow(out.fpl$rates), count)])
    }
    if (length(unique(out.fpl$rates$RGRm.hi))  == length(per.t$ril)) {
      per.t$RGRm.hi<-unique(out.fpl$rates$RGRm.hi)
    } else {
      count<-table(out.fpl$rates$times)[1]
      per.t$RGRm.hi<-length(out.fpl$rates$RGRm.hi[seq(1, nrow(out.fpl$rates), count)])
    }
    per.t<-as.data.frame(per.t)
    ril_final_rates.fpl<-rbind(ril_final_rates.fpl, per.t)
    # Get parameters
    parmas<-out.fpl$params
    parmas<-c(as.character(r),as.character(t),parmas)
    names(parmas)[1:2]<-c("genotype", "treatment")
    ril_final_parameters.fpl<-rbind(ril_final_parameters.fpl, parmas)
    # Get model fit summary
    fit_sum<-out.fpl$summary
    fit_sum<-c(as.character(r),as.character(t),fit_sum)
    names(fit_sum)[1:2]<-c("genotype", "treatment")
    ril_final_model_fit.fpl<-rbind(ril_final_model_fit.fpl, fit_sum)
  }
}

# Write to .csv
setwd("/Users/mfeldman/Desktop/setaria_height_paper/results/predicted_height/bellweather_ril")
write.csv(ril_final_model_fit.fpl, 'ril_final_model_fit.fpl.csv', row.names=F, quote=F)
write.csv(ril_final_parameters.fpl, 'ril_final_parameters.fpl.csv', row.names=F, quote=F)
write.csv(ril_final_rates.fpl, 'ril_final_rates.fpl.csv', row.names=F, quote=F)



###########################################################################################


###########
# Compare models
###########
setwd("/Users/mfeldman/Desktop/setaria_height_paper/results/predicted_height/bellweather_ril")

ril_final_model_fit.logis<-read.csv("ril_final_model_fit.logis.csv")
ril_final_model_fit.gomp<-read.csv("ril_final_model_fit.gomp.csv")
ril_final_model_fit.fpl<-read.csv("ril_final_model_fit.fpl.csv")
ril_final_parameters.logis<-read.csv("ril_final_parameters.logis.csv")
ril_final_parameters.gomp<-read.csv("ril_final_parameters.gomp.csv")
ril_final_parameters.fpl<-read.csv("ril_final_parameters.fpl.csv")
ril_final_rates.logis<-read.csv("ril_final_rates.logis.csv")
ril_final_rates.gomp<-read.csv("ril_final_rates.gomp.csv")
ril_final_rates.fpl<-read.csv("ril_final_rates.fpl.csv")


ril_final_model_fit.logis<-as.data.frame(ril_final_model_fit.logis)
colnames(ril_final_model_fit.logis)[c(3:5)]<-c("R2.logis", "AIC.logis", "RMSE.logis")

ril_final_model_fit.fpl<-as.data.frame(ril_final_model_fit.fpl)
colnames(ril_final_model_fit.fpl)[c(3:5)]<-c("R2.fpl", "AIC.fpl", "RMSE.fpl")

ril_final_model_fit.gomp<-as.data.frame(ril_final_model_fit.gomp)
colnames(ril_final_model_fit.gomp)[c(3:5)]<-c("R2.gomp", "AIC.gomp", "RMSE.gomp")

# Merge all 3 tables into a common table
ril_final_model_fit_table<-merge(ril_final_model_fit.logis, ril_final_model_fit.gomp, by=c('genotype', 'treatment'), all=T)
ril_final_model_fit_table<-merge(ril_final_model_fit_table, ril_final_model_fit.fpl, by=c('genotype', 'treatment'), all=T)


# Fetch the model with lowest AIC score
ril_fit_selection<-c()
for(i in 1:nrow(ril_final_model_fit_table)) {
  r<-ril_final_model_fit_table[i,c(4,7,10)]
  r<-as.numeric(as.character(unlist(r)))
  names(r)<-c('logis', 'gomp', 'fpl')
  c<-which.min(r)
  selected_fit<-c(as.character(ril_final_model_fit_table[i,1]), as.character(ril_final_model_fit_table[i,2]), names(c))
  ril_fit_selection<-rbind(ril_fit_selection, selected_fit)
}

rownames(ril_fit_selection)<-c(1:nrow(ril_fit_selection))
ril_fit_selection<-as.data.frame(ril_fit_selection)
ril_wet<-ril_fit_selection[ril_fit_selection$V2 == 'wet',]
ril_dry<-ril_fit_selection[ril_fit_selection$V2 == 'dry',]

table(ril_wet$V3)
table(ril_dry$V3)


# Select best model and write to a dataframe
ril_gomp.best<-ril_fit_selection[ril_fit_selection$V3 == 'gomp',c(1,2)]
ril_logis.best<-ril_fit_selection[ril_fit_selection$V3 == 'logis',c(1,2)]
ril_fpl.best<-ril_fit_selection[ril_fit_selection$V3 == 'fpl',c(1,2)]


# Extract value and rate estimates from gomp
ril_best.fit_rates<-c()
for(i in 1:nrow(ril_gomp.best)) {
  l<-ril_final_rates.gomp[as.character(ril_final_rates.gomp$ril) == as.character(ril_gomp.best[i,1]) & as.character(ril_final_rates.gomp$treatment) == as.character(ril_gomp.best[i,2]),]
  ril_best.fit_rates<-rbind(ril_best.fit_rates, l)
}

for(i in 1:nrow(ril_logis.best)) {
  l<-ril_final_rates.logis[as.character(ril_final_rates.logis$ril) == as.character(ril_logis.best[i,1]) & as.character(ril_final_rates.logis$treatment) == as.character(ril_logis.best[i,2]),]
  ril_best.fit_rates<-rbind(ril_best.fit_rates, l)
}

for(i in 1:nrow(ril_fpl.best)) {
  l<-ril_final_rates.fpl[as.character(ril_final_rates.fpl$ril) == as.character(ril_fpl.best[i,1]) & as.character(ril_final_rates.fpl$treatment) == as.character(ril_fpl.best[i,2]),]
  ril_best.fit_rates<-rbind(ril_best.fit_rates, l)
}


write.csv(ril_best.fit_rates, file=c("ril_best.fit_logistic_estimates_height.csv"), quote=F, row.names=F)

##################################################
# LETS MAKE PLOTS
##################################################
setwd( "/Users/mfeldman/Dropbox/setaria_height_paper/data/bellweather_ril")
ril.mean<-read.csv('ril_bellweather.mean.csv')
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
ril_best.fit_rates<-read.csv("ril_best.fit_logistic_estimates_height.csv")


rate_id<-sort(unique(as.character(ril_best.fit_rates$ril)))
growth_rate_report<-c()
pdf("ril_best.fit_logistic_estimates_height.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for(i in 1:length(rate_id)) {
  r<-rate_id[i]
  ril_set<-ril.mean[(ril.mean$genotype == r),]
  ril_rates<-ril_best.fit_rates[(ril_best.fit_rates$ril == r),]
  plant_ids<-unique(ril_set$plant_id)
  max.b<-max(max(ril_set$height), max(ril_rates$M.hi))
  set<-ril_set[ril_set$plant_id == plant_ids[1],]
  # Start making plots
  if (set[1,'treatment'] == 'wet') {l.color<-c("light blue")}
  if (set[1,'treatment'] == 'dry') {l.color<-c("gold")}
  plot(set$height~set$dap_i, type='l', xlim=c(8,33), ylim=c(0, max.b), col=l.color, xlab="Days after planting", ylab="height", main=r)
  if(length(plant_ids) >1) {
    for (j in 2:length(plant_ids)) {
      set<-ril_set[ril_set$plant_id == plant_ids[j],]
      if (set[1,'treatment'] == 'wet') {l.color<-c("light blue")}
      if (set[1,'treatment'] == 'dry') {l.color<-c("gold")}
      points(set$height~set$dap_i, type='l', xlim=c(8,33), ylim=c(0, max.b), col=l.color, xlab="Day", ylab="height")
    }
  }
  
  rate.w<-ril_rates[ril_rates$treatment == 'wet', ]
  if (nrow(rate.w) > 0) {
    l.color<-c("blue")
    p.color<-c("dark blue")
    max.rate.w<-max(rate.w$AGR)
    day.w<-rate.w[rate.w$AGR == max.rate.w, 'dap_i']
    max.val.w<-rate.w[rate.w$AGR == max.rate.w, 'M'] 
    lines(rate.w$M~rate.w$dap_i, lwd=2, col=l.color)
    points(day.w, max.val.w, cex=1.5, col=p.color, pch=18)
    #lines(rate.w$M.hi~rate.w$dap_i, lty=2, col=c('navy'))
    #lines(rate.w$M.lo~rate.w$dap_i, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$M.lo, rev(rate.w$M.hi)), col=COL.CI, lty=2)
  }
  
  rate.d<-ril_rates[ril_rates$treatment == 'dry', ]  
  if (nrow(rate.d) > 0) {
    l.color<-c("orange")
    p.color<-c("dark orange")
    max.rate.d<-max(rate.d$AGR)
    day.d<-rate.d[rate.d$AGR == max.rate.d, 'dap_i']
    max.val.d<-rate.d[rate.d$AGR == max.rate.d, 'M'] 
    lines(rate.d$M~rate.d$dap_i, lwd=2, col=l.color)
    points(day.d, max.val.d, cex=1.5, col=p.color, pch=18)
    #lines(rate.d$M.hi~rate.d$dap_i, lty=2, col=c('brown'))
    #lines(rate.d$M.lo~rate.d$dap_i, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$M.lo, rev(rate.d$M.hi)), col=COL.CI)
    
  }
  
  # treatment drought
  rate.d<-ril_rates[ril_rates$treatment == 'dry', ]  
  if (nrow(rate.d) > 0) {
    max.height.d<-max(rate.d$M)
    max.agr.d<-max(rate.d$AGR)
    max.agr.day.d<-rate.d[rate.d$AGR == max.agr.d, 'dap_i']
    max.agr.height.d<-rate.d[rate.d$AGR == max.agr.d, 'M']
    rgrt.agr.max.d<-rate.d[rate.d$AGR == max.agr.d, 'RGRt']
    max.rgrt.d<-max(rate.d$RGRt)
    max.rgrt.day.d<-rate.d[rate.d$RGRt == max.rgrt.d, 'dap_i']
    min.rgrt.d<-min(rate.d$RGRt)
    min.rgrt.day.d<-rate.d[rate.d$RGRt == min.rgrt.d, 'dap_i']
    rgrm.agr.max.d<-rate.d[rate.d$AGR == max.agr.d, 'RGRm']
    max.rgrm.d<-max(rate.d$RGRm)
    max.rgrm.day.d<-rate.d[rate.d$RGRm == max.rgrm.d, 'dap_i']
    min.rgrm.d<-min(rate.d$RGRm)
    min.rgrm.day.d<-rate.d[rate.d$RGRm == min.rgrm.d, 'dap_i']
  }
  # treatment well watered
  rate.w<-ril_rates[ril_rates$treatment == 'wet', ]
  if (nrow(rate.w) > 0) {
    max.height.w<-max(rate.w$M)
    max.agr.w<-max(rate.w$AGR)
    max.agr.day.w<-rate.w[rate.w$AGR == max.agr.w, 'dap_i']
    max.agr.height.w<-rate.w[rate.w$AGR == max.agr.w, 'M']
    rgrt.agr.max.w<-rate.w[rate.w$AGR == max.agr.w, 'RGRt']
    max.rgrt.w<-max(rate.w$RGRt)
    max.rgrt.day.w<-rate.w[rate.w$RGRt == max.rgrt.w, 'dap_i']
    min.rgrt.w<-min(rate.w$RGRt)
    min.rgrt.day.w<-rate.w[rate.w$RGRt == min.rgrt.w, 'dap_i']
    rgrm.agr.max.w<-rate.w[rate.w$AGR == max.agr.w, 'RGRm']
    max.rgrm.w<-max(rate.w$RGRm)
    max.rgrm.day.w<-rate.w[rate.w$RGRm == max.rgrm.w, 'dap_i']
    min.rgrm.w<-min(rate.w$RGRm)
    min.rgrm.day.w<-rate.w[rate.w$RGRm == min.rgrm.w, 'dap_i']
  }
  # Generate the report on a per/ril basis
  if (length(unique(ril_rates$treatment)) > 1) {
    ril_entry<-c(r, max.height.w, max.height.d, max.height.w - max.height.d, max.agr.w, max.agr.d, max.agr.w - max.agr.d, max.agr.day.w, max.agr.day.d, max.agr.day.w - max.agr.day.d, max.agr.height.w, max.agr.height.d, max.agr.height.w - max.agr.height.d)
    ril_entry<-c(ril_entry, rgrt.agr.max.w, rgrt.agr.max.d, rgrt.agr.max.w - rgrt.agr.max.d, max.rgrt.w, max.rgrt.d, max.rgrt.w - max.rgrt.d, min.rgrt.w, min.rgrt.d, min.rgrt.w - min.rgrt.d)
    ril_entry<-c(ril_entry, rgrm.agr.max.w, rgrm.agr.max.d, rgrm.agr.max.w - rgrm.agr.max.d, max.rgrm.w, max.rgrm.d, max.rgrm.w - max.rgrm.d, min.rgrm.w, min.rgrm.d, min.rgrm.w - min.rgrm.d)
    growth_rate_report<-rbind(growth_rate_report, ril_entry)
  }
  #colnames(growth_rate_report)<-c('genotype', 'max_height_wet', 'max_height_dry', 'max_height_diff','max_AGR_wet','max_AGR_dry', 'max_AGR_diff','max_AGR_day_wet','max_AGR_day_dry', 'max_AGR_day_diff', 'height_max_AGR_wet', 'height_max_AGR_dry', 'height_max_AGR_diff', )
  # Plot rates
  if (length(unique(ril_rates$treatment)) > 1) {
    max.r<-max(max(rate.w$AGR), max(rate.d$AGR))  
    plot(rate.w$AGR~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='AGR', ylim=c(0,max.r))
    lines(rate.d$AGR~rate.d$dap_i, col="orange")
    points(max.agr.day.w, max.agr.w, pch=18, col="dark blue")
    points(max.agr.day.d, max.agr.d, pch=18, col="dark orange")
    # Add confidence intervals
    #lines(rate.w$AGR.hi~rate.w$dap_i, lty=2, col=c('navy'))
    #lines(rate.w$AGR.lo~rate.w$dap_i, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$AGR.lo, rev(rate.w$AGR.hi)), col=COL.CI, lty=2)
    #lines(rate.d$AGR.hi~rate.d$dap_i, lty=2, col=c('brown'))
    #lines(rate.d$AGR.lo~rate.d$dap_i, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$AGR.lo, rev(rate.d$AGR.hi)), col=COL.CI, lty=2)
    
    max.RGRt<-max(max(rate.w$RGRt), max(rate.d$RGRt))
    min.RGRt<-min(min(rate.w$RGRt), min(rate.d$RGRt))
    plot(rate.w$RGRt~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='RGRt', ylim=c(min.RGRt,max.RGRt))
    lines(rate.d$RGRt~rate.d$dap_i, col="orange")
    points(max.agr.day.w, rgrt.agr.max.w, pch=18, col="dark blue")
    points(max.agr.day.d, rgrt.agr.max.d, pch=18, col="dark orange")
    # Add confidence intervals
    #lines(rate.w$RGRt.hi~rate.w$dap_i, lty=2, col=c('navy'))
    #lines(rate.w$RGRt.lo~rate.w$dap_i, lty=2, col=c('navy'))
    # polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$RGRt.lo, rev(rate.w$RGRt.hi)), col=COL.CI, lty=2)
    #lines(rate.d$RGRt.hi~rate.d$dap_i, lty=2, col=c('brown'))
    #lines(rate.d$RGRt.lo~rate.d$dap_i, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$RGRt.lo, rev(rate.d$RGRt.hi)), col=COL.CI, lty=2)
    
    max.RGRm<-max(max(rate.w$RGRm), max(rate.d$RGRm))
    min.RGRm<-min(min(rate.w$RGRm), min(rate.d$RGRm))
    plot(rate.w$RGRm~rate.w$M, type="l", col="blue", xlab='height', ylab='RGRm', ylim=c(min.RGRm,max.RGRm))
    lines(rate.d$RGRm~rate.d$M, col="orange")
    points(max.agr.height.w, rgrm.agr.max.w, pch=18, col="dark blue")
    points(max.agr.height.d, rgrm.agr.max.d, pch=18, col="dark orange")
    #lines(rate.w$RGRm.hi~rate.w$M, lty=2, col=c('navy'))
    #lines(rate.w$RGRm.lo~rate.w$M, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$M, rev(rate.w$M)), y=c(rate.w$RGRm.lo, rev(rate.w$RGRm.hi)), col=COL.CI, lty=2)
    #lines(rate.d$RGRm.hi~rate.d$M, lty=2, col=c('brown'))
    #lines(rate.d$RGRm.lo~rate.d$M, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$M, rev(rate.d$M)), y=c(rate.d$RGRm.lo, rev(rate.d$RGRm.hi)), col=COL.CI, lty=2)
  }
  
  if (length(unique(ril_rates$treatment)) == 1) {
    if (nrow(rate.w) > 0) {
      max.agr.w<-max(rate.w$AGR)
      max.agr.day.w<-rate.w[rate.w$AGR == max.agr.w, 'dap_i']
      plot(rate.w$AGR~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='AGR', ylim=c(0,max.agr.w))
      points(max.agr.day.w, max.agr.w, pch=18, col="dark blue")
      #lines(rate.w$AGR.hi~rate.w$dap_i, lty=2, col=c('navy'))
      #lines(rate.w$AGR.lo~rate.w$dap_i, lty=2, col=c('navy'))
      #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$AGR.lo, rev(rate.w$AGR.hi)), col=COL.CI, lty=2)
      
      max.RGRt<-max(rate.w$RGRt)
      min.RGRt<-min(rate.w$RGRt)
      plot(rate.w$RGRt~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='RGRt', ylim=c(min.RGRt,max.RGRt))
      points(max.agr.day.w, rgrt.agr.max.w, pch=18, col="dark blue")
      #lines(rate.w$RGRt.hi~rate.w$dap_i, lty=2, col=c('navy'))
      #lines(rate.w$RGRt.lo~rate.w$dap_i, lty=2, col=c('navy'))
      #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$RGRt.lo, rev(rate.w$RGRt.hi)), col=COL.CI, lty=2)
      
      max.RGRm<-max(rate.w$RGRm)
      min.RGRm<-min(rate.w$RGRm)
      plot(rate.w$RGRm~rate.w$M, type="l", col="blue", xlab='height', ylab='RGRm', ylim=c(min.RGRm,max.RGRm))
      points(max.agr.height.w, rgrm.agr.max.w, pch=18, col="dark blue")
      #lines(rate.w$RGRm.hi~rate.w$M, lty=2, col=c('navy'))
      #lines(rate.w$RGRm.lo~rate.w$M, lty=2, col=c('navy'))
      #polygon(x=c(rate.w$M, rev(rate.w$M)), y=c(rate.w$RGRm.lo, rev(rate.w$RGRm.hi)), col=COL.CI, lty=2)
      
    }
    if (nrow(rate.d) > 0) {
      max.agr.d<-max(rate.d$AGR)
      max.agr.day.d<-rate.d[rate.d$AGR == max.agr.d, 'dap_i']
      plot(rate.d$AGR~rate.d$dap_i, type="l", col="orange", xlab='Days after planting', ylab='AGR', ylim=c(0,max.agr.d))
      points(max.agr.day.d, max.agr.d, pch=18, col="dark orange")
      #lines(rate.d$AGR.hi~rate.d$dap_i, lty=2, col=c('brown'))
      #lines(rate.d$AGR.lo~rate.d$dap_i, lty=2, col=c('brown'))
      #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$AGR.lo, rev(rate.d$AGR.hi)), col=COL.CI, lty=2)
      
      max.RGRt<-max(rate.d$RGRt)
      min.RGRt<-min(rate.d$RGRt)
      plot(rate.d$RGRt~rate.d$dap_i, type="l", col="orange", xlab='Days after planting', ylab='RGRt', ylim=c(min.RGRt,max.RGRt))
      points(max.agr.day.d, rgrt.agr.max.d, pch=18, col="dark orange")
      #lines(rate.d$RGRt.hi~rate.d$dap_i, lty=2, col=c('brown'))
      #lines(rate.d$RGRt.lo~rate.d$dap_i, lty=2, col=c('brown'))
      #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$RGRt.lo, rev(rate.d$RGRt.hi)), col=COL.CI, lty=2)
      
      max.RGRm<-max(rate.d$RGRm)
      min.RGRm<-min(rate.d$RGRm)
      plot(rate.d$RGRm~rate.d$M, type="l", col="orange", xlab='height', ylab='RGRm', ylim=c(min.RGRm,max.RGRm))
      points(max.agr.height.d, rgrm.agr.max.d, pch=18, col="dark orange") 
      #lines(rate.d$RGRm.hi~rate.d$M, lty=2, col=c('brown'))
      #lines(rate.d$RGRm.lo~rate.d$M, lty=2, col=c('brown'))
      #polygon(x=c(rate.d$M, rev(rate.d$M)), y=c(rate.d$RGRm.lo, rev(rate.d$RGRm.hi)), col=COL.CI, lty=2)
    }
  }
}
dev.off()

###### NOT DONE #####
#colnames(dp1_best.fit_rates)[1]<-c("genotype")
#dp1_best.fit_rates<-dp1_best.fit_rates[,c(1:4,7)]
#write.csv(dp1_best.fit_rates, file="dp1_best.fit_logistic_estimates_height.csv", quote=F, row.names=F)




### Lets make some plots of the raw and rate data
setwd("/Users/mfeldman/Desktop/")

rate_id<-sort(unique(as.character(final_rates$ril)))
growth_rate_report<-c()
pdf("RIL_rates_gomp.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for(i in 1:length(rate_id)) {
  r<-rate_id[i]
  rilset<-ril.vis.shape[(ril.vis.shape$genotype == r),]
  rilrates<-final_rates[(final_rates$ril == r),]
  plant_ids<-unique(rilset$plant_id)
  max.b<-max(max(rilset$height_above_bound), max(rilrates$M.hi))
  set<-rilset[rilset$plant_id == plant_ids[1],]
  # Start making plots
  if (set[1,'treatment'] == 'wet') {l.color<-c("light blue")}
  if (set[1,'treatment'] == 'dry') {l.color<-c("gold")}
  plot(set$height_above_bound~set$dap_i, type='l', xlim=c(8,33), ylim=c(0, max.b), col=l.color, xlab="Days after planting", ylab="height_above_bound", main=r)
  if(length(plant_ids) >1) {
    for (j in 2:length(plant_ids)) {
      set<-rilset[rilset$plant_id == plant_ids[j],]
      if (set[1,'treatment'] == 'wet') {l.color<-c("light blue")}
      if (set[1,'treatment'] == 'dry') {l.color<-c("gold")}
      points(set$height_above_bound~set$dap_i, type='l', xlim=c(8,33), ylim=c(0, max.b), col=l.color, xlab="Day", ylab="height_above_bound")
    }
  }
  
  rate.w<-rilrates[rilrates$treamtment == 'wet', ]
  if (nrow(rate.w) > 0) {
    l.color<-c("blue")
    p.color<-c("dark blue")
    max.rate.w<-max(rate.w$AGR)
    day.w<-rate.w[rate.w$AGR == max.rate.w, 'dap_i']
    max.val.w<-rate.w[rate.w$AGR == max.rate.w, 'M'] 
    lines(rate.w$M~rate.w$dap_i, lwd=2, col=l.color)
    points(day.w, max.val.w, cex=1.5, col=p.color, pch=18)
    lines(rate.w$M.hi~rate.w$dap_i, lty=2, col=c('navy'))
    lines(rate.w$M.lo~rate.w$dap_i, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$M.lo, rev(rate.w$M.hi)), col=COL.CI, lty=2)
  }
  
  rate.d<-rilrates[rilrates$treamtment == 'dry', ]  
  if (nrow(rate.d) > 0) {
    l.color<-c("orange")
    p.color<-c("dark orange")
    max.rate.d<-max(rate.d$AGR)
    day.d<-rate.d[rate.d$AGR == max.rate.d, 'dap_i']
    max.val.d<-rate.d[rate.d$AGR == max.rate.d, 'M'] 
    lines(rate.d$M~rate.d$dap_i, lwd=2, col=l.color)
    points(day.d, max.val.d, cex=1.5, col=p.color, pch=18)
    lines(rate.d$M.hi~rate.d$dap_i, lty=2, col=c('brown'))
    lines(rate.d$M.lo~rate.d$dap_i, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$M.lo, rev(rate.d$M.hi)), col=COL.CI)
    
  }
  
  # treatment drought
  rate.d<-rilrates[rilrates$treamtment == 'dry', ]  
  if (nrow(rate.d) > 0) {
    max.height_above_bound.d<-max(rate.d$M)
    max.agr.d<-max(rate.d$AGR)
    max.agr.day.d<-rate.d[rate.d$AGR == max.agr.d, 'dap_i']
    max.agr.height_above_bound.d<-rate.d[rate.d$AGR == max.agr.d, 'M']
    rgrt.agr.max.d<-rate.d[rate.d$AGR == max.agr.d, 'RGRt']
    max.rgrt.d<-max(rate.d$RGRt)
    max.rgrt.day.d<-rate.d[rate.d$RGRt == max.rgrt.d, 'dap_i']
    min.rgrt.d<-min(rate.d$RGRt)
    min.rgrt.day.d<-rate.d[rate.d$RGRt == min.rgrt.d, 'dap_i']
    rgrm.agr.max.d<-rate.d[rate.d$AGR == max.agr.d, 'RGRm']
    max.rgrm.d<-max(rate.d$RGRm)
    max.rgrm.day.d<-rate.d[rate.d$RGRm == max.rgrm.d, 'dap_i']
    min.rgrm.d<-min(rate.d$RGRm)
    min.rgrm.day.d<-rate.d[rate.d$RGRm == min.rgrm.d, 'dap_i']
  }
  # treatment well watered
  rate.w<-rilrates[rilrates$treamtment == 'wet', ]
  if (nrow(rate.w) > 0) {
    max.height_above_bound.w<-max(rate.w$M)
    max.agr.w<-max(rate.w$AGR)
    max.agr.day.w<-rate.w[rate.w$AGR == max.agr.w, 'dap_i']
    max.agr.height_above_bound.w<-rate.w[rate.w$AGR == max.agr.w, 'M']
    rgrt.agr.max.w<-rate.w[rate.w$AGR == max.agr.w, 'RGRt']
    max.rgrt.w<-max(rate.w$RGRt)
    max.rgrt.day.w<-rate.w[rate.w$RGRt == max.rgrt.w, 'dap_i']
    min.rgrt.w<-min(rate.w$RGRt)
    min.rgrt.day.w<-rate.w[rate.w$RGRt == min.rgrt.w, 'dap_i']
    rgrm.agr.max.w<-rate.w[rate.w$AGR == max.agr.w, 'RGRm']
    max.rgrm.w<-max(rate.w$RGRm)
    max.rgrm.day.w<-rate.w[rate.w$RGRm == max.rgrm.w, 'dap_i']
    min.rgrm.w<-min(rate.w$RGRm)
    min.rgrm.day.w<-rate.w[rate.w$RGRm == min.rgrm.w, 'dap_i']
  }
  # Generate the report on a per/ril basis
  if (length(unique(rilrates$treamtment)) > 1) {
    ril_entry<-c(r, max.height_above_bound.w, max.height_above_bound.d, max.height_above_bound.w - max.height_above_bound.d, max.agr.w, max.agr.d, max.agr.w - max.agr.d, max.agr.day.w, max.agr.day.d, max.agr.day.w - max.agr.day.d, max.agr.height_above_bound.w, max.agr.height_above_bound.d, max.agr.height_above_bound.w - max.agr.height_above_bound.d)
    ril_entry<-c(ril_entry, rgrt.agr.max.w, rgrt.agr.max.d, rgrt.agr.max.w - rgrt.agr.max.d, max.rgrt.w, max.rgrt.d, max.rgrt.w - max.rgrt.d, min.rgrt.w, min.rgrt.d, min.rgrt.w - min.rgrt.d)
    ril_entry<-c(ril_entry, rgrm.agr.max.w, rgrm.agr.max.d, rgrm.agr.max.w - rgrm.agr.max.d, max.rgrm.w, max.rgrm.d, max.rgrm.w - max.rgrm.d, min.rgrm.w, min.rgrm.d, min.rgrm.w - min.rgrm.d)
    growth_rate_report<-rbind(growth_rate_report, ril_entry)
  }
  # Plot rates
  if (length(unique(rilrates$treamtment)) > 1) {
    max.r<-max(max(rate.w$AGR.hi), max(rate.d$AGR.hi))  
    plot(rate.w$AGR~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='AGR', ylim=c(0,max.r))
    lines(rate.d$AGR~rate.d$dap_i, col="orange")
    points(max.agr.day.w, max.agr.w, pch=18, col="dark blue")
    points(max.agr.day.d, max.agr.d, pch=18, col="dark orange")
    # Add confidence intervals
    lines(rate.w$AGR.hi~rate.w$dap_i, lty=2, col=c('navy'))
    lines(rate.w$AGR.lo~rate.w$dap_i, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$AGR.lo, rev(rate.w$AGR.hi)), col=COL.CI, lty=2)
    lines(rate.d$AGR.hi~rate.d$dap_i, lty=2, col=c('brown'))
    lines(rate.d$AGR.lo~rate.d$dap_i, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$AGR.lo, rev(rate.d$AGR.hi)), col=COL.CI, lty=2)
    
    max.RGRt<-max(max(rate.w$RGRt), max(rate.d$RGRt))
    min.RGRt<-min(min(rate.w$RGRt), min(rate.d$RGRt))
    plot(rate.w$RGRt~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='RGRt', ylim=c(min.RGRt,max.RGRt))
    lines(rate.d$RGRt~rate.d$dap_i, col="orange")
    points(max.agr.day.w, rgrt.agr.max.w, pch=18, col="dark blue")
    points(max.agr.day.d, rgrt.agr.max.d, pch=18, col="dark orange")
    # Add confidence intervals
    lines(rate.w$RGRt.hi~rate.w$dap_i, lty=2, col=c('navy'))
    lines(rate.w$RGRt.lo~rate.w$dap_i, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$RGRt.lo, rev(rate.w$RGRt.hi)), col=COL.CI, lty=2)
    lines(rate.d$RGRt.hi~rate.d$dap_i, lty=2, col=c('brown'))
    lines(rate.d$RGRt.lo~rate.d$dap_i, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$RGRt.lo, rev(rate.d$RGRt.hi)), col=COL.CI, lty=2)
    
    max.RGRm<-max(max(rate.w$RGRm), max(rate.d$RGRm))
    min.RGRm<-min(min(rate.w$RGRm), min(rate.d$RGRm))
    plot(rate.w$RGRm~rate.w$M, type="l", col="blue", xlab='height_above_bound', ylab='RGRm', ylim=c(min.RGRm,max.RGRm))
    lines(rate.d$RGRm~rate.d$M, col="orange")
    points(max.agr.height_above_bound.w, rgrm.agr.max.w, pch=18, col="dark blue")
    points(max.agr.height_above_bound.d, rgrm.agr.max.d, pch=18, col="dark orange")
    lines(rate.w$RGRm.hi~rate.w$M, lty=2, col=c('navy'))
    lines(rate.w$RGRm.lo~rate.w$M, lty=2, col=c('navy'))
    #polygon(x=c(rate.w$M, rev(rate.w$M)), y=c(rate.w$RGRm.lo, rev(rate.w$RGRm.hi)), col=COL.CI, lty=2)
    lines(rate.d$RGRm.hi~rate.d$M, lty=2, col=c('brown'))
    lines(rate.d$RGRm.lo~rate.d$M, lty=2, col=c('brown'))
    #polygon(x=c(rate.d$M, rev(rate.d$M)), y=c(rate.d$RGRm.lo, rev(rate.d$RGRm.hi)), col=COL.CI, lty=2)
  }
  
  if (length(unique(rilrates$treamtment)) == 1) {
    if (nrow(rate.w) > 0) {
      max.agr.w<-max(rate.w$AGR)
      max.agr.day.w<-rate.w[rate.w$AGR == max.agr.w, 'dap_i']
      plot(rate.w$AGR~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='AGR', ylim=c(0,max.agr.w))
      points(max.agr.day.w, max.agr.w, pch=18, col="dark blue")
      lines(rate.w$AGR.hi~rate.w$dap_i, lty=2, col=c('navy'))
      lines(rate.w$AGR.lo~rate.w$dap_i, lty=2, col=c('navy'))
      #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$AGR.lo, rev(rate.w$AGR.hi)), col=COL.CI, lty=2)
      
      max.RGRt<-max(rate.w$RGRt)
      min.RGRt<-min(rate.w$RGRt)
      plot(rate.w$RGRt~rate.w$dap_i, type="l", col="blue", xlab='Days after planting', ylab='RGRt', ylim=c(min.RGRt,max.RGRt))
      points(max.agr.day.w, rgrt.agr.max.w, pch=18, col="dark blue")
      lines(rate.w$RGRt.hi~rate.w$dap_i, lty=2, col=c('navy'))
      lines(rate.w$RGRt.lo~rate.w$dap_i, lty=2, col=c('navy'))
      #polygon(x=c(rate.w$dap_i, rev(rate.w$dap_i)), y=c(rate.w$RGRt.lo, rev(rate.w$RGRt.hi)), col=COL.CI, lty=2)
      
      max.RGRm<-max(rate.w$RGRm)
      min.RGRm<-min(rate.w$RGRm)
      plot(rate.w$RGRm~rate.w$M, type="l", col="blue", xlab='height_above_bound', ylab='RGRm', ylim=c(min.RGRm,max.RGRm))
      points(max.agr.height_above_bound.w, rgrm.agr.max.w, pch=18, col="dark blue")
      lines(rate.w$RGRm.hi~rate.w$M, lty=2, col=c('navy'))
      lines(rate.w$RGRm.lo~rate.w$M, lty=2, col=c('navy'))
      #polygon(x=c(rate.w$M, rev(rate.w$M)), y=c(rate.w$RGRm.lo, rev(rate.w$RGRm.hi)), col=COL.CI, lty=2)
      
    }
    if (nrow(rate.d) > 0) {
      max.agr.d<-max(rate.d$AGR)
      max.agr.day.d<-rate.d[rate.d$AGR == max.agr.d, 'dap_i']
      plot(rate.d$AGR~rate.d$dap_i, type="l", col="orange", xlab='Days after planting', ylab='AGR', ylim=c(0,max.agr.d))
      points(max.agr.day.d, max.agr.d, pch=18, col="dark orange")
      lines(rate.d$AGR.hi~rate.d$dap_i, lty=2, col=c('brown'))
      lines(rate.d$AGR.lo~rate.d$dap_i, lty=2, col=c('brown'))
      #polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$AGR.lo, rev(rate.d$AGR.hi)), col=COL.CI, lty=2)
      
      max.RGRt<-max(rate.d$RGRt)
      min.RGRt<-min(rate.d$RGRt)
      plot(rate.d$RGRt~rate.d$dap_i, type="l", col="orange", xlab='Days after planting', ylab='RGRt', ylim=c(min.RGRt,max.RGRt))
      points(max.agr.day.d, rgrt.agr.max.d, pch=18, col="dark orange")
      lines(rate.d$RGRt.hi~rate.d$dap_i, lty=2, col=c('brown'))
      lines(rate.d$RGRt.lo~rate.d$dap_i, lty=2, col=c('brown'))
     # polygon(x=c(rate.d$dap_i, rev(rate.d$dap_i)), y=c(rate.d$RGRt.lo, rev(rate.d$RGRt.hi)), col=COL.CI, lty=2)
      
      max.RGRm<-max(rate.d$RGRm)
      min.RGRm<-min(rate.d$RGRm)
      plot(rate.d$RGRm~rate.d$M, type="l", col="orange", xlab='height_above_bound', ylab='RGRm', ylim=c(min.RGRm,max.RGRm))
      points(max.agr.height_above_bound.d, rgrm.agr.max.d, pch=18, col="dark orange") 
      lines(rate.d$RGRm.hi~rate.d$M, lty=2, col=c('brown'))
      lines(rate.d$RGRm.lo~rate.d$M, lty=2, col=c('brown'))
     # polygon(x=c(rate.d$M, rev(rate.d$M)), y=c(rate.d$RGRm.lo, rev(rate.d$RGRm.hi)), col=COL.CI, lty=2)
    }
  }
}
dev.off()


##### Lets start with height_above_bound

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/predicted_height/bellweather_ril/height_above_bound")
height.l<-read.csv('ril_best.fit_logistic_estimates_height.csv')
height.l<-height.l[,c(1:4)]

days<-sort(unique(height.l$dap_i))

ril_height_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  data<-height.l[height.l$dap_i == day, ]
  colnames(data)[4]<-paste('height', day, sep="_")
  data<-data[,c(1,2,4)]
  if (d == 1) {
    ril_height_qtl<-data
  }
  
  if (d > 1) {
    ril_height_qtl<-merge(ril_height_qtl, data, by=c('ril', 'treatment'), all=T)
  }
}

# Add misc columns
ril_height_qtl$Obs<-c(1:nrow(ril_height_qtl))
ril_height_qtl$experiment<-rep("BP14", nrow(ril_height_qtl))
ril_height_qtl$year<-rep("2014", nrow(ril_height_qtl))
ril_height_qtl$plot<-rep("bellweater", nrow(ril_height_qtl))
ril_height_qtl$plot_id<-rep("bellweater", nrow(ril_height_qtl))
ril_height_qtl$sampling<-rep("bellweater", nrow(ril_height_qtl))

ril_height_qtl<-ril_height_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_height_qtl)[7]<-c("id")

write.csv(ril_height_qtl, file="ril_height_above_bound_qtl.csv", quote=F, row.names=F)




#######################################################################################
# Defined functions
#######################################################################################

output.mono.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
  coef <- coef(fit)
  params <- transform_param.mono(coef)
  K = params[1]; r = params[2]; M0 = params[3]
  fitted <- fit$fitted
  resid  <- fit$residual
  data <- data.frame(fitted = fitted, resid = resid)
  eq   <- bquote(paste(.(round(r, 4)) - .(round(K - M0, 2)) * e^{.(round(r, 4))*t}))
  mss <- sum((fitted - mean(fitted))^2)
  rss <- sum(resid^2)
  R2  <- mss/(mss + rss)
  rmse <- sqrt(rss)
  AIC <- AIC(fit)
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  rates = data.frame(
    times = times,
    M     = K-(exp(-r*times)*(K-M0)),
    AGR   = r*exp(-r*times)*(K-M0)) #dM/dt
  rates$RGRt <- rates$AGR/rates$M
  rates$RGRm <- (r*(K-rates$M))/rates$M
  if(LOG ==T){
    rates$RGRt <- rates$AGR
    rates$RGRm <- r*(K-rates$M)
    rates$AGR  <- rates$AGR*exp(rates$M)
  }
  if(CI == T){
    cov <- fit$varBeta
    y <- x <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
    x$K  <- y$Asym
    x$r  <- exp(y$lrc)
    x$M0 <- y$R0 #untransform best-fit parameters to K, r and M0
    M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
    for(i in 1:nrow(x)){
      x.i  <- x[i,]
      M[i,]    <-  x.i$K-exp(-x.i$r*times)*(x.i$K-x.i$M0)
      AGR[i,]  <-  x.i$r*exp(-x.i$r*times)*(x.i$K-x.i$M0)
      RGRt[i,] <- AGR[i,]/M[i,]
      RGRm[i,] <- x.i$r*(x.i$K - M[i,])/M[i,]
      if(LOG ==T){
        RGRt[i,] <- AGR[i,]
        RGRm[i,] <- (x.i$r*(x.i$K-M[i,]))
        AGR[i,]  <- AGR[i,]*exp(M[i,])
      }
    }
    CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
  } else {
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
  }
  return(out)
}


output.logis.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
  coef <- coef(fit)
  params <- transform_param.logis(coef)
  K = params[1]; r = params[2]; M0 = params[3]
  fitted <- fit$fitted
  resid  <- fit$residuals
  data <- data.frame(fitted = fitted, resid = resid)
  eq   <- bquote(paste(.(round(M0*r, 4)) /(.(round(M0, 3)) + .(round(K-M0, 2)) * e^{.(round(-r, 3))*t})))
  mss <- sum((fitted - mean(fitted))^2)
  rss <- sum(resid^2)
  R2  <- mss/(mss + rss)
  rmse <- sqrt(rss)
  AIC <- AIC(fit)
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  rates = data.frame(
    times = times,
    M    = (M0*K)/(M0+(K-M0)*exp(-r*times)),
    AGR  = (r*M0*K*(K-M0)*exp(-r*times))/(M0+(K-M0)*exp(-r*times))^2)
  rates$RGRt <- rates$AGR/rates$M
  rates$RGRm <- r*(1 - rates$M/K)
  if(LOG ==T){
    rates$RGRt <- rates$AGR
    rates$RGRm <- r*rates$M*(1-rates$M/K)
    rates$AGR  <- rates$AGR*exp(rates$M)
  }
  if(CI == T){
    cov   <- fit$varBeta
    x <- y <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
    x$K  <- y$Asym
    x$r  <- 1/y$xmid
    x$M0 <- y$Asym/(1 + exp(y$xmid/y$scal)) #untransform best-fit parameters to K, r and M0
    M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
    for(i in 1:nrow(x)){
      x.i  <- x[i,]
      M[i,]     <- (x.i$M0*x.i$K)/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))
      AGR[i,]   <- (x.i$r*x.i$M0*x.i$K*(x.i$K-x.i$M0)*exp(-x.i$r*times))/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))^2
      RGRt[i,] <- AGR[i,]/M[i,]
      RGRm[i,]  <-  x.i$r*(1 - M[i,]/x.i$K)
      if(LOG ==T){
        RGRt[i,] <- AGR[i,]
        RGRm[i,] <- x.i$r*M[i,]*(1 - M[i,]/x.i$K)
        AGR[i,]  <- AGR[i,]*exp(M[i,])
      }
    }
    CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
  } else {
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
  }
  return(out)
}

output.fpl.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
  params <- coef <- coef(fit)
  M0 = params[1]; K = params[2]; xmid = params[3]; r = params[4]
  fitted <- fit$fitted
  resid  <- fit$residuals
  data <- data.frame(fitted = fitted, resid = resid)
  eq   <- bquote(paste(.(round(M0, 4))+ (.(round(K-M0, 4)) /1+e^{(.(round(xmid, 3))*t)/.(round(r, 3))})))
  mss <- sum((fitted - mean(fitted))^2)
  rss <- sum(resid^2)
  R2  <- mss/(mss + rss)
  rmse <- sqrt(rss)
  AIC <- AIC(fit)
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  exp <- exp((xmid-times)/r)
  rates = data.frame(
    times = times,
    M    = M0+(K-M0)/(1+exp),
    AGR  = ((K-M0)*exp)/(r*(1+exp)^2)
  )
  rates$RGRt <- rates$AGR/rates$M
  rates$RGRm <- (M0-rates$M)*(K-rates$M)/(r*rates$M*(M0-K))
  if(LOG == T){
    rates$RGRt <- rates$AGR
    rates$RGRm <- (M0-rates$M)*(K-rates$M)/(r*(M0-K))
    rates$AGR  <- rates$AGR*exp(rates$M)
  }
  if(CI == T){
    cov   <- fit$varBeta
    x <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
    names(x) <- c("M0", "K", "xmid", "r")
    M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
    for(i in 1:nrow(x)){
      x.i  <- x[i,]
      exp <- exp((x.i$xmid-times)/x.i$r)
      M[i,]     <- x.i$M0 + (x.i$K-x.i$M0)/(1+exp)
      AGR[i,]   <- ((x.i$K-x.i$M0)*exp)/(x.i$r*(1+exp)^2)
      RGRt[i,] <- AGR[i,]/M[i,]
      RGRm[i,]  <- x.i$r*(1 - (M[i,]/x.i$K)^(1/x.i$xmid))
      #    	RGRm[i,]  <- ((x.i$M0-M[i,])*(x.i$K-M[i,]))/(M[i,]*x.i$r*(x.i$M0-x.i$K))
      if(LOG == T){
        RGRt[i,] <- AGR[i,]
        RGRm[i,] <- (x.i$M0-M[i,])*(x.i$K-M[i,])/(x.i$r*(x.i$M0-x.i$K))
        AGR[i,]  <- AGR[i,]*exp(M[i,]) 
      }
    }
    CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
  } else {
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
  }
  return(out)
}

output.gomp.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
  coef <- coef(fit)
  params <- transform_param.gomp(coef)
  K = params[1]; r = params[2]; M0 = params[3]
  fitted <- fit$fitted
  resid  <- fit$residual
  data <- data.frame(fitted = fitted, resid = resid)
  eq   <- bquote(paste(.(round(K, 2)) %*% .(round(M0 / K, 5)) ^ e^{.(round(r, 3))*t}))
  mss <- sum((fitted - mean(fitted))^2)
  rss <- sum(resid^2)
  R2  <- mss/(mss + rss)
  rmse <- sqrt(rss)
  AIC <- AIC(fit)
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  rates = data.frame(
    times = times,
    M    = K*((M0/K)^exp(-r*times)),
    AGR  = r*K*exp(-r*times)*log(K/M0)*(M0/K)^exp(-r*times))
  rates$RGRt <- rates$AGR/rates$M
  rates$RGRm <- r*log(K/rates$M)
  if(LOG == T){
    rates$RGRt <- rates$AGR
    rates$RGRm <- r*rates$M*log(K/rates$M)
    rates$AGR  <- rates$AGR*exp(rates$M)
  }
  if(CI == T){
    cov  <- fit$varBeta
    x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
    x$K  <- x$Asym
    x$M0 <- x$K/exp(x$b2)
    x$r  <- -log(x$b3)
    M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
    for(i in 1:nrow(x)){
      x.i  <- x[i,]
      M[i,]    <- x.i$K*((x.i$M0/x.i$K)^exp(-x.i$r*times))
      AGR[i,]  <- x.i$r*x.i$K*exp(-x.i$r*times)*log(x.i$K/x.i$M0)*(x.i$M0/x.i$K)^exp(-x.i$r*times)
      RGRt[i,] <- AGR[i,]/M[i,]
      RGRm[i,] <- x.i$r*log(x.i$K/M[i,])
      if(LOG ==T){
        RGRt[i,] <- AGR[i,]
        RGRm[i,] <- x.i$r*M[i,]*log(x.i$K/M[i,])
        AGR[i,]  <- AGR[i,]*exp(M[i,])
      }
    }
    CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
  } else {
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
  }
  return(out)
}


transform_param.logis <- function(coef){
  K = coef[1]
  r = 1/(coef[3])
  M0 =  K/(1 + exp(coef[2]/coef[3])) #untransform best-fit parameters to K, r and M0
  if(is.data.frame(K)){
    out <- cbind(K, r, M0)
  } else {
    out <- c(K, r, M0)
  }
  names(out) <- c("K", "r", "M0")
  return(out)
}

summarizer <- function(dat, alpha){
  n <- length(dat)
  quantiles <- c(alpha/2, 1-(alpha/2))
  CIs <- data.frame(matrix(NA, ncol(dat[[1]]), n*2))
  names(CIs) <- paste(rep(names(dat), each = 2), c("lo", "hi"), sep = ".")
  for(i in 1:n){
    CIs[,(2*i-1):(2*i)] <- t(apply(dat[[i]],    2, quantile, quantiles, na.rm = T))
  }
  return(CIs)
}


# The following three functions transform the parameters of the logistic, gompertz and monomolecular models (respectively), to put them onto a scale most easily comprable with the other models.
transform_param.logis <- function(coef){
  K = coef[1]
  r = 1/(coef[3])
  M0 =  K/(1 + exp(coef[2]/coef[3])) #untransform best-fit parameters to K, r and M0
  if(is.data.frame(K)){
    out <- cbind(K, r, M0)
  } else {
    out <- c(K, r, M0)
  }
  names(out) <- c("K", "r", "M0")
  return(out)
}

transform_param.gomp <- function(coef){
  K  <- coef[1]
  M0 <- K/exp(coef[2])
  r  <- -log(coef[3])
  out <- c(K, r, M0)
  names(out) <- c("K", "r", "M0")
  return(out)
}


transform_param.mono <- function(coef){
  K  <- coef[1]
  M0 <- coef[2] 
  r  <- exp(coef[3])
  out <- c(K, r, M0)
  names(out) <- c("K", "r", "M0")
  return(out)
}

setwd("/Users/mfeldman/Desktop")
pdf("BP14_A10_logistic_Replicates.pdf")
plot(set$height~set$dap_i, type='l', xlim=c(8,33), ylim=c(0, max.b), col=l.color, xlab="Days after planting", ylab="Height (cm)")


