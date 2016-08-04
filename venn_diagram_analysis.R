
################### 
library(ggplot2)
library(stringr)
library(VennDiagram)


################### 
library(ggplot2)
library(stringr)

##############################################################################################################################
# Define fxns
##############################################################################################################################

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


################### Lets make some fxns
remove_dup_qtl<-function(temp){
  all_qtl<-sort(table(temp$marker), decreasing=T)
  if (length(all_qtl) == 1) {
    treatments<-as.character(unique(temp$treatment))
    if(length(treatments) == 1) { 
      m.names<<-c(m.names, names(all_qtl)[1])
      # <<- means change the global variable (chr<<-max) changes the global variable chr to local variable max
      chr<<-c(chr,unique(temp[temp$marker == names(all_qtl)[1],'chr']))
      pos<<-c(pos,unique(temp[temp$marker == names(all_qtl)[1],'pos']))
      t<-as.character(unique(temp$treatment))
      condition<<-c(condition, t)
      qtl_count<<-c(qtl_count, 1)
      
      max.lod<<-c(max.lod, max(temp[temp$marker == names(all_qtl)[1],'lod'], na.rm = T))
      max.prop.var<<-c(max.prop.var, max(temp[temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
      max.fx<<-c(max.fx, max(temp[temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
      max.fx_se<<-c(max.fx_se, max(temp[temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
      max.L.CI_pos<<-c(max.L.CI_pos, max(temp[temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
      max.R.CI_pos<<-c(max.R.CI_pos, max(temp[temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
      
      med.lod<<-c(med.lod, median(temp[temp$marker == names(all_qtl)[1],'lod'], na.rm = T))
      med.prop.var<<-c(med.prop.var, median(temp[temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
      med.fx<<-c(med.fx, median(temp[temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
      med.fx_se<<-c(med.fx_se, median(temp[temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
      med.L.CI_pos<<-c(med.L.CI_pos, median(temp[temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
      med.R.CI_pos<<-c(med.R.CI_pos, median(temp[temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
      
      min.lod<<-c(min.lod, min(temp[temp$marker == names(all_qtl)[1],'lod'], na.rm = T))
      min.prop.var<<-c(min.prop.var, min(temp[temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
      min.fx<<-c(min.fx, min(temp[temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
      min.fx_se<<-c(min.fx_se, min(temp[temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
      min.L.CI_pos<<-c(min.L.CI_pos, min(temp[temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
      min.R.CI_pos<<-c(min.R.CI_pos, min(temp[temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
      
      print(chr) 
      print(pos)
      print(qtl_count)
      print(med.lod)
    } 
    if(length(treatments) > 1){
      for(t in treatments) {
        t.temp<-temp[temp$treatment == t,]
        m.names<<-c(m.names, names(all_qtl)[1])
        # <<- means change the global variable (chr<<-max) changes the global variable chr to local variable max
        chr<<-c(chr,unique(t.temp[t.temp$marker == names(all_qtl)[1],'chr']))
        pos<<-c(pos,unique(t.temp[t.temp$marker == names(all_qtl)[1],'pos']))
        condition<<-c(condition, t)
        qtl_count<<-c(qtl_count, 1)
        
        max.lod<<-c(max.lod, max(t.temp[t.temp$marker == names(all_qtl)[1],'lod'], na.rm=T))
        max.prop.var<<-c(max.prop.var, max(t.temp[t.temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
        max.fx<<-c(max.fx, max(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
        max.fx_se<<-c(max.fx_se, max(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
        max.L.CI_pos<<-c(max.L.CI_pos, max(t.temp[t.temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
        max.R.CI_pos<<-c(max.R.CI_pos, max(t.temp[t.temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
        
        med.lod<<-c(med.lod, median(t.temp[t.temp$marker == names(all_qtl)[1],'lod'], na.rm=T))
        med.prop.var<<-c(med.prop.var, median(t.temp[t.temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
        med.fx<<-c(med.fx, median(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
        med.fx_se<<-c(med.fx_se, median(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
        med.L.CI_pos<<-c(med.L.CI_pos, median(t.temp[t.temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
        med.R.CI_pos<<-c(med.R.CI_pos, median(t.temp[t.temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
        
        min.lod<<-c(min.lod, min(t.temp[t.temp$marker == names(all_qtl)[1],'lod'], na.rm=T))
        min.prop.var<<-c(min.prop.var, min(t.temp[t.temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
        min.fx<<-c(min.fx, min(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
        min.fx_se<<-c(min.fx_se, min(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
        min.L.CI_pos<<-c(min.L.CI_pos, min(t.temp[t.temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
        min.R.CI_pos<<-c(min.R.CI_pos, min(t.temp[t.temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
        
        print(chr) 
        print(pos)
        print(qtl_count)
        print(med.lod)
      }
    }
  }
  if (length(all_qtl) > 1) {
    
    name<-names(all_qtl)[1]
    tester<-temp[temp$marker == name,]
    treatments<-as.character(unique(tester$treatment))
    
    if(length(treatments) == 1) { 
      ave.pos<-mean(temp[temp$marker == name, 'pos'])
      m.name<-names(all_qtl)[1]
      cr<-unique(temp[temp$marker == names(all_qtl)[1],'chr'])
      po<-unique(temp[temp$marker == names(all_qtl)[1],'pos'])
      max.pos<-ave.pos+10
      min.pos<-ave.pos-10
      subset<-temp[temp$pos > min.pos & temp$pos < max.pos,]
      if(length(unique(subset$treatment)) == 1){
        m.names<<-c(m.names, m.name)
        chr<<-c(chr, cr)
        pos<<-c(pos, po)
        qtl_c<-nrow(subset)
        
        x.lod<-max(subset$lod, na.rm = T)
        x.prop.var<-max(subset$prop.var, na.rm = T)
        x.fx<-max(subset$additive.fx, na.rm = T)
        x.fx_se<-max(subset$additive.fx_se, na.rm = T)
        x.L.CI_pos<-max(subset$L.CI_pos, na.rm = T)
        x.R.CI_pos<-max(subset$R.CI_pos, na.rm = T)
        
        m.lod<-median(subset$lod, na.rm = T)
        m.prop.var<-median(subset$prop.var, na.rm = T)
        m.fx<-median(subset$additive.fx, na.rm = T)
        m.fx_se<-median(subset$additive.fx_se, na.rm = T)
        m.L.CI_pos<-median(subset$L.CI_pos, na.rm = T)
        m.R.CI_pos<-median(subset$R.CI_pos, na.rm = T)
        
        n.lod<-min(subset$lod, na.rm = T)
        n.prop.var<-min(subset$prop.var, na.rm = T)
        n.fx<-min(subset$additive.fx, na.rm = T)
        n.fx_se<-min(subset$additive.fx_se, na.rm = T)
        n.L.CI_pos<-min(subset$L.CI_pos, na.rm = T)
        n.R.CI_pos<-min(subset$R.CI_pos, na.rm = T)
        
        temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
        print(ave.pos) 
        print(chr) 
        print(pos)
        #print(collapsed_qtl)
        print(med.lod)
        condition<<-c(condition, treatments[1])
        qtl_count<<-c(qtl_count, qtl_c)
        
        max.lod<<-c(max.lod, x.lod)
        max.prop.var<<-c(max.prop.var, x.prop.var)
        max.fx<<-c(max.fx, x.fx)
        max.fx_se<<-c(max.fx_se, x.fx_se)
        max.L.CI_pos<<-c(max.L.CI_pos, x.L.CI_pos)
        max.R.CI_pos<<-c(max.R.CI_pos, x.R.CI_pos)
        
        med.lod<<-c(med.lod, m.lod)
        med.prop.var<<-c(med.prop.var, m.prop.var)
        med.fx<<-c(med.fx, m.fx)
        med.fx_se<<-c(med.fx_se, m.fx_se)
        med.L.CI_pos<<-c(med.L.CI_pos, m.L.CI_pos)
        med.R.CI_pos<<-c(med.R.CI_pos, m.R.CI_pos)
        
        min.lod<<-c(min.lod, n.lod)
        min.prop.var<<-c(min.prop.var, n.prop.var)
        min.fx<<-c(min.fx, n.fx)
        min.fx_se<<-c(min.fx_se, n.fx_se)
        min.L.CI_pos<<-c(min.L.CI_pos, n.L.CI_pos)
        min.R.CI_pos<<-c(min.R.CI_pos, n.R.CI_pos)
        
        remove_dup_qtl(temp)
      }
      if(length(unique(subset$treatment)) > 1){
        subset.ts<-unique(subset$treatment)
        for (t in subset.ts) {
          m.names<<-c(m.names, m.name)
          chr<<-c(chr, cr)
          pos<<-c(pos, po)
          t.subset<-subset[subset$treatment == t,]
          qtl_c<-nrow(t.subset)
          
          x.lod<-max(t.subset$lod, na.rm = T)
          x.prop.var<-max(t.subset$prop.var, na.rm = T)
          x.fx<-max(t.subset$additive.fx, na.rm = T)
          x.fx_se<-max(t.subset$additive.fx_se, na.rm = T)
          x.L.CI_pos<-max(t.subset$L.CI_pos, na.rm = T)
          x.R.CI_pos<-max(t.subset$R.CI_pos, na.rm = T)
          
          m.lod<-median(t.subset$lod, na.rm = T)
          m.prop.var<-median(t.subset$prop.var, na.rm = T)
          m.fx<-median(t.subset$additive.fx, na.rm = T)
          m.fx_se<-median(t.subset$additive.fx_se, na.rm = T)
          m.L.CI_pos<-median(t.subset$L.CI_pos, na.rm = T)
          m.R.CI_pos<-median(t.subset$R.CI_pos, na.rm = T)
          
          n.lod<-min(t.subset$lod, na.rm = T)
          n.prop.var<-min(t.subset$prop.var, na.rm = T)
          n.fx<-min(t.subset$additive.fx, na.rm = T)
          n.fx_se<-min(t.subset$additive.fx_se, na.rm = T)
          n.L.CI_pos<-min(t.subset$L.CI_pos, na.rm = T)
          n.R.CI_pos<-min(t.subset$R.CI_pos, na.rm = T)
          
          print(ave.pos) 
          print(chr) 
          print(pos)
          #print(collapsed_qtl)
          print(med.lod)
          
          condition<<-c(condition, t)
          qtl_count<<-c(qtl_count, qtl_c)
          
          max.lod<<-c(max.lod, x.lod)
          max.prop.var<<-c(max.prop.var, x.prop.var)
          max.fx<<-c(max.fx, x.fx)
          max.fx_se<<-c(max.fx_se, x.fx_se)
          max.L.CI_pos<<-c(max.L.CI_pos, x.L.CI_pos)
          max.R.CI_pos<<-c(max.R.CI_pos, x.R.CI_pos)
          
          med.lod<<-c(med.lod, m.lod)
          med.prop.var<<-c(med.prop.var, m.prop.var)
          med.fx<<-c(med.fx, m.fx)
          med.fx_se<<-c(med.fx_se, m.fx_se)
          med.L.CI_pos<<-c(med.L.CI_pos, m.L.CI_pos)
          med.R.CI_pos<<-c(med.R.CI_pos, m.R.CI_pos)
          
          min.lod<<-c(min.lod, n.lod)
          min.prop.var<<-c(min.prop.var, n.prop.var)
          min.fx<<-c(min.fx, n.fx)
          min.fx_se<<-c(min.fx_se, n.fx_se)
          min.L.CI_pos<<-c(min.L.CI_pos, n.L.CI_pos)
          min.R.CI_pos<<-c(min.R.CI_pos, n.R.CI_pos)
          
        }
        temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
        remove_dup_qtl(temp)
      }
    }
    if(length(treatments) > 1) {
      #for (t in treatments) {
      ave.pos<-mean(temp[temp$marker == name, 'pos'])
      m.name<-names(all_qtl)[1]
      cr<-unique(temp[temp$marker == names(all_qtl)[1],'chr'])
      po<-unique(temp[temp$marker == names(all_qtl)[1],'pos'])
      max.pos<-ave.pos+10
      min.pos<-ave.pos-10
      subset<-temp[temp$pos > min.pos & temp$pos < max.pos,]
      subset.ts<-unique(subset$treatment)
      for (t in subset.ts) {
        m.names<<-c(m.names, m.name)
        chr<<-c(chr, cr)
        pos<<-c(pos, po)
        t.subset<-subset[subset$treatment == t,]
        qtl_c<-nrow(t.subset)
        
        x.lod<-max(t.subset$lod, na.rm = T)
        x.prop.var<-max(t.subset$prop.var, na.rm = T)
        x.fx<-max(t.subset$additive.fx, na.rm = T)
        x.fx_se<-max(t.subset$additive.fx_se, na.rm = T)
        x.L.CI_pos<-max(t.subset$L.CI_pos, na.rm = T)
        x.R.CI_pos<-max(t.subset$R.CI_pos, na.rm = T)
        
        m.lod<-median(t.subset$lod, na.rm = T)
        m.prop.var<-median(t.subset$prop.var, na.rm = T)
        m.fx<-median(t.subset$additive.fx, na.rm = T)
        m.fx_se<-median(t.subset$additive.fx_se, na.rm = T)
        m.L.CI_pos<-median(t.subset$L.CI_pos, na.rm = T)
        m.R.CI_pos<-median(t.subset$R.CI_pos, na.rm = T)
        
        n.lod<-min(t.subset$lod, na.rm = T)
        n.prop.var<-min(t.subset$prop.var, na.rm = T)
        n.fx<-min(t.subset$additive.fx, na.rm = T)
        n.fx_se<-min(t.subset$additive.fx_se, na.rm = T)
        n.L.CI_pos<-min(t.subset$L.CI_pos, na.rm = T)
        n.R.CI_pos<-min(t.subset$R.CI_pos, na.rm = T)
        
        print(ave.pos) 
        print(chr) 
        print(pos)
        #print(collapsed_qtl)
        print(med.lod)
        
        condition<<-c(condition, t)
        qtl_count<<-c(qtl_count, qtl_c)
        
        max.lod<<-c(max.lod, x.lod)
        max.prop.var<<-c(max.prop.var, x.prop.var)
        max.fx<<-c(max.fx, x.fx)
        max.fx_se<<-c(max.fx_se, x.fx_se)
        max.L.CI_pos<<-c(max.L.CI_pos, x.L.CI_pos)
        max.R.CI_pos<<-c(max.R.CI_pos, x.R.CI_pos)
        
        med.lod<<-c(med.lod, m.lod)
        med.prop.var<<-c(med.prop.var, m.prop.var)
        med.fx<<-c(med.fx, m.fx)
        med.fx_se<<-c(med.fx_se, m.fx_se)
        med.L.CI_pos<<-c(med.L.CI_pos, m.L.CI_pos)
        med.R.CI_pos<<-c(med.R.CI_pos, m.R.CI_pos)
        
        min.lod<<-c(min.lod, n.lod)
        min.prop.var<<-c(min.prop.var, n.prop.var)
        min.fx<<-c(min.fx, n.fx)
        min.fx_se<<-c(min.fx_se, n.fx_se)
        min.L.CI_pos<<-c(min.L.CI_pos, n.L.CI_pos)
        min.R.CI_pos<<-c(min.R.CI_pos, n.R.CI_pos)
        
      }
      temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
      remove_dup_qtl(temp)
    }
  }
}


condense_qtl<-function(input){
  
  chrs<-sort(unique(input$chr))
  m.names<<-c()
  chr<<-c()
  pos<<-c()
  condition<<-c()
  qtl_count<<-c()
  med.lod<<-c()
  med.prop.var<<-c()
  med.fx<<-c()
  med.fx_se<<-c()
  med.L.CI_pos<<-c()
  med.R.CI_pos<<-c()
  
  max.lod<<-c()
  max.prop.var<<-c()
  max.fx<<-c()
  max.fx_se<<-c()
  max.L.CI_pos<<-c()
  max.R.CI_pos<<-c()
  
  min.lod<<-c()
  min.prop.var<<-c()
  min.fx<<-c()
  min.fx_se<<-c()
  min.L.CI_pos<<-c()
  min.R.CI_pos<<-c()
  
  for(ch in chrs) {
    temp<-input[input$chr == ch,]
    temp$marker<-as.character(temp$marker)
    remove_dup_qtl(temp)
  }
  
  input.collapsed<-as.data.frame(cbind(m.names, chr, pos, condition, qtl_count, max.lod, max.prop.var, max.fx, max.fx_se, max.L.CI_pos, max.R.CI_pos,med.lod, med.prop.var, med.fx, med.fx_se, med.L.CI_pos, med.R.CI_pos,min.lod, min.prop.var, min.fx, min.fx_se, min.L.CI_pos, min.R.CI_pos))        
  return(input.collapsed)
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



##############################################################################################################################
# Load in data
##############################################################################################################################
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/qtl/summary_tables")

# Load in data
DR13<-read.csv("DR13_height_blup_concatenated_summary_table.csv", header=F)
DR14<-read.csv("DR14_height_blup_concatenated_summary_table.csv", header=F)
DL13<-read.csv("dinneney_traits_concatenated_summary_table.csv", header=F)
BP14<-read.csv("BP14_ril_height_above_bound_concatenated_summary_table.csv", header=F)
DN13<-read.csv("DN13_height_blup_concatenated_summary_table.csv", header=F)
DN14<-read.csv("DN14_height_blup_concatenated_summary_table.csv", header=F)

# Add column names
colnames(DR13)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
colnames(DR14)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
colnames(DL13)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
colnames(BP14)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
colnames(DN13)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
colnames(DN14)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

# Remove uninteresting phenotypes from DL13
DL13<-DL13[DL13$trait == 'plant_height_DL13', ]

all.st<-rbind(DR13, DR14, DL13, BP14, DN13, DN14)

all.st_uni<-condense_qtl(all.st)
all.st_uni<-unify_marker(all.st)

raw.st_uni<-all.st_uni[all.st_uni$type == 'raw',]
comp.st_uni<-all.st_uni[all.st_uni$type == 'comp',]

# Lets start finding unique QTL in each of the experimental contrasts...
raw.st_uni$exp<-paste(raw.st_uni$exp, raw.st_uni$year, sep="")

# Break out by experiment
DR13_uni<-raw.st_uni[raw.st_uni$exp == 'DR13', ]
DR14_uni<-raw.st_uni[raw.st_uni$exp == 'DR14', ]
DL13_uni<-raw.st_uni[raw.st_uni$exp == 'DL13', ]
BP14_uni<-raw.st_uni[raw.st_uni$exp == 'BP14', ]
DN13_uni<-raw.st_uni[raw.st_uni$exp == 'DN13', ]
DN14_uni<-raw.st_uni[raw.st_uni$exp == 'DN14', ]

# Break out by treatment
DR13_uni.w<-unique(DR13_uni[DR13_uni$treatment == 'wet','marker'])
DR13_uni.d<-unique(DR13_uni[DR13_uni$treatment == 'dry','marker'])
DR14_uni.w<-unique(DR14_uni[DR14_uni$treatment == 'wet','marker'])
DR14_uni.d<-unique(DR14_uni[DR14_uni$treatment == 'dry','marker'])

DL13_uni.w<-unique(DL13_uni[DL13_uni$treatment == 'wet','marker'])
DL13_uni.d<-unique(DL13_uni[DL13_uni$treatment == 'dry','marker'])

BP14_uni.w<-unique(BP14_uni[BP14_uni$treatment == 'wet','marker'])
BP14_uni.d<-unique(BP14_uni[BP14_uni$treatment == 'dry','marker'])

DN13_uni.s<-unique(DN13_uni[DN13_uni$treatment == 'sparse','marker'])
DN13_uni.d<-unique(DN13_uni[DN13_uni$treatment == 'dense','marker'])
DN14_uni.s<-unique(DN14_uni[DN14_uni$treatment == 'sparse','marker'])
DN14_uni.d<-unique(DN14_uni[DN14_uni$treatment == 'dense','marker'])

# Get all unique QTL
DR.all_nr<-unique(c(DR13_uni.w, DR13_uni.d, DR14_uni.w, DR14_uni.d, BP14_uni.w, BP14_uni.d,DL13_uni.w, DL13_uni.d))

# For all QTL in DR13 which are the ones that do not overlap between treatment groups (Unique QTL among treatment groups)
DR13_nr<-setdiff(union(DR13_uni.d, DR13_uni.w), intersect(DR13_uni.w, DR13_uni.d))
DR14_nr<-setdiff(union(DR14_uni.d, DR14_uni.w), intersect(DR14_uni.w, DR14_uni.d))

BP14_nr<-setdiff(union(BP14_uni.d, BP14_uni.w), intersect(BP14_uni.w, BP14_uni.d))
DL13_nr<-setdiff(union(DL13_uni.d, DL13_uni.w), intersect(DL13_uni.w, DL13_uni.d))

# Same for density expriments
DN.all_nr<-unique(c(DN13_uni.d, DN13_uni.s, DN14_uni.d, DN14_uni.s))

DN13_nr<-setdiff(union(DN13_uni.s, DN13_uni.d), intersect(DN13_uni.s, DN13_uni.d))
DN14_nr<-setdiff(union(DN14_uni.s, DN14_uni.d), intersect(DN14_uni.s, DN14_uni.d))

# Make plot of all non-redundant QTL
grid.newpage()
draw.quintuple.venn(area1=length(DR.all_nr), 
                    area2=length(DR13_nr), 
                    area3=length(DR14_nr),
                    area4=length(BP14_nr), 
                    area5=length(DL13_nr),
                    n12=length(intersect(unique(DR.all_nr), unique(DR13_nr))),
                    n13=length(intersect(unique(DR.all_nr), unique(DR14_nr))),
                    n14=length(intersect(unique(DR.all_nr), unique(BP14_nr))),
                    n15=length(intersect(unique(DR.all_nr), unique(DL13_nr))),
                    n23=length(intersect(unique(DR13_nr), unique(DR14_nr))),
                    n24=length(intersect(unique(DR13_nr), unique(BP14_nr))),
                    n25=length(intersect(unique(DR13_nr), unique(DL13_nr))),
                    n34=length(intersect(unique(DR14_nr), unique(BP14_nr))),
                    n35=length(intersect(unique(DR14_nr), unique(DL13_nr))),
                    n45=length(intersect(unique(BP14_nr), unique(DL13_nr))),
                    n123=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), unique(DR14_nr)))),
                    n124=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), unique(BP14_nr)))),
                    n125=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), unique(DL13_nr)))),
                    n134=length(intersect(unique(DR.all_nr), intersect(unique(DR14_nr), unique(BP14_nr)))),
                    n135=length(intersect(unique(DR.all_nr), intersect(unique(DR14_nr), unique(DL13_nr)))),
                    n145=length(intersect(unique(DR.all_nr), intersect(unique(BP14_nr), unique(DL13_nr)))),
                    n234=length(intersect(unique(DR13_nr), intersect(unique(DR14_nr), unique(BP14_nr)))),
                    n235=length(intersect(unique(DR13_nr), intersect(unique(DR14_nr), unique(DL13_nr)))),
                    n245=length(intersect(unique(DR13_nr), intersect(unique(BP14_nr), unique(DL13_nr)))),
                    n345=length(intersect(unique(DR14_nr), intersect(unique(BP14_nr), unique(DL13_nr)))),
                    n1234=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), intersect(unique(DR14_nr), unique(BP14_nr))))),
                    n1235=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), intersect(unique(DR14_nr), unique(DL13_nr))))),
                    n1245=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), intersect(unique(BP14_nr), unique(DL13_nr))))),
                    n1345=length(intersect(unique(DR.all_nr), intersect(unique(DR14_nr), intersect(unique(BP14_nr), unique(DL13_nr))))),
                    n2345=length(intersect(unique(DR13_nr), intersect(unique(DR14_nr), intersect(unique(BP14_nr), unique(DL13_nr))))),
                    n12345=length(intersect(unique(DR.all_nr), intersect(unique(DR13_nr), intersect(unique(DR14_nr), intersect(unique(BP14_nr), unique(DL13_nr)))))),
                    category=c("All", "DR13", "DR14", "BP14", "DL13"),
                    fill=c("green","light blue", "dark blue", "orange", "red")
)

dev.off()


# Same figure without all QTL
grid.newpage()
draw.quad.venn(area1=length(DR13_nr), 
               area2=length(DR14_nr),
               area3=length(BP14_nr), 
               area4=length(DL13_nr),
               n12=length(intersect(unique(DR13_nr), unique(DR14_nr))),
               n13=length(intersect(unique(DR13_nr), unique(BP14_nr))),
               n14=length(intersect(unique(DR13_nr), unique(DL13_nr))),
               n23=length(intersect(unique(DR14_nr), unique(BP14_nr))),
               n24=length(intersect(unique(DR14_nr), unique(DL13_nr))),
               n34=length(intersect(unique(BP14_nr), unique(DL13_nr))),
               n123=length(intersect(unique(DR13_nr), intersect(unique(DR14_nr), unique(BP14_nr)))),
               n124=length(intersect(unique(DR13_nr), intersect(unique(DR14_nr), unique(DL13_nr)))),
               n134=length(intersect(unique(DR13_nr), intersect(unique(BP14_nr), unique(DL13_nr)))),
               n234=length(intersect(unique(DR14_nr), intersect(unique(BP14_nr), unique(DL13_nr)))),
               n1234=length(intersect(unique(DR13_nr), intersect(unique(DR14_nr), intersect(unique(BP14_nr), unique(DL13_nr))))),
               category=c("DR13", "DR14", "BP14", "DL13"),
               fill=c("green","light blue", "dark blue", "orange")
)


dev.off()


# Get QTL specific for wet or dry treatments: wet only (wo) or dry only (do)
DR13_wo<-setdiff(DR13_uni.w, DR13_uni.d)
DR13_do<-setdiff(DR13_uni.d, DR13_uni.w)

DR14_wo<-setdiff(DR14_uni.w, DR14_uni.d)
DR14_do<-setdiff(DR14_uni.d, DR14_uni.w)

BP14_wo<-setdiff(BP14_uni.w, BP14_uni.d)
BP14_do<-setdiff(BP14_uni.d, BP14_uni.w)

DL13_wo<-setdiff(DL13_uni.w, DL13_uni.d)
DL13_do<-setdiff(DL13_uni.d, DL13_uni.w)


# Get QTL specific for dense or sparse treatments: dense only (do) or sparse only (so)

DN13_so<-setdiff(DN13_uni.s, DN13_uni.d)
DN13_do<-setdiff(DN13_uni.d, DN13_uni.s)

DN14_so<-setdiff(DN14_uni.s, DN14_uni.d)
DN14_do<-setdiff(DN14_uni.d, DN14_uni.s)


##################################################################
# Make plot of overlapping QTL found in field drought experiments
# Figure 5
##################################################################

grid.newpage()
pdf("Figure_5a.pdf")
draw.quad.venn(area1=length(DR13_uni.w), 
               area2=length(DR13_uni.d),
               area3=length(DR14_uni.w), 
               area4=length(DR14_uni.d),
               n12=length(intersect(unique(DR13_uni.w), unique(DR13_uni.d))),
               n13=length(intersect(unique(DR13_uni.w), unique(DR14_uni.w))),
               n14=length(intersect(unique(DR13_uni.w), unique(DR14_uni.d))),
               n23=length(intersect(unique(DR13_uni.d), unique(DR14_uni.w))),
               n24=length(intersect(unique(DR13_uni.d), unique(DR14_uni.d))),
               n34=length(intersect(unique(DR14_uni.w), unique(DR14_uni.d))),
               n123=length(intersect(unique(DR13_uni.w), intersect(unique(DR13_uni.d), unique(DR14_uni.w)))),
               n124=length(intersect(unique(DR13_uni.w), intersect(unique(DR13_uni.d), unique(DR14_uni.d)))),
               n134=length(intersect(unique(DR13_uni.w), intersect(unique(DR14_uni.w), unique(DR14_uni.d)))),
               n234=length(intersect(unique(DR13_uni.d), intersect(unique(DR14_uni.w), unique(DR14_uni.d)))),
               n1234=length(intersect(unique(DR13_uni.w), intersect(unique(DR13_uni.d), intersect(unique(DR14_uni.w), unique(DR14_uni.d))))),
               category=c("2013 Wet", "2013 Dry", "2014 Wet", "2014 Dry"),
               fill=c("light blue","yellow", "dark blue", "orange"),
               cex=2,
               cat.cex=1.5
               
)


# Make plot of overlapping QTL found in field density experiments
grid.newpage()
pdf("Figure_5b.pdf")
draw.quad.venn(area1=length(DN13_uni.s), 
               area2=length(DN13_uni.d),
               area3=length(DN14_uni.s), 
               area4=length(DN14_uni.d),
               n12=length(intersect(unique(DN13_uni.s), unique(DN13_uni.d))),
               n13=length(intersect(unique(DN13_uni.s), unique(DN14_uni.s))),
               n14=length(intersect(unique(DN13_uni.s), unique(DN14_uni.d))),
               n23=length(intersect(unique(DN13_uni.d), unique(DN14_uni.s))),
               n24=length(intersect(unique(DN13_uni.d), unique(DN14_uni.d))),
               n34=length(intersect(unique(DN14_uni.s), unique(DN14_uni.d))),
               n123=length(intersect(unique(DN13_uni.s), intersect(unique(DN13_uni.d), unique(DN14_uni.s)))),
               n124=length(intersect(unique(DN13_uni.s), intersect(unique(DN13_uni.d), unique(DN14_uni.d)))),
               n134=length(intersect(unique(DN13_uni.s), intersect(unique(DN14_uni.s), unique(DN14_uni.d)))),
               n234=length(intersect(unique(DN13_uni.d), intersect(unique(DN14_uni.s), unique(DN14_uni.d)))),
               n1234=length(intersect(unique(DN13_uni.s), intersect(unique(DN13_uni.d), intersect(unique(DN14_uni.s), unique(DN14_uni.d))))),
               category=c("2013 Sparse", "2013 Dense", "2014 Sparse", "2014 Dense"),
               fill=c("light green","grey80", "dark green", "grey30"),
               cex=2,
               cat.cex=1.5
               
)


dev.off()



dr.field<-raw.st_uni[grep('DR.*', raw.st_uni$exp), ]
dr.cont<-raw.st_uni[grepl('BP14.*', raw.st_uni$exp) | grepl('DL13.*', raw.st_uni$exp),]
all.field<-raw.st_uni[grepl('DR.*', raw.st_uni$exp) | grepl('DN.*', raw.st_uni$exp),]

drf.marker<-unique(dr.field$marker)
drc.marker<-unique(dr.cont$marker)
f.marker<-unique(all.field$marker)


#dr.field.v.cont<-setdiff(union(drf.marker, drc.marker), intersect(drf.marker, drc.marker))
#all.field.v.cont<-setdiff(union(f.marker, drc.marker), intersect(f.marker, drc.marker))


# Controlled vs. Field Drought Only
grid.newpage()
draw.pairwise.venn(area1=length(drc.marker),
                   area2=length(drf.marker),
                   cross.area=length(intersect(drf.marker, drc.marker)),
                   c("Controlled", "Field Drought Experiment Only"),
                   fill=c("blue", "green")
)

dev.off()

# Lets find out what the overlapping and significant QTL are...
dr.field.v.cont.diff<-setdiff(union(drf.marker, drc.marker), intersect(drf.marker, drc.marker))
dr.field.v.cont.same<-intersect(drf.marker, drc.marker) 
dr.cont.only<-setdiff(dr.field.v.cont.diff, drf.marker)
dr.field.only<-setdiff(dr.field.v.cont.diff, drc.marker)

# Controlled vs. Field Drought
grid.newpage()
pdf("Figure_5c.pdf")
draw.pairwise.venn(area1=length(drc.marker),
                   area2=length(f.marker),
                   cross.area=length(intersect(f.marker, drc.marker)),
                   c("Controlled", "Field"),
                   fill=c("blue", "green"),
                   cex=2,
                   cat.cex=1.5
                   
)

dev.off()

all.field.v.cont.diff<-setdiff(union(f.marker, drc.marker), intersect(f.marker, drc.marker))
all.field.v.cont.same<-intersect(f.marker, drc.marker) 
all.cont.only<-setdiff(all.field.v.cont.diff, f.marker)
all.field.only<-setdiff(all.field.v.cont.diff, drc.marker)


##################################################################
# Prepare Supplemental Figure 6 (Figure_S6)
##################################################################

# Make venn diagrams of overlap between QTL detected in treatment block and experimental locations 
# within controlled environmental studies

grid.newpage()
pdf("Figure_S6a.pdf")
draw.quad.venn(area1=length(BP14_uni.w), 
               area2=length(BP14_uni.d),
               area3=length(DL13_uni.w), 
               area4=length(DL13_uni.d),
               n12=length(intersect(unique(BP14_uni.w), unique(BP14_uni.d))),
               n13=length(intersect(unique(BP14_uni.w), unique(DL13_uni.w))),
               n14=length(intersect(unique(BP14_uni.w), unique(DL13_uni.d))),
               n23=length(intersect(unique(BP14_uni.d), unique(DL13_uni.w))),
               n24=length(intersect(unique(BP14_uni.d), unique(DL13_uni.d))),
               n34=length(intersect(unique(DL13_uni.w), unique(DL13_uni.d))),
               n123=length(intersect(unique(BP14_uni.w), intersect(unique(BP14_uni.d), unique(DL13_uni.w)))),
               n124=length(intersect(unique(BP14_uni.w), intersect(unique(BP14_uni.d), unique(DL13_uni.d)))),
               n134=length(intersect(unique(BP14_uni.w), intersect(unique(DL13_uni.w), unique(DL13_uni.d)))),
               n234=length(intersect(unique(BP14_uni.d), intersect(unique(DL13_uni.w), unique(DL13_uni.d)))),
               n1234=length(intersect(unique(BP14_uni.w), intersect(unique(BP14_uni.d), intersect(unique(DL13_uni.w), unique(DL13_uni.d))))),
               category=c("Bellweather Wet", "Bellweather Dry", "Carnegie Wet", "Carnegie Dry"),
               fill=c("dark blue","orange", "light blue", "yellow"),
               cex=2,
               cat.cex=1.5
)


dev.off()


# WET ONLY
grid.newpage()
pdf("Figure_S6b.pdf")
draw.quad.venn(area1=length(DR13_wo), 
               area2=length(DR14_wo),
               area3=length(BP14_wo), 
               area4=length(DL13_wo),
               n12=length(intersect(unique(DR13_wo), unique(DR14_wo))),
               n13=length(intersect(unique(DR13_wo), unique(BP14_wo))),
               n14=length(intersect(unique(DR13_wo), unique(DL13_wo))),
               n23=length(intersect(unique(DR14_wo), unique(BP14_wo))),
               n24=length(intersect(unique(DR14_wo), unique(DL13_wo))),
               n34=length(intersect(unique(BP14_wo), unique(DL13_wo))),
               n123=length(intersect(unique(DR13_wo), intersect(unique(DR14_wo), unique(BP14_wo)))),
               n124=length(intersect(unique(DR13_wo), intersect(unique(DR14_wo), unique(DL13_wo)))),
               n134=length(intersect(unique(DR13_wo), intersect(unique(BP14_wo), unique(DL13_wo)))),
               n234=length(intersect(unique(DR14_wo), intersect(unique(BP14_wo), unique(DL13_wo)))),
               n1234=length(intersect(unique(DR13_wo), intersect(unique(DR14_wo), intersect(unique(BP14_wo), unique(DL13_wo))))),
               category=c("DR13 Wet", "DR14 Wet", "BP14 Wet", "DL13 Wet"),
               fill=c("green","light blue", "dark blue", "orange")
)

dev.off()

# DRY ONLY
grid.newpage()
pdf("Figure_S6c.pdf")
draw.quad.venn(area1=length(DR13_do), 
               area2=length(DR14_do),
               area3=length(BP14_do), 
               area4=length(DL13_do),
               n12=length(intersect(unique(DR13_do), unique(DR14_do))),
               n13=length(intersect(unique(DR13_do), unique(BP14_do))),
               n14=length(intersect(unique(DR13_do), unique(DL13_do))),
               n23=length(intersect(unique(DR14_do), unique(BP14_do))),
               n24=length(intersect(unique(DR14_do), unique(DL13_do))),
               n34=length(intersect(unique(BP14_do), unique(DL13_do))),
               n123=length(intersect(unique(DR13_do), intersect(unique(DR14_do), unique(BP14_do)))),
               n124=length(intersect(unique(DR13_do), intersect(unique(DR14_do), unique(DL13_do)))),
               n134=length(intersect(unique(DR13_do), intersect(unique(BP14_do), unique(DL13_do)))),
               n234=length(intersect(unique(DR14_do), intersect(unique(BP14_do), unique(DL13_do)))),
               n1234=length(intersect(unique(DR13_do), intersect(unique(DR14_do), intersect(unique(BP14_do), unique(DL13_do))))),
               category=c("DR13 Dry", "DR14 Dry", "BP14 Dry", "DL13 Dry"),
               fill=c("green","light blue", "dark blue", "orange")
)


dev.off()

# Sparse
grid.newpage()
pdf("Figure_S6d.pdf")
draw.pairwise.venn(area1=length(DN13_so),
                   area2=length(DN14_so),
                   cross.area=length(intersect(DN13_so, DN14_so)),
                   c("2013 Sparse", "2014 Sparse"),
                   fill=c("purple", "orange")
                   
)

dev.off()

# Dense
grid.newpage()
pdf("Figure_S6e.pdf")
draw.pairwise.venn(area1=length(DN13_do),
                   area2=length(DN14_do),
                   cross.area=length(intersect(DN13_do, DN14_do)),
                   c("2013 Dense", "2014 Dense"),
                   fill=c("purple", "orange")
                   
)

dev.off()


