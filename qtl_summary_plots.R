# Script for plotting the results of QTL studies of plant height

# Load libraries
library(ggplot2)
library(VennDiagram)

##############################################################################################################################
# Define fxns
##############################################################################################################################

#####################
###### This is a function to get unique qtl from a qtl summary table
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


##############################################################################################################################
# Load in data
##############################################################################################################################

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/qtl/summary_tables")

# Illinois field experiments
DN13.qtl<-read.csv("DN13_height_blup_concatenated_summary_table.csv", header=F)
colnames(DN13.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
DN14.qtl<-read.csv("DN14_height_blup_concatenated_summary_table.csv", header=F)
colnames(DN14.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
DR13.qtl<-read.csv("DR13_height_blup_concatenated_summary_table.csv", header=F)
colnames(DR13.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
DR14.qtl<-read.csv("DR14_height_blup_concatenated_summary_table.csv", header=F)
colnames(DR14.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

# Carnegie
DL13.qtl<-read.csv("dinneney_traits_concatenated_summary_table.csv", header=F)
colnames(DL13.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
DL13.qtl<-DL13.qtl[grep('height', DL13.qtl$trait),]
DL13.qtl$marker<-as.character(DL13.qtl$marker)
# Bellweather
BP14.qtl<-read.csv("BP14_ril_height_above_bound_concatenated_summary_table.csv", header=F)
colnames(BP14.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

# Lets put the comparison (treatment difference traits) in a seperate data.frame
DN13.diff.qtl<-DN13.qtl[DN13.qtl$type == "comp",]
DN13.qtl<-DN13.qtl[DN13.qtl$type == "raw",]
DN14.diff.qtl<-DN14.qtl[DN14.qtl$type == "comp",]
DN14.qtl<-DN14.qtl[DN14.qtl$type == "raw",]

DR13.diff.qtl<-DR13.qtl[DR13.qtl$type == "comp",]
DR13.qtl<-DR13.qtl[DR13.qtl$type == "raw",]
DR14.diff.qtl<-DR14.qtl[DR14.qtl$type == "comp",]
DR14.qtl<-DR14.qtl[DR14.qtl$type == "raw",]

BP14.diff.qtl<-BP14.qtl[BP14.qtl$type == "comp",]
BP14.qtl<-BP14.qtl[BP14.qtl$type == "raw",]

DL13.diff.qtl<-DL13.qtl[DL13.qtl$type == "comp",]
DL13.qtl<-DL13.qtl[DL13.qtl$type == "raw",]


# Combine Density
DN.qtl<-rbind(DN13.qtl, DN14.qtl)

# Combine Drought in field
DR_field.qtl<-rbind(DR13.qtl, DR14.qtl)

# Combine all drought
DR_all.qtl<-rbind(DR13.qtl, DR14.qtl, DL13.qtl, BP14.qtl)
# Make a column in this data.frame to specify field v. controlled environment
DR_all.qtl$env<-rep(NA, nrow(DR_all.qtl))

# Combine all field
field_all.qtl<-rbind(DR13.qtl, DR14.qtl, DN13.qtl, DN14.qtl)

# Combine all controlled 
DR_cont.qtl<-rbind(DL13.qtl, BP14.qtl)

# Combine all qtl
all.qtl<-rbind(DL13.qtl, BP14.qtl, DR13.qtl, DR14.qtl, DN13.qtl, DN14.qtl)
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
write.csv(all.qtl, file="Table_S9.csv", quote=F, row.names=F)

# Combine all diff qtl
all.diff.qtl<-rbind(DL13.diff.qtl, BP14.diff.qtl, DR13.diff.qtl, DR14.diff.qtl, DN13.diff.qtl, DN14.diff.qtl)

##############################################################################################################################
# Load in genetic map for plotting
##############################################################################################################################

# Load in genetic map
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/results/genetic_map/")
map<-read.table("GBS_map_A10xB100_v0.96.csv", header=T, sep=",")
chrs<-t(map[1,2:ncol(map)])
pos<-t(map[2,2:ncol(map)])


genome<-cbind(chrs, pos)
colnames(genome)<-c('chrs', 'pos')
genome<-as.data.frame(genome)
genome$chrs<-as.numeric(as.character(genome$chrs))
genome$pos<-as.numeric(as.character(genome$pos))

c1.max<-max(genome[genome$chrs == 1, 'pos'])
c2.max<-max(genome[genome$chrs == 2, 'pos'])
c3.max<-max(genome[genome$chrs == 3, 'pos'])
c4.max<-max(genome[genome$chrs == 4, 'pos'])
c5.max<-max(genome[genome$chrs == 5, 'pos'])
c6.max<-max(genome[genome$chrs == 6, 'pos'])
c7.max<-max(genome[genome$chrs == 7, 'pos'])
c8.max<-max(genome[genome$chrs == 8, 'pos'])
c9.max<-max(genome[genome$chrs == 9, 'pos'])

blank_data<-data.frame(chr=c('1','1','2','2','3','3','4','4','5','5','6','6','7','7','8','8','9','9'), x=c(0,c1.max,0,c2.max,0,c3.max,0,c4.max,0,c5.max,0,c6.max,0,c7.max,0,c8.max,0,c9.max), y=0)
all.qtl$chr<-factor(all.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-all.qtl$additive.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

all.qtl$plot.char<-plot.char
all.qtl$plot.char<-as.factor(all.qtl$plot.char)

all.qtl$group<-paste(all.qtl$exp, all.qtl$year, all.qtl$treatment, sep="_")

treatments<-as.character(all.qtl$treatment)
treatment.name<-unique(treatments)
plot.col<-c()
for(i in 1:length(treatments)){
  logical<-treatments[i] == treatment.name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

all.qtl$plot.col<-plot.col

##############################################################################################################################
# Make the plots of all QTL identified in different types of contrasts
##############################################################################################################################

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")
# This is Figure_3
pdf("Figure_3.pdf")
p<-ggplot() + geom_point(data = all.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + scale_fill_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))
dev.off()

# Lets just look at density
all_density<-all.qtl[all.qtl$exp == "DN",]

p<-ggplot() + geom_point(data = all_density, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + scale_fill_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

# Now just at all field studies
all_field<-all.qtl[all.qtl$exp != "BP" & all.qtl$exp !="DL",]

p<-ggplot() + geom_point(data = all_field, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + scale_fill_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

# Now just at all controlled environmental studies
just_controlled<-all.qtl[all.qtl$exp == "BP" | all.qtl$exp == "DL",]

p<-ggplot() + geom_point(data = just_controlled, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + scale_fill_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

##############################################################################################################################
# Collapse these summary tables into 10 cM radius
#############################################################################################################################

###########################################################################
# BP14
###########################################################################
BP14.collapsed<-condense_qtl(BP14.qtl)

###########################################################################
# DN13
###########################################################################
DN13.collapsed<-condense_qtl(DN13.qtl)

###########################################################################
# DN14
###########################################################################
DN14.collapsed<-condense_qtl(DN14.qtl)

###########################################################################
# DR13
###########################################################################
DR13.collapsed<-condense_qtl(DR13.qtl)

###########################################################################
# DR14
###########################################################################
DR14.collapsed<-condense_qtl(DR14.qtl)

###########################################################################
# DR_field
###########################################################################
DR_field.collapsed<-condense_qtl(DR_field.qtl)

###########################################################################
# All_qtl
###########################################################################
all.collapsed<-condense_qtl(all.qtl)

###########################################################################
# all.diff.qtl
###########################################################################
all.diff.collapsed<-condense_qtl(all.diff.qtl)

all.diff.qtl_exp<-all.diff.qtl
all.diff.qtl_exp$treatment<-paste(all.diff.qtl$treatment, paste(all.diff.qtl$exp, all.diff.qtl$year, sep=""), sep="_")
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
write.csv(all.diff.qtl_exp, file="Table_S10.csv", quote=F, row.names=F)

# Combine all QTL and diff QTL
all.qtl_w.diff_exp<-rbind(all.qtl[,c(1:20)], all.diff.qtl_exp)
all.qtl_w.diff_exp.collapsed<-condense_qtl(all.qtl_w.diff_exp)
#setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
#write.csv(all.qtl_w.diff_exp.collapsed, file="Table_S10.csv", quote=F, row.names=F)

###########################################################################
# DR_cont
###########################################################################
DR_cont.collapsed<-condense_qtl(DR_cont.qtl)

###########################################################################
# DN_cont
###########################################################################
DN_cont.collapsed<-condense_qtl(DN_cont.qtl)

###########################################################################
# Write out condensed QTL to table (Table 2)
###########################################################################
setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/main_figure_table")
table_2<-all.collapsed
colnames(table_2)<-c("Marker","Chromosome", "Position", "Treatment", "Number of times observed", 
                     "Maximum LOD score", "Maximum proportion of variance", "Maximum effect size", "Standard error of maximum effect size", "Maximum confidience interval left boarder", "Maximum confidience interval right boarder",
                     "Median LOD score", "Median proportion of variance", "Median effect size", "Standard error of median effect size", "Median confidience interval left boarder", "Nedian confidience interval right boarder",
                     "Minimum LOD score", "Minimum proportion of variance", "Minimum effect size", "Standard error of minimum effect size", "Minimum confidience interval left boarder", "Minimum confidience interval right boarder")

write.csv(table_2, file="Table_2.csv", quote=F, row.names=F)

###########################################################################
# Lets make some plots of the collapsed QTL
###########################################################################
all.collapsed.qtl<-all.collapsed

all.collapsed.qtl$chr<-factor(all.collapsed.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-all.collapsed.qtl$med.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

all.collapsed.qtl$plot.char<-plot.char
all.collapsed.qtl$plot.char<-as.factor(all.collapsed.qtl$plot.char)

#all.collapsed.qtl$group<-paste(all.collapsed.qtl$exp, all.collapsed.qtl$year, all.collapsed.qtl$treatment, sep="_")

treatments<-all.collapsed.qtl$condition
treatment.name<-unique(treatments)
plot.col<-c()
for(i in 1:length(treatments)){
  logical<-treatments[i] == treatment.name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

all.collapsed.qtl$plot.col<-plot.col

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
pdf("Figure_S5.pdf")
p<-ggplot() + geom_point(data = all.collapsed.qtl, aes(x = as.numeric(as.character(pos)), y = as.numeric(as.character(med.prop.var)), shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + scale_fill_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))
dev.off()


###########################################################################
# Write out QTL found in each experimental treatment block (for static QTL analysis)
###########################################################################
# Get a list of QTL/markers found in each experimental treatment block
marker_names<-names(table(all.collapsed$m.names)[table(all.collapsed$m.names) >3])
# Subset the data frame to include only these markers
all.env.qtl<-all.collapsed[all.collapsed$m.names %in% marker_names,]

# Convert to numeric
for(col in c(2,3,5:23)){
  all.env.qtl[,col]<-as.numeric(as.character(all.env.qtl[,col]))
}

# Get a summary and print to file for more detailed analysis
all.env.qtl.ag<-aggregate(all.env.qtl[,-c(1,4)], by=list(all.env.qtl$m.names), mean)
colnames(all.env.qtl.ag)[1]<-c("marker")
all.env.qtl.ag$med.prop.var

known_marker<-all.env.qtl.ag[,c(1,2,3,11)]
write.csv(known_marker, "markers_for_fixed_model.csv", quote=F, row.names=F)

###########################################################################
# Lets now plot the QTLs from the difference as well.
###########################################################################

all.qtl_w.diff<-rbind(all.qtl[,c(1:20)], all.diff.qtl)

all.qtl_w.diff.collapsed<-condense_qtl(all.qtl_w.diff)
all.qtl_w.diff.collapsed.qtl<-all.qtl_w.diff.collapsed

all.qtl_w.diff.collapsed.qtl$chr<-factor(all.qtl_w.diff.collapsed.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-all.qtl_w.diff.collapsed.qtl$med.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

all.qtl_w.diff.collapsed.qtl$plot.char<-plot.char
all.qtl_w.diff.collapsed.qtl[all.qtl_w.diff.collapsed.qtl$condition == 'diff', 'plot.char']<-c('20')
all.qtl_w.diff.collapsed.qtl$plot.char<-as.factor(all.qtl_w.diff.collapsed.qtl$plot.char)

#all.qtl_w.diff.collapsed.qtl$group<-paste(all.qtl_w.diff.collapsed.qtl$exp, all.qtl_w.diff.collapsed.qtl$year, all.qtl_w.diff.collapsed.qtl$treatment, sep="_")

treatments<-all.qtl_w.diff.collapsed.qtl$condition
treatment.name<-unique(treatments)
plot.col<-c()
for(i in 1:length(treatments)){
  logical<-treatments[i] == treatment.name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

all.qtl_w.diff.collapsed.qtl$plot.col<-plot.col

setwd("/Users/mfeldman/Dropbox/setaria_height_paper/manuscript/supplemental_figure_table")
pdf("Figure_S7.pdf")
p<-ggplot() + geom_point(data = all.qtl_w.diff.collapsed.qtl, aes(x = as.numeric(as.character(pos)), y = as.numeric(as.character(med.prop.var)), shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(20,24,25)) 
print(p + scale_color_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue", "5" = "purple")) + scale_fill_manual(values=c("3" = "gray39", "4" = "green", "1" = "orange", "2" = "blue", "5" = "purple")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))
dev.off()


