
## SUPPLEMENTARY FIGURE 4

################################################
##
## to fix according to how we update the manuscript 

##############################################

library(readr)
library(ggplot2)

patients=c('P8','P12','P2','P6','P10','P59','P41','P63','P7','P3','P1')
clusters<-c()
cells<-c()
events<-c()
stem<-c()

for (p in patients){
    f<-read_tsv(paste0('../',p,'/local_samples.tsv'))
    clusters<-append(clusters,length(unique(f$cluster)))
    events<-append(events,nrow(get(p)))
    cells<-append(cells,nrow(f))
    
    if (p=='P6'){
        stem<-append(stem,length(which(f$cluster=='3')))
    }
    else if(p=='P7'){
        stem<-append(stem,length(which(f$cluster=='1')))
    }
    else if(p=='P59'){
        stem<-append(stem,length(which(f$cluster=='3')))
    } 
    else if(p=='P63'){
        stem<-append(stem,length(which(f$cluster=='0')))
    } 
    else {
        stem<-append(stem,length(which(f$cluster=='2')))
    }
}

new_df<-data.frame(patients=patients,n_cells=cells, n_clusters=clusters,AS_events=events, `stem-like size`=stem)
write_tsv(new_df,'metadata_patients_stemlike.tsv')

read_ts()

# cells plot

l<-summary(lm(n_cells~AS_events,data=new_df))
slope<-round(l$coefficients['AS_events','Estimate'],3)
pval <-round(l$coefficients['AS_events', "Pr(>|t|)"],2)

png(paste0('figures/scatter_stemlike_patients_cells.png'),height = 2500, width = 2500, res = 600)

ggplot(new_df, aes(x=AS_events, y=n_cells)) + geom_point(size=3,shape=5,stroke=2) +
    geom_smooth(method = "lm",) + theme_bw() + labs(x='AS events', y='Cells')+ ylim(0,2000)+
    annotate("text", x = 350,y = 100, label = paste0('m: ',slope,'\n','adj-p: ',pval), color = "blue", size = 4,hjust=1)+
    
    theme( axis.text.x = element_text(color = "black",size=12),
           axis.text.y = element_text(color = "black",size=12),
           axis.title.x = element_text(size=14),
           axis.title.y = element_text(size=14))

dev.off()

# stem plot

l<-summary(lm(stem~AS_events,data=new_df))
slope<-round(l$coefficients['AS_events','Estimate'],3)
pval <-round(l$coefficients['AS_events', "Pr(>|t|)"],2)

png(paste0('figures/scatter_stemlike_patients_stem_size.png'),height = 2500, width = 2500, res = 600)

ggplot(new_df, aes(x=AS_events, y=stem)) + geom_point(size=3,shape=2,stroke=2) +
    geom_smooth(method = "lm",) + theme_bw() + labs(x='AS events', y='Stem-like size')+ ylim(0,500)+
    annotate("text", x = 350,y = 440, label = paste0('m: ',slope,'\n','adj-p: ',pval), color = "blue", size = 4,hjust=1)+
    
    theme( axis.text.x = element_text(color = "black",size=12),
           axis.text.y = element_text(color = "black",size=12),
           axis.title.x = element_text(size=14),
           axis.title.y = element_text(size=14))

dev.off()
