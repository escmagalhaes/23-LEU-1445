#updated by E.S. de Camargo Magalhaes on 31/01/2024 using R v4.3.2 "Eye Holes"

#####Create datatables and color object----

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(matlab))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(progenyClust))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RCy3))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(corrplot))

#Import patient dataset
setwd('/data') #customize folder for importing patient data file
rawtable<-read_excel("dataset_old.xlsx") 

#Dataframe with protein data only
rppa_data<-as.data.frame(rawtable[c(-2:-83)])
row.names(rppa_data)<-rppa_data$id
rppa_data<-rppa_data[,c(-1)]
rppa_data<-(t(rppa_data))
rppa_data<-as.data.frame(rppa_data)

#Dataframe with clinical data only
clin_data<-as.data.frame(rawtable[,c(1:83)])
str(clin_data)
selected<-colnames(clin_data[,c(3:15,25:ncol(clin_data))])
clin_data[selected]<-lapply(clin_data[selected],factor)
str(clin_data)

#Create color object
mycolors1<-c('red2','dodgerblue','darkgoldenrod','blue3','purple3','green3','darkorange2','navajowhite','palegreen','lightsalmon')

#####Create VH and CC patient subsets to identify prognostic protrins (PS1)----

##Protein list for VH patients
vh<- tibble::rownames_to_column(as.data.frame(t(rppa_data)), "id")
vh$id<-as.integer(vh$id)
vh<-merge(clin_data,vh,by='id')
vh<-vh[vh$hma=='yes' & vh$venetoclax=='yes' & !is.na(vh$hma),]
rownames(vh)<-vh$id 
dim(vh)

#Several quantile splits for all proteins regarding survival
vh_surv<-vh %>% relocate((ncol(clin_data)+1):ncol(vh), everything())
pval_ps1<-data.frame(matrix(NA,nrow=nrow(rppa_data),ncol=6))
pval_ps1[[1]]<-colnames(vh_surv[1:nrow(rppa_data)])
colnames(pval_ps1)<-c('ptn_names',paste0('surv_',c(paste0('Q',c(2:6)))))

for (i in 1:nrow(rppa_data)) {
  for (j in 2:6) {
    vh_surv$ntile<-factor(ntile(vh_surv[[i]],j))
    pval_ps1[[i,j]]<-surv_pvalue(survfit(Surv(surv_time, status) ~ vh_surv$ntile, data=vh_surv))$pval.txt
  }}

row.names(pval_ps1)<-pval_ps1$ptn_names
pval_ps1<-pval_ps1[,c(-1)]

for (i in 1:ncol(pval_ps1)) {
  pval_ps1[[i]]<-gsub("p = ", "", pval_ps1[[i]])
  pval_ps1[[i]]<-gsub("p < ", "", pval_ps1[[i]])
}

pval_ps1<-mutate_all(pval_ps1, function(x) as.numeric(as.character(x)))

#Create protein a list of VH protein subsets by choosing proteins according to the p-value
vh_ps1_list<-list(as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q2<0.01,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q3<0.01,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q4<0.01,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q5<0.01,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q6<0.01,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q2<0.05,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q3<0.05,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q4<0.05,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q5<0.05,]))]))
                  ,as.data.frame(t(vh[,c(rownames(pval_ps1[pval_ps1$surv_Q6<0.05,]))]))
)

for (i in 1:length(vh_ps1_list)){print(dim(vh_ps1_list[[i]]))}

# Same Protein lists, but for CC patients
cc<- tibble::rownames_to_column(as.data.frame(t(rppa_data)), "id")
cc$id<-as.integer(cc$id)
cc<-merge(clin_data,cc, by='id')
cc<-cc[cc$hma=='no' & cc$venetoclax=='no' & !is.na(cc$hma) |
         cc$hma=='yes' & cc$venetoclax=='yes' & cc$arac=='yes'& !is.na(cc$hma),]
#Include ARAC+VH patients
rownames(cc)<-cc$id 
dim(cc)

cc_ps1_list<-list(as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q2<0.01,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q3<0.01,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q4<0.01,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q5<0.01,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q6<0.01,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q2<0.05,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q3<0.05,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q4<0.05,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q5<0.05,]))]))
                  ,as.data.frame(t(cc[,c(rownames(pval_ps1[pval_ps1$surv_Q6<0.05,]))]))
)


##PS1 protein groups analysis----

##The following section section was designed for exploratory data analysis##
##Output is very large##
### DO NOT RUN ###

#Define distance functions based on correlation
#row_dist_fun <- function(x) as.dist((1-cor(t(x),method='pearson')))

#VH patients
#for (i in 1:length(vh_ps1_list)){
#  skip_to_next <- FALSE
#  tryCatch({
#    setwd('/results')#custom destination folder for file generation
#    pdf(file=paste(c('ps1_vh_gp'),as.character(i),'.pdf',sep=''))

#    #Clustering with progenyclust()
#    pg_vh<-progenyClust(t(vh_ps1_list[[i]]),FUNclust=hclust.progenyClust,method='both',score.invert = F,
#                            ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
#   vh_clusters<-cbind(t(vh_ps1_list[[i]]),as.data.frame(pg_vh$cluster))
#   vh_clusters<-t(vh_clusters)

#    #Create object with clusters and clinical features
#    vh_clin<-vh[,c(-(ncol(clin_data)+1):-ncol(vh))]
#    vh_clin<-as.data.frame(cbind(vh_clin,as.data.frame(pg_vh$cluster)))

#    #Create object with clusters, protein expression and clinical features and filter AraC patients
#    vh_clusters<-tibble::rownames_to_column(as.data.frame(t(vh_ps1_list[[i]])),'id')
#    vh_clusters$id<-as.integer(vh_clusters$id)
#    vh_clusters<-merge(vh_clusters,vh_clin,by='id')
#    rownames(vh_clusters)<-vh_clusters$id
#    vh_clusters<-vh_clusters[vh$hma=='yes' & !is.na(vh$hma) & vh$venetoclax=='yes' & vh$arac=='no'
#                             ,c(2:(ncol(t(vh_ps1_list[[i]]))+1),(ncol(vh_clusters)-ncol(pg_vh$cluster)+1):ncol(vh_clusters))]
#    vh_clusters<-t(vh_clusters)
#    vh_clin<-vh_clin[vh$hma=='yes' & !is.na(vh$hma) & vh$venetoclax=='yes' & vh$arac=='no',]

#    #Loop for plotting survival plots by cluster
#    for (j in (ncol(vh_clin)-ncol(pg_vh$cluster)+1):ncol(vh_clin)) {
#      vh_clin[[j]]<-factor(vh_clin[[j]])
#      levels(vh_clin[[j]])<-paste0('C', c(1:nlevels(vh_clin[[j]])))
#      km_surv_fit_vh_avg<-survfit(Surv(surv_time,status)~1,data=vh_clin)
#      km_surv_fit_vh_cluster<-survfit(Surv(surv_time,status)~vh_clin[[j]],data=vh_clin)
#      vh_surv_fit<-list(Overall = km_surv_fit_vh_avg, Cluster = km_surv_fit_vh_cluster)
#      print(ggsurvplot_combine(vh_surv_fit,data=vh_clin,pal=c('black',mycolors1[1:nlevels(vh_clin[[j]])])
#                               ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(vh_clin[[j]]))))
#                               ,title= paste(c("Survival Clusters VH patients set "),as.character(i),sep='')
#                               ,xlab='Time (years)',xlim = c(0,10),break.x.by=1,conf.int=F
#                               ,font.main = 20,font.x = 15,font.y = 15,ylab='Cumulative Survival probability',pval=F
#                               ,risk.table=T,tables.col = "strata",risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      
#    }

#    #Loop for plotting remission plots by cluster
#    for (j in (ncol(vh_clin)-ncol(pg_vh$cluster)+1):ncol(vh_clin)) {
#      tryCatch({
#        vh_rem<-vh_clin[!is.na(vh_clin$relapse),]
#        km_rem_fit_vh_avg<-survfit(Surv(rem_time, relapse)~1,data=vh_rem)
#        km_rem_fit_vh_cluster<-survfit(Surv(rem_time,relapse)~vh_rem[[j]],data=vh_rem)
#        vh_rem_fit<-list(Overall = km_rem_fit_vh_avg, Cluster = km_rem_fit_vh_cluster)
#        print(ggsurvplot_combine(vh_rem_fit,data=vh_clin,pal=c('black',mycolors1[1:nlevels(vh_rem[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(vh_rem[[j]]))))
#                                 ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)'
#                                 ,title= paste(c("Remission Clusters VH patients set "),as.character(i),sep='')
#                                 ,font.main = 20,font.x = 15,font.y = 15,ylab='Complete Remission probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }

#    #Loop for plotting heatmaps by cluster
#    for (j in (nrow(vh_clusters)-ncol(pg_vh$cluster)+1):nrow(vh_clusters)) {
#      vh_clusters<-vh_clusters[,order(vh_clusters[j,])]
#      vh_anno<-as.data.frame(factor(vh_clusters[j,]))
#      colnames(vh_anno)<-'cluster'
#      levels(vh_anno$cluster)<-paste0('C',c(1:nlevels(vh_anno$cluster)))
#      order_by_cl_vh<-vh_clusters[1:nrow(vh_ps1_list[[i]]),]
#      vh_colors<-mycolors1[1:nlevels(vh_anno$cluster)]
#      names(vh_colors)<-unique(vh_anno$cluster)
#      vh_colors<-list(cluster=vh_colors)

#      if(min(vh_ps1_list[[i]])>=-2 & max(vh_ps1_list[[i]])<=2){
#        vh_break=seq(-2.4,2.4,0.4)
#      }else if(min(vh_ps1_list[[i]])>=-2){
#        vh_break=c(seq(-2.4,2,0.4),max(vh_ps1_list[[i]]))
#      }else if(max(vh_ps1_list[[i]])<=2){
#        vh_break=c(min(vh_ps1_list[[i]]),seq(-2,2.4,0.4))
#      }else{
#        vh_break=c(min(vh_ps1_list[[i]]),seq(-2,2,0.4),max(vh_ps1_list[[i]]))
#      }

#      pheatmap(order_by_cl_vh,annotation_colors=vh_colors,main=paste(c("Protein Clusters VH patients set "),as.character(i),sep='')
#               ,annotation_col=vh_anno,cluster_rows=T,cluster_cols=F,clustering_method ='ward.D2',treeheight_row=0
#               ,col=jet.colors(length(vh_break)-1),breaks=vh_break,scale="none",show_colnames=F,fontsize_row=8)
#    }
#    dev.off()
#    dev.off()
#  }, error = function(e) {skip_to_next <<- TRUE})
#  if(skip_to_next) { next }
#}


#CC patients
#for (i in 1:length(cc_ps1_list)){
#  skip_to_next <- FALSE
#  tryCatch({
#    setwd('/results')#custom destination folder for file generation
#    pdf(file=paste(c('ps1_cc_gp'),as.character(i),'.pdf',sep=''))
#    
#    #Clustering with progenyclust()
#    pg_cc<-progenyClust(t(cc_ps1_list[[i]]),FUNclust=hclust.progenyClust,method='both',score.invert = F,
#                               ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
#    cc_clusters<-cbind(t(cc_ps1_list[[i]]),as.data.frame(pg_cc$cluster))
#    cc_clusters<-t(cc_clusters)
# 
#    #Create object with clusters and clinical features
#    cc_clin<-cc[,c(-(ncol(clin_data)+1):-ncol(cc))]
#    cc_clin<-as.data.frame(cbind(cc_clin,as.data.frame(pg_cc$cluster)))

#    #Create object with clusters, protein expression and clinical features and filter patients
#    cc_clusters<-tibble::rownames_to_column(as.data.frame(t(cc_ps1_list[[i]])),'id')
#    cc_clusters$id<-as.integer(cc_clusters$id)
#    cc_clusters<-merge(cc_clusters,cc_clin,by='id')
#    rownames(cc_clusters)<-cc_clusters$id
#    cc_clusters<-cc_clusters[cc_clusters$hma=='no' & !is.na(cc_clusters$hma) 
#                                           & cc_clusters$venetoclax=='no' 
#                                           & cc_clusters$arac=='yes' 
#                                           | cc_clusters$hma=='yes' & cc_clusters$venetoclax=='yes'
#                                           & !is.na(cc_clusters$hma) & cc_clusters$arac=='yes'
#                                           ,c(2:(ncol(t(cc_ps1_list[[i]]))+1),(ncol(cc_clusters)-ncol(pg_cc$cluster)+1):ncol(cc_clusters))]
#    cc_clusters<-t(cc_clusters)
#    cc_clin<-cc_clin[cc_clin$hma=='no' & !is.na(cc_clin$hma) 
#                                   & cc_clin$venetoclax=='no' 
#                                   & cc_clin$arac=='yes' 
#                                   | cc_clin$hma=='yes' & cc_clin$venetoclax=='yes'
#                                   & !is.na(cc_clin$hma) & cc_clin$arac=='yes',]

#    #Loop for plotting survival plots by cluster
#    for (j in (ncol(cc_clin)-ncol(pg_cc$cluster)+1):ncol(cc_clin)) {
#      cc_clin[[j]]<-factor(cc_clin[[j]])
#      levels(cc_clin[[j]])<-paste0('C', c(1:nlevels(cc_clin[[j]])))
#      km_surv_fit_cc_avg<-survfit(Surv(surv_time,status)~1,data=cc_clin)
#      km_surv_fit_cc_cluster<-survfit(Surv(surv_time,status)~cc_clin[[j]],data=cc_clin)
#      cc_surv_fit<-list(Overall = km_surv_fit_cc_avg, Cluster = km_surv_fit_cc_cluster)
#      print(ggsurvplot_combine(cc_surv_fit,data=cc_clin,pal=c('black',mycolors1[1:nlevels(cc_clin[[j]])])
#                               ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(cc_clin[[j]]))))
#                               ,title= paste(c("Survival Clusters CC patients set "),as.character(i),sep='')
#                               ,xlab='Time (years)',xlim = c(0,10),break.x.by=1,conf.int=F
#                               ,font.main = 20,font.x = 15,font.y = 15,ylab='Cumulative Survival probability',pval=T
#                               ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                               ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#    }
#    
#    #Loop for plotting remission plots by cluster
#    for (j in (ncol(cc_clin)-ncol(pg_cc$cluster)+1):ncol(cc_clin)) {
#      tryCatch({
#        cc_rem<-cc_clin[!is.na(cc_clin$relapse),]
#        km_rem_fit_cc_avg<-survfit(Surv(rem_time, relapse)~1,data=cc_rem)
#        km_rem_fit_cc_cluster<-survfit(Surv(rem_time,relapse)~cc_rem[[j]],data=cc_rem)
#        cc_rem_fit<-list(Overall = km_rem_fit_cc_avg, Cluster = km_rem_fit_cc_cluster)
#        print(ggsurvplot_combine(cc_rem_fit,data=cc_rem,pal=c('black',mycolors1[1:nlevels(cc_rem[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(cc_rem[[j]]))))
#                                 ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)'
#                                 ,title= paste(c("Remission Clusters CC patients set "),as.character(i),sep='')
#                                 ,font.main = 20,font.x = 15,font.y = 15,ylab='Complete Remission probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }

#    #Loop for plotting heatmaps by cluster
#    for (j in (nrow(cc_clusters)-ncol(pg_cc$cluster)+1):nrow(cc_clusters)) {
#      cc_clusters<-cc_clusters[,order(cc_clusters[j,])]
#      cc_anno<-as.data.frame(factor(cc_clusters[j,]))
#      colnames(cc_anno)<-'cluster'
#      levels(cc_anno$cluster)<-paste0('C',c(1:nlevels(cc_anno$cluster)))
#      order_by_cl_cc<-as.data.frame(cc_clusters[1:(nrow(cc_ps1_list[[i]])-2),])
#      order_by_cl_cc<-mutate_all(order_by_cl_cc, function(x) as.numeric(as.character(x)))
#      cc_colors<-mycolors1[1:nlevels(cc_anno$cluster)]
#      names(cc_colors)<-unique(cc_anno$cluster)
#      cc_colors<-list(cluster=cc_colors)

#      if(min(cc_ps1_list[[i]])>=-2 & max(cc_ps1_list[[i]])<=2){
#        cc_break=seq(-2.4,2.4,0.4)
#      }else if(min(cc_ps1_list[[i]])>=-2){
#        cc_break=c(seq(-2.4,2,0.4),max(cc_ps1_list[[i]]))
#      }else if(max(cc_ps1_list[[i]])<=2){
#        cc_break=c(min(cc_ps1_list[[i]]),seq(-2,2.4,0.4))
#      }else{
#        cc_break=c(min(cc_ps1_list[[i]]),seq(-2,2,0.4),max(cc_ps1_list[[i]]))
#      }

#      pheatmap(order_by_cl_cc,annotation_colors=cc_colors,main=paste(c("Protein Clusters CC patients set "),as.character(i),sep='')
#               ,annotation_col=cc_anno,cluster_rows=T,cluster_cols=F,clustering_method ='ward.D2',treeheight_row=0
#               ,col=jet.colors(length(cc_break)-1),breaks=cc_break,scale="none",show_colnames=F,fontsize_row=8)
#      
#    }
#    dev.off()
#    dev.off()
#  }, error = function(e) {skip_to_next <<- TRUE})
#  if(skip_to_next) { next }
#}

####### END OF SECTION ###

###VH_PS1_gp6x3clusters----

#Clustering with VH protein group 6 with 3 clusters (defined according to the exploratory data analysis)
pg_vh_ps1<-progenyClust(t(vh_ps1_list[[6]]),FUNclust=hclust.progenyClust,method='gap',ncluster=2:10,size=10
                        ,iteration=100,repeats=10,nrandom=10)
vh_ps1_clusters<-as.data.frame(pg_vh_ps1$cluster)
vh_ps1_clusters$C3<-factor(vh_ps1_clusters$C3)
levels(vh_ps1_clusters$C3)<-c('C2','C1','C3')
vh_ps1_clusters$C3<-relevel(vh_ps1_clusters$C3,'C3')
vh_ps1_clusters$C3<-relevel(vh_ps1_clusters$C3,'C2')
vh_ps1_clusters$C3<-relevel(vh_ps1_clusters$C3,'C1')


#Create objects with clusters and clinical features
vh_ps1_clin<-vh[,c(-(ncol(clin_data)+1):-ncol(vh))]
vh_ps1_clin<-as.data.frame(cbind(vh_ps1_clin,as.data.frame(vh_ps1_clusters$C3)))
colnames(vh_ps1_clin)<-c(colnames(vh_ps1_clin[1:ncol(clin_data)]),'cluster')
dim(vh_ps1_clin)

#Create object with clusters, protein expression and clinical features to filter AraC patients
vh_ps1_ptn<-tibble::rownames_to_column(as.data.frame(t(vh_ps1_list[[6]])),'id')
vh_ps1_ptn$id<-as.integer(vh_ps1_ptn$id)
vh_ps1_ptn<-merge(vh_ps1_ptn,vh_ps1_clin,by='id')
rownames(vh_ps1_ptn)<-vh_ps1_ptn$id
vh_ps1_ptn<-vh_ps1_ptn[vh_ps1_ptn$hma=='yes'& !is.na(vh_ps1_ptn$hma) & vh_ps1_ptn$venetoclax=='yes' & vh_ps1_ptn$arac=='no'
                       ,c(2:(ncol(t(vh_ps1_list[[6]]))+1),(ncol(vh_ps1_ptn)))]
dim(vh_ps1_ptn)

#Filter out AraC patients from Clinical data
vh_ps1_clin<-vh_ps1_clin[vh_ps1_clin$hma=='yes' & !is.na(vh_ps1_clin$hma) & vh_ps1_clin$venetoclax=='yes'& vh_ps1_clin$arac=='no',]
dim(vh_ps1_clin)

###CC_PS1_gp6x3clusters----

#Clustering with the same protein list (VH protein group 6) with 3 clusters
pg_cc_ps1<-progenyClust(t(cc_ps1_list[[6]]),FUNclust=hclust.progenyClust,method='gap',
                        ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
cc_ps1_clusters<-as.data.frame(pg_cc_ps1$cluster)
cc_ps1_clusters$C3<-factor(cc_ps1_clusters$C3)
levels(cc_ps1_clusters$C3)<-c('C2','C3','C1')
cc_ps1_clusters$C3<-relevel(cc_ps1_clusters$C3,'C3')
cc_ps1_clusters$C3<-relevel(cc_ps1_clusters$C3,'C2')
cc_ps1_clusters$C3<-relevel(cc_ps1_clusters$C3,'C1')

#Create objects with clusters and clinical features
cc_ps1_clin<-cc[,c(-(ncol(clin_data)+1):-ncol(cc))]
cc_ps1_clin<-as.data.frame(cbind(cc_ps1_clin,as.data.frame(cc_ps1_clusters$C3)))
colnames(cc_ps1_clin)<-c(colnames(cc_ps1_clin[1:ncol(clin_data)]),'cluster')
dim(cc_ps1_clin)

#Create object with clusters, protein expression and clinical features and filter AraC patients
cc_ps1_ptn<-tibble::rownames_to_column(as.data.frame(t(cc_ps1_list[[6]])),'id')
cc_ps1_ptn$id<-as.integer(cc_ps1_ptn$id)
cc_ps1_ptn<-merge(cc_ps1_ptn,cc_ps1_clin,by='id')
rownames(cc_ps1_ptn)<-cc_ps1_ptn$id
cc_ps1_ptn<-cc_ps1_ptn[cc_ps1_ptn$hma=='no' & !is.na(cc_ps1_ptn$hma) & cc_ps1_ptn$venetoclax=='no' & cc_ps1_ptn$arac=='yes'|
                         cc_ps1_ptn$hma=='yes' & cc_ps1_ptn$venetoclax=='yes' & !is.na(cc_ps1_ptn$hma) & cc_ps1_ptn$arac=='yes'
                       ,c(2:(ncol(t(cc_ps1_list[[6]]))+1),ncol(cc_ps1_ptn))]
dim(cc_ps1_ptn)

#Filter AraC patients in clinical data
cc_ps1_clin<-cc_ps1_clin[cc_ps1_clin$hma=='no' & !is.na(cc_ps1_clin$hma) & cc_ps1_clin$venetoclax=='no'  & cc_ps1_clin$arac=='yes'|
                           cc_ps1_clin$hma=='yes' & cc_ps1_clin$venetoclax=='yes'& !is.na(cc_ps1_clin$hma) & cc_ps1_clin$arac=='yes',]
dim(cc_ps1_clin)

####Combine VH and CC data for PS1----

##Combine clinical dataframes
ps1_clin_join<-vh_ps1_clin
ps1_clin_join$trt<-rep('VH',nrow(ps1_clin_join))
ps1_cc_clin_join<-cc_ps1_clin
ps1_cc_clin_join$trt<-rep('CC',nrow(ps1_cc_clin_join))
ps1_clin_join<-as.data.frame(rbind(ps1_clin_join,ps1_cc_clin_join))
ps1_clin_join$cluster2<-interaction(ps1_clin_join$cluster,ps1_clin_join$trt)
levels(ps1_clin_join$cluster2)<-c('C1 CC','C2 CC','C3 CC','C1 VH','C2 VH','C3 VH')
ps1_clin_join$cluster2<-relevel(ps1_clin_join$cluster2,'C3 CC')
ps1_clin_join$cluster2<-relevel(ps1_clin_join$cluster2,'C2 CC')
ps1_clin_join$cluster2<-relevel(ps1_clin_join$cluster2,'C1 CC')
ps1_clin_join$cluster2<-relevel(ps1_clin_join$cluster2,'C3 VH')
ps1_clin_join$cluster2<-relevel(ps1_clin_join$cluster2,'C2 VH')
ps1_clin_join$cluster2<-relevel(ps1_clin_join$cluster2,'C1 VH')
dim(ps1_clin_join)

#Merge ptn dataframes
ps1_ptn_join<-vh_ps1_ptn
levels(ps1_ptn_join$cluster)<-c('PS1 VH C1','PS1 VH C2','PS1 VH C3')
cc_ptn_join<-cc_ps1_ptn
levels(cc_ptn_join$cluster)<-c('PS1 CC C1','PS1 CC C2','PS1 CC C3')
ps1_ptn_join<-as.data.frame(rbind(ps1_ptn_join,cc_ptn_join))

#####Create Cluster C2 patient subsets to identify prognostic proteins (PS2)----

##Protein list by CC patients
cc_c2<- tibble::rownames_to_column(as.data.frame(t(rppa_data)), "id")
cc_c2$id<-as.integer(cc_c2$id)
cc_c2<-merge(cc_ps1_clin,cc_c2, by='id')
cc_c2<-cc_c2[cc_c2$cluster=='C2',]
rownames(cc_c2)<-cc_c2$id 
dim(cc_c2)

#Several quantile splits for all proteins regarding survival
cc_c2_surv<-cc_c2 %>% relocate((ncol(clin_data)+2):ncol(cc_c2), everything())
pval_ps2<-data.frame(matrix(NA,nrow=nrow(rppa_data),ncol=6))
pval_ps2[[1]]<-colnames(cc_c2_surv[1:nrow(rppa_data)])
colnames(pval_ps2)<-c('ptn_names',paste0('surv_',c(paste0('Q',c(2:6)))))

for (i in 1:nrow(rppa_data)) {
  for (j in 2:6) {
    cc_c2_surv$ntile<-factor(ntile(cc_c2_surv[[i]],j))
    pval_ps2[[i,j]]<-surv_pvalue(survfit(Surv(surv_time,status)~cc_c2_surv$ntile,data=cc_c2_surv))$pval.txt
  }}

row.names(pval_ps2)<-pval_ps2$ptn_names
pval_ps2<-pval_ps2[,c(-1)]

for (i in 1:ncol(pval_ps2)) {
  pval_ps2[[i]]<-gsub("p = ", "", pval_ps2[[i]])
  pval_ps2[[i]]<-gsub("p < ", "", pval_ps2[[i]])
}

pval_ps2<-mutate_all(pval_ps2, function(x) as.numeric(as.character(x)))

#Create protein a list of CC subsets by choosing proteins according to the p-value
cc_c2_list<-list(as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q2<0.01,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q3<0.01,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q4<0.01,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q5<0.01,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q6<0.01,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q2<0.05,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q3<0.05,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q4<0.05,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q5<0.05,]))]))
                 ,as.data.frame(t(cc_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q6<0.05,]))]))
)

for (i in 1:length(cc_c2_list)){print(dim(cc_c2_list[[i]]))}

#Create protein lists for VH patients based on CC subsets
vh_c2<- tibble::rownames_to_column(as.data.frame(t(rppa_data)), "id")
vh_c2$id<-as.integer(vh_c2$id)
vh_c2<-merge(vh_ps1_clin,vh_c2, by='id')
vh_c2<-vh_c2[vh_c2$cluster=='C2',]
rownames(vh_c2)<-vh_c2$id 
dim(vh_c2)

vh_c2_list<-list(as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q2<0.01,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q3<0.01,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q4<0.01,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q5<0.01,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q6<0.01,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q2<0.05,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q3<0.05,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q4<0.05,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q5<0.05,]))]))
                 ,as.data.frame(t(vh_c2[,c(rownames(pval_ps2[pval_ps2$surv_Q6<0.05,]))]))
)

##PS2 protein groups analysis----

##The following section section was designed for exploratory data analysis##
##Output is very large##
#### DO NOT RUN ###

#VH_C2 patients
#for (i in 1:length(vh_c2_list)){
#  skip_to_next <- FALSE
#  tryCatch({
#    setwd('/results')#custom destination folder for file generation
#    pdf(file=paste(c('ps2_vh_gp'),as.character(i),'.pdf',sep=''))
#    
#    #Clustering with progenyclust()
#    pg_vh_c2<-progenyClust(t(vh_c2_list[[i]]),FUNclust=hclust.progenyClust,method='both',score.invert = F,
#                                   ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
#    vh_c2_clusters<-cbind(t(vh_c2_list[[i]]),as.data.frame(pg_vh_c2$cluster))
#    vh_c2_clusters<-t(vh_c2_clusters)
#    d_row_vh_c2<-row_dist_fun(vh_c2_list[[i]])
#    
#    #Create object with clusters and clinical features
#    vh_c2_clin<-vh_c2[,c(-(ncol(clin_data)+1):-ncol(vh_c2))]
#    vh_c2_clin<-as.data.frame(cbind(vh_c2_clin,as.data.frame(pg_vh_c2$cluster)))
#    
#    #Loop for plotting survival plots by cluster
#    for (j in (ncol(vh_c2_clin)-ncol(pg_vh_c2$cluster)+1):ncol(vh_c2_clin)) {
#      tryCatch({
#        vh_c2_clin[[j]]<-factor(vh_c2_clin[[j]])
#        levels(vh_c2_clin[[j]])<-paste0('C', c(1:nlevels(vh_c2_clin[[j]])))
#        km_surv_fit_vh_c2_avg<-survfit(Surv(surv_time,status)~1,data=vh_c2_clin)
#        km_surv_fit_vh_c2_cluster<-survfit(Surv(surv_time,status)~vh_c2_clin[[j]],data=vh_c2_clin)
#        vh_c2_surv_fit<-list(Overall = km_surv_fit_vh_c2_avg, Cluster = km_surv_fit_vh_c2_cluster)
#        print(ggsurvplot_combine(vh_c2_surv_fit,data=vh_c2_clin,pal=c('black',mycolors1[1:nlevels(vh_c2_clin[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(vh_c2_clin[[j]]))))
#                                 ,title= paste(c("Survival PS2 VH gp "),as.character(i),sep='')
#                                 ,xlab='Time (years)',xlim = c(0,9),break.x.by=1,conf.int=F
#                                 ,font.main = 16,font.x = 15,font.y = 15,ylab='Cumulative Survival probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }
#    #Loop for plotting remission plots by cluster
#    for (j in (ncol(vh_c2_clin)-ncol(pg_vh_c2$cluster)+1):ncol(vh_c2_clin)) {
#      tryCatch({
#        vh_c2_rem<-vh_c2_clin[!is.na(vh_c2_clin$relapse),]
#        km_rem_fit_vh_c2_avg<-survfit(Surv(rem_time, relapse)~1,data=vh_c2_rem)
#        km_rem_fit_vh_c2_cluster<-survfit(Surv(rem_time,relapse)~vh_c2_rem[[j]],data=vh_c2_rem)
#        vh_c2_rem_fit<-list(Overall = km_rem_fit_vh_c2_avg, Cluster = km_rem_fit_vh_c2_cluster)
#        print(ggsurvplot_combine(vh_c2_rem_fit,data=vh_c2_rem,pal=c('black',mycolors1[1:nlevels(vh_c2_rem[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(vh_c2_rem[[j]]))))
#                                 ,xlim = c(0,9),break.x.by=1,conf.int=F,xlab='Time (years)'
#                                 ,title= paste(c("Remission PS2 VH gp"),as.character(i),sep='')
#                                 ,font.main = 16,font.x = 15,font.y = 15,ylab='Complete Remission probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }
#    
#    #Loop for plotting heatmaps by cluster
#    for (j in (nrow(vh_c2_clusters)-ncol(pg_vh_c2$cluster)+1):nrow(vh_c2_clusters)) {
#      vh_c2_clusters<-vh_c2_clusters[,order(vh_c2_clusters[j,])]
#      vh_c2_anno<-as.data.frame(factor(vh_c2_clusters[j,]))
#      colnames(vh_c2_anno)<-'cluster'
#      levels(vh_c2_anno$cluster)<-paste0('C',c(1:nlevels(vh_c2_anno$cluster)))
#      order_by_cl_vh_c2<-vh_c2_clusters[1:nrow(vh_c2_list[[i]]),]
#      vh_c2_colors<-mycolors1[1:nlevels(vh_c2_anno$cluster)]
#      names(vh_c2_colors)<-unique(vh_c2_anno$cluster)
#      vh_c2_colors<-list(cluster=vh_c2_colors)
#      
#      if(min(vh_c2_list[[i]])>=-2 & max(vh_c2_list[[i]])<=2){
#        vh_c2_break=seq(-2.4,2.4,0.4)
#      }else if(min(vh_c2_list[[i]])>=-2){
#        vh_c2_break=c(seq(-2.4,2,0.4),max(vh_c2_list[[i]]))
#      }else if(max(vh_c2_list[[i]])<=2){
#        vh_c2_break=c(min(vh_c2_list[[i]]),seq(-2,2.4,0.4))
#      }else{
#        vh_c2_break=c(min(vh_c2_list[[i]]),seq(-2,2,0.4),max(vh_c2_list[[i]]))
#      }
#      
#      pheatmap(order_by_cl_vh_c2,annotation_colors=vh_c2_colors,main=paste(c("PS2 VH Clusters gp"),as.character(i),sep='')
#               ,annotation_col=vh_c2_anno,cluster_rows=T,cluster_cols=F,clustering_method ='ward.D2'
#               ,clustering_distance_rows=d_row_vh_c2,col=jet.colors(length(vh_c2_break)-1)
#               ,breaks=vh_c2_break,scale="none",show_colnames=F,fontsize_row=8)
#    }
#    dev.off()
#    dev.off()
#  }, error = function(e) {skip_to_next <<- TRUE})
#  if(skip_to_next) { next }
#}

#CC_C2 patients
#for (i in 1:length(cc_c2_list)){
#  skip_to_next <- FALSE
#  tryCatch({
#    setwd('/results')#custom destination folder for file generation
#    pdf(file=paste(c('ps2_cc_gp'),as.character(i),'.pdf',sep=''))
#    
#    #Clustering with progenyclust()
#    pg_cc_c2<-progenyClust(t(cc_c2_list[[i]]),FUNclust=hclust.progenyClust,method='both',score.invert = F,
#                                         ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
#    cc_c2_clusters<-cbind(t(cc_c2_list[[i]]),as.data.frame(pg_cc_c2$cluster))
#    cc_c2_clusters<-t(cc_c2_clusters)
#    d_row_cc_c2<-row_dist_fun(cc_c2_list[[i]])
#    
#    #Create object with clusters and clinical features
#    cc_c2_clin<-cc_c2[,c(-(ncol(clin_data)+1):-ncol(cc_c2))]
#    cc_c2_clin<-as.data.frame(cbind(cc_c2_clin,as.data.frame(pg_cc_c2$cluster)))
#   
#    #Loop for plotting survival plots by cluster
#    for (j in (ncol(cc_c2_clin)-ncol(pg_cc_c2$cluster)+1):ncol(cc_c2_clin)) {
#      cc_c2_clin[[j]]<-factor(cc_c2_clin[[j]])
#      levels(cc_c2_clin[[j]])<-paste0('C', c(1:nlevels(cc_c2_clin[[j]])))
#      km_surv_fit_cc_c2_avg<-survfit(Surv(surv_time,status)~1,data=cc_c2_clin)
#      km_surv_fit_cc_c2_cluster<-survfit(Surv(surv_time,status)~cc_c2_clin[[j]],data=cc_c2_clin)
#      cc_c2_surv_fit<-list(Overall = km_surv_fit_cc_c2_avg, Cluster = km_surv_fit_cc_c2_cluster)
#      print(ggsurvplot_combine(cc_c2_surv_fit,data=cc_c2_clin,pal=c('black',mycolors1[1:nlevels(cc_c2_clin[[j]])])
#                               ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(cc_c2_clin[[j]]))))
#                               ,title= paste(c("Survival PS2 CC patients gp "),as.character(i),sep='')
#                               ,xlab='Time (years)',xlim = c(0,9),break.x.by=1,conf.int=F
#                               ,font.main = 16,font.x = 15,font.y = 15,ylab='Cumulative Survival probability',pval=T
#                               ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                               ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#    }
#    #Loop for plotting remission plots by cluster
#    for (j in (ncol(cc_c2_clin)-ncol(pg_cc_c2$cluster)+1):ncol(cc_c2_clin)) {
#      tryCatch({
#        cc_c2_rem<-cc_c2_clin[!is.na(cc_c2_clin$relapse),]
#        km_rem_fit_cc_c2_avg<-survfit(Surv(rem_time, relapse)~1,data=cc_c2_rem)
#        km_rem_fit_cc_c2_cluster<-survfit(Surv(rem_time,relapse)~cc_c2_rem[[j]],data=cc_c2_rem)
#        cc_c2_rem_fit<-list(Overall = km_rem_fit_cc_c2_avg, Cluster = km_rem_fit_cc_c2_cluster)
#        print(ggsurvplot_combine(cc_c2_rem_fit,data=cc_c2_rem,pal=c('black',mycolors1[1:nlevels(cc_c2_rem[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(cc_c2_rem[[j]]))))
#                                 ,xlim = c(0,9),break.x.by=1,conf.int=F,xlab='Time (years)'
#                                 ,title= paste(c("Remission PS2 CC gp"),as.character(i),sep='')
#                                 ,font.main = 16,font.x = 15,font.y = 15,ylab='Complete Remission probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }
#    
#    #Loop for plotting heatmaps by cluster
#    for (j in (nrow(cc_c2_clusters)-ncol(pg_cc_c2$cluster)+1):nrow(cc_c2_clusters)) {
#      cc_c2_clusters<-cc_c2_clusters[,order(cc_c2_clusters[j,])]
#      cc_c2_anno<-as.data.frame(factor(cc_c2_clusters[j,]))
#      colnames(cc_c2_anno)<-'cluster'
#      levels(cc_c2_anno$cluster)<-paste0('C',c(1:nlevels(cc_c2_anno$cluster)))
#      order_by_cl_cc_c2<-cc_c2_clusters[1:nrow(cc_c2_list[[i]]),]
#      cc_c2_colors<-mycolors1[1:nlevels(cc_c2_anno$cluster)]
#      names(cc_c2_colors)<-unique(cc_c2_anno$cluster)
#      cc_c2_colors<-list(cluster=cc_c2_colors)
#      
#      if(min(cc_c2_list[[i]])>=-2 & max(cc_c2_list[[i]])<=2){
#        cc_c2_break=seq(-2.4,2.4,0.4)
#      }else if(min(cc_c2_list[[i]])>=-2){
#        cc_c2_break=c(seq(-2.4,2,0.4),max(cc_c2_list[[i]]))
#      }else if(max(cc_c2_list[[i]])<=2){
#        cc_c2_break=c(min(cc_c2_list[[i]]),seq(-2,2.4,0.4))
#      }else{
#        cc_c2_break=c(min(cc_c2_list[[i]]),seq(-2,2,0.4),max(cc_c2_list[[i]]))
#      }
#      
#      pheatmap(order_by_cl_cc_c2,annotation_colors=cc_c2_colors,main=paste(c("PS2 CC Clusters gp"),as.character(i),sep='')
#               ,annotation_col=cc_c2_anno,cluster_rows=T,cluster_cols=F,clustering_method ='ward.D2'
#               ,clustering_distance_rows=d_row_cc_c2,col=jet.colors(length(cc_c2_break)-1)
#               ,breaks=cc_c2_break,scale="none",show_colnames=F,fontsize_row=8)
#    }
#    dev.off()
#    dev.off()
#  }, error = function(e) {skip_to_next <<- TRUE})
#  if(skip_to_next) { next }
#}

####### END OF SECTION ###

###VH_PS2_gp2x2clusters----

#Clustering with VH protein group 2 with 2 clusters (defined according to the exploratory data analysis)
pg_vh_ps2<-progenyClust(t(vh_c2_list[[2]]),FUNclust=hclust.progenyClust,method='gap',ncluster=2:10
                        ,size=10,iteration=100,repeats=10,nrandom=10)
vh_ps2_clusters<-as.data.frame(pg_vh_ps2$cluster)
vh_ps2_clusters$C2<-factor(vh_ps2_clusters$C2)
levels(vh_ps2_clusters$C2)<-c('C1','C2')
vh_ps2_clusters$C2<-relevel(vh_ps2_clusters$C2,'C2')
vh_ps2_clusters$C2<-relevel(vh_ps2_clusters$C2,'C1')

#Create objects with clusters and clinical features
vh_ps2_clin<-vh_c2[,c(-(ncol(clin_data)+1):-ncol(vh_c2))]
vh_ps2_clin<-as.data.frame(cbind(vh_ps2_clin,as.data.frame(vh_ps2_clusters$C2)))
colnames(vh_ps2_clin)<-c(colnames(vh_ps2_clin[1:ncol(clin_data)]),'cluster')
dim(vh_ps2_clin)

#Create objects with clusters and protein data
vh_ps2_ptn<-cbind(t(vh_c2_list[[2]]),as.data.frame(vh_ps2_clusters$C2))
colnames(vh_ps2_ptn)<-c(colnames(t(vh_c2_list[[2]])),'cluster')
dim(vh_ps2_ptn)

###CC_PS2_gp2x2clusters----

#Clustering with CC protein group 2 with 2 clusters (defined according to the exploratory data analysis)
pg_cc_ps2<-progenyClust(t(cc_c2_list[[2]]),FUNclust=hclust.progenyClust,method='gap',
                        ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
cc_ps2_clusters<-as.data.frame(pg_cc_ps2$cluster)
cc_ps2_clusters$C2<-factor(cc_ps2_clusters$C2)
levels(cc_ps2_clusters$C2)<-c('C2','C1')
cc_ps2_clusters$C2<-relevel(cc_ps2_clusters$C2,'C2')
cc_ps2_clusters$C2<-relevel(cc_ps2_clusters$C2,'C1')

#Create objects with clusters and clinical features
cc_ps2_clin<-cc_c2[,c(-(ncol(clin_data)+1):-ncol(cc_c2))]
cc_ps2_clin<-as.data.frame(cbind(cc_ps2_clin,as.data.frame(cc_ps2_clusters$C2)))
colnames(cc_ps2_clin)<-c(colnames(cc_ps2_clin[1:ncol(clin_data)]),'cluster')
dim(cc_ps2_clin)

#Create objects with clusters and protein data
cc_ps2_ptn<-cbind(t(cc_c2_list[[2]]),as.data.frame(cc_ps2_clusters$C2))
colnames(cc_ps2_ptn)<-c(colnames(t(cc_c2_list[[2]])),'cluster')
dim(cc_ps2_ptn)


####Combine VH and CC data for PS2----

##Combine clinical dataframes
ps2_clin_join<-vh_ps2_clin
ps2_clin_join$trt<-rep('VH',nrow(ps2_clin_join))
ps2_cc_clin_join<-cc_ps2_clin
ps2_cc_clin_join$trt<-rep('CC',nrow(ps2_cc_clin_join))
ps2_clin_join<-as.data.frame(rbind(ps2_clin_join,ps2_cc_clin_join))
ps2_clin_join<-ps2_clin_join
ps2_clin_join$cluster2<-interaction(ps2_clin_join$cluster,ps2_clin_join$trt)
levels(ps2_clin_join$cluster2)<-c('C1 CC','C2 CC','C1 VH','C2 VH')
ps2_clin_join$cluster2<-relevel(ps2_clin_join$cluster2,'C2 CC')
ps2_clin_join$cluster2<-relevel(ps2_clin_join$cluster2,'C1 CC')
ps2_clin_join$cluster2<-relevel(ps2_clin_join$cluster2,'C2 VH')
ps2_clin_join$cluster2<-relevel(ps2_clin_join$cluster2,'C1 VH')
dim(ps2_clin_join)

#Merge ptn dataframes
ps2_ptn_join<-vh_ps2_ptn
levels(ps2_ptn_join$cluster)<-c('PS2 VH C1','PS2 VH C2')
ps2_cc_ptn_join<-cc_ps2_ptn
levels(ps2_cc_ptn_join$cluster)<-c('PS2 CC C1','PS2 CC C2')
ps2_ptn_join<-as.data.frame(rbind(ps2_ptn_join,ps2_cc_ptn_join))

###Create Cluster C3 patient subsets toto identify prognostic proteins (PS3)----

##Protein list by CC patients
cc_c3<- tibble::rownames_to_column(as.data.frame(t(rppa_data)), "id")
cc_c3$id<-as.integer(cc_c3$id)
cc_c3<-merge(cc_ps1_clin,cc_c3, by='id')
cc_c3<-cc_c3[cc_c3$cluster=='C3',]
rownames(cc_c3)<-cc_c3$id 
dim(cc_c3)

#Several quantile splits for all proteins regarding survival
cc_c3_surv<-cc_c3 %>% relocate((ncol(clin_data)+2):ncol(cc_c3), everything())
pval_ps3<-data.frame(matrix(NA,nrow=nrow(rppa_data),ncol=6))
pval_ps3[[1]]<-colnames(cc_c3_surv[1:nrow(rppa_data)])
colnames(pval_ps3)<-c('ptn_names',paste0('surv_',c(paste0('Q',c(2:6)))))

for (i in 1:nrow(rppa_data)) {
  for (j in 2:6) {
    cc_c3_surv$ntile<-factor(ntile(cc_c3_surv[[i]],j))
    pval_ps3[[i,j]]<-surv_pvalue(survfit(Surv(surv_time, status) ~ cc_c3_surv$ntile, data=cc_c3_surv))$pval.txt
  }}

row.names(pval_ps3)<-pval_ps3$ptn_names
pval_ps3<-pval_ps3[,c(-1)]

for (i in 1:ncol(pval_ps3)) {
  pval_ps3[[i]]<-gsub("p = ", "", pval_ps3[[i]])
  pval_ps3[[i]]<-gsub("p < ", "", pval_ps3[[i]])
}

pval_ps3<-mutate_all(pval_ps3, function(x) as.numeric(as.character(x)))

#Create protein a list of CC subsets by choosing proteins according to the p-value
cc_c3_list<-list(as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q2<0.01,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q3<0.01,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q4<0.01,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q5<0.01,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q6<0.01,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q2<0.05,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q3<0.05,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q4<0.05,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q5<0.05,]))]))
                 ,as.data.frame(t(cc_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q6<0.05,]))]))
)

for (i in 1:length(cc_c3_list)){print(dim(cc_c3_list[[i]]))}

#Create protein lists for VH patients based on CC subsets
vh_c3<- tibble::rownames_to_column(as.data.frame(t(rppa_data)), "id")
vh_c3$id<-as.integer(vh_c3$id)
vh_c3<-merge(vh_ps1_clin,vh_c3, by='id')
vh_c3<-vh_c3[vh_c3$cluster=='C3',]
rownames(vh_c3)<-vh_c3$id 
dim(vh_c3)

vh_c3_list<-list(as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q2<0.01,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q3<0.01,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q4<0.01,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q5<0.01,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q6<0.01,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q2<0.05,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q3<0.05,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q4<0.05,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q5<0.05,]))]))
                 ,as.data.frame(t(vh_c3[,c(rownames(pval_ps3[pval_ps3$surv_Q6<0.05,]))]))
)

##PS3 protein groups analysis----

##The following section section was designed for exploratory data analysis##
##Output is very large##
### DO NOT RUN ###

#VH patients
#for (i in 1:length(vh_c3_list)){
#  skip_to_next <- FALSE
#  tryCatch({
#    setwd('/results')#custom destination folder for file generation
#    pdf(file=paste(c('ps3_vh_gp'),as.character(i),'.pdf',sep=''))
#    
#    #Clustering with progenyclust()
#    pg_vh_c3<-progenyClust(t(vh_c3_list[[i]]),FUNclust=hclust.progenyClust,method='both',score.invert = F,
#                                   ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
#    vh_c3_clusters<-cbind(t(vh_c3_list[[i]]),as.data.frame(pg_vh_c3$cluster))
#    vh_c3_clusters<-t(vh_c3_clusters)
#    d_row_vh_c3<-row_dist_fun(vh_c3_list[[i]])
#    
#    #Create object with clusters and clinical features
#    vh_c3_clin<-vh_c3[,c(-(ncol(clin_data)+1):-ncol(vh_c3))]
#    vh_c3_clin<-as.data.frame(cbind(vh_c3_clin,as.data.frame(pg_vh_c3$cluster)))
#    
#    #Loop for plotting survival plots by cluster
#    for (j in (ncol(vh_c3_clin)-ncol(pg_vh_c3$cluster)+1):ncol(vh_c3_clin)) {
#      tryCatch({
#        vh_c3_clin[[j]]<-factor(vh_c3_clin[[j]])
#        levels(vh_c3_clin[[j]])<-paste0('C', c(1:nlevels(vh_c3_clin[[j]])))
#        km_surv_fit_vh_c3_avg<-survfit(Surv(surv_time,status)~1,data=vh_c3_clin)
#        km_surv_fit_vh_c3_cluster<-survfit(Surv(surv_time,status)~vh_c3_clin[[j]],data=vh_c3_clin)
#        vh_c3_surv_fit<-list(Overall = km_surv_fit_vh_c3_avg, Cluster = km_surv_fit_vh_c3_cluster)
#        print(ggsurvplot_combine(vh_c3_surv_fit,data=vh_c3_clin,pal=c('black',mycolors1[1:nlevels(vh_c3_clin[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(vh_c3_clin[[j]]))))
#                                 ,title= paste(c("Survival PS2 VH gp "),as.character(i),sep='')
#                                 ,xlab='Time (years)',xlim = c(0,9),break.x.by=1,conf.int=F
#                                 ,font.main = 16,font.x = 15,font.y = 15,ylab='Cumulative Survival probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }
#    #Loop for plotting remission plots by cluster
#    for (j in (ncol(vh_c3_clin)-ncol(pg_vh_c3$cluster)+1):ncol(vh_c3_clin)) {
#      tryCatch({
#        vh_c3_rem<-vh_c3_clin[!is.na(vh_c3_clin$relapse),]
#        km_rem_fit_vh_c3_avg<-survfit(Surv(rem_time, relapse)~1,data=vh_c3_rem)
#        km_rem_fit_vh_c3_cluster<-survfit(Surv(rem_time,relapse)~vh_c3_rem[[j]],data=vh_c3_rem)
#        vh_c3_rem_fit<-list(Overall = km_rem_fit_vh_c3_avg, Cluster = km_rem_fit_vh_c3_cluster)
#        print(ggsurvplot_combine(vh_c3_rem_fit,data=vh_c3_rem,pal=c('black',mycolors1[1:nlevels(vh_c3_rem[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(vh_c3_rem[[j]]))))
#                                 ,xlim = c(0,9),break.x.by=1,conf.int=F,xlab='Time (years)'
#                                 ,title= paste(c("Remission PS2 VH gp"),as.character(i),sep='')
#                                 ,font.main = 16,font.x = 15,font.y = 15,ylab='Complete Remission probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }
#    
#    #Loop for plotting heatmaps by cluster
#    for (j in (nrow(vh_c3_clusters)-ncol(pg_vh_c3$cluster)+1):nrow(vh_c3_clusters)) {
#      vh_c3_clusters<-vh_c3_clusters[,order(vh_c3_clusters[j,])]
#      vh_c3_anno<-as.data.frame(factor(vh_c3_clusters[j,]))
#      colnames(vh_c3_anno)<-'cluster'
#      levels(vh_c3_anno$cluster)<-paste0('C',c(1:nlevels(vh_c3_anno$cluster)))
#      order_by_cl_vh_c3<-vh_c3_clusters[1:nrow(vh_c3_list[[i]]),]
#      vh_c3_colors<-mycolors1[1:nlevels(vh_c3_anno$cluster)]
#      names(vh_c3_colors)<-unique(vh_c3_anno$cluster)
#      vh_c3_colors<-list(cluster=vh_c3_colors)
#      
#      if(min(vh_c3_list[[i]])>=-2 & max(vh_c3_list[[i]])<=2){
#        vh_c3_break=seq(-2.4,2.4,0.4)
#      }else if(min(vh_c3_list[[i]])>=-2){
#        vh_c3_break=c(seq(-2.4,2,0.4),max(vh_c3_list[[i]]))
#      }else if(max(vh_c3_list[[i]])<=2){
#        vh_c3_break=c(min(vh_c3_list[[i]]),seq(-2,2.4,0.4))
#      }else{
#        vh_c3_break=c(min(vh_c3_list[[i]]),seq(-2,2,0.4),max(vh_c3_list[[i]]))
#      }
#      
#      pheatmap(order_by_cl_vh_c3,annotation_colors=vh_c3_colors,main=paste(c("PS2 VH Clusters gp"),as.character(i),sep='')
#               ,annotation_col=vh_c3_anno,cluster_rows=T,cluster_cols=F,clustering_method ='ward.D2'
#               ,clustering_distance_rows=d_row_vh_c3,col=jet.colors(length(vh_c3_break)-1)
#               ,breaks=vh_c3_break,scale="none",show_colnames=F,fontsize_row=8)
#    }
#    dev.off()
#    dev.off()
#  }, error = function(e) {skip_to_next <<- TRUE})
#  if(skip_to_next) { next }
#}
#
##CC patients
#for (i in 1:length(cc_c3_list)){
#  skip_to_next <- FALSE
#  tryCatch({
#    setwd('/results')#custom destination folder for file generation
#    pdf(file=paste(c('ps3_cc_gp'),as.character(i),'.pdf',sep=''))
#    
#    #Clustering with progenyclust()
#    pg_cc_c3<-progenyClust(t(cc_c3_list[[i]]),FUNclust=hclust.progenyClust,method='both',score.invert = F,
#                                      ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
#    cc_c3_clusters<-cbind(t(cc_c3_list[[i]]),as.data.frame(pg_cc_c3$cluster))
#    cc_c3_clusters<-t(cc_c3_clusters)
#    d_row_cc_c3<-row_dist_fun(cc_c3_list[[i]])
#    
#    #Create object with clusters and clinical features
#    cc_c3_clin<-cc_c3[,c(-(ncol(clin_data)+1):-ncol(cc_c3))]
#    cc_c3_clin<-as.data.frame(cbind(cc_c3_clin,as.data.frame(pg_cc_c3$cluster)))
#    
#    #Loop for plotting survival plots by cluster
#    for (j in (ncol(cc_c3_clin)-ncol(pg_cc_c3$cluster)+1):ncol(cc_c3_clin)) {
#      cc_c3_clin[[j]]<-factor(cc_c3_clin[[j]])
#      levels(cc_c3_clin[[j]])<-paste0('C', c(1:nlevels(cc_c3_clin[[j]])))
#      km_surv_fit_cc_c3_avg<-survfit(Surv(surv_time,status)~1,data=cc_c3_clin)
#      km_surv_fit_cc_c3_cluster<-survfit(Surv(surv_time,status)~cc_c3_clin[[j]],data=cc_c3_clin)
#      cc_c3_surv_fit<-list(Overall = km_surv_fit_cc_c3_avg, Cluster = km_surv_fit_cc_c3_cluster)
#      print(ggsurvplot_combine(cc_c3_surv_fit,data=cc_c3_clin,pal=c('black',mycolors1[1:nlevels(cc_c3_clin[[j]])])
#                               ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(cc_c3_clin[[j]]))))
#                               ,title= paste(c("Survival PS2 CC patients gp "),as.character(i),sep='')
#                               ,xlab='Time (years)',xlim = c(0,9),break.x.by=1,conf.int=F
#                               ,font.main = 16,font.x = 15,font.y = 15,ylab='Cumulative Survival probability',pval=T
#                               ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                               ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#    }
#    #Loop for plotting remission plots by cluster
#    for (j in (ncol(cc_c3_clin)-ncol(pg_cc_c3$cluster)+1):ncol(cc_c3_clin)) {
#      tryCatch({
#        cc_c3_rem<-cc_c3_clin[!is.na(cc_c3_clin$relapse),]
#        km_rem_fit_cc_c3_avg<-survfit(Surv(rem_time, relapse)~1,data=cc_c3_rem)
#        km_rem_fit_cc_c3_cluster<-survfit(Surv(rem_time,relapse)~cc_c3_rem[[j]],data=cc_c3_rem)
#        cc_c3_rem_fit<-list(Overall = km_rem_fit_cc_c3_avg, Cluster = km_rem_fit_cc_c3_cluster)
#        print(ggsurvplot_combine(cc_c3_rem_fit,data=cc_c3_rem,pal=c('black',mycolors1[1:nlevels(cc_c3_rem[[j]])])
#                                 ,legend.title=" ",legend.labs=c('Overall',paste0('Cluster',c(1:nlevels(cc_c3_rem[[j]]))))
#                                 ,xlim = c(0,9),break.x.by=1,conf.int=F,xlab='Time (years)'
#                                 ,title= paste(c("Remission PS2 CC gp"),as.character(i),sep='')
#                                 ,font.main = 16,font.x = 15,font.y = 15,ylab='Complete Remission probability',pval=T
#                                 ,pval.method.coord=c(6,0.73),pval.coord=c(6,0.65),risk.table=T,tables.col = "strata"
#                                 ,risk.table.height = 0.25,ggtheme = theme_bw(),size = 0.75))
#      }, error = function(e) {skip_to_next <<- TRUE})
#      if(skip_to_next) { next }
#    }
#    
#    #Loop for plotting heatmaps by cluster
#    for (j in (nrow(cc_c3_clusters)-ncol(pg_cc_c3$cluster)+1):nrow(cc_c3_clusters)) {
#      cc_c3_clusters<-cc_c3_clusters[,order(cc_c3_clusters[j,])]
#      cc_c3_anno<-as.data.frame(factor(cc_c3_clusters[j,]))
#      colnames(cc_c3_anno)<-'cluster'
#      levels(cc_c3_anno$cluster)<-paste0('C',c(1:nlevels(cc_c3_anno$cluster)))
#      order_by_cl_cc_c3<-cc_c3_clusters[1:nrow(cc_c3_list[[i]]),]
#      cc_c3_colors<-mycolors1[1:nlevels(cc_c3_anno$cluster)]
#      names(cc_c3_colors)<-unique(cc_c3_anno$cluster)
#      cc_c3_colors<-list(cluster=cc_c3_colors)
#      
#      if(min(cc_c3_list[[i]])>=-2 & max(cc_c3_list[[i]])<=2){
#        cc_c3_break=seq(-2.4,2.4,0.4)
#      }else if(min(cc_c3_list[[i]])>=-2){
#        cc_c3_break=c(seq(-2.4,2,0.4),max(cc_c3_list[[i]]))
#      }else if(max(cc_c3_list[[i]])<=2){
#        cc_c3_break=c(min(cc_c3_list[[i]]),seq(-2,2.4,0.4))
#      }else{
#        cc_c3_break=c(min(cc_c3_list[[i]]),seq(-2,2,0.4),max(cc_c3_list[[i]]))
#      }
#      
#      pheatmap(order_by_cl_cc_c3,annotation_colors=cc_c3_colors,main=paste(c("PS2 CC Clusters gp"),as.character(i),sep='')
#               ,annotation_col=cc_c3_anno,cluster_rows=T,cluster_cols=F,clustering_method ='ward.D2'
#               ,clustering_distance_rows=d_row_cc_c3,col=jet.colors(length(cc_c3_break)-1)
#               ,breaks=cc_c3_break,scale="none",show_colnames=F,fontsize_row=8)
#    }
#    dev.off()
#    dev.off()
#  }, error = function(e) {skip_to_next <<- TRUE})
#  if(skip_to_next) { next }
#}

### END OF SECTION ###

####VH_PS3_gp4x2clusters----

#Clustering VH protein group 4 with 2 clusters (defined according to the exploratory data analysis)
pg_vh_ps3<-progenyClust(t(vh_c3_list[[4]]),FUNclust=hclust.progenyClust,method='gap',ncluster=2:10
                        ,size=10,iteration=100,repeats=10,nrandom=10)
vh_ps3_clusters<-as.data.frame(pg_vh_ps3$cluster)
vh_ps3_clusters$C2<-factor(vh_ps3_clusters$C2)
levels(vh_ps3_clusters$C2)<-c('C1','C2')
vh_ps3_clusters$C2<-relevel(vh_ps3_clusters$C2,'C2')
vh_ps3_clusters$C2<-relevel(vh_ps3_clusters$C2,'C1')

#Create objects with clusters and clinical features
vh_ps3_clin<-vh_c3[,c(-(ncol(clin_data)+1):-ncol(vh_c3))]
vh_ps3_clin<-as.data.frame(cbind(vh_ps3_clin,as.data.frame(vh_ps3_clusters$C2)))
colnames(vh_ps3_clin)<-c(colnames(vh_ps3_clin[1:ncol(clin_data)]),'cluster')
dim(vh_ps3_clin)

#Create objects with clusters and protein data
vh_ps3_ptn<-cbind(t(vh_c3_list[[4]]),as.data.frame(vh_ps3_clusters$C2))
colnames(vh_ps3_ptn)<-c(colnames(t(vh_c3_list[[4]])),'cluster')
dim(vh_ps3_ptn)


####CC_PS3_gp4x2clusters----

#Clustering CC protein group 4 with 2 clusters (defined according to the exploratory data analysis)
pg_cc_ps3<-progenyClust(t(cc_c3_list[[4]]),FUNclust=hclust.progenyClust,method='gap',
                        ncluster=2:10,size=10,iteration=100,repeats=10,nrandom=10)
cc_ps3_clusters<-as.data.frame(pg_cc_ps3$cluster)
cc_ps3_clusters$C2<-factor(cc_ps3_clusters$C2)
levels(cc_ps3_clusters$C2)<-c('C2','C1')
cc_ps3_clusters$C2<-relevel(cc_ps3_clusters$C2,'C2')
cc_ps3_clusters$C2<-relevel(cc_ps3_clusters$C2,'C1')

#Create objects with clusters and clinical features
cc_ps3_clin<-cc_c3[,c(-(ncol(clin_data)+1):-ncol(cc_c3))]
cc_ps3_clin<-as.data.frame(cbind(cc_ps3_clin,as.data.frame(cc_ps3_clusters$C2)))
colnames(cc_ps3_clin)<-c(colnames(cc_ps3_clin[1:ncol(clin_data)]),'cluster')
dim(cc_ps3_clin)

#Create objects with clusters and protein data
cc_ps3_ptn<-cbind(t(cc_c3_list[[4]]),as.data.frame(cc_ps3_clusters$C2))
colnames(cc_ps3_ptn)<-c(colnames(t(cc_c3_list[[4]])),'cluster')
dim(cc_ps3_ptn)

####Combine VH and CC PS3 dataframes----

##Combine KM plots
ps3_clin_join<-vh_ps3_clin
ps3_clin_join$trt<-rep('VH',nrow(ps3_clin_join))
ps3_cc_clin_join<-cc_ps3_clin
ps3_cc_clin_join$trt<-rep('CC',nrow(ps3_cc_clin_join))
ps3_clin_join<-as.data.frame(rbind(ps3_clin_join,ps3_cc_clin_join))
ps3_clin_join<-ps3_clin_join
ps3_clin_join$cluster2<-interaction(ps3_clin_join$cluster,ps3_clin_join$trt)
levels(ps3_clin_join$cluster2)<-c('C1 CC','C2 CC','C1 VH','C2 VH')
ps3_clin_join$cluster2<-relevel(ps3_clin_join$cluster2,'C2 CC')
ps3_clin_join$cluster2<-relevel(ps3_clin_join$cluster2,'C1 CC')
ps3_clin_join$cluster2<-relevel(ps3_clin_join$cluster2,'C2 VH')
ps3_clin_join$cluster2<-relevel(ps3_clin_join$cluster2,'C1 VH')
dim(ps3_clin_join)

#Merge ptn dataframes
ps3_ptn_join<-vh_ps3_ptn
levels(ps3_ptn_join$cluster)<-c('PS3 VH C1','PS3 VH C2')
ps3_cc_ptn_join<-cc_ps3_ptn
levels(ps3_cc_ptn_join$cluster)<-c('PS3 CC C1','PS3 CC C2')
ps3_ptn_join<-as.data.frame(rbind(ps3_ptn_join,ps3_cc_ptn_join))


####Combine Analysis of PS1,2 and 3----

##Combine clinical dataframes
clin_merge1<-ps1_clin_join[ps1_clin_join$cluster=='C1',]
clin_merge1$cluster<-droplevels(clin_merge1$cluster)
levels(clin_merge1$cluster)<-c('C1')
clin_merge2<-ps2_clin_join
levels(clin_merge2$cluster)<-c('C2','C3')
clin_merge3<-ps3_clin_join
levels(clin_merge3$cluster)<-c('C4','C5')
clin_merge<-as.data.frame(rbind(clin_merge1,clin_merge2,clin_merge3))
clin_merge$cluster2<-interaction(clin_merge$cluster,clin_merge$trt)
levels(clin_merge$cluster2)<-c('C1 CC','C2 CC','C3 CC','C4 CC','C5 CC','C1 VH','C2 VH','C3 VH','C4 VH','C5 VH')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C5 CC')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C4 CC')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C3 CC')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C2 CC')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C1 CC')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C5 VH')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C4 VH')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C3 VH')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C2 VH')
clin_merge$cluster2<-relevel(clin_merge$cluster2,'C1 VH')
clin_merge$cluster<-relevel(clin_merge$cluster,'C5')
clin_merge$cluster<-relevel(clin_merge$cluster,'C4')
clin_merge$cluster<-relevel(clin_merge$cluster,'C3')
clin_merge$cluster<-relevel(clin_merge$cluster,'C2')
clin_merge$cluster<-relevel(clin_merge$cluster,'C1')
dim(clin_merge)

#Combine PS1,2 and 3 proteins in a single dataset
names_ps1<-colnames(ps1_ptn_join[,-ncol(ps1_ptn_join)])
names_ps2<-colnames(ps2_ptn_join[,-ncol(ps2_ptn_join)])
names_ps3<-colnames(ps3_ptn_join[,-ncol(ps3_ptn_join)])
ps_names<-process_region_data(Venn(list(names_ps1,names_ps2,names_ps3)))$item
ptn_names<-c(ps_names[[1]],ps_names[[2]],ps_names[[3]],ps_names[[4]],ps_names[[5]],ps_names[[6]],ps_names[[7]])
vh_merge<-vh[vh$hma=='yes'& !is.na(vh$hma) & vh$venetoclax=='yes' & vh$arac=='no',c(-2:-ncol(clin_data))]
cc_merge<-cc[cc$hma=='no' & !is.na(cc$hma) & cc$venetoclax=='no' & cc$arac=='yes'|
               cc$hma=='yes' & cc$venetoclax=='yes' & !is.na(cc$hma) & cc$arac=='yes',c(-2:-ncol(clin_data))]
ptn_merge<-rbind(vh_merge,cc_merge)
merge_cl<-clin_merge[,c('id','cluster')]
ptn_merge<-merge(ptn_merge,merge_cl,by='id')
dim(ptn_merge)
ptn_merge<-ptn_merge[,c('id',ptn_names,'cluster')]
rownames(ptn_merge)<-ptn_merge$id
ptn_merge$id<-as.numeric(ptn_merge$id)
ptn_merge<-ptn_merge[,-1]
dim(ptn_merge)

####Apply PS1,2 and 3 clusters to the updated dataset----
#Create objects with cluster and treatment identification by PS set
ps1_clusters<-as.data.frame(ps1_clin_join[,c('id','cluster','trt','cluster2')])
rownames(ps1_clusters)<-NULL
ps2_clusters<-as.data.frame(ps2_clin_join[,c('id','cluster','trt','cluster2')])
rownames(ps2_clusters)<-NULL
ps3_clusters<-as.data.frame(ps3_clin_join[,c('id','cluster','trt','cluster2')])
rownames(ps3_clusters)<-NULL
ps_clusters<-as.data.frame(clin_merge[,c('id','cluster','trt','cluster2')])
rownames(ps_clusters)<-NULL

#Import updated dataset
setwd('/data') #customize folder for importing patient data file
rawtable2<-read_excel("dataset_updated.xlsx")

#Create and merge dataframe with clinical data only with cluster+treatment data
clin_data2<-as.data.frame(rawtable2[,c(1:83)])
str(clin_data2)
selected<-colnames(clin_data2[,c(3:15,25:ncol(clin_data2))])
clin_data2[selected]<-lapply(clin_data2[selected],factor)
str(clin_data2)
ps1_clin_join2<-merge(clin_data2,ps1_clusters,by='id')
ps2_clin_join2<-merge(clin_data2,ps2_clusters,by='id')
ps3_clin_join2<-merge(clin_data2,ps3_clusters,by='id')
clin_merge2<-merge(clin_data2,ps_clusters,by='id')

#Create and merge dataframe with protein data only with cluster+treatment data
rppa_data2<-as.data.frame(rawtable2[c(-2:-83)])
row.names(rppa_data2)<-rppa_data2$id
ps1_ptn_join2<-merge(rppa_data2,ps1_clusters,by='id')
ps1_ptn_join2$id<-as.integer(ps1_ptn_join2$id)
rownames(ps1_ptn_join2)<-ps1_ptn_join2$id 
ps2_ptn_join2<-merge(rppa_data2,ps2_clusters,by='id')
ps2_ptn_join2$id<-as.integer(ps2_ptn_join2$id)
rownames(ps2_ptn_join2)<-ps2_ptn_join2$id 
ps3_ptn_join2<-merge(rppa_data2,ps3_clusters,by='id')
ps3_ptn_join2$id<-as.integer(ps3_ptn_join2$id)
rownames(ps3_ptn_join2)<-ps3_ptn_join2$id 

merge_ptn<-merge(rppa_data2,ps_clusters,by='id')
merge_ptn$id<-as.integer(merge_ptn$id)
rownames(merge_ptn)<-merge_ptn$id 

####Generate Figures and Tables for manuscript----

####Supplementary Tables S2 and S3: Materials and Methods data from Protein Selector Set###
#Create dataframe with p-values for all proteins and for selected proteins in PS1, 2 and 3
pval_all_ps1<- tibble::rownames_to_column(as.data.frame(pval_ps1), "Protein")
colnames(pval_all_ps1)<-c('Protein Name','Median Split PS1','Tertiles PS1','Quartiles PS1','Quintiles PS1','Sextiles PS1')
pval_all_ps2<- tibble::rownames_to_column(as.data.frame(pval_ps2), "Protein")
colnames(pval_all_ps2)<-c('Protein Name','Median Split PS2','Tertiles PS2','Quartiles PS2','Quintiles PS2','Sextiles PS2')
pval_all_ps3<- tibble::rownames_to_column(as.data.frame(pval_ps3), "Protein")
colnames(pval_all_ps3)<-c('Protein Name','Median Split PS3','Tertiles PS3','Quartiles PS3','Quintiles PS3','Sextiles PS3')
pval_all_table<-merge(pval_all_ps1,pval_all_ps2,by='Protein Name')
pval_all_table<-merge(pval_all_table,pval_all_ps3,by='Protein Name')

#Create dataframe with pvalues from selected PS proteins
ps1_pval_tab<-tibble::rownames_to_column(pval_ps1, "protein")
ps1_pval_tab<-as.data.frame(ps1_pval_tab[ps1_pval_tab$surv_Q2<0.05,c(1:2)])
colnames(ps1_pval_tab)<-c('PS1 proteins','p_value')
ps1_pval_tab <- ps1_pval_tab[order(ps1_pval_tab$p_value),]
rownames(ps1_pval_tab)<-NULL

ps2_pval_tab<-tibble::rownames_to_column(pval_ps2, "protein")
ps2_pval_tab<-as.data.frame(ps2_pval_tab[ps2_pval_tab$surv_Q3<0.01,c(1,3)])
colnames(ps2_pval_tab)<-c('PS2 proteins','p_value')
ps2_pval_tab <- ps2_pval_tab[order(ps2_pval_tab$p_value),]
rownames(ps2_pval_tab)<-NULL

ps3_pval_tab<-tibble::rownames_to_column(pval_ps3, "protein")
ps3_pval_tab<-as.data.frame(ps3_pval_tab[ps3_pval_tab$surv_Q5<0.01,c(1,3)])
colnames(ps3_pval_tab)<-c('PS3 proteins','p_value')
ps3_pval_tab <- ps3_pval_tab[order(ps3_pval_tab$p_value),]
rownames(ps3_pval_tab)<-NULL

#Export data to excel files
sup_s2<-createWorkbook()
sh_sup_s2<-createSheet(sup_s2," ")
cell_style<-CellStyle(sup_s2)+Font(sup_s2,isBold=TRUE)+Border()
addDataFrame(data.frame("Supplementary Table S2. P-values of comparisons between quantiles for each 411 proteins according to the quantile group"=double(),check.names=FALSE),sheet=sh_sup_s2,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style)
addDataFrame(pval_all_table,sheet=sh_sup_s2,startColumn=1,startRow=2,row.names=FALSE,colnamesStyle=cell_style)
setwd('/results')#custom destination folder for file generation
saveWorkbook(sup_s2,"Supplementary_table_S2.xlsx")

sup_s3<-createWorkbook()
sh_sup_s3<-createSheet(sup_s3," ")
cell_style<-CellStyle(sup_s3)+Font(sup_s3,isBold=TRUE)+Border()
addDataFrame(data.frame("Supplementary Table S3. P-values of selected proteins according to each protein selector set"=double(),check.names=FALSE),sheet=sh_sup_s3,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style)
addDataFrame(ps1_pval_tab,sheet=sh_sup_s3,startColumn=1,startRow=2,row.names=FALSE,colnamesStyle=cell_style)
addDataFrame(ps2_pval_tab,sheet=sh_sup_s3,startColumn=3,startRow=2,row.names=FALSE,colnamesStyle=cell_style)
addDataFrame(ps3_pval_tab,sheet=sh_sup_s3,startColumn=5,startRow=2,row.names=FALSE,colnamesStyle=cell_style)
setwd('/results')#custom destination folder for file generation
saveWorkbook(sup_s3,"Supplementary_table_S3.xlsx")

####Figure 1: PS1 analyis###
#PS1 heatmap (Figure 1A)#
ps1_join_cl<-ps1_ptn_join2[,c('CD4','CDKN1B','CDKN1B.pS10','ETS1','LCK','UGT1A1','SPARC','TNFRSF4'
                             ,'SGK1','SOX17','DLST','SOX2','EZR','JMJD6','CASP3.cle','BIRC2'
                             ,'PPARA','TYRO3','NDUFB4','ASS1','RHEB','AKT3','AKR1C3','HDAC3','SCD'
                             ,'WIPI2','RAF1.pS338','MKNK1','STAT3.pY705','PKM2'
                             ,'MAPK14.pT180_Y182','EZH2','WEE1','CDK1_2_3.pT14','KIT','EIF2AK2'
                             ,'ASH2L','SMARCB1','MEN1','BRD4','CDKN1B.pT198','ARID1A','HSF1.pS326'
                             ,'EIF4E.pS209','EIF4EBP1.pS65','MAPK14','TSC2','CHEK2','CASP3'
                             ,'EP300','HSF1','LYN','STAT3.pS727','NOTCH3','CD74','cluster')]

#Merge clusters dataframe with clinical data for annotation
ps1_join_cl<- tibble::rownames_to_column(ps1_join_cl, "id")
ps1_join_clin_subset<-ps1_clin_join2[,c('id','trt','aml_gp','cyto_risk','complex_kar','mut_tp53','mut_dnmt3'
                                       ,'mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras','mut_ptpn11')]

ps1_join_cl<-merge(ps1_join_cl,ps1_join_clin_subset,by='id')
ps1_join_cl$cyto_risk<-factor(ps1_join_cl$cyto_risk,levels=c('favorable','intermediary','unfavorable',NA)
                              ,labels = c('favorable','intermediary','unfavorable','N/A'), exclude = NULL)
selected_var<-c('complex_kar','mut_tp53','mut_dnmt3','mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras','mut_ptpn11')
ps1_join_cl[selected_var]<-lapply(ps1_join_cl[selected_var], function(x) factor(x,levels=c('no','yes',NA)
                                                                                ,labels = c('no','yes','N/A'), exclude = NULL))

#Reorder clustering object
ps1_join_cl<-ps1_join_cl[order(ps1_join_cl[,'cluster']),]

#Create annotation object
ps1_join_ann<-as.data.frame(ps1_join_cl[,c('id','cluster','trt','aml_gp','cyto_risk','complex_kar','mut_tp53','mut_dnmt3'
                                           ,'mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras'
                                           ,'mut_ptpn11')])

row.names(ps1_join_ann)<-c(1:nrow(ps1_join_ann))
ps1_join_ann<-as.data.frame(ps1_join_ann[,-1])
ps1_join_ann<-ps1_join_ann %>% select('cluster','trt',everything())
colnames(ps1_join_ann)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk','Complex Karyotype','TP53 Mutation','DNMT3 Mutation'
                          ,'IDH Mutation','SRSF2 Mutation','FLT3 Mutation','ASXL1 Mutation','RAS Mutation'
                          ,'PTPN11 Mutation')

#Clean clustering object
row.names(ps1_join_cl)<-c(1:nrow(ps1_join_cl))
ps1_join_cl<-as.data.frame(ps1_join_cl[,c(c(2:(nrow(vh_ps1_list[[6]])+1)))])
ps1_join_cl<-as.data.frame(t(ps1_join_cl))
ps1_join_cl<-mutate_all(ps1_join_cl, function(x) as.numeric(as.character(x)))

#Create color object
ps1_join_cl_colors<-list()
ps1_join_cl_colors[[1]]<-mycolors1[1:nlevels(ps1_join_ann$'Cluster')]
names(ps1_join_cl_colors[[1]])<-unique(ps1_join_ann$'Cluster')
ps1_join_cl_colors[[2]]<-c(pal_jco()(2))
names(ps1_join_cl_colors[[2]])<-unique(ps1_join_ann$'Treatment')
ps1_join_cl_colors[[3]]<-c(pal_jco()(3))
names(ps1_join_cl_colors[[3]])<-unique(ps1_join_ann$'AML Group')
ps1_join_cl_colors[[4]]<-c(mycolors1[1],pal_jco()(2)[2:1],pal_jco()(3)[3])
names(ps1_join_cl_colors[[4]])<-c('favorable','intermediary','unfavorable','N/A')
for (i in 5:ncol(ps1_join_ann)){
  ps1_join_cl_colors[[i]]<-c(pal_jco()(3))
  names(ps1_join_cl_colors[[i]])<-c('yes','no','N/A')
}
names(ps1_join_cl_colors)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk','Complex Karyotype'
                             ,'TP53 Mutation','DNMT3 Mutation','IDH Mutation','SRSF2 Mutation'
                             ,'FLT3 Mutation','ASXL1 Mutation','RAS Mutation','PTPN11 Mutation')

#Normalize sample signal
if(min(ps1_join_cl)>=-2 & max(ps1_join_cl)<=2){
  ps1_join_break=seq(-2.4,2.4,0.4)
}else if(min(ps1_join_cl)>=-2){
  ps1_join_break=c(seq(-2.4,2,0.4),max(ps1_join_cl))
}else if(max(ps1_join_cl)<=2){
  ps1_join_break=c(min(ps1_join_cl),seq(-2,2.4,0.4))
}else{
  ps1_join_break=c(min(ps1_join_cl),seq(-2,2,0.4),max(ps1_join_cl))
}

ps1_join_ht<-pheatmap(ps1_join_cl,annotation_colors=ps1_join_cl_colors
                      ,main="Protein Selector 1 (PS1) Patients (N=419)"
                      ,annotation_col=ps1_join_ann,cluster_rows=F,cluster_cols=F
                      ,clustering_method ='ward.D2',clustering_distance_rows=d_row_ps1_join
                      ,fontsize=12,col=jet.colors(length(ps1_join_break)-1),breaks=ps1_join_break
                      ,scale="none",show_colnames=F,fontsize_row=12)

setwd('/results')#custom destination folder for file generation
ggsave('Figure_1A.pdf',ps1_join_ht, height=8.27,width=11.69,units='in')
#Legends were adjusted using adobe illustrator

##Survival ps1 (figure 1B)##
sv_fit_ps1_join2<-surv_fit(Surv(surv_time,status)~trt+cluster,data=ps1_clin_join2)
sv_plot_ps1_join2<-ggsurvplot(sv_fit_ps1_join2,data=ps1_clin_join2,legend.labs=c('Overall','C1 CC','C2 CC','C3 CC','C1 VH','C2 VH','C3 VH')
                              ,pal=c('black',rep(mycolors1[1:nlevels(ps1_clin_join2$cluster)],2))
                              ,legend.title=" ",linetype=c('solid',rep('twodash',3),rep('solid',3))
                              ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Overall Survival PS1',xlab='Time (years)'
                              ,font.main = 16,font.x = 14,font.y = 14,ylab='Cumulative Survival probability'#,risk.table.height = 0.25,tables.col = "strata"
                              ,risk.table=F,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                              ,size=0.5,fontsize=3,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)
sv_plot_ps1_join2$plot<-sv_plot_ps1_join2$plot+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

pval_table_ps1_sv<-tableGrob(symnum(pairwise_survdiff(Surv(surv_time,status)~cluster2,data=ps1_clin_join2)$p.value
                                    ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                    ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                      core=list(fg_params=list(cex=1,fontface=2)
                                                ,bg_params=list(fill='white',col='black'))
                                      ,colhead=list(fg_params=list(col='white',fontface=2)
                                                    ,bg_params=list(fill=c(rep(mycolors1[1:3],2)),col=NA))
                                      ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                    ,bg_params = list(fill=c('white',c(mycolors1[2:3],mycolors1[1:3])),col=NA)))
                             ,widths=unit(c(rep(1,5)),c(rep('cm',5))))

sv_ps1_new<-ggarrange(sv_plot_ps1_join2$plot,pval_table_ps1_sv,ncol=1,nrow=2,heights=c(0.7,0.3))

setwd('/results')#custom destination folder for file generation
ggsave('Figure_1B.pdf',sv_ps1_new, height=11.69,width=8.27,units='in')

#PS2 heatmap (Figure 1C)#
ps2_join_cl<-ps2_ptn_join2[,c('JAG1','G6PD','NOTCH1.cle','ADM'
                             ,'STK4','ITGAL','CBX7','BABAM1.pS29','SGK3','GSKA_B','SPI1'
                             ,'PAX6','NOTCH1','EIF4G1','HSP90AA1_B1','SMAD2','FASN'
                             ,'FOXM1','H3K27Me3','ASNS'
                             ,'GATA1','WEE1.pS642','BCL2','PDCD4','LMNB1'
                             ,'FOS','MECOM','CDH1','BMI1','H3K27Ac','NLN','IRS1','SUZ12'
                             ,'IGFBP2','cluster')]

#Merge clusters dataframe with clinical data for annotation
ps2_join_cl<- tibble::rownames_to_column(ps2_join_cl, "id")
ps2_join_clin_subset<-ps2_clin_join2[,c('id','trt','aml_gp','cyto_risk','complex_kar','mut_tp53','mut_dnmt3'
                                       ,'mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras','mut_ptpn11')]

ps2_join_cl<-merge(ps2_join_cl,ps2_join_clin_subset,by='id')
ps2_join_cl$cyto_risk<-factor(ps2_join_cl$cyto_risk,levels=c('favorable','intermediary','unfavorable',NA)
                              ,labels = c('favorable','intermediary','unfavorable','N/A'), exclude = NULL)
selected_var<-c('complex_kar','mut_tp53','mut_dnmt3','mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras','mut_ptpn11')
ps2_join_cl[selected_var]<-lapply(ps2_join_cl[selected_var], function(x) factor(x,levels=c('no','yes',NA)
                                                                                ,labels = c('no','yes','N/A'), exclude = NULL))

#Reorder clustering object
ps2_join_cl<-ps2_join_cl[order(ps2_join_cl[,'cluster']),]

#Create annotation object
ps2_join_ann<-as.data.frame(ps2_join_cl[,c('id','cluster','trt','aml_gp','cyto_risk','complex_kar','mut_tp53','mut_dnmt3'
                                           ,'mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras'
                                           ,'mut_ptpn11')])
row.names(ps2_join_ann)<-c(1:nrow(ps2_join_ann))
ps2_join_ann<-as.data.frame(ps2_join_ann[,-1])
ps2_join_ann<-ps2_join_ann %>% select('cluster','trt',everything())
colnames(ps2_join_ann)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk','Complex Karyotype','TP53 Mutation','DNMT3 Mutation'
                          ,'IDH Mutation','SRSF2 Mutation','FLT3 Mutation','ASXL1 Mutation','RAS Mutation'
                          ,'PTPN11 Mutation')

#Clean clustering object
row.names(ps2_join_cl)<-c(1:nrow(ps2_join_cl))
ps2_join_cl<-as.data.frame(ps2_join_cl[,c(c(2:(nrow(cc_c2_list[[2]])+1)))])
ps2_join_cl<-as.data.frame(t(ps2_join_cl))
ps2_join_cl<-mutate_all(ps2_join_cl, function(x) as.numeric(as.character(x)))

#Create color object
ps2_join_cl_colors<-list()
ps2_join_cl_colors[[1]]<-mycolors1[4:5]
names(ps2_join_cl_colors[[1]])<-unique(ps2_join_ann$'Cluster')
ps2_join_cl_colors[[2]]<-c(pal_jco()(2)[2:1])
names(ps2_join_cl_colors[[2]])<-unique(ps2_join_ann$'Treatment')
ps2_join_cl_colors[[3]]<-c(pal_jco()(3))
names(ps2_join_cl_colors[[3]])<-unique(ps2_join_ann$'AML Group')
ps2_join_cl_colors[[4]]<-c(mycolors1[1],pal_jco()(2)[2:1],pal_jco()(3)[3])
names(ps2_join_cl_colors[[4]])<-c('favorable','intermediary','unfavorable','N/A')
for (i in 5:ncol(ps2_join_ann)){
  ps2_join_cl_colors[[i]]<-c(pal_jco()(3))
  names(ps2_join_cl_colors[[i]])<-c('yes','no','N/A')
}
names(ps2_join_cl_colors)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk','Complex Karyotype'
                             ,'TP53 Mutation','DNMT3 Mutation','IDH Mutation','SRSF2 Mutation'
                             ,'FLT3 Mutation','ASXL1 Mutation','RAS Mutation','PTPN11 Mutation')

#Normalize sample signal
if(min(ps2_join_cl)>=-2 & max(ps2_join_cl)<=2){
  ps2_join_break=seq(-2.4,2.4,0.4)
}else if(min(ps2_join_cl)>=-2){
  ps2_join_break=c(seq(-2.4,2,0.4),max(ps2_join_cl))
}else if(max(ps2_join_cl)<=2){
  ps2_join_break=c(min(ps2_join_cl),seq(-2,2.4,0.4))
}else{
  ps2_join_break=c(min(ps2_join_cl),seq(-2,2,0.4),max(ps2_join_cl))
}

ps2_join_ht<-pheatmap(ps2_join_cl,annotation_colors=ps2_join_cl_colors
                      ,main="Protein Selector 2 (PS2) C2 Patients (N=182)"
                      ,annotation_col=ps2_join_ann,cluster_rows=F,cluster_cols=F
                      ,clustering_method ='ward.D2',clustering_distance_rows=d_row_ps2_join
                      ,fontsize=12,col=jet.colors(length(ps2_join_break)-1),breaks=ps2_join_break
                      ,scale="none",show_colnames=F,fontsize_row=12)

setwd('/results')#custom destination folder for file generation
ggsave('Figure_1C.pdf',ps2_join_ht, height=8.27,width=11.69,units='in')
#Legends were adjusted using adobe illustrator

#PS3 heatmap (Figure 1D)#
#Remove treatment order
join_ps3_cl<-ps3_ptn_join2[,c('GAB2','H3K27Me3','HK2','TUBA1A'
                             ,'EIF4E','BRAF.pS445','SMARCA2','MTOR','EIF4G1','ARID1A','ADM'
                             ,'SPI1','PDL1','PIM2','PEA15','MAP1LC3A_B','STAT1','INPP4B'
                             ,'XPF','HSPB1.pS82','WEE1.pS642','EIF2AK2'
                             ,'TP53BP1','PTEN','LMNB1','EEF2K','AURORA_A_B_C.pT288_232_198'
                             ,'HSF1.pS326','RAB11','cluster')]

#Merge clusters dataframe with clinical data for annotation
join_ps3_cl<- tibble::rownames_to_column(join_ps3_cl, "id")
join_ps3_clin_subset<-ps3_clin_join2[,c('id','trt','aml_gp','cyto_risk','complex_kar','mut_tp53','mut_dnmt3'
                                       ,'mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras','mut_ptpn11')]

join_ps3_cl<-merge(join_ps3_cl,join_ps3_clin_subset,by='id')
join_ps3_cl$cyto_risk<-factor(join_ps3_cl$cyto_risk,levels=c('favorable','intermediary','unfavorable',NA)
                              ,labels = c('favorable','intermediary','unfavorable','N/A'), exclude = NULL)
selected_var<-c('complex_kar','mut_tp53','mut_dnmt3','mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras','mut_ptpn11')
join_ps3_cl[selected_var]<-lapply(join_ps3_cl[selected_var], function(x) factor(x,levels=c('no','yes',NA)
                                                                                ,labels = c('no','yes','N/A'), exclude = NULL))

#Reorder clustering object
join_ps3_cl<-join_ps3_cl[order(join_ps3_cl[,'cluster']),]

#Create annotation object
join_ps3_ann<-as.data.frame(join_ps3_cl[,c('id','cluster','trt','aml_gp','cyto_risk','complex_kar','mut_tp53','mut_dnmt3'
                                           ,'mut_idh','mut_srsf2','mut_flt3','mut_aslx1','mut_ras'
                                           ,'mut_ptpn11')])
row.names(join_ps3_ann)<-c(1:nrow(join_ps3_ann))
join_ps3_ann<-as.data.frame(join_ps3_ann[,-1])
join_ps3_ann<-join_ps3_ann %>% select('cluster','trt',everything())
colnames(join_ps3_ann)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk','Complex Karyotype','TP53 Mutation','DNMT3 Mutation'
                          ,'IDH Mutation','SRSF2 Mutation','FLT3 Mutation','ASXL1 Mutation','RAS Mutation'
                          ,'PTPN11 Mutation')

#Clean clustering object
row.names(join_ps3_cl)<-c(1:nrow(join_ps3_cl))
join_ps3_cl<-as.data.frame(join_ps3_cl[,c(c(2:(nrow(vh_c3_list[[4]])+1)))])
join_ps3_cl<-as.data.frame(t(join_ps3_cl))
join_ps3_cl<-mutate_all(join_ps3_cl, function(x) as.numeric(as.character(x)))

#Create color object
join_ps3_cl_colors<-list()
join_ps3_cl_colors[[1]]<-mycolors1[6:7]
names(join_ps3_cl_colors[[1]])<-unique(join_ps3_ann$'Cluster')
join_ps3_cl_colors[[2]]<-c(pal_jco()(2)[1:2])
names(join_ps3_cl_colors[[2]])<-unique(join_ps3_ann$'Treatment')
join_ps3_cl_colors[[3]]<-c(pal_jco()(3))
names(join_ps3_cl_colors[[3]])<-unique(join_ps3_ann$'AML Group')
join_ps3_cl_colors[[4]]<-c(mycolors1[1],pal_jco()(2)[2:1],pal_jco()(3)[3])
names(join_ps3_cl_colors[[4]])<-c('favorable','intermediary','unfavorable','N/A')
for (i in 5:ncol(join_ps3_ann)){
  join_ps3_cl_colors[[i]]<-c(pal_jco()(3))
  names(join_ps3_cl_colors[[i]])<-c('yes','no','N/A')
}
names(join_ps3_cl_colors)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk','Complex Karyotype'
                             ,'TP53 Mutation','DNMT3 Mutation','IDH Mutation','SRSF2 Mutation'
                             ,'FLT3 Mutation','ASXL1 Mutation','RAS Mutation','PTPN11 Mutation')

#Normalize sample signal
if(min(join_ps3_cl)>=-2 & max(join_ps3_cl)<=2){
  join_ps3_break=seq(-2.4,2.4,0.4)
}else if(min(join_ps3_cl)>=-2){
  join_ps3_break=c(seq(-2.4,2,0.4),max(join_ps3_cl))
}else if(max(join_ps3_cl)<=2){
  join_ps3_break=c(min(join_ps3_cl),seq(-2,2.4,0.4))
}else{
  join_ps3_break=c(min(join_ps3_cl),seq(-2,2,0.4),max(join_ps3_cl))
}

join_ps3_ht<-pheatmap(join_ps3_cl,annotation_colors=join_ps3_cl_colors
                      ,main="Protein Selector 3 (PS3) Patients (N=146)"
                      ,annotation_col=join_ps3_ann,cluster_rows=F,cluster_cols=F
                      ,clustering_method ='ward.D2',clustering_distance_rows=d_row_join_ps3
                      ,fontsize=12,col=jet.colors(length(join_ps3_break)-1),breaks=join_ps3_break
                      ,scale="none",show_colnames=F,fontsize_row=12)

setwd('/results')#custom destination folder for file generation
ggsave('Figure_1D.pdf',join_ps3_ht, height=8.27,width=11.69,units='in')
#Legends were adjusted using adobe illustrator

##Survival PS2 (figure 1E)##
sv_fit_ps2_join2<-surv_fit(Surv(surv_time,status)~trt+cluster,data=ps2_clin_join2)
sv_plot_ps2_join2<-ggsurvplot(sv_fit_ps2_join2,data=ps2_clin_join2,legend.labs=c('Overall','C1 CC','C2 CC','C1 VH','C2 VH')
                              ,pal=c('black',rep(mycolors1[4:5],2))
                              ,legend.title=" ",linetype=c('solid',rep('twodash',2),rep('solid',2))
                              ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Overall Survival PS2',xlab='Time (years)'
                              ,font.main = 16,font.x = 14,font.y = 14,ylab='Cumulative Survival probability'#,risk.table.height = 0.25,tables.col = "strata"
                              ,risk.table=F,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                              ,size=0.5,fontsize=3,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)
sv_plot_ps2_join2$plot<-sv_plot_ps2_join2$plot+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

pval_table_ps2_sv<-tableGrob(symnum(pairwise_survdiff(Surv(surv_time,status)~cluster2,data=ps2_clin_join2)$p.value
                                    ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                    ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                      core=list(fg_params=list(cex=1,fontface=2)
                                                ,bg_params=list(fill='white',col='black'))
                                      ,colhead=list(fg_params=list(col='white',fontface=2)
                                                    ,bg_params=list(fill=c(rep(mycolors1[4:5],2)),col=NA))
                                      ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                    ,bg_params = list(fill=c('white',c(mycolors1[5],mycolors1[4:5])),col=NA)))
                             ,widths=unit(c(rep(2.8,3)),c(rep('cm',3))))

sv_ps2_new<-ggarrange(sv_plot_ps2_join2$plot,pval_table_ps2_sv,ncol=1,nrow=2,heights=c(0.7,0.3))

setwd('/results')#custom destination folder for file generation
ggsave('Figure_1E.pdf',sv_ps2_new, height=11.69,width=8.27,units='in')

##Survival PS3 (figure 1F)##
sv_fit_ps3_join2<-surv_fit(Surv(surv_time,status)~trt+cluster,data=ps3_clin_join2)
sv_plot_ps3_join2<-ggsurvplot(sv_fit_ps3_join2,data=ps3_clin_join2,legend.labs=c('Overall','C1 CC','C2 CC','C1 VH','C2 VH')
                              ,pal=c('black',rep(mycolors1[6:7],2))
                              ,legend.title=" ",linetype=c('solid',rep('twodash',2),rep('solid',2))
                              ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Overall Survival PS3',xlab='Time (years)'
                              ,font.main = 16,font.x = 14,font.y = 14,ylab='Cumulative Survival probability'#,risk.table.height = 0.25,tables.col = "strata"
                              ,risk.table=F,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                              ,size=0.5,fontsize=3,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)
sv_plot_ps3_join2$plot<-sv_plot_ps3_join2$plot+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

pval_table_ps3_sv<-tableGrob(symnum(pairwise_survdiff(Surv(surv_time,status)~cluster2,data=ps3_clin_join2)$p.value
                                    ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                    ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                      core=list(fg_params=list(cex=1,fontface=2)
                                                ,bg_params=list(fill='white',col='black'))
                                      ,colhead=list(fg_params=list(col='white',fontface=2)
                                                    ,bg_params=list(fill=c(rep(mycolors1[6:7],2)),col=NA))
                                      ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                    ,bg_params = list(fill=c('white',mycolors1[7],mycolors1[6:7]),col=NA)))
                             ,widths=unit(c(rep(2.8,3)),c(rep('cm',3))))

sv_ps3_new<-ggarrange(sv_plot_ps3_join2$plot,pval_table_ps3_sv,ncol=1,nrow=2,heights=c(0.7,0.3))

setwd('/results')#custom destination folder for file generation
ggsave('Figure_1F.pdf',sv_ps3_new, height=11.69,width=8.27,units='in')

####Figure 2: Combined analysis of PS1, 2 and 3 ###
#Heatmap of PS1, 2 and 3 combined (Figure 2A)#
#Set order of protein list
merge_cl<-merge_ptn[,c('PKM2','NOTCH1.cle','CD4','CDKN1B','CDKN1B.pS10','ETS1','LCK','ADM'
                        ,'UGT1A1','JAG1','CBX7','STAT1','PEA15','STK4','EZR','JMJD6'
                        ,'CASP3.cle','PAX6','INPP4B','G6PD','SPARC','TNFRSF4'
                        ,'SGK1','SOX17','DLST','SOX2','PDL1','PIM2','NOTCH1','BIRC2','PPARA','TYRO3','NDUFB4','ASS1','RHEB'
                        ,'AKT3','AKR1C3','HDAC3','SCD','WIPI2','EIF4G1','HSP90AA1_B1','FASN','FOXM1'
                        ,'H3K27Me3','ASNS','RAF1.pS338','MKNK1','STAT3.pY705','GAB2','TUBA1A','EIF4E'
                        ,'ITGAL','GSKA_B','MAP1LC3A_B','SPI1','MAPK14.pT180_Y182','HSPB1.pS82'
                        ,'XPF','EZH2','GATA1','WEE1','TP53BP1','CDK1_2_3.pT14','WEE1.pS642','EIF2AK2'
                        ,'PTEN','EEF2K','PDCD4','LMNB1','KIT','BCL2','ASH2L','SMARCB1','MEN1','BRD4'
                        ,'CDKN1B.pT198','NLN','HSF1.pS326','CD74','ARID1A','SMARCA2','H3K27Ac'
                        ,'EIF4E.pS209','EIF4EBP1.pS65','MECOM','FOS','CDH1','SMAD2','SGK3','HK2'
                        ,'BABAM1.pS29','MAPK14','TSC2','CHEK2','CASP3','EP300','HSF1','LYN'
                        ,'STAT3.pS727','NOTCH3','BMI1','IRS1','SUZ12','IGFBP2','MTOR'
                        ,'AURORA_A_B_C.pT288_232_198','RAB11','BRAF.pS445'
                        ,'cluster')]

#Merge clusters dataframe with clinical data for annotation
merge_cl<- tibble::rownames_to_column(merge_cl, "id")
clin_merge_subset<-clin_merge2[,c('id','trt','aml_gp','cyto_risk','mut_cebpa','mut_gata2','mut_flt3','mut_npm1')]

merge_cl<-merge(merge_cl,clin_merge_subset,by='id')
merge_cl$cyto_risk<-factor(merge_cl$cyto_risk,levels=c('favorable','intermediary','unfavorable',NA)
                           ,labels = c('favorable','intermediary','unfavorable','N/A'), exclude = NULL)
selected_var<-c('mut_cebpa','mut_gata2','mut_flt3','mut_npm1')
merge_cl[selected_var]<-lapply(merge_cl[selected_var], function(x) factor(x,levels=c('no','yes',NA)
                                                                          ,labels = c('no','yes','N/A'), exclude = NULL))

#Reorder clustering object
merge_cl<-merge_cl[order(merge_cl[,'cluster']),]

#Create column annotation object
merge_ann<-as.data.frame(merge_cl[,c('id','cluster','trt','aml_gp','cyto_risk','mut_cebpa','mut_gata2','mut_flt3','mut_npm1')])
row.names(merge_ann)<-c(1:nrow(merge_ann))
merge_ann<-as.data.frame(merge_ann[,-1])
colnames(merge_ann)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk'
                       ,'CEBPA Mutation','GATA2 Mutation','FLT3 Mutation','NPM1 Mutation')

#Clean clustering object
row.names(merge_cl)<-c(1:nrow(merge_cl))
merge_cl<-as.data.frame(merge_cl[,2:(length(ptn_names)+1)])
merge_cl<-as.data.frame(t(merge_cl))
merge_cl<-mutate_all(merge_cl, function(x) as.numeric(as.character(x)))

#Create row annotation object
merge_anno_row<-as.data.frame(ptn_names)
merge_anno_row$ps<-factor(c(rep('PS1',52),rep('PS2',28),rep('PS3',20)
                            ,rep('PS1 and PS2',3),rep('PS2 and PS3',6))
                          ,levels=c('PS1','PS2','PS3','PS1 and PS2','PS2 and PS3'))
colnames(merge_anno_row)<-c('names','ps')
merge_anno_row<-merge_anno_row[order(match(merge_anno_row$names,rownames(merge_cl))),]
merge_anno_row<-as.data.frame(merge_anno_row[,-1])
colnames(merge_anno_row)<-c('Protein Selector')
rownames(merge_anno_row)<-rownames(merge_cl)

#Create color object
merge_cl_colors<-list()
merge_cl_colors[[1]]<-mycolors1[c(1,4:7)]
names(merge_cl_colors[[1]])<-unique(merge_ann$'Cluster')
merge_cl_colors[[2]]<-c(pal_jco()(2))
names(merge_cl_colors[[2]])<-unique(merge_ann$'Treatment')
merge_cl_colors[[3]]<-c(pal_jco()(3))
names(merge_cl_colors[[3]])<-unique(merge_ann$'AML Group')
merge_cl_colors[[4]]<-c(mycolors1[1],pal_jco()(2)[2:1],pal_jco()(3)[3])
names(merge_cl_colors[[4]])<-c('favorable','intermediary','unfavorable','N/A')
for (i in 5:ncol(merge_ann)){
  merge_cl_colors[[i]]<-c(pal_jco()(3))
  names(merge_cl_colors[[i]])<-c('yes','no','N/A')
}
merge_cl_colors[[9]]<-c(pal_jco()(2),mycolors1[c(1,5,6)])
names(merge_cl_colors[[9]])<-c('PS1','PS2','PS3','PS1 and PS2','PS2 and PS3')

names(merge_cl_colors)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk'
                          ,'CEBPA Mutation','GATA2 Mutation','FLT3 Mutation','NPM1 Mutation'
                          ,'Protein Selector')

#Normalize sample signal
if(min(merge_cl)>=-2 & max(merge_cl)<=2){
  merge_break=seq(-2.4,2.4,0.4)
}else if(min(merge_cl)>=-2){
  merge_break=c(seq(-2.4,2,0.4),max(merge_cl))
}else if(max(merge_cl)<=2){
  merge_break=c(min(merge_cl),seq(-2,2.4,0.4))
}else{
  merge_break=c(min(merge_cl),seq(-2,2,0.4),max(merge_cl))
}

merge_ht<-pheatmap(merge_cl,annotation_colors=merge_cl_colors
                   ,main="Protein Selector 1, 2 and 3 combined (N=419)"
                   ,annotation_row=merge_anno_row,cluster_rows=F,cluster_cols=F
                   ,clustering_method ='ward.D2',annotation_col=merge_ann,fontsize=11
                   ,col=jet.colors(length(merge_break)-1),breaks=merge_break,scale="none"
                   ,show_colnames=F,fontsize_row=10,border_color=F,annotation_names_row=F)

setwd('/results')#custom destination folder for file generation
ggsave('Figure_2A.pdf',merge_ht, height=8.27,width=11.69,units='in')
#Legends were adjusted using adobe illustrator

##Survival PS1,2 and 3 combined (figure 2B)##
sv_fit_merge2<-surv_fit(Surv(surv_time,status)~trt+cluster,data=clin_merge2)
sv_plot_merge2<-ggsurvplot(sv_fit_merge2,data=clin_merge2,legend.labs=c('Overall','C1 CC','C2 CC','C3 CC','C4 CC','C5 CC','C1 VH','C2 VH','C3 VH','C4 VH','C5 VH')
                           ,pal=c('black',rep(mycolors1[c(1,4:7)],2))
                           ,legend.title=" ",linetype=c('solid',rep('twodash',5),rep('solid',5))
                           ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Overall Survival PS1,2 and 3',xlab='Time (years)'
                           ,font.main = 16,font.x = 14,font.y = 14,ylab='Cumulative Survival probability'#,risk.table.height = 0.25,tables.col = "strata"
                           ,risk.table=F,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                           ,size=0.5,fontsize=3,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)
sv_plot_merge2$plot<-sv_plot_merge2$plot+
  theme(legend.text=element_text(size = 10),legend.key.size=unit(0.5, "cm")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

pval_table_merged_sv<-tableGrob(symnum(pairwise_survdiff(Surv(surv_time,status)~cluster2,data=clin_merge2)$p.value
                                       ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                       ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                         core=list(fg_params=list(cex=1,fontface=2)
                                                   ,bg_params=list(fill='white',col='black'))
                                         ,colhead=list(fg_params=list(col='white',fontface=2)
                                                       ,bg_params=list(fill=c(rep(mycolors1[c(1,4:7)],2)),col=NA))
                                         ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                       ,bg_params = list(fill=c('white',c(mycolors1[c(4:7)],mycolors1[c(1,4:7)])),col=NA)))
                                ,widths=unit(c(rep(0.001,9)),c(rep('cm',9))))

sv_merge3<-ggarrange(sv_plot_merge2$plot,pval_table_merged_sv,ncol=1,nrow=2,heights=c(0.7,0.3))

setwd('/results')#custom destination folder for file generation
ggsave('Figure_2B.pdf',sv_merge3, height=11.69,width=8.27,units='in')

##Top Correlations between PS proteins (Figure 2C)###
setwd('/results')#custom destination folder for file generation
pdf(file='Figure_2C.pdf',height=11.69,width=8.27)
cor_ps_full<-cor(merge_ptn[,c(ptn_names)],method='pearson') #use previously created objects to create matrix
diag(cor_ps_full)<-0 #set diagonal values to 0, intead of 1
sel_cor<-apply(abs(cor_ps_full)>=0.6,1,any) #select rows with top correlations (coeff>0.6) 
cor_ps<-cor_ps_full[sel_cor,sel_cor] #select rows and columns by name of top correlations
dim(cor_ps)
test_cor_ps<-cor.mtest(merge_ptn[,row.names(cor_ps)],conf.level = 0.95)
cor_plot_ps<-corrplot(cor_ps,order='hclust',hclust.method='ward.D2',method='color'
                      ,col=rev(COL2('RdYlBu')),tl.srt=45,cl.pos='b',diag=F
                      ,p.mat=test_cor_ps$p,tl.col='black',type='lower',mar=c(0,0,1,0)
                      ,sig.level=c(0.001,0.01,0.05),pch.cex=0.5,insig='label_sig',cl.ratio=0.1
                      ,tl.cex=0.85,pch.col='grey20',cl.cex=0.85,addgrid.col=NA 
                      ,title="Top correlation between PS1, 2 and 3 proteins")
dev.off()
dev.off()

####Supplementary Figure 1: Remission analysis for each Proteins Selector Set###
##Remission PS1 (Supplementary figure 1A)##
ps1_rem_join2<-ps1_clin_join2[!is.na(ps1_clin_join2$relapse),]
rem_fit_ps1_join2<-surv_fit(Surv(rem_time,relapse)~trt+cluster,data=ps1_rem_join2)
rem_plot_ps1_join2<-ggsurvplot(rem_fit_ps1_join2,data=ps1_rem_join2,legend.labs=c('Overall','C1 CC','C2 CC','C3 CC','C1 VH','C2 VH','C3 VH')
                               ,pal=c('black',mycolors1[1:nlevels(ps1_clin_join2$cluster)],'red3','dodgerblue1','darkgoldenrod3')
                               ,legend.title=" ",linetype=c('solid',rep('twodash',3),rep('solid',3))
                               ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Remission Duration PS1',xlab='Time (years)'
                               ,font.main = 16,font.x = 14,font.y = 14,ylab='Complete Remission probability',risk.table.height = 0.25,tables.col = "strata"
                               ,risk.table=T,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                               ,size=0.5,fontsize=4,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)

rem_plot_ps1_join2$plot<-rem_plot_ps1_join2$plot+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

rem_plot_ps1_join2$table<-rem_plot_ps1_join2$table+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,legend.position='none',plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text=element_text(size=12))

pval_table_ps1_rem<-tableGrob(symnum(pairwise_survdiff(Surv(rem_time,relapse)~cluster2,data=ps1_rem_join2)$p.value
                                     ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                     ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                       core=list(fg_params=list(cex=1,fontface=2)
                                                 ,bg_params=list(fill='white',col='black'))
                                       ,colhead=list(fg_params=list(col='white',fontface=2)
                                                     ,bg_params=list(fill=c(rep(mycolors1[1:3],2)),col=NA))
                                       ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                     ,bg_params = list(fill=c('white',c(mycolors1[2:3],mycolors1[1:3])),col=NA)))
                              ,widths=unit(c(rep(2.5,5)),c(rep('cm',5))))

rem_ps1_new<-ggarrange(rem_plot_ps1_join2$plot,rem_plot_ps1_join2$table,pval_table_ps1_rem,ncol=1,nrow=3,heights=c(0.7,0.3,0.25))

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S1A.pdf',rem_ps1_new, height=11.69,width=8.27,units='in')

##Remission PS2 (Supplementary figure S1B)##
ps2_rem_join2<-ps2_clin_join2[!is.na(ps2_clin_join2$relapse),]
rem_fit_ps2_join2<-surv_fit(Surv(rem_time,relapse)~trt+cluster,data=ps2_rem_join2)
rem_plot_ps2_join2<-ggsurvplot(rem_fit_ps2_join2,data=ps2_rem_join2,legend.labs=c('Overall','C1 CC','C2 CC','C1 VH','C2 VH')
                               ,pal=c('black',mycolors1[4:5],'blue2','purple2')
                               ,legend.title=" ",linetype=c('solid',rep('twodash',2),rep('solid',2))
                               ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Remission Duration PS2',xlab='Time (years)'
                               ,font.main = 16,font.x = 14,font.y = 14,ylab='Complete Remission probability',risk.table.height = 0.25,tables.col = "strata"
                               ,risk.table=T,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                               ,size=0.5,fontsize=4,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)

rem_plot_ps2_join2$plot<-rem_plot_ps2_join2$plot+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

rem_plot_ps2_join2$table<-rem_plot_ps2_join2$table+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,legend.position='none',plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text=element_text(size=12))


pval_table_ps2_rem<-tableGrob(symnum(pairwise_survdiff(Surv(rem_time,relapse)~cluster2,data=ps2_rem_join2)$p.value
                                     ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                     ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                       core=list(fg_params=list(cex=1,fontface=2)
                                                 ,bg_params=list(fill='white',col='black'))
                                       ,colhead=list(fg_params=list(col='white',fontface=2)
                                                     ,bg_params=list(fill=c(rep(mycolors1[4:5],2)),col=NA))
                                       ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                     ,bg_params=list(fill=c('white',c(mycolors1[5],mycolors1[4:5])),col=NA)))
                              ,widths=unit(c(rep(4.2,3)),c(rep('cm',3))))

rem_ps2_new<-ggarrange(rem_plot_ps2_join2$plot,rem_plot_ps2_join2$table,pval_table_ps2_rem,ncol=1,nrow=3,heights=c(0.7,0.3,0.25))

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S1B.pdf',rem_ps2_new, height=11.69,width=8.27,units='in')

##Remission PS3 (Supplementary figure S1C)##
ps3_rem_join2<-ps3_clin_join2[!is.na(ps3_clin_join2$relapse),]
rem_fit_ps3_join2<-surv_fit(Surv(rem_time,relapse)~trt+cluster,data=ps3_rem_join2)
rem_plot_ps3_join2<-ggsurvplot(rem_fit_ps3_join2,data=ps3_rem_join2,legend.labs=c('Overall','C1 CC','C2 CC','C1 VH','C2 VH')
                               ,pal=c('black',mycolors1[6:7],'green2','darkorange1')
                               ,legend.title=" ",linetype=c('solid',rep('twodash',2),rep('solid',2))
                               ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Remission Duration PS3',xlab='Time (years)'
                               ,font.main = 16,font.x = 14,font.y = 14,ylab='Complete Remission probability',risk.table.height = 0.25,tables.col = "strata"
                               ,risk.table=T,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                               ,size=0.5,fontsize=4,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)

rem_plot_ps3_join2$plot<-rem_plot_ps3_join2$plot+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

rem_plot_ps3_join2$table<-rem_plot_ps3_join2$table+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,legend.position='none',plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text=element_text(size=12))


pval_table_ps3_rem<-tableGrob(symnum(pairwise_survdiff(Surv(rem_time,relapse)~cluster2,data=ps3_rem_join2)$p.value
                                     ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                     ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                       core=list(fg_params=list(cex=1,fontface=2)
                                                 ,bg_params=list(fill='white',col='black'))
                                       ,colhead=list(fg_params=list(col='white',fontface=2)
                                                     ,bg_params=list(fill=c(rep(mycolors1[6:7],2)),col=NA))
                                       ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                     ,bg_params=list(fill=c('white',mycolors1[7],mycolors1[6:7]),col=NA)))
                              ,widths=unit(c(rep(4.2,3)),c(rep('cm',3))))

rem_ps3_new<-ggarrange(rem_plot_ps3_join2$plot,rem_plot_ps3_join2$table,pval_table_ps3_rem,ncol=1,nrow=3,heights=c(0.7,0.3,0.25))

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S1C.pdf',rem_ps3_new, height=11.69,width=8.27,units='in')

##Remission PS1, 2 and 3 combined (Supplementary figure S1D)##
rem_data_merge2<-clin_merge2[!is.na(clin_merge2$relapse),]
rem_fit_merge2<-surv_fit(Surv(rem_time,relapse)~trt+cluster,data=rem_data_merge2)
rem_plot_merge2<-ggsurvplot(rem_fit_merge2,data=rem_data_merge2,legend.labs=c('Overall','C1 CC','C2 CC','C3 CC','C4 CC','C5 CC','C1 VH','C2 VH','C3 VH','C4 VH','C5 VH')
                            ,pal=c('black',mycolors1[c(1,4:7)],'red3','blue2','purple2','green2','darkorange1')
                            ,legend.title=" ",linetype=c('solid',rep('twodash',5),rep('solid',5))
                            ,xlim = c(0,10),break.x.by=1,conf.int=F,title= 'Remission Duration PS1,2 and 3',xlab='Time (years)'
                            ,font.main = 16,font.x = 14,font.y = 14,ylab='Complete Remission probability',risk.table.height = 0.25,tables.col = "strata"
                            ,risk.table=T,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05)
                            ,size=0.5,fontsize=4,add.all=T,pval=T,pval.coord = c(6,0.9),pval.size=6,font.tickslab = 12,censor.size=3)

rem_plot_merge2$plot<-rem_plot_merge2$plot+
  theme(legend.text=element_text(size = 10),legend.key.size=unit(0.5, "cm")
        ,plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

rem_plot_merge2$table<-rem_plot_merge2$table+
  theme(legend.text=element_text(size = 12),legend.key.size=unit(1.25, "line")
        ,legend.position='none',plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text=element_text(size=12))


pval_table_merged_rem<-tableGrob(symnum(pairwise_survdiff(Surv(rem_time,relapse)~cluster2,data=rem_data_merge2)$p.value
                                        ,cutpoints=c(0,0.001,0.01,0.05, 1),symbols=c("p<0.001","p<0.01","p<0.05","ns ")
                                        ,abbr.colnames=F,na='-'),theme=ttheme_minimal(
                                          core=list(fg_params=list(cex=1,fontface=2)
                                                    ,bg_params=list(fill='white',col='black'))
                                          ,colhead=list(fg_params=list(col='white',fontface=2)
                                                        ,bg_params=list(fill=c(rep(mycolors1[c(1,4:7)],2)),col=NA))
                                          ,rowhead=list(fg_params=list(col='white',fontface=2)
                                                        ,bg_params = list(fill=c('white',c(mycolors1[c(4:7)],mycolors1[c(1,4:7)])),col=NA)))
                                 ,widths=unit(c(rep(0.001,9)),c(rep('cm',9))))

rem_merge3<-ggarrange(rem_plot_merge2$plot,rem_plot_merge2$table,pval_table_merged_rem,ncol=1,nrow=3,heights=c(0.7,0.3,0.25))

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S1D.pdf',rem_merge3, height=11.69,width=8.27,units='in')

####Supplementary figure S2: Survival and Remission by cluster###
#Adjust color object
mycolors2<-c('red2','blue3','purple3','green3','darkorange2')
mycolors3<-c('red3','blue2','purple2','green2','darkorange1')
sv_cl_list<-list()
rem_cl_list<-list()

for (i in 1:5){
  clin_cl<-clin_merge2[clin_merge2$cluster==paste(c("C"),as.character(i),sep=''),]
  clin_cl$cluster<-droplevels(clin_cl$cluster)
  levels(clin_cl$cluster)<-c(paste(c("C"),as.character(i),sep=''))
  clin_cl$cluster2<-droplevels(clin_cl$cluster2)
  
  ##Survival##
  km_surv_fit_clin_cl<-surv_fit(Surv(surv_time,status)~trt+cluster,data=clin_cl)
  clin_cl_plot<-ggsurvplot(km_surv_fit_clin_cl,data=clin_cl,pal=c('black',mycolors2[i],mycolors3[i])
                           ,legend.labs=c('Overall',paste(c("CC C"),as.character(i),sep=''),paste(c("VH C"),as.character(i),sep=''))
                           ,xlim = c(0,10),break.x.by=1,conf.int=F,title=paste(c("OS Cluster"),as.character(i),sep=''),xlab='Time (years)',legend.title=" "
                           ,font.main = 8,font.x = 7,font.y = 7,ylab='Cumulative Survival probability',linetype=c('solid','twodash','solid')
                           ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,risk.table=T,tables.col = "strata",add.all=T
                           ,font.tickslab = 7,fontsize =2.5,censor.size=2,pval=T,pval.coord = c(4,0.95),pval.size=3,risk.table.height = 0.35)
  clin_cl_plot$plot<-clin_cl_plot$plot+
    theme(legend.text=element_text(size=7),legend.key.size=unit(0.5,"cm"))
  clin_cl_plot$table<-clin_cl_plot$table+
    theme(legend.text=element_text(size=7),legend.key.size=unit(0.5,"cm")
          ,legend.position='none',axis.text=element_text(size=7),axis.title.x=element_text(size=7)
          ,plot.title=element_text(size=7))
   sv_cl_list[[i]]<-clin_cl_plot
}

for (i in 1:5){
  rem_cl<-clin_merge2[clin_merge2$cluster==paste(c("C"),as.character(i),sep='') & !is.na(clin_merge2$relapse),]
  rem_cl$cluster<-droplevels(rem_cl$cluster)
  levels(rem_cl$cluster)<-c(paste(c("C"),as.character(i),sep=''))
  
  ##Remission##
  km_rem_fit_rem_cl<-surv_fit(Surv(rem_time,relapse)~trt+cluster,data=rem_cl)
  rem_cl_rem_plot<-ggsurvplot(km_rem_fit_rem_cl,data=rem_cl,pal=c('black',mycolors2[i],mycolors3[i])
                              ,legend.labs=c('Overall',paste(c("CC C"),as.character(i),sep=''),paste(c("VH C"),as.character(i),sep=''))
                              ,xlim = c(0,10),break.x.by=1,conf.int=F,title=paste(c("RD Cluster"),as.character(i),sep=''),xlab='Time (years)',legend.title=" "
                              ,font.main = 8,font.x = 7,font.y = 7,ylab='Complete Remission probability',linetype=c('solid','twodash','solid')
                              ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,risk.table=T,tables.col = "strata",add.all=T
                              ,font.tickslab = 7,fontsize=2.5,censor.size=2,pval=T,pval.coord = c(4,0.95),pval.size=3,risk.table.height = 0.35)
  rem_cl_rem_plot$plot<-rem_cl_rem_plot$plot+
    theme(legend.text=element_text(size = 7), legend.key.size = unit(0.5, "cm"))
  rem_cl_rem_plot$table<-rem_cl_rem_plot$table+
    theme(legend.text=element_text(size=7),legend.key.size=unit(0.5,"cm")
          ,legend.position='none',axis.text=element_text(size=7),axis.title.x=element_text(size=7)
          ,plot.title=element_text(size=7))
  rem_cl_list[[i]]<-rem_cl_rem_plot
}

#Make KM plots for cluster 5 removing favorable cytogenetics cases
clin_cl_c5<-clin_merge2[clin_merge2$cluster=="C5",]
clin_cl_c5$cluster<-droplevels(clin_cl_c5$cluster)
levels(clin_cl_c5$cluster)<-"C5"
clin_cl_c5$cluster2<-droplevels(clin_cl_c5$cluster2)
dim(clin_cl_c5)
clin_cl_c5<-clin_cl_c5[!(clin_cl_c5$cyto_risk=='favorable'),]
dim(clin_cl_c5)
rem_cl_c5<-clin_cl_c5[!is.na(clin_cl_c5$relapse),]

##Survival##
km_surv_fit_clin_cl_c5<-surv_fit(Surv(surv_time,status)~trt+cluster,data=clin_cl_c5)
clin_cl_c5_plot<-ggsurvplot(km_surv_fit_clin_cl_c5,data=clin_cl_c5,pal=c('black',mycolors2[5],mycolors3[5])
                            ,legend.labs=c('Overall',"CC C5","VH C5")
                            ,xlim = c(0,10),break.x.by=1,conf.int=F,title="OS C5 without Fav. Cyto",xlab='Time (years)',legend.title=" "
                            ,font.main = 8,font.x = 7,font.y = 6,ylab='Cumulative Survival probability',linetype=c('solid','twodash','solid')
                            ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,risk.table=T,tables.col = "strata",add.all=T
                            ,font.tickslab = 7,fontsize =2.5,censor.size=2,pval=T,pval.coord = c(4,0.95),pval.size=3,risk.table.height = 0.35)
clin_cl_c5_plot$plot<-clin_cl_c5_plot$plot+
  theme(legend.text=element_text(size=7),legend.key.size=unit(0.5, "cm"))
clin_cl_c5_plot$table<-clin_cl_c5_plot$table+
  theme(legend.text=element_text(size=7),legend.key.size=unit(0.5,"cm")
        ,legend.position='none',axis.text=element_text(size=7),axis.title.x=element_text(size=7)
        ,plot.title=element_text(size=7))
sv_cl_list[[6]]<-clin_cl_c5_plot

##Remission##
km_rem_fit_rem_cl_c5<-surv_fit(Surv(rem_time,relapse)~trt+cluster,data=rem_cl_c5)
rem_cl_c5_rem_plot<-ggsurvplot(km_rem_fit_rem_cl_c5,data=rem_cl_c5,pal=c('black',mycolors2[5],mycolors3[5])
                               ,legend.labs=c('Overall',"CC C5","VH C5")
                               ,xlim = c(0,10),break.x.by=1,conf.int=F,title="RD C5 without Fav. Cyto",xlab='Time (years)',legend.title=" "
                               ,font.main = 8,font.x = 7,font.y = 7,ylab='Complete Remission probability',linetype=c('solid','twodash','solid')
                               ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,risk.table=T,tables.col = "strata",add.all=T
                               ,font.tickslab = 7,fontsize=2.5,censor.size=2,pval=T,pval.coord = c(4,0.95),pval.size=3,risk.table.height = 0.35)
rem_cl_c5_rem_plot$plot<-rem_cl_c5_rem_plot$plot+
  theme(legend.text=element_text(size = 7), legend.key.size = unit(0.5, "cm"))
rem_cl_c5_rem_plot$table<-rem_cl_c5_rem_plot$table+
  theme(legend.text=element_text(size=7),legend.key.size=unit(0.5,"cm")
        ,legend.position='none',axis.text=element_text(size=7),axis.title.x=element_text(size=7)
        ,plot.title=element_text(size=7))
rem_cl_list[[6]]<-rem_cl_c5_rem_plot

km_cl_plots<-arrange_ggsurvplots(list(
  sv_cl_list[[1]],sv_cl_list[[3]],sv_cl_list[[5]],rem_cl_list[[1]],rem_cl_list[[3]],rem_cl_list[[5]]
  ,sv_cl_list[[2]],sv_cl_list[[4]],sv_cl_list[[6]],rem_cl_list[[2]],rem_cl_list[[4]],rem_cl_list[[6]])
  ,nrow=3,ncol=4,print=F)

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S2.pdf',km_cl_plots, width=11.69,height=8.27,units='in')

####Supplementary Figure S3: Full correlation plot and protein networks###
##Full correlation plot (Supplementary Figure S3A)###
setwd('/results')#custom destination folder for file generation
pdf(file='Supplementary_Figure_S3A.pdf',height=11.69,width=8.27)
cor_ps_full<-cor(merge_ptn[,c(ptn_names)],method='pearson')
test_cor_ps_full<-cor.mtest(merge_ptn[,row.names(cor_ps_full)],conf.level = 0.95)
cor_plot_ps_full<-corrplot(cor_ps_full,order='hclust',hclust.method='ward.D2',method='color'
                           ,col=rev(COL2('RdYlBu')),tl.srt=45,cl.pos='b',diag=F
                           ,p.mat=test_cor_ps_full$p,tl.col='black',type='lower',mar=c(0,0,1,0)
                           ,sig.level=c(0.001,0.01,0.05),pch.cex=0.25,insig='label_sig',cl.ratio=0.1
                           ,tl.cex=0.4,pch.col='grey20',cl.cex=0.6,addgrid.col=NA 
                           ,title="Correlations between PS1, 2 and 3 proteins")
dev.off()
dev.off()

###Suplementary Table S4: Correlation coefficients and p-values for all comparisons
sup_s4<-createWorkbook()
sh1_sup_s4<-createSheet(sup_s4,"Coefficients")
sh2_sup_s4<-createSheet(sup_s4,"p-values")
cell_style2<-CellStyle(sup_s4)+Font(sup_s4,isBold=TRUE)+Border()
addDataFrame(data.frame("Supplementary Table S4. Matrix of correlation coefficients p-values for each comparison"=double()
                        ,check.names=FALSE),sheet=sh1_sup_s4,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(data.frame("Supplementary Table S4. Matrix of p-values for each comparison"=double()
                        ,check.names=FALSE),sheet=sh2_sup_s4,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(cor_ps_full,sheet=sh1_sup_s4,startColumn=1,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(test_cor_ps_full$p,sheet=sh2_sup_s4,startColumn=1,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
setwd('/results')#custom destination folder for file generation
saveWorkbook(sup_s4,"Supplementary_table_S4.xlsx")

######Generate Protein Networks in Cytoscape (Supplementary Figure S3B-G)t###
#Adjust protein expression dataset of PS1,2 and 3 combined for Cytoscape
ps_all_ptn<-merge_ptn[,c('PKM2','NOTCH1.cle','CD4','CDKN1B','CDKN1B.pS10','ETS1','LCK','ADM'
                         ,'UGT1A1','JAG1','CBX7','STAT1','PEA15','STK4','EZR','JMJD6'
                         ,'CASP3.cle','PAX6','INPP4B','G6PD','SPARC','TNFRSF4'
                         ,'SGK1','SOX17','DLST','SOX2','PDL1','PIM2','NOTCH1','BIRC2','PPARA','TYRO3','NDUFB4','ASS1','RHEB'
                         ,'AKT3','AKR1C3','HDAC3','SCD','WIPI2','EIF4G1','HSP90AA1_B1','FASN','FOXM1'
                         ,'H3K27Me3','ASNS','RAF1.pS338','MKNK1','STAT3.pY705','GAB2','TUBA1A','EIF4E'
                         ,'ITGAL','GSKA_B','MAP1LC3A_B','SPI1','MAPK14.pT180_Y182','HSPB1.pS82'
                         ,'XPF','EZH2','GATA1','WEE1','TP53BP1','CDK1_2_3.pT14','WEE1.pS642','EIF2AK2'
                         ,'PTEN','EEF2K','PDCD4','LMNB1','KIT','BCL2','ASH2L','SMARCB1','MEN1','BRD4'
                         ,'CDKN1B.pT198','NLN','HSF1.pS326','CD74','ARID1A','SMARCA2','H3K27Ac'
                         ,'EIF4E.pS209','EIF4EBP1.pS65','MECOM','FOS','CDH1','SMAD2','SGK3','HK2'
                         ,'BABAM1.pS29','MAPK14','TSC2','CHEK2','CASP3','EP300','HSF1','LYN'
                         ,'STAT3.pS727','NOTCH3','BMI1','IRS1','SUZ12','IGFBP2','MTOR'
                         ,'AURORA_A_B_C.pT288_232_198','RAB11','BRAF.pS445'
                         ,'cluster')]
ps_all<-ps_all_ptn[,sort(colnames(ps_all_ptn))]
ps_all<-ps_all %>% relocate('cluster', .after=everything())

#Calculate mean expression by cluster
ps_all_mean<-data.frame(matrix(NA,nrow=6,ncol=(ncol(ps_all)-1)))
colnames(ps_all_mean)<-colnames(ps_all)[1:(ncol(ps_all)-1)]
rownames(ps_all_mean)<-c(paste0('mean_ps_all_c',c(1:5)),'mean_ps_all')

for (i in 1:nrow(ps_all_mean)){
  ps_all_mean[i,]<-apply(ps_all[ps_all$cluster==paste('C',as.character(i),sep=''),1:(ncol(ps_all)-1)],2,function(x){mean(x)})
}
ps_all_mean[6,]<-apply(ps_all[,1:(ncol(ps_all)-1)],2,function(x){mean(x)})
ps_all_mean<-as.data.frame(t(ps_all_mean))

#Merge dataframe to identify the PS of every protein
ps_id<-as.data.frame(c(rep(1,52),rep(2,28),rep(3,20),rep(4,3),rep(5,6)))
#Keep the same order for ID the PS set, but change it to numbers for further analysis
rownames(ps_id)<-ptn_names
colnames(ps_id)<-c('ps_id')
ps_all_mean<-merge(ps_all_mean,ps_id,by=0)
row.names(ps_all_mean)<-ps_all_mean[,1]
ps_all_mean<-ps_all_mean[,-1]

#Adjust protein names for importing dataframes into StringApp 
ps_all_mean_adj<-ps_all_mean[!row.names(ps_all_mean) %in% c('H3K27Ac','H3K27Me3'),]
#Remove H3 PTM since STRING_APP does not accept those 
ps_all_mean_adj$rppa_ptn<-row.names(ps_all_mean_adj)
row.names(ps_all_mean_adj)<-gsub("CDK1_2_3.pT14","CDK1",row.names(ps_all_mean_adj))
ps_all_mean_adj$rppa_ptn<-gsub("CDK1_2_3.pT14","CDK1.pT14",ps_all_mean_adj$rppa_ptn)
row.names(ps_all_mean_adj)<-gsub("AURORA_A_B_C.pT288_232_198","AURKA",row.names(ps_all_mean_adj))
ps_all_mean_adj$rppa_ptn<-gsub("AURORA_A_B_C.pT288_232_198","AURORA_A.pT288",ps_all_mean_adj$rppa_ptn)
row.names(ps_all_mean_adj)<-gsub("HSP90AA1_B1","HSP90AA1",row.names(ps_all_mean_adj))
ps_all_mean_adj$rppa_ptn<-gsub("HSP90AA1_B1","HSP90AA1",ps_all_mean_adj$rppa_ptn)
row.names(ps_all_mean_adj)<-gsub("MAP1LC3A_B","MAP1LC3A",row.names(ps_all_mean_adj))
ps_all_mean_adj$rppa_ptn<-gsub("MAP1LC3A_B","MAP1LC3A",ps_all_mean_adj$rppa_ptn)
row.names(ps_all_mean_adj)<-gsub("GSKA_B","GSK3A",row.names(ps_all_mean_adj))
ps_all_mean_adj$rppa_ptn<-gsub("GSKA_B","GSKA",ps_all_mean_adj$rppa_ptn)
row.names(ps_all_mean_adj)<-gsub("RAB11","RAB11A",row.names(ps_all_mean_adj))
ps_all_mean_adj$rppa_ptn<-gsub("RAB11","RAB11A",ps_all_mean_adj$rppa_ptn)
ps_all_mean_adj_subset<-ps_all_mean_adj[c('AURKA','CDK1','HSP90AA1','MAP1LC3A','GSK3A','RAB11A'),] 
ps_all_mean_adj_subset1<-ps_all_mean_adj_subset
rownames(ps_all_mean_adj_subset1)<-c('AURKB','CDK2','HSP90AB1','MAP1LC3B','GSK3B','RAB11B')
ps_all_mean_adj_subset1$rppa_ptn<-c('AURORA_B.pT232','CDK2.pT14','HSP90AB1','MAP1LC3B','GSKB','RAB11B')
ps_all_mean_adj_subset2<-ps_all_mean_adj_subset[c('AURKA','CDK1'),]
rownames(ps_all_mean_adj_subset2)<-c('AURKC','CDK3')
ps_all_mean_adj_subset2$rppa_ptn<-c('AURORA_C.pT198','CDK3.pT14')
ps_all_mean_adj<-rbind(ps_all_mean_adj,ps_all_mean_adj_subset1,ps_all_mean_adj_subset2)
ps_all_mean_adj<-ps_all_mean_adj[order(row.names(ps_all_mean_adj)),]
dim(ps_all_mean_adj)

#Adjust names of proteins that only have PTM and others to match STRING algorithm
row.names(ps_all_mean_adj)<-gsub("BABAM1.pS29","BABAM1",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("BRAF.pS445","BRAF",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("RAF1.pS338","RAF1",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("EIF4EBP1.pS65","EIF4EBP1",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("HSPB1.pS82","HSPB1",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("RAF1.pS338","RAF1",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("STAT3.pS727","STAT3",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("PKM2","PKM",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("PDL1","CD274",row.names(ps_all_mean_adj))
row.names(ps_all_mean_adj)<-gsub("XPF","ERCC4",row.names(ps_all_mean_adj))

#Network 1 only has the total proteins
ps_all_mean1<-ps_all_mean_adj[!row.names(ps_all_mean_adj) %in% c("CASP3.cle","CDKN1B.pS10","CDKN1B.pT198"
                                                                 ,"EIF4E.pS209","HSF1.pS326","MAPK14.pT180_Y182"
                                                                 ,'NOTCH1.cle','STAT3.pY705','WEE1.pS642'),]

#Network 2 only has the PTM proteins and CDKN1B.pS10
ps_all_mean2<-ps_all_mean_adj[!row.names(ps_all_mean_adj) %in% c("CASP3","CDKN1B.pT198","CDKN1B"
                                                                 ,"EIF4E","HSF1","MAPK14","NOTCH1","STAT3","WEE1"),]
row.names(ps_all_mean2)<-row.names(ps_all_mean1)

#Network 3 only has the PTM proteins and CDKN1B.pT198
ps_all_mean3<-ps_all_mean_adj[!row.names(ps_all_mean_adj) %in% c("CASP3","CDKN1B.pS10","CDKN1B"
                                                                 ,"EIF4E","HSF1","MAPK14","NOTCH1"
                                                                 ,"STAT3","WEE1"),]
row.names(ps_all_mean3)<-row.names(ps_all_mean1)
ps_all_mean_list<-list(ps_all_mean1,ps_all_mean2,ps_all_mean3)

##Open Cytoscape software v3.10.1##
##Connect Cytoscape
cytoscapePing()
cytoscapeVersionInfo()

#Install necessary cytoscape apps
#install app 'Legend Creator' mannually
installApp('stringApp')
installApp('enhancedGraphics')

##Create Protein networks in Cytoscape##
#It is necessary to create 3 separate networks and merge afterwards to account for all ptns#
#Remove all PTM naming and repetead protein names to import for STRING analysis#
#cat(row.names(ps_all_mean1),sep=',') #this gives the list of protein already formatted for input
for (i in 1:3){
  strg_cmd<-'string protein query query="ADM,AKR1C3,AKT3,ARID1A,ASH2L,ASNS,ASS1,AURKA,AURKB,AURKC,BABAM1,BCL2,BIRC2,BMI1,BRAF,BRD4,CASP3,CBX7,CD4
,CD74,CDH1,CDK1,CDK2,CDK3,CDKN1B,CHEK2,DLST,EEF2K,EIF2AK2,EIF4E,EIF4EBP1,EIF4G1,EP300,ETS1,EZH2,EZR,FASN
,FOS,FOXM1,G6PD,GAB2,GATA1,GSK3A,GSK3B,HDAC3,HK2,HSF1,HSP90AA1,HSP90AB1,HSPB1,IGFBP2,INPP4B,IRS1,ITGAL
,JAG1,JMJD6,KIT,LCK,LMNB1,LYN,MAP1LC3A,MAP1LC3B,MAPK14,MECOM,MEN1,MKNK1,MTOR,NDUFB4,NLN,NOTCH1,NOTCH3
,PAX6,PDCD4,CD274,PEA15,PIM2,PKM,PPARA,PTEN,RAB11A,RAB11B,RAF1,RHEB,SCD,SGK1,SGK3,SMAD2,SMARCA2,SMARCB1
,SOX17,SOX2,SPARC,SPI1,STAT1,STAT3,STK4,SUZ12,TNFRSF4,TP53BP1,TSC2,TUBA1A,TYRO3,UGT1A1,WEE1,WIPI2,ERCC4"
  species="Homo sapiens" limit=0 cutoff=0.4'
  commandsRun(strg_cmd)
  
  #Load mean expression values of each cluster
  loadTableData(ps_all_mean_list[[i]],table.key.column="display name") 
  
  #Rename network
  renameNetwork(noquote(paste('ps_all_net',as.character(i),sep='')))
}

########################Merge the 3 networks in cytoscape environment###
#Tools>Merge>Networks;Union; Select c('ps_all_net1','ps_all_net2','ps_all_net3'); Advanced Options>Matching Columns='rppa_ptn'
setCurrentNetwork(network="Merged Network")

##Adjust Visual style
mystyle<-"dataStyle"
createVisualStyle(mystyle)
setVisualStyle(mystyle)

layoutNetwork('force-directed defaultSpringLength=1000 defaultSpringCoefficient=0.000003')
setNodeShapeDefault("ellipse",mystyle)
setNodeLabelMapping('rppa_ptn',mystyle)
setNodeColorDefault("#AAAAAA",mystyle)
setNodeSizeDefault(250,mystyle)
setNodeFontSizeDefault(100,mystyle)
setEdgeLineWidthDefault(3,mystyle)
setNodeBorderWidthDefault(40,mystyle)
#setNodeBorderColorDefault('#FD5903',mystyle)
#setNodeFontFaceDefault('SansSerif,bold,100',mystyle)

#Set Border color according to Protein Selector Set
ps_id_val<-c(1:5)
bord_col<-c('#1973BB','#EDC01C',mycolors1[c(1,5,6)])
setNodeBorderColorMapping('ps_id',ps_id_val,bord_col,style.name=mystyle,mapping.type='d')

##Visualize expression data (color inside nodes)
#Normalize sample signal
if(min(ps_all_mean)>=-2 & max(ps_all_mean)<=2){
  ps_all_node_breaks=seq(-2.4,2.4,0.4)
}else if(min(ps_all_mean)>=-2){
  ps_all_node_breaks=c(seq(-2.4,2,0.4),max(ps_all_mean))
}else if(max(ps_all_mean)<=2){
  ps_all_node_breaks=c(min(ps_all_mean),seq(-2,2.4,0.4))
}else{
  ps_all_node_breaks=c(min(ps_all_mean),seq(-2,2,0.4),max(ps_all_mean))
}

#Create color object
ps_all_node_colors<-jet.colors(length(ps_all_node_breaks))
setNodeColorMapping('mean_ps_all',ps_all_node_breaks,ps_all_node_colors,style.name=mystyle)

#Create Subnetworks
#Top correlated proteins network
top_cor_ptn<-c(row.names(cor_ps)) #get names from correlation matrix object
setCurrentNetwork(network="Merged Network")
createSubnetwork(top_cor_ptn,subnetwork.name='top_cor_ptn',nodes.by.col="rppa_ptn")

#Histone modifiers network
epi<-c('ASH2L','BRD4','EP300','EZH2','HDAC3','JMJD6','MEN1','SMARCB1','BMI1','CBX7','MECOM'
       ,'SUZ12','SMARCA2','ARID1A')
setCurrentNetwork(network="Merged Network")
createSubnetwork(epi,subnetwork.name='epi',nodes.by.col="rppa_ptn")

#Cell cycle and DDR network
cc_ddr<-c('CDK1.pT14','CDK2.pT14','CDK3.pT14','CDKN1B','CDKN1B.pS10','CDKN1B.pT198','CHEK2'
          ,'SPARC','UGT1A1','WEE1','BABAM1.pS29','FOXM1','MECOM','AURORA_A.pT288','AURORA_B.pT232'
          ,'AURORA_C.pT198','MAP1LC3A','MAP1LC3B','TP53BP1','XPF','WEE1.pS642')
setCurrentNetwork(network="Merged Network")
createSubnetwork(cc_ddr,subnetwork.name='cc_ddr',nodes.by.col="rppa_ptn")

#Ribosomal and transcription network
rb_trns<-c('EIF4E.pS209','EIF4EBP1.pS65','ETS1','SOX17','WIPI2','FOS','PAX6','PDCD4','EEF2K'
           ,'EIF4E','ARID1A','EIF2AK2','EIF4G1','SPI1')
setCurrentNetwork(network="Merged Network")
createSubnetwork(rb_trns,subnetwork.name='rb_trns',nodes.by.col="rppa_ptn")

#Metabolism network
metab<-c('EP300','HK2','FASN','SCD','AKR1C3','DLST','G6PD','PKM2','PPARA','GSKA','GSKB'
         ,'NDUFB4','ASNS','ASS1','ADM')
setCurrentNetwork(network="Merged Network")
createSubnetwork(metab,subnetwork.name='metab',nodes.by.col="rppa_ptn")

##################Adjust the space between nodes in all Networks in Cytoscape environment to your own preference###
#################Add legends in Cytoscape environment with Legend Creator and export figures###

setwd('/results') #custom destination folder for file generation
setCurrentNetwork(network="Merged Network")
exportImage(filename='Supplementary_Figure_S3B',type='SVG')
setCurrentNetwork(network="top_cor_ptn")
exportImage(filename='Supplementary_Figure_S3C',type='SVG')
setCurrentNetwork(network="epi")
exportImage(filename='Supplementary_Figure_S3D',type='SVG')
setCurrentNetwork(network="cc_ddr")
exportImage(filename='Supplementary_Figure_S3E',type='SVG')
setCurrentNetwork(network="rb_trns")
exportImage(filename='Supplementary_Figure_S3F',type='SVG')
setCurrentNetwork(network="metab")
exportImage(filename='Supplementary_Figure_S3G',type='SVG')

####Table1 and Supplementary table S6: Clinical and molecular tables### 
#Compute function to calculate p-values, using fisher exact test#
fisher.test.simulate.p.values <- function(data, variable, by, ...) {
  result <- list()
  test_results <- stats::fisher.test(data[[variable]],data[[by]],simulate.p.value=T,B=10000)
  result$p <- test_results$p.value
  result$test <- test_results$method
  result
}

#Expanded Clinical datatables (Supplementary Table S6)#
merge_tab<-clin_merge2 %>%
  select(cluster,age,age_gp,gender,race,wbc,blasts,hgb,plt,b2m,aml_gp,cyto_risk,cyto_prog,complex_kar,diploid,del5_5q,del7_7q
         ,t_8_21,t_9_11,t_6_9,t_15_17,t_11q23,inv16,del12,trisomy6,trisomy8,non_t_9_11_t_11q23,t_9_22_Ph_pos
         ,mut_aslx1,mut_cebpa,mut_dnmt3,mut_ezh2,mut_flt3,mut_flt3_d835,mut_flt3_itd,mut_gata2,mut_idh,mut_idh,mut_idh1
         ,mut_idh2,mut_jak2,mut_kit,mut_kmt2a,mut_mll,mut_npm1,mut_ptpn11,mut_ras,mut_runx1,mut_srsf2,mut_tet2,mut_tp53
         ,mut_u2af,mut_wt1)%>%
  tbl_summary(by=cluster,statistic=list(all_continuous()~"{mean}({sd})",all_categorical()~"{n}/{N}({p}%)")
              ,percent='col',digits=all_continuous()~1,type=all_categorical()~"categorical",missing='ifany',missing_text='N/A'
              ,label=list(cluster~"Cluster",age~'Age (years)',age_gp~'Age group (years)',gender~"Gender",race~"Race",aml_gp~"AML group"
                          ,wbc~'White Blood Cell count (K/uL)',blasts~'Blasts (%)',hgb~'Hemoglobin (g/dl)'
                          ,plt~'Platelets (K/uL)',b2m~'Serum B2M (ug/mL)'
                          ,cyto_prog~'Cytogenetic Prognosis',cyto_risk~'Cytogenetic Risk'
                          ,complex_kar~'Complex Karyotype',diploid~'Diploid Karyotype',t_9_22_Ph_pos~'t(9;22)Ph+'
                          ,del5_5q~'-5/5q-',del7_7q~'-7/7q-',t_8_21~'t(8;21)',t_15_17~'t(15;17)'
                          ,t_9_11~'t(9;11)',t_6_9~'t(6;9)',t_11q23~'t(11q23)',inv16~'Inv 16',del12~'Del 12'
                          ,trisomy6~'Trisomy 6',trisomy8~'Trisomy 8',non_t_9_11_t_11q23~'non-t(9;11)/t(11q23)'
                          ,mut_aslx1~"ASLX1 Mut",mut_cebpa~"CEBPA Mut",mut_dnmt3~"DNMT3 Mut",mut_ezh2~'EZH2 Mut'
                          ,mut_flt3~'FLT3 Mut',mut_flt3_d835~"FLT3_D835 Mut",mut_flt3_itd~"FLT3_ITD Mut",mut_gata2~"GATA2 Mut"
                          ,mut_idh~"IDH Mut",mut_idh1~"IDH1 Mut",mut_idh2~"IDH2 Mut",mut_jak2~"JAK2 Mut",mut_kit~"KIT Mut"
                          ,mut_kmt2a~"KMT2A Mut",mut_mll~"MLL Mut",mut_npm1~"NPM1 Mut",mut_ptpn11~"PTPN11 Mut"
                          ,mut_ras~"RAS Mut",mut_runx1~"RUNX1 Mut",mut_srsf2~"SRSF2 Mut",mut_tet2~"TET2 Mut"
                          ,mut_tp53~"TP53 Mut",mut_u2af~"U2AF Mut",mut_wt1~"WT1 Mut")) %>%
  add_p(test = list(all_categorical() ~ "fisher.test.simulate.p.values"),pvalue_fun = ~style_pvalue(.x, digits = 2))%>%
  bold_p()

merge_all<-clin_merge2 %>%
  select(age,age_gp,gender,race,wbc,blasts,hgb,plt,b2m,aml_gp,cyto_risk,cyto_prog,complex_kar,diploid,del5_5q,del7_7q
         ,t_8_21,t_9_11,t_6_9,t_15_17,t_11q23,inv16,del12,trisomy6,trisomy8,non_t_9_11_t_11q23,t_9_22_Ph_pos
         ,mut_aslx1,mut_cebpa,mut_dnmt3,mut_ezh2,mut_flt3,mut_flt3_d835,mut_flt3_itd,mut_gata2,mut_idh,mut_idh,mut_idh1
         ,mut_idh2,mut_jak2,mut_kit,mut_kmt2a,mut_mll,mut_npm1,mut_ptpn11,mut_ras,mut_runx1,mut_srsf2,mut_tet2,mut_tp53
         ,mut_u2af,mut_wt1)%>%
  tbl_summary(statistic=list(all_continuous()~"{mean}({sd})",all_categorical()~"{n}/{N}({p}%)")
              ,digits=all_continuous()~1,type=all_categorical()~"categorical",missing='ifany',missing_text='N/A'
              ,label=list(age~'Age (years)',age_gp~'Age group (years)',gender~"Gender",race~"Race",aml_gp~"AML group"
                          ,wbc~'White Blood Cell count (K/uL)',blasts~'Blasts (%)',hgb~'Hemoglobin (g/dl)'
                          ,plt~'Platelets (K/uL)',b2m~'Serum B2M (ug/mL)'
                          ,cyto_prog~'Cytogenetic Prognosis',cyto_risk~'Cytogenetic Risk'
                          ,complex_kar~'Complex Karyotype',diploid~'Diploid Karyotype',t_9_22_Ph_pos~'t(9;22)Ph+'
                          ,del5_5q~'-5/5q-',del7_7q~'-7/7q-',t_8_21~'t(8;21)',t_15_17~'t(15;17)'
                          ,t_9_11~'t(9;11)',t_6_9~'t(6;9)',t_11q23~'t(11q23)',inv16~'Inv 16',del12~'Del 12'
                          ,trisomy6~'Trisomy 6',trisomy8~'Trisomy 8',non_t_9_11_t_11q23~'non-t(9;11)/t(11q23)'
                          ,mut_aslx1~"ASLX1 Mut",mut_cebpa~"CEBPA Mut",mut_dnmt3~"DNMT3 Mut",mut_ezh2~'EZH2 Mut'
                          ,mut_flt3~'FLT3 Mut',mut_flt3_d835~"FLT3_D835 Mut",mut_flt3_itd~"FLT3_ITD Mut",mut_gata2~"GATA2 Mut"
                          ,mut_idh~"IDH Mut",mut_idh1~"IDH1 Mut",mut_idh2~"IDH2 Mut",mut_jak2~"JAK2 Mut",mut_kit~"KIT Mut"
                          ,mut_kmt2a~"KMT2A Mut",mut_mll~"MLL Mut",mut_npm1~"NPM1 Mut",mut_ptpn11~"PTPN11 Mut"
                          ,mut_ras~"RAS Mut",mut_runx1~"RUNX1 Mut",mut_srsf2~"SRSF2 Mut",mut_tet2~"TET2 Mut"
                          ,mut_tp53~"TP53 Mut",mut_u2af~"U2AF Mut",mut_wt1~"WT1 Mut"))

exp_tab<-tbl_merge(list(merge_tab, merge_all)) %>%
  bold_labels()%>%
  modify_header(label~'**Variable**')%>%
  modify_caption("**Supplementary Table S6.** Expanded Demographic, Clinical and Molecular Characteristics")%>%
  modify_spanning_header(c(paste0('stat_', c(1:nlevels(clin_merge2$cluster)),'_1'))~'**Patient Cluster**','stat_0_2'~'**Overall**','p.value_1'~'.')

setwd('/results')#custom destination folder for file generation
gtsave(data=as_gt(exp_tab),filename="Supplementary_table_S6.html")
#Open the saved .html file with excel and then save as .xlsx

#Only significant variables Clinical datatables (Table1)#
merge_tab2<-clin_merge2 %>%
  select(cluster,age,wbc,blasts,plt,aml_gp,cyto_risk,complex_kar,del5_5q,del7_7q,inv16
         ,mut_aslx1,mut_cebpa,mut_dnmt3,mut_ezh2,mut_flt3,mut_flt3_d835,mut_flt3_itd,mut_npm1,mut_tp53
  )%>%
  tbl_summary(by=cluster,statistic=list(all_continuous()~"{mean}({sd})",all_categorical()~"{p}%")
              ,percent='col',digits=all_continuous()~1
              ,type=all_dichotomous()~'dichotomous',missing='no',missing_text='N/A'
              ,value=list(aml_gp~'secondary',cyto_risk~'unfavorable')
              ,label=list(cluster~"Cluster"
                          ,age~'Age (years)',wbc~'White Blood Cell count (K/uL)'
                          ,aml_gp~"Secondary AML"
                          ,blasts~'Blasts (%)',plt~'Platelets (K/uL)',cyto_risk~'Unfavorable Cytogenetics'
                          ,complex_kar~'Complex Karyotype'
                          ,del5_5q~'-5/5q-',del7_7q~'-7/7q-',inv16~'Inv 16'
                          ,mut_aslx1~"ASLX1 Mutation",mut_cebpa~"CEBPA Mutation",mut_dnmt3~"DNMT3 Mutation",mut_ezh2~'EZH2 Mutation'
                          ,mut_flt3~'FLT3 Mutation',mut_flt3_d835~"FLT3 D835 Mutation",mut_flt3_itd~"FLT3 ITD Mutation"
                          ,mut_npm1~"NPM1 Mutation",mut_tp53~"TP53 Mutation"
              )) %>%
  add_p(test = list(all_categorical() ~ "fisher.test.simulate.p.values"),pvalue_fun = ~style_pvalue(.x, digits = 2))%>%
  bold_p()

merge_all2<-clin_merge2 %>%
  select(age,wbc,blasts,plt,aml_gp,cyto_risk,complex_kar,del5_5q,del7_7q,inv16
         ,mut_aslx1,mut_cebpa,mut_dnmt3,mut_ezh2,mut_flt3,mut_flt3_d835,mut_flt3_itd,mut_npm1,mut_tp53)%>%
  tbl_summary(statistic=list(all_continuous()~"{mean}({sd})",all_categorical()~"{p}%")
              ,digits=all_continuous()~1
              ,type=all_dichotomous()~'dichotomous',missing='no',missing_text='N/A'
              ,value=list(aml_gp~'secondary',cyto_risk~'unfavorable')
              ,label=list(age~'Age (years)',wbc~'White Blood Cell count (K/uL)',aml_gp~"Secondary AML"
                          ,blasts~'Blasts (%)',plt~'Platelets (K/uL)',cyto_risk~'Unfavorable Cytogenetics'
                          ,complex_kar~'Complex Karyotype',del5_5q~'-5/5q-',del7_7q~'-7/7q-',inv16~'Inv 16'
                          ,mut_aslx1~"ASLX1 Mutation",mut_cebpa~"CEBPA Mutation",mut_dnmt3~"DNMT3 Mutation",mut_ezh2~'EZH2 Mutation'
                          ,mut_flt3~'FLT3 Mutation',mut_flt3_d835~"FLT3 D835 Mutation",mut_flt3_itd~"FLT3 ITD Mutation"
                          ,mut_npm1~"NPM1 Mutation",mut_tp53~"TP53 Mutation"))

sig_tab<-tbl_merge(list(merge_all2,merge_tab2)) %>%
  bold_labels()%>%
  modify_header(label~'**Variable**')%>%
  modify_caption("**Table 1.** Significant Demographic, Clinical and Molecular Characteristics")%>%
  modify_spanning_header('stat_0_1'~'**Overall**',c(paste0('stat_',c(1:nlevels(clin_merge2$cluster)),'_2'))~'**Patient Cluster**','p.value_2'~'.')

setwd('/results')#custom destination folder for file generation
gtsave(data=as_gt(sig_tab),filename="Table1.html")
#Open the saved .html and then save as .pdf

####Supplementary Figure S4: OS and RD curves split 5 clusters by varibles###
var_sp2<-clin_merge2
var_sp2$male<-factor(ifelse(var_sp2$gender=='male',c('yes'),c('no')))
var_sp2$female<-factor(ifelse(var_sp2$gender=='female',c('yes'),c('no')))
var_sp2$primary_aml<-factor(ifelse(var_sp2$aml_gp=='primary',c('yes'),c('no')))
var_sp2$secondary_aml<-factor(ifelse(var_sp2$aml_gp=='secondary',c('yes'),c('no')))
var_sp2$cyto_favorable<-factor(ifelse(var_sp2$cyto_risk=='favorable',c('yes'),c('no')))
var_sp2$cyto_intermed<-factor(ifelse(var_sp2$cyto_risk=='intermediary',c('yes'),c('no')))
var_sp2$cyto_unfavorable<-factor(ifelse(var_sp2$cyto_risk=='unfavorable',c('yes'),c('no')))
var_sp2$C1<-factor(ifelse(var_sp2$cluster=='C1',c('yes'),c('no')))
var_sp2$C2<-factor(ifelse(var_sp2$cluster=='C2',c('yes'),c('no')))
var_sp2$C3<-factor(ifelse(var_sp2$cluster=='C3',c('yes'),c('no')))
var_sp2$C4<-factor(ifelse(var_sp2$cluster=='C4',c('yes'),c('no')))
var_sp2$C5<-factor(ifelse(var_sp2$cluster=='C5',c('yes'),c('no')))
##Create several categorical variables by splitting levels into 'yes' or 'no' variables

var_sp2<-var_sp2 %>% relocate('mut_aslx1','mut_cebpa','mut_dnmt3','mut_ezh2','mut_gata2'
                              ,'mut_idh','mut_idh1','mut_idh2','mut_flt3','mut_flt3_itd','mut_flt3_d835'
                              ,'mut_jak2','mut_kit','mut_kmt2a','mut_mll','mut_npm1','mut_ptpn11','mut_ras'
                              ,'mut_runx1','mut_srsf2','mut_tet2','mut_tp53','mut_u2af','mut_wt1'
                              ,'age_40_minus','age_41_55','age_56_70','age_70_plus','male','female'
                              ,'race_white','race_black','race_hisp','race_asian','race_other'
                              ,'primary_aml','secondary_aml','cyto_favorable','cyto_intermed','cyto_unfavorable'
                              ,'complex_kar','diploid','del5_5q','del7_7q','inv16'
                              ,'del12','t_8_21','t_11q23','t_9_11','t_6_9','trisomy8'
                              ,'trisomy6','non_t_9_11_t_11q23','t_9_22_Ph_pos'
                              ,'C1','C2','C3','cluster', everything())

#Remove variables that do not have all cluster levels
var_sp2_sv<-var_sp2[,c('mut_aslx1','mut_cebpa','mut_dnmt3','mut_gata2'
                       ,'mut_idh','mut_idh1','mut_idh2','mut_flt3','mut_flt3_itd'
                       ,'mut_jak2','mut_kit','mut_mll','mut_npm1','mut_ptpn11','mut_ras'
                       ,'mut_runx1','mut_srsf2','mut_tet2','mut_tp53','mut_wt1'
                       ,'age_40_minus','age_41_55','age_56_70','age_70_plus','male','female'
                       ,'race_white','race_black','race_hisp','race_asian'
                       ,'primary_aml','secondary_aml','cyto_favorable','cyto_intermed','cyto_unfavorable'
                       ,'complex_kar','diploid','del5_5q','del7_7q','inv16'
                       ,'del12','t_8_21','t_11q23','trisomy8'
                       ,'status','surv_time','relapse','rem_time','cluster')]

colnames(var_sp2_sv)<-c('ASXL1 Mutation','CEBPA Mutation','DNMT3 Mutation','GATA2 Mutation'
                        ,'IDH Mutation','IDH1 Mutation','IDH2 Mutation','FLT3 Mutation'
                        ,'FLT3_ITD Mutation','JAK2 Mutation','KIT Mutation','MLL Mutation'
                        ,'NPM1 Mutation','PTPN11 Mutation','RAS Mutation','RUNX1 Mutation'
                        ,'SRSF2 Mutation','TET2 Mutation','TP53 Mutation','WT1 Mutation'
                        ,'Age < 40','Age 41-55','Age 56-70','Age > 70','Males','Females'
                        ,'White Race','Black Race','Hispanic Race','Asian Race'
                        ,'Primary AML','Secondary AML','Fav. Cyto Risk','Intermed. Cyto Risk'
                        ,'Unfav. Cyto Risk','Complex Karyotype','Diploid Karyotype'
                        ,'-5/5q-','-7/7q-','Inv16','Del12','t(8;21)','t(11q23)'
                        ,'Trisomy 8','status','surv_time','relapse','rem_time','cluster')

#Remove variables that do not have all cluster levels 
var_sp2_rem<-var_sp2[,c('mut_dnmt3','mut_idh','mut_idh1','mut_idh2','mut_flt3','mut_flt3_itd'
                        ,'mut_npm1','mut_ptpn11','mut_ras'
                        ,'mut_srsf2','mut_tet2','mut_tp53'
                        ,'age_40_minus','age_41_55','age_56_70','age_70_plus','male','female'
                        ,'race_white','race_black','race_hisp'
                        ,'primary_aml','secondary_aml','cyto_favorable','cyto_intermed','cyto_unfavorable'
                        ,'complex_kar','diploid','del5_5q','del7_7q','inv16'
                        ,'t_8_21','t_11q23','trisomy8'
                        ,'status','surv_time','relapse','rem_time','cluster')]

colnames(var_sp2_rem)<-c('DNMT3 Mutation','IDH Mutation','IDH1 Mutation','IDH2 Mutation'
                         ,'FLT3 Mutation','FLT3_ITD Mutation'
                         ,'NPM1 Mutation','PTPN11 Mutation','RAS Mutation'
                         ,'SRSF2 Mutation','TET2 Mutation','TP53 Mutation'
                         ,'Age < 40','Age 41-55','Age 56-70','Age > 70','Males','Females'
                         ,'White Race','Black Race','Hispanic Race'
                         ,'Primary AML','Secondary AML','Fav. Cyto Risk','Intermed. Cyto Risk'
                         ,'Unfav. Cyto Risk','Complex Karyotype','Diploid Karyotype'
                         ,'-5/5q-','-7/7q-','Inv16','t(8;21)','t(11q23)'
                         ,'Trisomy 8','status','surv_time','relapse','rem_time','cluster')

#Survival
sv_var_list<-list()
for (i in 1:44){
  df_var_sv<-var_sp2_sv[var_sp2_sv[[i]]=='yes' & !is.na(var_sp2_sv[[i]]),]
  km_sv_fit_var<-survfit(Surv(surv_time,status)~df_var_sv$cluster,data=df_var_sv)
  plot_var_sv<-ggsurvplot(km_sv_fit_var,data=df_var_sv,pal=mycolors2
                          ,legend.title=" ",legend.labs=paste0('C',c(1:nlevels(df_var_sv$cluster)))
                          ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)'
                          ,title=paste(c("OS "),as.character(colnames(df_var_sv[i])),sep=' ')
                          ,font.main = 12,font.x = 10,font.y = 10,ylab='Survival probability',pval=T
                          ,pval.size=4,pval.coord=c(7.5,0.95),risk.table=T,risk.table.height = 0.3,tables.col = "strata"
                          ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,fontsize =2.8,censor.size=2)
  plot_var_sv$plot<-plot_var_sv$plot+
    theme(legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))
  plot_var_sv$table<-plot_var_sv$table+
    theme(legend.text=element_text(size=10),legend.key.size=unit(0.5,"cm"),axis.text=element_text(size=8)
          ,legend.position='none',axis.title.x=element_text(size=10),plot.title=element_text(size=10))
  sv_var_list[[i]]<-plot_var_sv
}
names(sv_var_list)<-c(colnames(var_sp2_sv[1:44]))

#Remission
rem_var_list<-list()
for (i in 1:34){
  df_var_rem<-var_sp2_rem[var_sp2_rem[[i]]=='yes' & !is.na(var_sp2_rem[[i]]),]
  km_rem_fit_var<-survfit(Surv(rem_time,relapse)~df_var_rem$cluster,data=df_var_rem)
  plot_var_rem<-ggsurvplot(km_rem_fit_var,data=df_var_rem,pal=mycolors2
                           ,legend.title=" ",legend.labs=paste0('C',c(1:nlevels(df_var_rem$cluster)))
                           ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)'
                           ,title=paste(c("RD "),as.character(colnames(df_var_rem[i])),sep=' ')
                           ,font.main = 12,font.x = 10,font.y = 10,ylab='Remission probability',pval=T
                           ,pval.size=4,pval.coord=c(7.5,0.95),risk.table=T,risk.table.height = 0.3,tables.col = "strata"
                           ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,fontsize =2.8,censor.size=2)
  plot_var_rem$plot<-plot_var_rem$plot+theme(legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))
  plot_var_rem$table<-plot_var_rem$table+
    theme(legend.text=element_text(size=10),legend.key.size=unit(0.5,"cm"),axis.text=element_text(size=8)
          ,legend.position='none',axis.title.x=element_text(size=10),plot.title=element_text(size=10))
  rem_var_list[[i]]<-plot_var_rem
}
names(rem_var_list)<-c(colnames(var_sp2_rem[1:34]))

var_split_plots<-arrange_ggsurvplots(list(
  sv_var_list[[21]],sv_var_list[[22]]
  ,sv_var_list[[23]],sv_var_list[[24]],sv_var_list[[25]],sv_var_list[[26]],sv_var_list[[27]],sv_var_list[[28]]
  ,sv_var_list[[29]],sv_var_list[[30]],sv_var_list[[31]],sv_var_list[[32]],sv_var_list[[33]],sv_var_list[[34]]
  ,sv_var_list[[35]],sv_var_list[[36]],sv_var_list[[37]],sv_var_list[[38]],sv_var_list[[39]],sv_var_list[[40]]
  ,sv_var_list[[41]],sv_var_list[[42]],sv_var_list[[43]],sv_var_list[[44]],sv_var_list[[1]],sv_var_list[[2]]
  ,sv_var_list[[3]],sv_var_list[[4]],sv_var_list[[5]],sv_var_list[[6]],sv_var_list[[7]],sv_var_list[[8]]
  ,sv_var_list[[9]],sv_var_list[[10]],sv_var_list[[11]],sv_var_list[[12]],sv_var_list[[13]],sv_var_list[[14]]
  ,sv_var_list[[15]],sv_var_list[[16]],sv_var_list[[17]],sv_var_list[[18]],sv_var_list[[19]],sv_var_list[[20]]
  
  ,rem_var_list[[14]],rem_var_list[[15]],rem_var_list[[16]],rem_var_list[[17]],rem_var_list[[18]],rem_var_list[[19]]
  ,rem_var_list[[20]],rem_var_list[[21]],rem_var_list[[22]],rem_var_list[[23]],rem_var_list[[24]],rem_var_list[[25]]
  ,rem_var_list[[26]],rem_var_list[[27]],rem_var_list[[28]],rem_var_list[[29]],rem_var_list[[30]],rem_var_list[[31]]
  ,rem_var_list[[32]],rem_var_list[[33]],rem_var_list[[34]],rem_var_list[[1]],rem_var_list[[2]],rem_var_list[[3]]
  ,rem_var_list[[4]],rem_var_list[[5]],rem_var_list[[6]],rem_var_list[[7]],rem_var_list[[8]],rem_var_list[[9]]
  ,rem_var_list[[10]],rem_var_list[[11]],rem_var_list[[12]])
  ,nrow=2,ncol=3,print=F)

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S4.pdf',var_split_plots, width=11.69,height=8.27,units='in')

####Supplementary Figure S5: OS and RD curves split 5 clusters x 2 trt (VH and CC) by variables###
var_sp3<-clin_merge2
var_sp3$age_70_minus<-factor(ifelse(var_sp3$age<=70,'yes','no'))
var_sp3$male<-factor(ifelse(var_sp3$gender=='male',c('yes'),c('no')))
var_sp3$female<-factor(ifelse(var_sp3$gender=='female',c('yes'),c('no')))
var_sp3$primary_aml<-factor(ifelse(var_sp3$aml_gp=='primary',c('yes'),c('no')))
var_sp3$secondary_aml<-factor(ifelse(var_sp3$aml_gp=='secondary',c('yes'),c('no')))
var_sp3$cyto_favorable<-factor(ifelse(var_sp3$cyto_risk=='favorable',c('yes'),c('no')))
var_sp3$cyto_intermed<-factor(ifelse(var_sp3$cyto_risk=='intermediary',c('yes'),c('no')))
var_sp3$cyto_unfavorable<-factor(ifelse(var_sp3$cyto_risk=='unfavorable',c('yes'),c('no')))

##Create several categorical variables by splitting levels into 'yes' or 'no' variables

var_sp3<-var_sp3 %>% relocate('mut_aslx1','mut_cebpa','mut_dnmt3','mut_ezh2','mut_gata2'
                              ,'mut_idh','mut_idh1','mut_idh2','mut_flt3','mut_flt3_itd','mut_flt3_d835'
                              ,'mut_jak2','mut_kit','mut_kmt2a','mut_mll','mut_npm1','mut_ptpn11','mut_ras'
                              ,'mut_runx1','mut_srsf2','mut_tet2','mut_tp53','mut_u2af','mut_wt1'
                              ,'age_70_minus','age_70_plus','male','female'
                              ,'race_white','race_black','race_hisp','race_asian','race_other'
                              ,'primary_aml','secondary_aml','cyto_favorable','cyto_intermed'
                              ,'cyto_unfavorable','complex_kar','diploid','del5_5q','del7_7q','inv16'
                              ,'del12','t_8_21','t_11q23','t_9_11','t_6_9','trisomy8'
                              ,'trisomy6','non_t_9_11_t_11q23','t_9_22_Ph_pos'
                              ,'trt','cluster', everything())

#Remove variables that do not have all cluster levels
var_sp3_sv<-var_sp3[,c('mut_dnmt3','mut_flt3','mut_flt3_itd'
                       ,'mut_runx1','mut_tet2','mut_tp53'
                       ,'age_70_minus','age_70_plus','male','female'
                       ,'race_white','primary_aml','secondary_aml'
                       ,'cyto_intermed','cyto_unfavorable','complex_kar','diploid','del5_5q'
                       ,'status','surv_time','relapse','rem_time','trt','cluster')]

colnames(var_sp3_sv)<-c('DNMT3 Mutation','FLT3 Mutation','FLT3_ITD Mutation'
                        ,'RUNX1 Mutation','TET2 Mutation','TP53 Mutation'
                        ,'Age < 70','Age > 70','Males','Females','White Race','Primary AML'
                        ,'Secondary AML','Intermed. Cyto Risk'
                        ,'Unfav. Cyto Risk','Complex Karyotype','Diploid Karyotype'
                        ,'-5/5q-','status','surv_time','relapse','rem_time','trt','cluster')

#Remove variables that do not have all cluster levels 
var_sp3_rem<-var_sp3[,c('mut_dnmt3','mut_flt3'
                        ,'mut_tet2','mut_tp53'
                        ,'age_70_minus','age_70_plus','male','female'
                        ,'race_white','primary_aml','secondary_aml'
                        ,'cyto_intermed','cyto_unfavorable','complex_kar','del5_5q'
                        ,'status','surv_time','relapse','rem_time','trt','cluster')]

colnames(var_sp3_rem)<-c('DNMT3 Mutation','FLT3 Mutation'
                         ,'TET2 Mutation','TP53 Mutation'
                         ,'Age < 70','Age > 70','Males','Females','White Race','Primary AML'
                         ,'Secondary AML','Intermed. Cyto Risk'
                         ,'Unfav. Cyto Risk','Complex Karyotype'
                         ,'-5/5q-','status','surv_time','relapse','rem_time','trt','cluster')

#Survival
sv_var_list2<-list()
for (i in 1:18){
  df_var_sv2<-var_sp3_sv[var_sp3_sv[[i]]=='yes' & !is.na(var_sp3_sv[[i]]),]
  km_sv_fit_var2<-survfit(Surv(surv_time,status)~trt+cluster,data=df_var_sv2)
  plot_var_sv2<-ggsurvplot(km_sv_fit_var2,data=df_var_sv2,pal=c(mycolors2,mycolors3)
                           ,legend.title=" ",legend.labs=c('C1 CC','C2 CC','C3 CC','C4 CC','C5 CC','C1 VH','C2 VH','C3 VH','C4 VH','C5 VH')
                           ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)',linetype=c(rep('twodash',5),rep('solid',5))
                           ,title=paste(c("OS "),as.character(colnames(df_var_sv2[i])),sep=' ')
                           ,font.main = 12,font.x = 10,font.y = 10,ylab='Survival probability',pval=T
                           ,pval.size=3,pval.coord=c(7.5,0.95),risk.table=T,risk.table.height = 0.4,tables.col = "strata"
                           ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,fontsize =2.5,censor.size=2)
  plot_var_sv2$plot<-plot_var_sv2$plot+theme(legend.text=element_text(size=8),legend.key.size = unit(0.4, "cm"))
  plot_var_sv2$table<-plot_var_sv2$table+
    theme(legend.text=element_text(size=8),legend.key.size=unit(0.4,"cm"),axis.text=element_text(size=8)
          ,legend.position='none',axis.title.x=element_text(size=10),plot.title=element_text(size=10))
  sv_var_list2[[i]]<-plot_var_sv2
}
names(sv_var_list2)<-c(colnames(var_sp3_sv[1:18]))

#Make a separate plot for favorable cytogenetics
fav_cyto<-var_sp3[var_sp3$cyto_favorable=='yes' & !is.na(var_sp3$cyto_favorable),]
km_surv_fit_fav_cyto<-survfit(Surv(surv_time,status)~trt+cluster,data=fav_cyto)
plot_fav_cyto<-ggsurvplot(km_surv_fit_fav_cyto,data=fav_cyto,pal=mycolors2
                          ,legend.title=" ",legend.labs=c('C1 CC','C2 CC','C3 CC','C4 CC','C5 CC')
                          ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)',linetype=rep('twodash',5)
                          ,title="OS Fav. Cyto Risk CC"
                          ,font.main = 12,font.x = 10,font.y = 10,ylab='Survival probability',pval=T
                          ,pval.size=3,pval.coord=c(7.5,0.95),risk.table=T,risk.table.height = 0.4,tables.col = "strata"
                          ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,fontsize =2.5,censor.size=2)
plot_fav_cyto$plot<-plot_fav_cyto$plot+theme(legend.text=element_text(size=8), legend.key.size = unit(0.4, "cm"))
plot_fav_cyto$table<-plot_fav_cyto$table+
  theme(legend.text=element_text(size=8),legend.key.size=unit(0.4,"cm"),axis.text=element_text(size=8)
        ,legend.position='none',axis.title.x=element_text(size=10),plot.title=element_text(size=10))
sv_var_list2[[19]]<-plot_fav_cyto

#Remission
rem_var_list2<-list()
for (i in 1:15){
  df_var_rem2<-var_sp3_rem[var_sp3_rem[[i]]=='yes' & !is.na(var_sp3_rem[[i]]),]
  km_rem_fit_var2<-survfit(Surv(rem_time,relapse)~trt+cluster,data=df_var_rem2)
  plot_var_rem2<-ggsurvplot(km_rem_fit_var2,data=df_var_rem2,pal=c(mycolors2,mycolors3)
                            ,legend.title=" ",legend.labs=c('C1 CC','C2 CC','C3 CC','C4 CC','C5 CC','C1 VH','C2 VH','C3 VH','C4 VH','C5 VH')
                            ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)',linetype=c(rep('twodash',5),rep('solid',5))
                            ,title=paste(c("RD "),as.character(colnames(df_var_rem2[i])),sep=' ')
                            ,font.main = 12,font.x = 10,font.y = 10,ylab='Remission probability',pval=T
                            ,pval.size=3,pval.coord=c(7.5,0.95),risk.table=T,risk.table.height = 0.4,tables.col = "strata"
                            ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,fontsize =2.5,censor.size=2)
  plot_var_rem2$plot<-plot_var_rem2$plot+theme(legend.text=element_text(size=8), legend.key.size = unit(0.4, "cm"))
  plot_var_rem2$table<-plot_var_rem2$table+
    theme(legend.text=element_text(size=8),legend.key.size=unit(0.4,"cm"),axis.text=element_text(size=8)
          ,legend.position='none',axis.title.x=element_text(size=10),plot.title=element_text(size=10))
  rem_var_list2[[i]]<-plot_var_rem2
}
names(rem_var_list2)<-c(colnames(var_sp3_rem[1:15]))

#Make a separate plot for favorable cytogenetics
rem_fav_cyto<-var_sp3[var_sp3$cyto_favorable=='yes' & !is.na(var_sp3$cyto_favorable),]
km_rem_fit_fav_cyto<-survfit(Surv(rem_time,relapse)~trt+cluster,data=rem_fav_cyto)
plot_rem_fav_cyto<-ggsurvplot(km_rem_fit_fav_cyto,data=rem_fav_cyto,pal=mycolors2
                              ,legend.title=" ",legend.labs=c('C1 CC','C2 CC','C3 CC','C4 CC','C5 CC')
                              ,xlim = c(0,10),break.x.by=1,conf.int=F,xlab='Time (years)',linetype=rep('twodash',5)
                              ,title="RD Fav. Cyto Risk CC"
                              ,font.main = 12,font.x = 10,font.y = 10,ylab='Remission probability',pval=T
                              ,pval.size=3,pval.coord=c(7.5,0.95),risk.table=T,risk.table.height = 0.4,tables.col = "strata"
                              ,ggtheme=theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05),size = 0.5,fontsize =2.5,censor.size=2)
plot_rem_fav_cyto$plot<-plot_rem_fav_cyto$plot+theme(legend.text=element_text(size=8), legend.key.size = unit(0.4, "cm"))
plot_rem_fav_cyto$table<-plot_rem_fav_cyto$table+
  theme(legend.text=element_text(size=8),legend.key.size=unit(0.4,"cm"),axis.text=element_text(size=8)
        ,legend.position='none',axis.title.x=element_text(size=10),plot.title=element_text(size=10))
rem_var_list2[[16]]<-plot_rem_fav_cyto

var_split_plots2<-arrange_ggsurvplots(list(
  sv_var_list2[[7]],sv_var_list2[[8]],sv_var_list2[[9]],sv_var_list2[[10]],sv_var_list2[[11]],sv_var_list2[[12]]
  ,sv_var_list2[[13]],sv_var_list2[[19]],sv_var_list2[[14]],sv_var_list2[[15]],sv_var_list2[[16]],sv_var_list2[[17]]
  ,sv_var_list2[[18]],sv_var_list2[[1]],sv_var_list2[[2]],sv_var_list2[[3]],sv_var_list2[[4]],sv_var_list2[[5]]
  ,sv_var_list2[[6]]
  
  ,rem_var_list2[[6]],rem_var_list2[[7]],rem_var_list2[[8]],rem_var_list2[[9]] ,rem_var_list2[[10]],rem_var_list2[[11]]
  ,rem_var_list2[[16]],rem_var_list2[[12]],rem_var_list2[[13]],rem_var_list2[[14]],rem_var_list2[[15]]
  ,rem_var_list2[[1]],rem_var_list2[[2]],rem_var_list2[[3]],rem_var_list2[[4]])
  ,nrow=2,ncol=3, print=F)

setwd('/results')#custom destination folder for file generation
ggsave('Supplementary_figure_S5.pdf',var_split_plots2, width=11.69,height=8.27,units='in')

####Table 2 and Supplementary Table S7: Cox proportional hazards analysis###
##Cox proportional hazards analysis grouped with new levels##
#Grouped by cluster and treatment (3 groups)
clin_merge2$cluster3<-factor(ifelse(clin_merge2$cluster2=='C1 VH'|
                                      clin_merge2$cluster2=='C2 CC'|
                                      clin_merge2$cluster2=='C4 CC','Group1'
                                    ,ifelse(clin_merge2$cluster2=='C1 CC'|
                                              clin_merge2$cluster2=='C2 VH'|
                                              clin_merge2$cluster2=='C3 CC','Group2','Group3'))
                             ,levels=c('Group1','Group2','Group3'))

#Expanded Univariate analysis (Supplementary Table S7)#
hma_vt_grouped_surv_uv2<-tbl_uvregression(
  clin_merge2[c('surv_time','status','cluster3','age','gender','race_white','race_black','race_hisp'
                ,'race_asian','race_other','wbc','blasts','hgb','plt','b2m','aml_gp'
                ,'complex_kar','diploid','del5_5q','del7_7q','t_8_21','t_9_11','t_6_9','t_11q23'
                ,'inv16','del12','trisomy8','trisomy6','non_t_9_11_t_11q23','t_9_22_Ph_pos'
                ,'mut_aslx1','mut_cebpa','mut_dnmt3','mut_ezh2','mut_flt3','mut_flt3_d835','mut_flt3_itd'
                ,'mut_gata2','mut_idh','mut_idh1','mut_idh2','mut_jak2','mut_kit','mut_kmt2a'
                ,'mut_mll','mut_npm1','mut_ptpn11','mut_ras','mut_runx1','mut_srsf2','mut_tet2','mut_tp53'
                ,'mut_u2af','mut_wt1'
                
  )]
  ,label=list(cluster3~"Cluster",age~'Age (years)',gender~"Gender",aml_gp~"AML group"
              ,race_white~'White (race)',race_black~'Black (race)',race_hisp~'Hispanic (race)'
              ,race_asian~'Asian (race)',race_other~'Other (race)'
              ,wbc~'White Blood Cell count (K/uL)',blasts~'Blasts (%)',hgb~'Hemoglobin (g/dl)'
              ,plt~'Platelets (K/uL)',b2m~'Serum B2M (ug/mL)'
              ,complex_kar~'Complex Karyotype',diploid~'Diploid Karyotype',t_9_22_Ph_pos~'t(9;22)Ph+'
              ,del5_5q~'-5/5q-',del7_7q~'-7/7q-',t_8_21~'t(8;21)'
              ,t_9_11~'t(9;11)',t_6_9~'t(6;9)',t_11q23~'t(11q23)',inv16~'Inv 16',del12~'Del 12'
              ,trisomy6~'Trisomy 6',trisomy8~'Trisomy 8',non_t_9_11_t_11q23~'non-t(9;11)/t(11q23)'
              ,mut_aslx1~"ASLX1 Mut",mut_cebpa~"CEBPA Mut",mut_dnmt3~"DNMT3 Mut",mut_ezh2~'EZH2 Mut'
              ,mut_flt3~'FLT3 Mut',mut_flt3_d835~"FLT3_D835 Mut",mut_flt3_itd~"FLT3_ITD Mut",mut_gata2~"GATA2 Mut"
              ,mut_idh~"IDH Mut",mut_idh1~"IDH1 Mut",mut_idh2~"IDH2 Mut",mut_jak2~"JAK2 Mut",mut_kit~"KIT Mut"
              ,mut_kmt2a~"KMT2A Mut",mut_mll~"MLL Mut",mut_npm1~"NPM1 Mut",mut_ptpn11~"PTPN11 Mut"
              ,mut_ras~"RAS Mut",mut_runx1~"RUNX1 Mut",mut_srsf2~"SRSF2 Mut",mut_tet2~"TET2 Mut"
              ,mut_tp53~"TP53 Mut",mut_u2af~"U2AF Mut",mut_wt1~"WT1 Mut"
  )
  ,add_estimate_to_reference_rows=T
  ,method=coxph,y=Surv(surv_time,status),exponentiate=T,hide_n=T,pvalue_fun=~style_pvalue(.x,digits=2))%>%bold_p()

hma_vt_grouped_rem_uv2<-tbl_uvregression(
  clin_merge2[c('rem_time','relapse','cluster3','age','gender','race_white','race_black','race_hisp'
                ,'race_asian','race_other','wbc','blasts','hgb','plt','b2m','aml_gp'
                ,'complex_kar','diploid','del5_5q','del7_7q','t_8_21','t_9_11','t_6_9','t_11q23'
                ,'inv16','del12','trisomy8','trisomy6','non_t_9_11_t_11q23','t_9_22_Ph_pos'
                ,'mut_aslx1','mut_cebpa','mut_dnmt3','mut_ezh2','mut_flt3','mut_flt3_d835','mut_flt3_itd'
                ,'mut_gata2','mut_idh','mut_idh1','mut_idh2','mut_jak2','mut_kit','mut_kmt2a'
                ,'mut_mll','mut_npm1','mut_ptpn11','mut_ras','mut_runx1','mut_srsf2','mut_tet2','mut_tp53'
                ,'mut_u2af','mut_wt1'
  )]
  ,label=list(cluster3~"Cluster",age~'Age (years)',gender~"Gender",aml_gp~"AML group"
              ,race_white~'White (race)',race_black~'Black (race)',race_hisp~'Hispanic (race)'
              ,race_asian~'Asian (race)',race_other~'Other (race)'
              ,wbc~'White Blood Cell count (K/uL)',blasts~'Blasts (%)',hgb~'Hemoglobin (g/dl)'
              ,plt~'Platelets (K/uL)',b2m~'Serum B2M (ug/mL)'
              ,complex_kar~'Complex Karyotype',diploid~'Diploid Karyotype',t_9_22_Ph_pos~'t(9;22)Ph+'
              ,del5_5q~'-5/5q-',del7_7q~'-7/7q-',t_8_21~'t(8;21)'
              ,t_9_11~'t(9;11)',t_6_9~'t(6;9)',t_11q23~'t(11q23)',inv16~'Inv 16',del12~'Del 12'
              ,trisomy6~'Trisomy 6',trisomy8~'Trisomy 8',non_t_9_11_t_11q23~'non-t(9;11)/t(11q23)'
              ,mut_aslx1~"ASLX1 Mut",mut_cebpa~"CEBPA Mut",mut_dnmt3~"DNMT3 Mut",mut_ezh2~'EZH2 Mut'
              ,mut_flt3~'FLT3 Mut',mut_flt3_d835~"FLT3_D835 Mut",mut_flt3_itd~"FLT3_ITD Mut",mut_gata2~"GATA2 Mut"
              ,mut_idh~"IDH Mut",mut_idh1~"IDH1 Mut",mut_idh2~"IDH2 Mut",mut_jak2~"JAK2 Mut",mut_kit~"KIT Mut"
              ,mut_kmt2a~"KMT2A Mut",mut_mll~"MLL Mut",mut_npm1~"NPM1 Mut",mut_ptpn11~"PTPN11 Mut"
              ,mut_ras~"RAS Mut",mut_runx1~"RUNX1 Mut",mut_srsf2~"SRSF2 Mut",mut_tet2~"TET2 Mut"
              ,mut_tp53~"TP53 Mut",mut_u2af~"U2AF Mut",mut_wt1~"WT1 Mut"
  )
  ,add_estimate_to_reference_rows=T
  ,method=coxph,y=Surv(rem_time,relapse),exponentiate=T,hide_n=T,pvalue_fun=~style_pvalue(.x,digits=2))%>%bold_p()

exp_uv<-tbl_merge(list(hma_vt_grouped_surv_uv2,hma_vt_grouped_rem_uv2)
          ,tab_spanner=c("**Univariate OS<br>(N=419)**","**Univariate CRD<br>(N=274)**"))%>%
  bold_labels()%>%
  modify_header(label~"**Variable**")%>%
  modify_caption("**Supplementary Table S7.** Expanded Univariate Cox Proportional Hazards of 
  Overall Survival (OS) and Complete Remission Duration (CRD)")

setwd('/results')#custom destination folder for file generation
gtsave(data=as_gt(exp_uv),filename="Supplementary_table_S7.html")
#Open the saved .html file with excel and then save as .xlsx

#UV and MV analysis (Table 2)#
#Significant Univariate analysis#
merge_sv_uv3<-tbl_uvregression(
  clin_merge2[c('surv_time','status','cluster3','age','race_white'
                ,'race_asian','race_black','aml_gp','blasts','hgb','b2m'
                ,'complex_kar','del5_5q','del7_7q','t_8_21'
                ,'inv16','del12'
                ,'mut_aslx1','mut_cebpa','mut_flt3','mut_flt3_d835','mut_flt3_itd'
                ,'mut_idh2','mut_jak2','mut_mll','mut_npm1','mut_ptpn11','mut_runx1','mut_tp53'
  )]
  ,label=list(cluster3~'Cluster membership',age~'Age (years)',aml_gp~"AML group"
              ,race_white~'White (race)',race_asian~'Asian (race)',race_black~'Black (race)'
              ,blasts~'Blasts (%)',hgb~'Hemoglobin (g/dL)',b2m~'Serum B2M (mug/mL)'
              ,complex_kar~'Complex Karyotype',del5_5q~'-5/5q-',del7_7q~'-7/7q-'
              ,t_8_21~'t(8;21)',inv16~'Inv 16',del12~'Del 12'
              ,mut_aslx1~"ASLX1 Mutation",mut_cebpa~"CEBPA Mutation",mut_flt3~'FLT3 Mutation'
              ,mut_flt3_d835~"FLT3 D835 Mutation",mut_flt3_itd~"FLT3 ITD Mutation"
              ,mut_idh2~"IDH2 Mutation",mut_jak2~"JAK2 Mutation",mut_mll~"MLL Mutation"
              ,mut_npm1~"NPM1 Mutation",mut_ptpn11~"PTPN11 Mutation",mut_runx1~"RUNX1 Mutation"
              ,mut_tp53~"TP53 Mutation")
  ,show_single_row=c('race_white'
                     ,'race_asian','race_black','aml_gp','blasts','hgb','b2m'
                     ,'complex_kar','del5_5q','del7_7q','t_8_21'
                     ,'inv16','del12'
                     ,'mut_aslx1','mut_cebpa','mut_flt3','mut_flt3_d835','mut_flt3_itd'
                     ,'mut_idh2','mut_jak2','mut_mll','mut_npm1','mut_ptpn11','mut_runx1','mut_tp53')
  ,add_estimate_to_reference_rows=T
  ,method=coxph,y=Surv(surv_time,status),exponentiate=T,hide_n=T,pvalue_fun=~style_pvalue(.x,digits=2))%>%bold_p()

merge_rem_uv3<-tbl_uvregression(
  clin_merge2[c('rem_time','relapse','cluster3','age','race_white'
                ,'race_asian','race_black','aml_gp','blasts','hgb','b2m'
                ,'complex_kar','del5_5q','del7_7q','t_8_21'
                ,'inv16','del12'
                ,'mut_aslx1','mut_cebpa','mut_flt3','mut_flt3_d835','mut_flt3_itd'
                ,'mut_idh2','mut_jak2','mut_mll','mut_npm1','mut_ptpn11','mut_runx1','mut_tp53'
  )]
  ,label=list(cluster3~'Cluster membership',age~'Age (years)',aml_gp~"AML group"
              ,race_white~'White (race)',race_asian~'Asian (race)',race_black~'Black (race)'
              ,blasts~'Blasts (%)',hgb~'Hemoglobin (g/dL)',b2m~'Serum B2M (ug/mL)'
              ,complex_kar~'Complex Karyotype',del5_5q~'-5/5q-',del7_7q~'-7/7q-'
              ,t_8_21~'t(8;21)',inv16~'Inv 16',del12~'Del 12'
              ,mut_aslx1~"ASLX1 Mutation",mut_cebpa~"CEBPA Mutation",mut_flt3~'FLT3 Mutation'
              ,mut_flt3_d835~"FLT3 D835 Mutation",mut_flt3_itd~"FLT3 ITD Mutation"
              ,mut_idh2~"IDH2 Mutation",mut_jak2~"JAK2 Mutation",mut_mll~"MLL Mutation"
              ,mut_npm1~"NPM1 Mutation",mut_ptpn11~"PTPN11 Mutation",mut_runx1~"RUNX1 Mutation"
              ,mut_tp53~"TP53 Mutation")
  ,show_single_row=c('race_white'
                     ,'race_asian','race_black','aml_gp','blasts','hgb','b2m'
                     ,'complex_kar','del5_5q','del7_7q','t_8_21'
                     ,'inv16','del12'
                     ,'mut_aslx1','mut_cebpa','mut_flt3','mut_flt3_d835','mut_flt3_itd'
                     ,'mut_idh2','mut_jak2','mut_mll','mut_npm1','mut_ptpn11','mut_runx1','mut_tp53')
  ,add_estimate_to_reference_rows=T
  ,method=coxph,y=Surv(rem_time,relapse),exponentiate=T,hide_n=T,pvalue_fun=~style_pvalue(.x,digits=2))%>%bold_p()

#Multivariate analysis#
merge_sv_mv3<-coxph(Surv(surv_time,status)~cluster3+age+race_asian+race_white
                    +aml_gp+complex_kar+del12+del5_5q+del7_7q+inv16+t_8_21
                    +mut_aslx1+mut_cebpa+mut_flt3+mut_idh2+mut_jak2+mut_npm1
                    +mut_ptpn11+mut_tp53,data = clin_merge2) %>%
  tbl_regression(label=list(cluster3~'Cluster membership',age~'Age (years)',aml_gp~'AML group'
                            ,complex_kar~"Complex Karyotype",race_white~'White (race)',race_asian~'Asian (race)'
                            ,del7_7q~'-7/7q-',t_8_21~'t(8;21)',inv16~'Inv 16',del12~'Del 12',del5_5q~'-5/5q-'
                            ,mut_flt3~'FLT3 Mutation',mut_aslx1~"ASLX1 Mutation",mut_tp53~"TP53 Mutation",mut_cebpa~"CEBPA Mutation"
                            ,mut_jak2~"JAK2 Mutation",mut_npm1~"NPM1 Mutation",mut_ptpn11~"PTPN11 Mutation",mut_idh2~"IDH2 Mutation")
                 ,show_single_row=c(race_asian,race_white,aml_gp,complex_kar,del12,del5_5q,del7_7q
                                    ,inv16,t_8_21,mut_aslx1,mut_cebpa,mut_flt3,mut_idh2,mut_jak2
                                    ,mut_npm1,mut_ptpn11,mut_tp53)
                 ,add_estimate_to_reference_rows=T,exponentiate=T,pvalue_fun=~style_pvalue(.x,digits=2))%>% bold_p()

merge_rem_mv3<-coxph(Surv(rem_time,relapse)~cluster3+age+race_black+aml_gp+complex_kar
                     +inv16+del5_5q+mut_flt3+mut_tp53
                     ,data = clin_merge2) %>%
  tbl_regression(label=list(cluster3~'Cluster membership',age~'Age (years)',aml_gp~'AML group'
                            ,race_black~'Black (race)',complex_kar~"Complex Karyotype",inv16~'Inv 16'
                            ,del5_5q~'-5/5q-',mut_flt3~'FLT3 Mutation',mut_tp53~"TP53 Mutation")
                 ,show_single_row=c(race_black,aml_gp,complex_kar,inv16,del5_5q,mut_flt3,mut_tp53)
                 ,add_estimate_to_reference_rows=T,exponentiate=T,pvalue_fun=~style_pvalue(.x,digits=2))%>% bold_p()

sig_uv_mv<-tbl_merge(list(merge_sv_uv3,merge_sv_mv3,merge_rem_uv3,merge_rem_mv3)
          ,tab_spanner=c("**Univariate OS<br>(N=419)**","**Multivariate OS<br>(N=419)**"
                         ,"**Univariate CRD<br>(N=274)**","**Multivariate CRD<br>(N=274)**"))%>%
  bold_labels()%>%
  modify_header(label~"**Variable**")%>%
  modify_caption("**Table 2.** Significant Univariate and Multivariate Cox Proportional Hazards of 
  Overall Survival (OS) and Complete Remission Duration (CRD)")

setwd('/results')#custom destination folder for file generation
gtsave(data=as_gt(sig_uv_mv),filename="Table2.html")
#Open the saved .html file and then save as .pdf

####Figure 3: Protein classifier Set  analysis###
#Export data for AI analysis
merge_ptn2<-merge_ptn
merge_ptn2$right_trt<-factor(ifelse(merge_ptn2$cluster=='C1','VH'
                                    ,ifelse(merge_ptn2$cluster=='C2','CC'
                                            ,ifelse(merge_ptn2$cluster=='C3','CC'
                                                    ,ifelse(merge_ptn2$cluster=='C4','CC','undefined'))))
                             ,levels=c('CC','VH','undefined'))
merge_ptn2<-merge_ptn2 %>% select('id','cluster','trt','cluster2','right_trt', everything())
colnames(merge_ptn2)[1:5]<-c('ID','Cluster','Treatment','Cluster2','Correct Treatment')

ptn_ai<-createWorkbook()
ptn_ai_sh<-createSheet(ptn_ai," ")
cell_style2<-CellStyle(ptn_ai)+Font(ptn_ai,isBold=TRUE)+Border()
addDataFrame(merge_ptn2,sheet=ptn_ai_sh,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style2)
saveWorkbook(ptn_ai,"input_ai_analysis.xlsx")

#Heatmap for Protein Classifier (Figure 3B)#
#Protein names imported from AI analysis results
merge_cl<-merge_ptn[,c('TGM2','NOTCH1.cle','NFE2L2','SOX2','NDUFB4','DUSP4','SMAD2.pS245_250_255'
                        ,'EZH2','ASH2L','SPI1','MAPK14.pT180_Y182','EIF4E.pS209','EIF4EBP1.pS65'
                        ,'RAD51','cluster')]

#Merge clusters dataframe with clinical data for annotation
merge_cl<- tibble::rownames_to_column(merge_cl, "id")
clin_merge_subset<-clin_merge2[,c('id','trt','aml_gp','cyto_risk','mut_cebpa','mut_gata2','mut_flt3','mut_npm1')]

merge_cl<-merge(merge_cl,clin_merge_subset,by='id')
merge_cl$cyto_risk<-factor(merge_cl$cyto_risk,levels=c('favorable','intermediary','unfavorable',NA)
                           ,labels = c('favorable','intermediary','unfavorable','N/A'), exclude = NULL)
selected_var<-c('mut_cebpa','mut_gata2','mut_flt3','mut_npm1')
merge_cl[selected_var]<-lapply(merge_cl[selected_var], function(x) factor(x,levels=c('no','yes',NA)
                                                                          ,labels = c('no','yes','N/A'), exclude = NULL))

#Reorder clustering object
merge_cl<-merge_cl[order(merge_cl[,'cluster']),]

#Create column annotation object
merge_ann<-as.data.frame(merge_cl[,c('id','cluster','trt','aml_gp','cyto_risk','mut_cebpa','mut_gata2','mut_flt3','mut_npm1')])
row.names(merge_ann)<-c(1:nrow(merge_ann))
merge_ann<-as.data.frame(merge_ann[,-1])
colnames(merge_ann)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk'
                       ,'CEBPA Mutation','GATA2 Mutation','FLT3 Mutation','NPM1 Mutation')

#Clean clustering object
row.names(merge_cl)<-c(1:nrow(merge_cl))
merge_cl<-as.data.frame(merge_cl[,2:15])
merge_cl<-as.data.frame(t(merge_cl))
merge_cl<-mutate_all(merge_cl, function(x) as.numeric(as.character(x)))

#Create color object
merge_cl_colors<-list()
merge_cl_colors[[1]]<-mycolors1[c(1,4:7)]
names(merge_cl_colors[[1]])<-unique(merge_ann$'Cluster')
merge_cl_colors[[2]]<-c(pal_jco()(2))
names(merge_cl_colors[[2]])<-unique(merge_ann$'Treatment')
merge_cl_colors[[3]]<-c(pal_jco()(3))
names(merge_cl_colors[[3]])<-unique(merge_ann$'AML Group')
merge_cl_colors[[4]]<-c(mycolors1[1],pal_jco()(2)[2:1],pal_jco()(3)[3])
names(merge_cl_colors[[4]])<-c('favorable','intermediary','unfavorable','N/A')
for (i in 5:ncol(merge_ann)){
  merge_cl_colors[[i]]<-c(pal_jco()(3))
  names(merge_cl_colors[[i]])<-c('yes','no','N/A')
}

names(merge_cl_colors)<-c('Cluster','Treatment','AML Group','Cytogenetic Risk'
                          ,'CEBPA Mutation','GATA2 Mutation','FLT3 Mutation','NPM1 Mutation')

#Normalize sample signal
if(min(merge_cl)>=-2 & max(merge_cl)<=2){
  merge_break=seq(-2.4,2.4,0.4)
}else if(min(merge_cl)>=-2){
  merge_break=c(seq(-2.4,2,0.4),max(merge_cl))
}else if(max(merge_cl)<=2){
  merge_break=c(min(merge_cl),seq(-2,2.4,0.4))
}else{
  merge_break=c(min(merge_cl),seq(-2,2,0.4),max(merge_cl))
}

ai_ht<-pheatmap(merge_cl,annotation_colors=merge_cl_colors,main="Protein Classifier system (N=14) by clusters",cluster_rows=F,treeheight_row=0, cluster_cols=F,clustering_method ='ward.D2',annotation_col=merge_ann,fontsize=11,col=jet.colors(length(merge_break)-1),breaks=merge_break,scale="none",show_colnames=F,fontsize_row=10,border_color=F,annotation_names_row=F)

setwd('/results')#custom destination folder for file generation
ggsave('Figure_3B.pdf',ai_ht, height=8.27,width=11.69,units='in')
#Legends were adjusted using adobe illustrator

####Figure 4: Differential expression analysis###
##Define Defferentially expressed proteins (DEP) for every cluster##
#Calculate p-values, adjusted p-values, and Log2-Fold-changes of each cluster vs all the others
merge_ptn3<-merge_ptn[,c(-1,-(ncol(merge_ptn)-1):-ncol(merge_ptn))]
dep_df<-data.frame(matrix(NA,nrow=(ncol(merge_ptn3)-1),ncol=15))
rownames(dep_df)<-colnames(merge_ptn3[1:(ncol(merge_ptn3)-1)])
colnames(dep_df)<-c(paste0('p_val_C',c(1:nlevels(merge_ptn3$cluster)))
                    ,paste0('padj_C',c(1:nlevels(merge_ptn3$cluster)))
                    ,paste0('LFC_C',c(1:nlevels(merge_ptn3$cluster))))
for (i in 1:(nlevels(merge_ptn3$cluster))) {
  for (j in 1:(ncol(merge_ptn3)-1)) {
    ptn_df<-merge_ptn3
    row.names(ptn_df)<-NULL
    ptn_df$cluster<-factor(ifelse(ptn_df$cluster==paste('C',as.character(i),sep='')
                                  ,paste('C',as.character(i),sep=''),'rest')
                           ,levels=c(paste('C',as.character(i),sep=''),'rest'))
    dep_df[[j,i]]<-wilcox_test(as.formula(paste(colnames(ptn_df[j]),'cluster',sep="~")),data=ptn_df)$p
    dep_df[[j,i+10]]<-(mean(ptn_df[ptn_df$cluster==paste('C',as.character(i),sep=''),j]))-(mean(ptn_df[ptn_df$cluster=='rest',j]))
  }
  dep_df[i+5]<-p.adjust(dep_df[[i]],"fdr",n=length(dep_df[[i]]))
}

#Determine up- and down-regulated DE proteins in every cluster
up_list<-list()
down_list<-list()
for (i in 1:(nlevels(merge_ptn3$cluster))) {
  up_ptn<-row.names(dep_df[dep_df[[i+5]]<(0.05) & dep_df[[i+10]]>(0.5),])
  down_ptn<-row.names(dep_df[dep_df[[i+5]]<(0.05) & dep_df[[i+10]]<(-0.5),])
  up_list[[i]]<-up_ptn
  down_list[[i]]<-down_ptn
}
names(up_list)<-c(paste0('up-regulated C',c(1:(nlevels(merge_ptn3$cluster)))))
names(down_list)<-c(paste0('down-regulated C',c(1:(nlevels(merge_ptn3$cluster)))))

dep_list<-list()
for (i in 1:(nlevels(merge_ptn3$cluster))) {
  dep_list[[i]]<-c(up_list[[i]],down_list[[i]])
}
names(dep_list)<-c(paste0('DE proteins C',c(1:(nlevels(merge_ptn3$cluster)))))

##Supplementary table S9: p-values, adjusted p-values and Log2 Fold change values for each comparison
dep_df_adj<-tibble::rownames_to_column(dep_df,"Protein")
colnames(dep_df_adj)<-c("Protein",paste0('p-value C',c(1:nlevels(merge_ptn3$cluster)),' vs rest')
                        ,paste0('adj. p-value C',c(1:nlevels(merge_ptn3$cluster)),' vs rest')
                        ,paste0('Log2-Fold-change C',c(1:nlevels(merge_ptn3$cluster)),' vs rest'))

sup_s9<-createWorkbook()
sh1_sup_s9<-createSheet(sup_s9," ")
cell_style2<-CellStyle(sup_s9)+Font(sup_s9,isBold=TRUE)+Border()
addDataFrame(data.frame("Supplementary Table S9. p-values, adjusted p-values and Log2 Fold change values for each comparison"=double()
                        ,check.names=FALSE),sheet=sh1_sup_s9,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(dep_df_adj,sheet=sh1_sup_s9,startColumn=1,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
setwd('/results')#custom destination folder for file generation
saveWorkbook(sup_s9,"Supplementary_table_S9.xlsx")

##Supplementary table S10: up- and down-regulated proteins in each cluster
sup_s10<-createWorkbook()
sh1_sup_s10<-createSheet(sup_s10," ")
cell_style2<-CellStyle(sup_s10)+Font(sup_s10,isBold=TRUE)+Border()
addDataFrame(data.frame("Supplementary Table S10. Differential expression proteins and their directionality (up- or down-regulated) for each cluster"=double()
                        ,check.names=FALSE),sheet=sh1_sup_s10,startColumn=1,startRow=1,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(dep_list[1],sheet=sh1_sup_s10,startColumn=1,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(up_list[1],sheet=sh1_sup_s10,startColumn=2,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(down_list[1],sheet=sh1_sup_s10,startColumn=3,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(dep_list[2],sheet=sh1_sup_s10,startColumn=4,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(up_list[2],sheet=sh1_sup_s10,startColumn=5,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(down_list[2],sheet=sh1_sup_s10,startColumn=6,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(dep_list[3],sheet=sh1_sup_s10,startColumn=7,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(up_list[3],sheet=sh1_sup_s10,startColumn=8,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(down_list[3],sheet=sh1_sup_s10,startColumn=9,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(dep_list[4],sheet=sh1_sup_s10,startColumn=10,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(up_list[4],sheet=sh1_sup_s10,startColumn=11,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(down_list[4],sheet=sh1_sup_s10,startColumn=12,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(dep_list[5],sheet=sh1_sup_s10,startColumn=13,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(up_list[5],sheet=sh1_sup_s10,startColumn=14,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
addDataFrame(down_list[5],sheet=sh1_sup_s10,startColumn=15,startRow=2,row.names=FALSE,colnamesStyle=cell_style2)
setwd('/results')#custom destination folder for file generation
saveWorkbook(sup_s10,"Supplementary_table_S10.xlsx")

#Figure 4A: Log2-fold-change by cluster of DE proteins of cluster C5 
#Heatmap with DEP C5 LFC for all clusters
#Prepare dataframe
lfc_dep<-dep_df[c('ZAP70','CD4','NOTCH1.cle','ERN1','FZR1','CDKN2A','LEF1','ETS1','S100A4'
                  ,'CDKN1B','VIM','HSPB1.pS82','EIF4E.pS209','HSF1.pS326','RPS6.pS240_244'
                  ,'MAPK14.pT180_Y182','PRMT1','KIT','BIRC5','CHEK1','CCNB1','EIF2AK2'
                  ,'CDK1_2_3.pT14','RB1.pS807_811'),11:15]
colnames(lfc_dep)<-c(1:5)

#Create Annotation object
lfc_dep_ann<-as.data.frame(paste0('C',c(1:5),' vs rest'))
colnames(lfc_dep_ann)<-'LFC by cluster'

#Create color object
lfc_dep_colors<-list()
lfc_dep_colors[[1]]<-mycolors1[c(1,4:7)]
names(lfc_dep_colors[[1]])<-paste0('C',c(1:5),' vs rest')
names(lfc_dep_colors)<-'LFC by cluster'

#Normalize sample signal
if(min(lfc_dep)>=-2 & max(lfc_dep)<=2){
  lfc_dep_break=seq(-2.4,2.4,0.4)
}else if(min(lfc_dep)>=-2){
  lfc_dep_break=c(seq(-2.4,2,0.4),max(lfc_dep))
}else if(max(lfc_dep)<=2){
  lfc_dep_break=c(min(lfc_dep),seq(-2,2.4,0.4))
}else{
  lfc_dep_break=c(min(lfc_dep),seq(-2,2,0.4),max(lfc_dep))
}

lfc_dep_ht<-pheatmap(lfc_dep,main="Log2−fold change (LFC) of Cluster C5 DE Proteins"
                     ,fontsize_row=10,show_colnames=F,annotation_col=lfc_dep_ann,annotation_colors=lfc_dep_colors
                     ,cluster_rows=F,cluster_cols=F,fontsize=12,treeheight_row=0,annotation_names_row=F
                     ,col=jet.colors(length(lfc_dep_break)-1),breaks=lfc_dep_break,scale="none")

setwd('/results')#custom destination folder for file generation
ggsave('Figure_4A.pdf',lfc_dep_ht,height=11.69,width=8.27,units='in')

##Figure 4B: DE Proteins C5 Protein Networks with mean expression and up- and down-regulated identified###
#Create DEP C5 dataset, calculate mean expression and identify up and down proteins
dep_c5_ptn<-as.data.frame(apply(merge_ptn3[merge_ptn3$cluster=='C5',c(up_list[[5]],down_list[[5]])],2,function(x){mean(x)}))
up_down<-c(rep(1,13),rep(2,11))
dep_c5_ptn<-as.data.frame(cbind(dep_c5_ptn,up_down))
dep_c5_ptn<-dep_c5_ptn[sort(row.names(dep_c5_ptn)),]
colnames(dep_c5_ptn)<-c('mean_exp','status')

#Adjust protein names for importing dataframes into StringApp 
dep_c5_ptn_adj<-dep_c5_ptn
dep_c5_ptn_adj$rppa_ptn<-row.names(dep_c5_ptn_adj)
row.names(dep_c5_ptn_adj)<-gsub("CDK1_2_3.pT14","CDK1",row.names(dep_c5_ptn_adj))
dep_c5_ptn_adj$rppa_ptn<-gsub("CDK1_2_3.pT14","CDK1.pT14",dep_c5_ptn_adj$rppa_ptn)
row.names(dep_c5_ptn_adj)<-gsub("EIF4E.pS209","EIF4E",row.names(dep_c5_ptn_adj))
row.names(dep_c5_ptn_adj)<-gsub("HSF1.pS326","HSF1",row.names(dep_c5_ptn_adj))
row.names(dep_c5_ptn_adj)<-gsub("HSPB1.pS82","HSPB1",row.names(dep_c5_ptn_adj))
row.names(dep_c5_ptn_adj)<-gsub("MAPK14.pT180_Y182","MAPK14",row.names(dep_c5_ptn_adj))
row.names(dep_c5_ptn_adj)<-gsub("NOTCH1.cle","NOTCH1",row.names(dep_c5_ptn_adj))
row.names(dep_c5_ptn_adj)<-gsub("RB1.pS807_811","RB1",row.names(dep_c5_ptn_adj))
row.names(dep_c5_ptn_adj)<-gsub("RPS6.pS240_244","RPS6",row.names(dep_c5_ptn_adj))
dep_c5_ptn_adj_subset<-dep_c5_ptn_adj['CDK1',] 
dep_c5_ptn_adj_subset1<-dep_c5_ptn_adj_subset
rownames(dep_c5_ptn_adj_subset1)<-c('CDK2')
dep_c5_ptn_adj_subset1$rppa_ptn<-c('CDK2.pT14')
dep_c5_ptn_adj_subset2<-dep_c5_ptn_adj_subset
rownames(dep_c5_ptn_adj_subset2)<-c('CDK3')
dep_c5_ptn_adj_subset2$rppa_ptn<-c('CDK3.pT14')
dep_c5_ptn_adj<-rbind(dep_c5_ptn_adj,dep_c5_ptn_adj_subset1,dep_c5_ptn_adj_subset2)
dep_c5_ptn_adj<-dep_c5_ptn_adj[order(row.names(dep_c5_ptn_adj)),]
dim(dep_c5_ptn_adj)

#######Open Cytoscape software v3.10.1###
##Connect Cytoscape
cytoscapePing()
cytoscapeVersionInfo()

#Close previous cytoscape session if open
closeSession(save.before.closing=F)

##Create Protein network in Cytoscape##
#Remove all PTM naming and repetead protein names to import for STRING analysis#
#cat(row.names(dep_c5_ptn_adj),sep=',') #this gives the list of protein already formatted for input
strg_cmd<-'string protein query query="BIRC5,CCNB1,CD4,CDK1,CDK2,CDK3,CDKN1B,CDKN2A,CHEK1
,EIF2AK2,EIF4E,ERN1,ETS1,FZR1,HSF1,HSPB1,KIT,LEF1,MAPK14,NOTCH1,PRMT1,RB1,RPS6,S100A4,VIM,ZAP70" 
  species="Homo sapiens" limit=0 cutoff=0.4'
commandsRun(strg_cmd)

#Load mean expression values of each cluster
loadTableData(dep_c5_ptn_adj,table.key.column="display name") 

#Rename network
renameNetwork(noquote('dep_c5_net'))

##Adjust Visual style
mystyle<-"dataStyle"
createVisualStyle(mystyle)
setVisualStyle(mystyle)

layoutNetwork('force-directed defaultSpringLength=1000 defaultSpringCoefficient=0.000003')
setNodeShapeDefault("ellipse",mystyle)
setNodeLabelMapping('rppa_ptn',mystyle)
setNodeColorDefault("#AAAAAA",mystyle)
setNodeSizeDefault(250,mystyle)
setNodeFontSizeDefault(100,mystyle)
setEdgeLineWidthDefault(3,mystyle)
setNodeBorderWidthDefault(40,mystyle)

#Set Border color according to Protein Selector Set
up_down_val<-c(1,2)
bord_col<-c('firebrick1','deepskyblue2')
setNodeBorderColorMapping('status',up_down_val,bord_col,style.name=mystyle,mapping.type='d')

##Visualize expression data (color inside nodes)
#Normalize sample signal
if(min(dep_c5_ptn)>=-2 & max(dep_c5_ptn)<=2){
  dep_c5_ptn_node_breaks=seq(-2.4,2.4,0.4)
}else if(min(dep_c5_ptn)>=-2){
  dep_c5_ptn_node_breaks=c(seq(-2.4,2,0.4),max(dep_c5_ptn))
}else if(max(dep_c5_ptn)<=2){
  dep_c5_ptn_node_breaks=c(min(dep_c5_ptn),seq(-2,2.4,0.4))
}else{
  dep_c5_ptn_node_breaks=c(min(dep_c5_ptn),seq(-2,2,0.4),max(dep_c5_ptn))
}

#Create color object
dep_c5_ptn_node_colors<-jet.colors(length(dep_c5_ptn_node_breaks))

#Attach expression data to nodes in each network, add legend and export
setNodeColorMapping('mean_exp',dep_c5_ptn_node_breaks,dep_c5_ptn_node_colors,style.name=mystyle)

##################Adjust the space between nodes in all Networks in Cytoscape environment to your own preference###
#################Add legends in Cytoscape environment with Legend Creator and export figures###
setwd('/results') #custom destination folder for file generation
exportImage(filename='Figure_4B',type='SVG')

###Enriched processes of DE proteins for every cluster (Figure 4C)###
#Export DE protein names to Enrichr and analyze: https://maayanlab.cloud/Enrichr/
setwd('/data') #adjust world directory setting
go_top<-read_excel("Supplementary_table_S11.xlsx",col_names=c('path','overlap','pval','p_adj','odds','score','genes')) 
go_top<-go_top[c(3:22),] #select top 20 paths and remove table headings
go_top$path<-gsub("WP\\w+ *", "", go_top$path)
go_top[,3:(ncol(go_top)-1)]<-mutate_all(go_top[,3:(ncol(go_top)-1)],function(x) as.numeric(as.character(x)))
go_top$path<-factor(go_top$path,levels=unique(go_top$path[order(go_top$score,decreasing=T)]))
enrich_plot<-ggplot(go_top,aes(x=path,y=score,fill=score))+ggtitle('Top Enriched Biological Processes DEP C5')+
  geom_bar(aes(fct_rev(path),score,fill=score),color='black',position=position_dodge(width=0.7),stat="identity",linewidth=0.1)+
  labs(y='Combined Score')+labs(x='')+scale_fill_continuous(limits=c(-30000,6000))+
  theme(aspect.ratio = 1/1,legend.position='',legend.title=element_blank(),axis.title.y=element_blank()
        ,axis.text.x=element_text(size=6,face="bold",color='black',margin=margin(t=15,r=0,b=0,l=0)),axis.text.y=element_text(size=6,face="bold",color='black')
        ,plot.title = element_text(size=8,hjust =0,face='bold',color='black'),axis.title.x=element_text(size=8,face="bold")
        ,panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()
        ,axis.line=element_line(linewidth=0.25,linetype="solid",colour="black"),panel.border=element_blank()
        ,axis.ticks=element_line(linewidth=0.25),axis.ticks.length=unit(0.025, "cm"))+guides(fill=guide_legend(byrow = TRUE))+coord_flip()

setwd('/results') #custom destination folder for file generation
ggsave('Figure_4C.pdf',enrich_plot,width=10,height=6,units='in')

####Supplementary Figure S6: Volcano plots for observing DEP directionality###
#Generate volcano plots
dep_vol_list<-list()
for (i in 1:(nlevels(merge_ptn3$cluster))){
  key_col<-ifelse(dep_df[[i+5]]<(0.05) & dep_df[[i+10]]<(-0.5),'deepskyblue2'
                  ,ifelse(dep_df[[i+5]]<(0.05) & dep_df[[i+10]]>(0.5),'firebrick1','gray50'))
  names(key_col)[key_col=='firebrick1']<-'Upregulated'
  names(key_col)[key_col=='gray50']<-'NS'
  names(key_col)[key_col=='deepskyblue2']<-'Downregulated'
  
  plot_vol<-EnhancedVolcano(dep_df,lab=row.names(dep_df),cutoffLineType='twodash',cutoffLineWidth=0.25
                            ,widthConnectors=0.15,pointSize=0.5,labSize=1.5,colAlpha=1,colCustom=key_col,caption=NULL,drawConnectors=T
                            ,legendPosition='bottom',legendLabSize=8,legendIconSize=2,FCcutoff=0.5,pCutoff=0.05,subtitle=NULL
                            ,arrowheads=F,border='partial',borderWidth=0.25,borderColour='black',xlim=c(-2.5,2.5),ylim=c(0,40)
                            ,title=paste('DE Proteins',as.character(gsub('LFC_','',colnames(dep_df[i+10]))))
                            ,x=colnames(dep_df[i+10]),y=colnames(dep_df[i]),boxedLabels=T,axisLabSize=10
                            ,titleLabSize=12,gridlines.major=F, gridlines.minor=F
  )+theme(axis.ticks=element_line(linewidth=0.25),axis.ticks.length = unit(0.1, "cm"))
  dep_vol_list[[i]]<-plot_vol
}

#Adjust y-axis of graphs
dep_vol_list[[2]]<-dep_vol_list[[2]]+coord_cartesian(ylim=c(0,10))
dep_vol_list[[3]]<-dep_vol_list[[3]]+coord_cartesian(ylim=c(0,10))
dep_vol_list[[4]]<-dep_vol_list[[4]]+coord_cartesian(ylim=c(0,35))
dep_vol_list[[5]]<-dep_vol_list[[5]]+coord_cartesian(ylim=c(0,20))

setwd('/results')#custom destination folder for file generation
ggsave(file='Supplementary_figure_S6.pdf',marrangeGrob(grobs=dep_vol_list,ncol=2,nrow=3,top=NULL),width=8.27,height=11.69,units='in')
