library(readxl)
library(dplyr)
library(ggplot2)
library(assignPOP) #for assignments tests
library(inlmisc) #for more Paul Tol color options

#Colors: Paul Tol color scheme
scale_taxon<-as.character(GetColors(3,scheme="muted"))

#data
data.folder<-"D:/Dropbox/Documents/research/mammals/porcupines/short_ms/functional_mandible_analysis/"

data.raw<-read_excel(paste(data.folder,"Porcupine modern BASIC DATA cleaned.xlsx",sep="")  )                   
data.fossils<-read_excel(paste(data.folder,"fossil_mandible_functional_traits.xlsx",sep=""))                     

# functions ------
#geometric mean function from StackOverflow. Thank you @paul-mcmurdie, @ben-bolker, @Gregor.
gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

# clean data ------
data.spp<-data.raw
data.spp$SPECIMEN<-gsub("(Erethizon dorsatum).*","\\1",data.spp$SPECIMEN,perl=TRUE)
data.spp$SPECIMEN<-gsub("(Coendou bicolor).*","\\1",data.spp$SPECIMEN,perl=TRUE)
data.spp$SPECIMEN<-gsub("(Coendou).*(mexicanus).*","\\1 \\2",data.spp$SPECIMEN,perl=TRUE)
data.spp$SPECIMEN<-gsub("(Coendou).*(insidiosus).*","\\1 \\2",data.spp$SPECIMEN,perl=TRUE)
data.spp$`CATALOG NUMBER`<-gsub("(.*) (.*)","\\1\\2",data.spp$`CATALOG NUMBER`,perl=TRUE)

data.erethizon<-data.spp #sub in all spp with cleaned names, look for residuals

# #or, use 1 spp
# #largest sample size comes from E. dorsatum, so use that
# data.erethizon<-grep("dorsatum",data.raw$SPECIMEN,perl=TRUE) %>% data.raw[.,]

#remove the "-" for NA characters and change to numeric
data.erethizon.na <- apply(data.erethizon, 2,  function(x) gsub("^-",NA,x)) %>% as.data.frame
data.erethizon.na$AGE <-as.numeric(data.erethizon.na$AGE)

data.erethizon.na[,c(6:40)] <-apply(data.erethizon.na[,c(6:40)],2, function(x) as.numeric(as.character(x)))

# summary statistics incisor ------
data.erethizon.clean<-dplyr::select(as.data.frame(data.erethizon.na),"CATALOG NUMBER","SPECIMEN","SEX",
                                    "AGE","ALVEOLAR LENGTH (mm)","WIDTH m1","LENGTH m1","WIDTH m2",
                                    "LENGTH m2","WIDTH m3","LENGTH m3","AP DIAMETER INCISOR (mm)",
                                    "WIDTH INCISOR (mm)","DEPTH OF MANDIBLE AT p4 (mm)",
                                    "DEPTH OF MANDIBLE AT m2 (mm)","LENGTH OF DIASTEMA (mm)",
                                    "LENGTH OF SYMPHYSIS (mm)","ANGLE BETWEEN  BUCCAL","m1-m2 AND INCISOR LINGUAL",
                                    "DISTANCE FROM m3 TO FORAMEN (mm)","LONGITUDINAL AXIS OF MANDIBULAR TOOTH ROW")


data.erethizon.clean$incisor.prop<-data.erethizon.clean$`AP DIAMETER INCISOR (mm)`/data.erethizon.clean$`WIDTH INCISOR (mm)` 
data.erethizon.clean$genus<-"Coendou"
data.erethizon.clean$genus[which(data.erethizon.clean$SPECIMEN=="Erethizon dorsatum")]<-"Erethizon"
data.erethizon.clean %>% group_by(genus) %>% 
  summarise(avg=mean(incisor.prop, na.rm=T),std.dev=sd(incisor.prop,na.rm=T),min=min(incisor.prop,na.rm=T),
            max=max(incisor.prop,na.rm=T)) %>% write.csv(.,"incisor_proportion_summary.csv")

data.erethizon.clean %>% group_by(genus) %>% 
  summarise(avg=mean(`ANGLE BETWEEN  BUCCAL`, na.rm=T),std.dev=sd(`ANGLE BETWEEN  BUCCAL`,na.rm=T),
            min=min(`ANGLE BETWEEN  BUCCAL`,na.rm=T),max=max(`ANGLE BETWEEN  BUCCAL`,na.rm=T)) %>% 
  write.csv(.,"procumbency_proportion_summary.csv")

# # explore ontogeny --------
# # It appeared to us that the porcupine mandible continues to grow throughout life,
# #meaning some traits will vary according to age. 
# #Can they all be thrown into the same statistical basket? Can your statistics detail growth patterns?
# 
# # #leave out p4/dp4: variable presence will prevent comparisons over ontogeny
# # data.erethizon.clean<-dplyr::select(as.data.frame(data.erethizon.na),"CATALOG NUMBER","SPECIMEN","SEX",
# #                                     "AGE","ALVEOLAR LENGTH (mm)","WIDTH m1","LENGTH m1","WIDTH m2",
# #                                     "LENGTH m2","WIDTH m3","LENGTH m3","AP DIAMETER INCISOR (mm)",
# #                                     "WIDTH INCISOR (mm)","DEPTH OF MANDIBLE AT p4 (mm)",
# #                                     "DEPTH OF MANDIBLE AT m2 (mm)","LENGTH OF DIASTEMA (mm)",
# #                                     "LENGTH OF SYMPHYSIS (mm)","ANGLE BETWEEN  BUCCAL","m1-m2 AND INCISOR LINGUAL",
# #                                     "DISTANCE FROM m3 TO FORAMEN (mm)","LONGITUDINAL AXIS OF MANDIBULAR TOOTH ROW")
# 
# # #molars don't change in size over ontogeny: what you would expect since they don't grow.
# # ggplot(data=data.erethizon.clean,aes(x=AGE,y=`WIDTH m1`,color=SPECIMEN)) + geom_point()
# # ggplot(data=data.erethizon.clean,aes(x=AGE,y=`LENGTH m1`,color=SPECIMEN)) + geom_point()
# # ggplot(data=data.erethizon.clean,aes(x=AGE,y=`WIDTH m2`,color=SPECIMEN)) + geom_point()
# # ggplot(data=data.erethizon.clean,aes(x=AGE,y=`LENGTH m2`,color=SPECIMEN)) + geom_point()
# # 
# # #other measurements appear to be more size & age related
# ggplot(data=data.erethizon.clean,aes(x=AGE,y=`ALVEOLAR LENGTH (mm)`,color=SPECIMEN)) + geom_point()
# # ggplot(data=data.erethizon.clean,aes(x=AGE,y=`LENGTH OF DIASTEMA (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.clean,aes(x=AGE,y=`ANGLE BETWEEN  BUCCAL`,color=SPECIMEN)) + geom_point()
# # ggplot(data=data.erethizon.clean,aes(x=AGE,y=`DEPTH OF MANDIBLE AT m2 (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.clean,aes(x=AGE,y=`WIDTH INCISOR (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.clean,aes(x=AGE,y=`AP DIAMETER INCISOR (mm)`,color=SPECIMEN)) + geom_point()

# size correction ------
#create a smaller dataset without molar sizes, things not to recovered from important fossils
data.erethizon.size<-dplyr::select(as.data.frame(data.erethizon.na),"CATALOG NUMBER","SPECIMEN","SEX",
                                   "AGE","m1-m2 AND INCISOR LINGUAL","ANGLE BETWEEN  BUCCAL",
                                   "AP DIAMETER INCISOR (mm)","WIDTH INCISOR (mm)",
                                   "ALVEOLAR LENGTH (mm)")
#filter incomplete linear measuremetn rows: affects geometric  transform
data.erethizon.size.complete<-filter(data.erethizon.size,
                            complete.cases(data.erethizon.size[,c(7:ncol(data.erethizon.size))]))

#Calculate the geometric mean for each specimen (row) and divide all variables by the geometric mean
#of linear measurements (ignore angle measurements: will be scaled. Also, don't seem to have same)
geo.mean<-apply(data.erethizon.size.complete[,c(7:ncol(data.erethizon.size.complete))],
                1,function(x) gm_mean(x %>% unlist))
erethizon.geo<-data.erethizon.size.complete[,c(7:ncol(data.erethizon.size.complete))]/geo.mean
data.erethizon.resize<-data.erethizon.size.complete
data.erethizon.resize[,c(7:ncol(data.erethizon.size.complete))]<-erethizon.geo

# #check.
# ggplot(data=data.erethizon.resize,aes(x=AGE,y=`ALVEOLAR LENGTH (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.resize,aes(x=AGE,y=`ANGLE BETWEEN  BUCCAL`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.resize,aes(x=AGE,y=`WIDTH INCISOR (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.resize,aes(x=AGE,y=`AP DIAMETER INCISOR (mm)`,color=SPECIMEN)) + geom_point()

# # check dimorphism -------
# ggplot(data=data.erethizon.size,aes(x=SEX,y=`AP DIAMETER INCISOR (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.size,aes(x=SEX,y=`DEPTH OF MANDIBLE AT m2 (mm)`,color=SPECIMEN)) + geom_point()
# ggplot(data=data.erethizon.size,aes(x=SEX,y=`ALVEOLAR LENGTH (mm)`,color=SPECIMEN)) + geom_point()
# #at first glance looks okay
# data.sex<-data.erethizon.size %>% filter(.,!is.na(SEX),SPECIMEN=="Erethizon dorsatum",AGE!=1)
# t.test(`ALVEOLAR LENGTH (mm)`~SEX,data=data.erethizon.size)
# t.test(`DEPTH OF MANDIBLE AT m2 (mm)`~SEX,data=data.erethizon.size)

# put all data together ------------
data.erethizon.alltime<-rbind(data.erethizon.size,
                              data.fossils[which(data.fossils$`CATALOG NUMBER`=="UF.223810"),])
# can one procumbency measure stand for another? -------
# ggplot(data=data.erethizon.size,aes(x=`m1-m2 AND INCISOR LINGUAL`,y=`ANGLE BETWEEN  BUCCAL`,
#                                      color=SPECIMEN)) + geom_point(size=3)
# procumbent1<-lm(`m1-m2 AND INCISOR LINGUAL`~`ANGLE BETWEEN  BUCCAL`,data=data.erethizon.size)
# # procumbent2<-lm(log(`m1-m2 AND INCISOR LINGUAL`)~log(`ANGLE BETWEEN  BUCCAL`),data=data.erethizon.clean)
# # procumbent3<-lm(`m1-m2 AND INCISOR LINGUAL`~`ANGLE BETWEEN  BUCCAL`,
# #  data=data.erethizon.clean[grep("dorsatum",data.erethizon.clean$SPECIMEN,perl=TRUE),])
# summary(procumbent1)
# # summary(procumbent3) #unlogged angles are better, full dataset better.
# 
# #predict missing Idaho measurement
# procumbent.predict<-lm(`ANGLE BETWEEN  BUCCAL`~`m1-m2 AND INCISOR LINGUAL`,data=data.erethizon.alltime)
# idaho<-which(data.erethizon.alltime$SPECIMEN=="Fossil: Idaho")
# data.erethizon.alltime$`ANGLE BETWEEN  BUCCAL`[idaho]<-data.erethizon.alltime$`m1-m2 AND INCISOR LINGUAL`[idaho]*procumbent.predict$coefficients[2] + procumbent.predict$coefficients[1]


# #check
# ggplot(data=data.erethizon.alltime,aes(x=`m1-m2 AND INCISOR LINGUAL`,y=`ANGLE BETWEEN  BUCCAL`,
#                                     color=SPECIMEN)) + geom_point(size=3)

# prep for PCA -----
#It looks like age category 1 still retains some residual differences from other age categories
#These are juveniles, in which the mandibular tooth row is incompletely erupted [< 1 yr, Roze 2009 p.49;
#mid-winter year 1 Sutton 1972 p.59]
#Then, lingual procumbency seems to be the limiting factor. Remove that column
#then remove any specimens with missing data
data.complete<-data.erethizon.alltime %>% #filter(.,AGE!=1) %>%
  dplyr::select(.,-c(`m1-m2 AND INCISOR LINGUAL`,SEX, AGE)) %>%
  filter(complete.cases(.))
tail(data.complete)
dim(data.complete)
colnames(data.complete)

# Standardize & transform measurements for PCA. (Granatosky et al. 2014)
#Mosimann amd James '79?
# Note that you don't have to correct angular measurements
# Calculate the geometric mean for each specimen (row) and divide all variables by the geometric mean
geo.mean<-apply(data.complete[,c(4:ncol(data.complete))],
                 1,function(x) gm_mean(x %>% unlist))
erethizon.geo<-data.complete[,c(4:ncol(data.complete))]/geo.mean
  
# After geometric correction, log-transform the data in order to make sure that measurements are normally
#distributed. After size correcting and log-transforming, you are ready to go with your PCA
#log() default is natural logarithm
erethizon.geo<-log(erethizon.geo)
data.corrected<-cbind(data.complete[,c(1:3)],erethizon.geo)

# PCA --------
#perform principal components analysis
#Use correlation (scale.=TRUE) and not covariance (scale.=FALSE), don't overweight larger variables like angles
PCA.por<-prcomp(data.corrected[,c(3:ncol(data.corrected))],scale.=TRUE)
PCA.data<-cbind(data.corrected,PCA.por$x)

PCA.function<-prcomp(data.corrected[,c(3:5)],scale.=TRUE)
PCA.data<-cbind(data.corrected,PCA.function$x)

#find % of variation explained by each PC
percents<-summary(PCA.function)$importance[2,]*100

# #make figure (initial)
# # cairo_pdf(paste("example","PCA.pdf",sep="_"),width=3.3,height=2.2)
# ggplot(data=PCA.data, aes(x=PC1,y=PC2))+
#   geom_point(aes(color=SPECIMEN),alpha=0.9,size=2) +
#   # scale_shape_manual(values=c(21,22,23,24,25)) +
#   # scale_fill_manual(values=color4porcupines) +
#   xlab(paste("PC 1 (",percents[1],"%)",sep="")) +
#   ylab(paste("PC 2 (",percents[2],"%)",sep="")) +
#   theme_classic() +
#   theme(text=element_text(size=8), legend.text = element_text(face="italic"),
#         legend.background = element_rect(colour = 'black', size=0.5))
# # dev.off()

#loadings
PCA.function$rotation
# PCA.por$rotation[,2] %>% abs %>% order(.,decreasing=TRUE) %>%
#   .[c(1:5)] %>% PCA.por$rotation[.,1]

# CVA? Assignments? --------
PCA.data$ID<-"Coendou"
PCA.data$ID[which(PCA.data$SPECIMEN=="Erethizon dorsatum")]<-"Erethizon"

#If only using poyeri, then run the line below
PCA.data$ID[nrow(PCA.data)]<-"unknown"

#If using all fossils, then run the line below
# PCA.data$ID[c((nrow(PCA.data)-nrow(data.fossils)+1):nrow(PCA.data))]<-"unknown"

PCA.data$climate<-"warm"
PCA.data$climate[which(PCA.data$SPECIMEN=="Erethizon dorsatum"|
                         PCA.data$SPECIMEN=="Echinoprocta rufescens")]<-"cold"
PCA.data$climate[c((nrow(PCA.data)-nrow(data.fossils)+1):nrow(PCA.data))]<-"unknown"

known_sample<-PCA.data[which(PCA.data$ID!="unknown"),c(1,3:5,10)]
unknown_sample<-PCA.data[which(PCA.data$ID=="unknown"),c(1,3:5,10)]

#ASSIGNMENT
# assign.MC(known_sample, train.inds=c(0.1,0.25,0.5),iterations=30, model="lda",
#           dir=paste(data.folder,"Result-folder/",sep=""),scaled=TRUE)
# assign.kfold(known_sample, k.fold=c(3, 4, 5), model="lda",
#              dir=paste(data.folder,"Result-folder2/",sep=""),scaled=TRUE)

accuMC <- accuracy.MC(dir = paste(data.folder,"Result-folder/",sep="")) #Use this function for Monte-Carlo cross-validation results
accuKF <- accuracy.kfold(dir = paste(data.folder,"Result-folder2/",sep="")) #Use this function for K-fold cross-validation results

# accuracy.plot(accuKF, pop=c("all", "Coendou", "Erethizon")) +
#   ylim(0, 1) #+ #Set y limit between 0 and 1
#   # annotate("segment",x=0.4,xend=3.6,y=0.5,yend=0.5,colour="red",size=1) #Add a red horizontal line at y = 0.5 (null assignment rate for 3 populations)
#
# membership.plot(dir = paste(data.folder,"Result-folder2/",sep=""))
# assign.matrix( dir= paste(data.folder,"Result-folder2/",sep=""), train.inds=c(10))

# #ASSIGNMENT
# assign.X(x1=known_sample, x2=unknown_sample, dir=paste(data.folder,"Result-folder3/",sep=""), model="svm")
classify.fossils<-read.table(paste(data.folder,"Result-folder3/AssignmentResult.txt",sep=""),header=TRUE)

PCA.data$classification<-PCA.data$ID

#If using only poyeri, then run the line below
PCA.data$ID[nrow(PCA.data)]<-PCA.data$SPECIMEN[nrow(PCA.data)]

# #If using all fossils, then run the line below
# PCA.data$ID[c((nrow(PCA.data)-nrow(data.fossils)+1):nrow(PCA.data))]<-PCA.data$SPECIMEN[c((nrow(PCA.data)-nrow(data.fossils)+1):nrow(PCA.data))]

classify.max<-apply(classify.fossils[,c(3:4)],1, max)
for(i in 1:nrow(classify.fossils)){
  PCA.data$classification[i+nrow(known_sample)]<-c("Coendou","Erethizon")[which(classify.fossils[i,c(3:4)]==classify.max[i])]
}

PCA.data.plot<-filter(PCA.data,ID!="Fossil: Cumberland Cave")

# pdf(paste(data.folder,"../Fig2_mandible_function_PCA.pdf",sep=""),width=3.425,height=4,pointsize=9,pagecentre=FALSE) #default units inches
cairo_pdf(paste(data.folder,"Fig2_mandible_function_PCA_v03.pdf",sep=""),width=3.425,height=4)
ggplot(data=PCA.data.plot, aes(x=PC1,y=PC2))+
  geom_point(aes(color=ID,shape=classification),alpha=0.9,size=2) +
  # scale_shape_manual(values=c(21,22,23,24,25)) +
  scale_color_manual(values=scale_taxon) +
  xlab(paste("PC 1 (",round(percents[1],1),"%)",sep="")) +
  ylab(paste("PC 2 (",round(percents[2],1),"%)",sep="")) +
  theme_classic() + #geom_text(aes(label=`CATALOG NUMBER`)) +
  theme(text=element_text(size=8), legend.text = element_text(face="italic"),
        legend.position="bottom",
        #legend.background = element_rect(colour = 'black', size=0.5),
        legend.box = "vertical")
dev.off()

pick.specimens<-which(PCA.data.plot$`CATALOG NUMBER` %in% c("UF13306","UF33211"))
ggplot(data=PCA.data.plot, aes(x=PC1,y=PC2))+
  geom_point(color="black",aes(shape=classification),alpha=0.9,size=2) +
  geom_point(data=PCA.data.plot[pick.specimens,c(7,8)],color="red")+
  xlab(paste("PC 1 (",round(percents[1],1),"%)",sep="")) +
  ylab(paste("PC 2 (",round(percents[2],1),"%)",sep="")) +
  theme_classic() + #geom_text(aes(label=`CATALOG NUMBER`)) +
  theme(text=element_text(size=8), legend.text = element_text(face="italic"),
        legend.position="bottom",
        #legend.background = element_rect(colour = 'black', size=0.5),
        legend.box = "vertical")
