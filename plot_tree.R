#phylogenetic trees for short porcupine paper

# libraries, dependencies ---------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
#try just a tidyverse version:
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)

# set paths ------
# setwd("D:/Dropbox/Documents/research/mammals/porcupines/phylogeny/")
setwd("C:/Users/nsvit/Dropbox/Documents/research/mammals/porcupines/phylogeny/")

parsimony.consensus.tree<-read.nexus("07_parsimony_bb/erethizontidae_2023_05_23_paup_con.tre")
bayesian.strict.tree<-read.mrbayes("09_combined_MrBayes_upmat/allcompat/erethizontidae_total_matrix.nex.con.tre")
SI.tip.tree<-read.mrbayes("10_morph_MrBayes_upmat_tip/allcompat/erethizontidae_matrix.nex.con.tre")
SI.morph.tree<-read.mrbayes("11_morph_MrBayes_upmat/allcompat/erethizontidae_matrix.nex.con.tre")

# tree edits --------
#Define which taxa should be highlighted on the tree (based on sampling)
focal.taxa<-c("Erethizon_poyeri") 

# main plot: complete dataset, Bayesian tip-dated -----
#need to substitute "_" for " "
bayesian.strict.tree@phylo$tip.label<-gsub("(.*)_(.*)","\\1 \\2",bayesian.strict.tree@phylo$tip.label)

#need to add BPP as node labels, improve color scheme
bayesian.strict.tree@data$prob<-as.numeric(bayesian.strict.tree@data$prob)

#need to fix timescale label. Solution from Guangchuang Yu
#https://github.com/YuLab-SMU/ggtree/issues/87

#tip ranges would be nice. Use the geom_range option below

backward.tree<-ggplot(bayesian.strict.tree,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
  theme_tree2() + #hexpand(0.02,direction=1) + #scale_x_continuous(labels = abs)
  geom_range(range='height_0.95HPD',alpha=0.4,size=2) +
  geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2) #+
# scale_color_continuous(low="gray",high="black") #+
#improve the color scheme
time.scale<-seq(from=-8,to=0,by=2)
revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
ggsave("bayesian.strict.tree.all.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 

getwd()
# supplementary figures, other trees ------
#need to substitute "_" for " "
SI.morph.tree@phylo$tip.label<-gsub("(.*)_(.*)","\\1 \\2",SI.morph.tree@phylo$tip.label)
SI.tip.tree@phylo$tip.label<-gsub("(.*)_(.*)","\\1 \\2",SI.tip.tree@phylo$tip.label)

#need to add BPP as node labels, improve color scheme
SI.morph.tree@data$prob<-as.numeric(SI.morph.tree@data$prob)
SI.tip.tree@data$prob<-as.numeric(SI.tip.tree@data$prob)

#need to fix timescale label. Solution from Guangchuang Yu
#https://github.com/YuLab-SMU/ggtree/issues/87
#tip ranges would be nice. Use the geom_range option below
backward.tree<-ggplot(SI.morph.tree,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
  theme_tree2() + 
  # geom_range(range='height_0.95HPD',alpha=0.4,size=2) +
  scale_x_continuous(labels = abs) +
  geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2)  #+
time.scale<-seq(from=-10,to=0,by=2)
revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
ggsave("SI_Bayes_morph_only.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 

backward.tree<-ggplot(SI.tip.tree,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
  theme_tree2() + 
  # geom_range(range='height_0.95HPD',alpha=0.4,size=2) +
  scale_x_continuous(labels = abs) +
  geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2)  #+
time.scale<-seq(from=-10,to=0,by=2)
revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
# scale_color_continuous(low="gray",high="black") #improve the color scheme
ggsave("SI_Bayes_tip.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 




# Bayesian tree ------
# plot(parsimony.consensus.tree[[1]])
plot(bayesian.mcc.tree)

tree<-bayesian.mcc.tree

#Define which taxa should be highlighted on the tree (preferably single clade)
clade.erethizon<-c("Erethizon_poyeri","Erethizon_kleini","Erethizon_dorsatum")

# #Find which tree edges define the clade
# branches.erethizon<-which.edge(tree,clade.erethizon)
# 
# branchwidth<-c(rep(1.5,length(1:(branches.erethizon[1]-1))),rep(3,length(branches.erethizon)),
#                rep(1.5,(max(tree$edge)-branches.erethizon[length(branches.erethizon)])))
labelmod<-rep(3,length(tree$tip.label))
labelmod[which(tree$tip.label=="Erethizon_poyeri")]<-4

#plot millions ofyears on y axis so peeps can see time calibration. Need to invert
porcupine.node.ages<-(node.depth.edgelength(tree)-max(node.depth.edgelength(tree)))*-1
porcupine.tip.dates<-porcupine.node.ages[c(1:17)]

pdf(file="Fig3_Tree.pdf",width=3.425,height=4,pointsize=9,pagecentre=FALSE) #default units inches
par(mar=c(3,1,0,0))
plot(tree,cex=0.8,font=labelmod, label.offset=0.2,no.margin=FALSE,
           root.edge=FALSE)
axisPhylo()
#add Bayesian posterior probability/bootstrap vals at nodes to show support
nodelabels(data.raw$bpp,data.raw$node,cex=.7, adj=c(1.2,1.4), bg="white", frame="none") 
dev.off()


# parsimony for supplementary information -----------

pdf(file="Fig_SI_ParsimonyTree.pdf",width=3.425,height=4,pointsize=9,pagecentre=FALSE) #default units inches
par(mar=c(3,1,0,0))
plot(parsimony.consensus.tree[[1]])
dev.off()


png(file="Fig_SI_ParsimonyTree.png",width=3.425,height=4,units="in", res=300) #default units inches
par(mar=c(3,1,0,0))
plot(parsimony.consensus.tree[[1]])
dev.off()

