#### 1. setup ####
#  load packages
library(BiocManager)
#BiocManager::install(c("phyloseq", "vegan", "igraph", "dplyr", "tidyr", "radiant.data", "tibble", "reshape2", "Hmisc", "gridExtra"))
library(stringr)
library(ggplot2)
library(phyloseq)
library(vegan)
library(igraph)
library(dplyr)
library(tidyr)
library(radiant.data)
library(tibble)
library(reshape2)
library(Hmisc)
library(gridExtra)
#  set wd
setwd("~/Desktop/RMBL/RMBL_gradients/")


#### 2a. wrangle data - separate AMF and PCB data ####
# #  load metadata
# md<-read.delim("~/Desktop/RMBL/RMBL_gradients/RMBL_gradients_2018_master_data.csv", sep=",", header = TRUE)
# md <- md[is.na(md$metaG_SID)==FALSE & md$SampleID!="C_27_wk_3",]
# row.names(md)<-md$SampleID
# 
# #  load p data
# p_all<-read.delim("sandbox/07032023_pcyc/p_all.csv",sep=",") # gff
# p_all <- p_all[p_all$sampleID!="C_27_wk_3",]
# p_all_fetab<-read.delim("sandbox/07032023_pcyc/p_all_fetab.csv",sep=",") # fetab
# p_all_fetab <- p_all_fetab[,colnames(p_all_fetab)!="C_27_wk_3"]
# p_all_fun <-read.delim("~/Desktop/dbs/phosphorus/p_annotations.csv",sep=",") # annotations
# row.names(p_all_fun)<-p_all_fun$product # not ordered the same as fetab'
# p_all[p_all$taxa=="unknown",]$taxa <- "Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned"
# 
# #  make bacterial species data
# p_all_fetab_tax<-read.delim("sandbox/07032023_pcyc/p_all_fetab_tax.csv",sep=",")
# p_all_fetab_tax <- p_all_fetab_tax[,colnames(p_all_fetab_tax)!="C_27_wk_3"]
# row.names(p_all_fetab_tax) <- str_replace(row.names(p_all_fetab_tax), "unknown","Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned")
# 
# tax_key<-row.names(p_all_fetab_tax)
# 
# p_all_tax_tab <- data.frame(matrix(nrow=3718,ncol=8))
# for(i in 1:3718){
#   this_tax<- str_split_1(tax_key[i], pattern=";")
#   for(j in 1:8){
#     p_all_tax_tab[i,j]<-this_tax[j]
#   }
# }
# row.names(p_all_tax_tab)<-tax_key
# colnames(p_all_tax_tab)<- c("Kingdom", "Phylum","Class", "Order","Family","Genus","Species","Strain")
# p_all_tax_tab<-na.omit(p_all_tax_tab)
# 
# 
# #  add species data to p_all
# p_all[,37:44]<-NA # initialize
# for(i in tax_key){p_all[p_all$taxa==i,37:44] <- p_all_tax_tab[i,]}
# p_all$Phylum <- NULL
# colnames(p_all)<-c(colnames(p_all)[1:35], colnames(p_all_tax_tab))
# rm(i,j,this_tax) # clean env
# 
# #  load in AM data
# load("~/Desktop/RMBL/RMBL_gradients/AMF/curation_output/Phylo_SSU_TSS.RData")
# Phylo_SSU_TSS <- subset_samples(Phylo_SSU_TSS, SampleID%in%colnames(p_all_fetab)) #  subset to just samples with metaG data:
# Phylo_SSU_TSS<-filter_taxa(Phylo_SSU_TSS, function(x) sum(x > 0) > (1), TRUE) #  filter out 0-read taxa:
# sample_names(Phylo_SSU_TSS) <- sample_data(Phylo_SSU_TSS)$SampleID #  match samplenames
# #Phylo_SSU_TSS_sg<-tax_glom(Phylo_SSU_TSS, "Species") #  merge by species (464 vs 85)
# Phylo_SSU_TSS_sg<-Phylo_SSU_TSS
# 
# #  extract otu and taxa tables:
# am <- as.data.frame(otu_table(Phylo_SSU_TSS_sg))
# am_fun <- as.data.frame(tax_table(Phylo_SSU_TSS_sg))
# 
# #  trim p fetab
# p_all_fetab <- p_all_fetab[,colnames(p_all_fetab)%in%colnames(am)]
# p_all_fetab_tax <- p_all_fetab_tax[,colnames(p_all_fetab_tax)%in%colnames(am)]
# 
# load("~/Desktop/RMBL/RMBL_gradients/P-cyc manuscript/network_import.RData")
# 
# #  66 samples with 464 AMF ASVs, 198 PCB products, and 3717 PCB taxa,  22427 observations of PCB genes

#### 2b. wrangle data - merge everything ####
# 
# # order and stack into 1 dataframe:
# p_all_fetab<-p_all_fetab[,sort((colnames(p_all_fetab)))]
# #p_all_fetab_tax<-p_all_fetab_tax[,sort((colnames(p_all_fetab_tax)))]
# am<-am[,sort((colnames(am)))]
# #network_tab<-rbind(am,p_all_fetab)
# 
# #  outline key functions
# ecs <- unique(p_all[p_all$EC_rank=="Hydrolases" & p_all$p_fate=="Mineralized from organic matter",]$product)
# ecs <- c(ecs, "myo-inositol-1(or 4)-monophosphatase")
# lms <- unique(p_all[p_all$pathway=="Lipid metabolism",]$product)
# lms<- lms[!(lms%in%ecs)]
# ps <- unique(p_all[p_all$pathway=="Phosphate Solubilization",]$product)
# rpo<-c("DNA-directed RNA polymerase subunit beta'","DNA-directed RNA polymerase subunit beta")
# key<-c(ecs,lms,ps)
# 
# #subset to key functions and phyla
# #p_all_key<-p_all[p_all$product%in%key & p_all$Phylum%in%c("Acidobacteria", "Proteobacteria", "Actinobacteria"),] #  147
# #p_all_key<-p_all_key[p_all$product%in%key,] #  373
# p_all_key<-p_all
# 
# #  bin products by taxonomic level to add more complexity to network
# #p_all_key$prod_phy <- paste(p_all_key$Genus, p_all_key$product, sep="_") # Genus = 3524
# p_all_key$prod_phy <- paste(p_all_key$Phylum, p_all_key$product, sep="_") # Phylum = 635
# #length(unique(p_all_key$prod_phy))
# 
# #  construct feature table for product x taxa
# p_all_fetab.p.g<-data.frame(matrix(nrow=length(unique(p_all_key$prod_phy)), ncol=66))
# for(i in 1:length(unique(p_all_key$prod_phy))){
#   for(j in 1:66){
#     prod<-unique(p_all_key$prod_phy)[i]
#     sample<-colnames(p_all_fetab)[j]
#     p_all_fetab.p.g[i,j]<-sum(p_all_key[p_all_key$prod_phy==prod & p_all_key$sampleID==sample,]$depth_adj)
#   }
# }
# colnames(p_all_fetab.p.g)<-colnames(p_all_fetab)
# row.names(p_all_fetab.p.g)<-unique(p_all_key$prod_phy)
# rm(i,j,prod,sample)
# 
# # construct annotation table
# key_prod <- data.frame(matrix(nrow=length(unique(p_all_key$prod_phy)), ncol=7))
# colnames(key_prod)<-c("Phylum","product","Richness","Pathway","Key_cat","EC","Mean_abun")
# row.names(key_prod)<-row.names(p_all_fetab.p.g)
# 
# # manually edit tax_tab for consistency
# p_all_tax_tab[p_all_tax_tab$Phylum=="Balneolaeota",]$Phylum<-"Bacteroidetes"
# 
# for(i in 1:length(unique(p_all_key$prod_phy))){
#   this_pp<-row.names(key_prod)[i]
#   key_prod[i,1]<-str_split_1(this_pp,"_")[1] # for phylum
#   #key_prod[i,1]<-unique(p_all_tax_tab[p_all_tax_tab$Genus==(str_split_1(this_pp,"_")[1]),]$Phylum) #  for genus
#   key_prod[i,2]<-str_split_1(this_pp,"_")[2]
#   key_prod[i,3]<-length(unique(p_all_key[p_all_key$prod_phy==this_pp,]$taxa))
#   key_prod[i,4]<-unique(p_all_key[p_all_key$prod_phy==this_pp,]$pathway)
#   key_prod[i,5]<-unique(p_all_key[p_all_key$prod_phy==this_pp,]$pathway) # need to fix this one
#   key_prod[i,6]<-unique(p_all_key[p_all_key$prod_phy==this_pp,]$EC_rank)
#   key_prod[i,7]<-mean(p_all_key[p_all_key$prod_phy==this_pp,]$depth_adj)
# }
# for(i in 1:length(unique(p_all_key$prod_phy))){
#   tp<-key_prod[i,2]
#   if(tp%in%ecs){key_prod[i,5] <- "Extracellular phosphatases"}
#   if(tp%in%lms){key_prod[i,5] <- "Phospholipid turnover"}
#   if(tp%in%ps){key_prod[i,5] <- "Phosphate solubilization"}
#   if(!(tp%in%key)){key_prod[i,5] <- "Other function"}
# }
# key_prod<-key_prod[,c(5,2,3,4,1,6,7)]
# # colnames(key_prod)<-colnames(am_fun)
# #
# p_all_fetab.p.g <- decostand(p_all_fetab.p.g, method = "total", MARGIN=2) #TSS scale it
# am <- decostand(am, method = "total", MARGIN=2) #TSS scale it
# network_tab<-rbind(am,p_all_fetab.p.g)
# network_tax<-rbind(am_fun,key_prod)
# #  change all other bacterial names to "Other Bacterial Phyla"
# network_tax[network_tax$Kingdom!="Fungi" & !(network_tax$Family%in%c("Acidobacteria", "Proteobacteria", "Actinobacteria")),]$Family <- "Other Bacterial Phyla"
# 
# #  add amf framework
# edapho <- c("Gigasporaceae", "Diversisporaceae")
# rhizo <- c("Glomeraceae", "Claroideoglomeraceae", "Paraglomeraceae")
# ans <- c("Archaeosporaceae", "Ambisporaceae", "Pacisporaceae" , "Acaulosporaceae")
# 
# network_tax$fg <- network_tax$Kingdom
# network_tax[network_tax$Family%in%edapho,]$fg <- "Edaphophilic"
# network_tax[network_tax$Family%in%rhizo,]$fg <- "Rhizophilic"
# network_tax[network_tax$Family%in%ans,]$fg <- "Ancestral"
# network_tax[network_tax$fg=="Fungi",]$fg <- "Unassigned"
# 
# network_tax$ft <- network_tax$Family
# network_tax[network_tax$Family%in%edapho,]$ft <- "Edaphophilic"
# network_tax[network_tax$Family%in%rhizo,]$ft <- "Rhizophilic"
# network_tax[network_tax$Family%in%ans,]$ft <- "Ancestral"

# load("~/Desktop/RMBL/RMBL_gradients/sandbox/network_genus.RData")
load("~/Desktop/RMBL/RMBL_gradients/sandbox/network_phylum.RData")


#### 3. set param ####
p<-0.88
p2<-(0-0.60)

ec<-"darkgrey"

## color attributes 
colors<-c("Edaphophilic"="#cafb6f",
           "Rhizophilic"="#71a3ce",
           "Ancestral"="#9acd9b",
          "Unassigned"="grey40",
           
           "Phospholipid turnover"="#ff0060",
           "Extracellular phosphatases"="#fb40fc",
           "Phosphate solubilization"="#8f30e3",
           "Other function"="grey",
           
           "Actinobacteria"="#d3914f", 
           "Acidobacteria"="#cd5c5c",
           "Proteobacteria"="#c5b74f",
           "Other Bacterial Phyla"="grey")


title_legend<-c("Edaphophilic AMF",  "Rhizophilic AMF", "Ancestral AMF", "Unassigned",
                 "Phospholipid turnover","Extracellular phosphatases","Phosphate solubilization", "Other Function", 
                 "Actinobacteria","Proteobacteria", "Acidobacteria", "Other Bacterial Phyla")

#### 4. initialize quant table ####
net_quant <- as.data.frame(matrix(ncol=43,nrow=18))
all_neg<-as.numeric()
all_pos<-as.numeric()
colnames(net_quant) <- c("Plot","Gradient","Elevation_est", 
                         "sign_cor_ratio", "cohes_pos", "cohes_neg", "mod", "numE", "numV",
                         "A","mean_pos", "mean_neg", "num_pos", "num_neg",
                         "pH", "PO4", 
                         "lms.V", "ps.V", "ecp.V", "edapho.V", "rhizo.V", "ans.V", 
                         "FF.E", "FB.E", "BB.E", 
                         "FF.E.pos", "FB.E.pos", "BB.E.pos", 
                         "FF.E.neg", "FB.E.neg", "BB.E.neg",
                         "lms.E.neg", "ps.E.neg", "ecp.E.neg", "edapho.E.neg", "rhizo.E.neg", "ans.E.neg",
                         "lms.E.pos", "ps.E.pos", "ecp.E.pos", "edapho.E.pos", "rhizo.E.pos", "ans.E.pos")
net_quant[,1:3]<-unique(md[,2:4])

#  add pH and PO4
for(i in 1:18){net_quant[i,15] <- mean(md[md$Plot==net_quant[i,1],]$pH)}
for(i in 1:18){net_quant[i,16] <- mean(md[md$Plot==net_quant[i,1],]$phosphate.g.g.1, na.rm=TRUE)}

#### 5. quantify  networks ####
for(i in 1:18){
print(i)
#  subset and calculate correlations
this_wk<-md[md$Plot==net_quant[i,1],]$SampleID
this_nw<-network_tab[,colnames(network_tab)%in%this_wk]  #  subset network to one mountain x elevation
net_correl <- cor(t(this_nw),method="kendall") #  calculate kendall rank correlations
net_correl[is.na(net_correl)] <- 0
net_correl[p > net_correl & net_correl > p2]<- 0  # apply thresholds
#net_correl[p > abs(net_correl)]<- 0  # apply thresholds

##  generate network file and trim nodes
set.seed(444) # <3 <3 <3
net_work <- graph_from_adjacency_matrix(net_correl,mode="upper",weighted=TRUE, diag=FALSE) #  create network file
net_work <- delete_vertices(net_work,degree(net_work)==0) #  remove nodes without edges

##  label vertices -- fungal families & p gene metabolic pathways
net_work_used <- data.frame(V(net_work)$name) #  extract list of species represented in the network
net_work_used$fg <- network_tax[row.names(network_tax)%in%net_work_used$V.net_work..name,]$fg # pcb function and am guild
net_work_used$ft <- network_tax[row.names(network_tax)%in%net_work_used$V.net_work..name,]$ft # pcb taxa and am guild
net_work_used[is.na(net_work_used$fg)==TRUE,2] <- "Unassigned" #  replace NA with unassigned
net_work_used[is.na(net_work_used$ft)==TRUE,3] <- "Unassigned"
V(net_work)$fg = net_work_used$fg #label vertex: pcb function and am guild
V(net_work)$ft = net_work_used$ft #label vertex: pcb taxa and am guild

# ## 3d. create grouping clusters
# net_work_g = net_work #  copy graph
# E(net_work_g)$weight = 1 #  assign weight to edges
# for(j in unique(V(net_work)$Kingdom)){ #  for each Family
#   GroupV = which(V(net_work)$Kingdom == j) #  pull indeces of vertexes in group i
#   net_work_g = add_edges(net_work_g, #  add edges to the graph - add weight of 5 to #list of all possible combinations
#                          combn(GroupV, 2),attr=list(weight=6))
# }
# LO = layout_with_fr(net_work_g) #create a layout based on G_Grouped

##  create grouping clusters
net_work_g = net_work #  copy graph
E(net_work_g)$weight = 1 #  assign base weight to edges
cluster_27<-as.data.frame(membership(cluster_louvain(net_work_g))) #  cluster with louvain algorithm
V(net_work_g)$cluster <-cluster_27$x #  add clusters to vertices
for(j in unique(V(net_work_g)$cluster)){
  GroupV = which(V(net_work_g)$cluster == j) #  pull indeces of vertexes in group i
  net_work_g = add_edges(net_work_g, combn(GroupV, 2),attr=list(weight=6))
} #  add weight of 6 to each cluster group
LO = layout_with_fr(net_work_g) #create a layout based on grouped network
rm(net_work_g)

##  quantify
cors <- melt(net_correl) #  melt correlation network
cors2<-cors[cors$value!=0  & is.na(cors$value)==FALSE & cors$value !=1 & cors$value != (0-1),] #  just pull out significant cors
cors2$Var1.g <- NA
cors2$Var2.g <- NA
for(j in 1:dim(cors2)[1]){
  cors2[j,4]<-network_tax[row.names(network_tax)==as.character(cors2[j,1]),8]
  cors2[j,5]<-network_tax[row.names(network_tax)==as.character(cors2[j,2]),8]
}
cors2$Var1.k <- cors2$Var1.g
cors2$Var2.k <- cors2$Var2.g
cors2$Var1.k<-str_replace_all(cors2$Var1.k, c("Ancestral" = "AM Fungi", "Unassigned"= "AM Fungi", "Rhizophilic"= "AM Fungi", "Edaphophilic"= "AM Fungi", "Other function" = "P Cycling Bacteria", "Phospholipid turnover" = "P Cycling Bacteria", "Extracellular phosphatases" = "P Cycling Bacteria", "Phosphate solubilization" = "P Cycling Bacteria"))
cors2$Var2.k<-str_replace_all(cors2$Var2.k,  c("Ancestral" = "AM Fungi", "Unassigned"= "AM Fungi", "Rhizophilic"= "AM Fungi", "Edaphophilic"= "AM Fungi", "Other function" = "P Cycling Bacteria", "Phospholipid turnover" = "P Cycling Bacteria", "Extracellular phosphatases" = "P Cycling Bacteria", "Phosphate solubilization" = "P Cycling Bacteria"))
cors2$edge_type <- paste(cors2$Var1.k, cors2$Var2.k, sep="-")
values<-cors2$value

#  cohesion 
#net_quant[i,5] <- cohesion(net_work)
#  negative to positive ratio:
net_quant[i,4] <- length(values[values<0])/length(values[values>0])
all_neg <- c(all_neg, values[values<0])
all_pos <- c(all_pos, values[values>0])
#  number of edges (connections)
net_quant[i,8] <- ecount(net_work)
#  number of vertices (points)
net_quant[i,9] <- vcount(net_work)
#net_quant[i,9] <- vcount(net_work)
#  modularity
E(net_work)$weight <- E(net_work)$weight  + 1 #  add 1 to make them positive
E(net_work)$weight[is.na(E(net_work)$weight)]<-0 #  remove NAs
clusters<-cluster_louvain(net_work, weights = E(net_work)$weight, resolution = 1)
net_quant[i,7] <- modularity(clusters)
#  assortativity
V(net_work)$fg <- factor(V(net_work)$fg) #  vertex value must be factor or int for assortativity calculation
net_quant[i,10] <- assortativity_nominal(net_work, types=V(net_work)$fg)
#net_quant[i,10] <- assortativity_degree(net_work)

#  mean positive cor value
net_quant[i,11] <- mean(values[values>0])
net_quant[i,12] <- mean(values[values<0])
net_quant[i,13] <- length(values[values>0])
net_quant[i,14] <- length(values[values<0])
#  number of vertices per guild/functional group
Vs <- as.character(unique(c(cors2$Var1, cors2$Var2)))
net_quant[i,17] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Phospholipid turnover",])])
net_quant[i,18] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Phosphate solubilization",])])
net_quant[i,19] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Extracellular phosphatases",])])
net_quant[i,20] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Edaphophilic",])])
net_quant[i,21] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Rhizophilic",])])
net_quant[i,22] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Ancestral",])])
#  number of F-F, F-B, or B-B edges
net_quant[i,23]<-dim(cors2[cors2$edge_type == "AM Fungi-AM Fungi",])[1]
net_quant[i,24]<-dim(cors2[cors2$edge_type == "AM Fungi-P Cycling Bacteria" | cors2$edge_type == "P Cycling Bacteria-AM Fungi",])[1]
net_quant[i,25]<-dim(cors2[cors2$edge_type == "P Cycling Bacteria-P Cycling Bacteria",])[1]
#  number of F-F, F-B, or B-B positive edges
net_quant[i,26]<-dim(cors2[cors2$edge_type == "AM Fungi-AM Fungi" & cors2$value > 0,])[1]
net_quant[i,27]<-dim(cors2[cors2$value > 0 & cors2$edge_type == "AM Fungi-P Cycling Bacteria" | cors2$edge_type == "P Cycling Bacteria-AM Fungi",])[1]
net_quant[i,28]<-dim(cors2[cors2$value > 0 & cors2$edge_type == "P Cycling Bacteria-P Cycling Bacteria",])[1]
#  number of F-F, F-B, or B-B negative edges
net_quant[i,29]<-dim(cors2[cors2$edge_type == "AM Fungi-AM Fungi" & cors2$value < 0,])[1]
net_quant[i,30]<-dim(cors2[cors2$value < 0 & cors2$edge_type == "AM Fungi-P Cycling Bacteria" | cors2$edge_type == "P Cycling Bacteria-AM Fungi",])[1]
net_quant[i,31]<-dim(cors2[cors2$value < 0 & cors2$edge_type == "P Cycling Bacteria-P Cycling Bacteria",])[1]
#  number of negative edges per guild/group
Vs <- as.character(unique(c(cors2[cors2$value <0,]$Var1, cors2[cors2$value <0,]$Var2)))
net_quant[i,32] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Phospholipid turnover",])])
net_quant[i,33] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Phosphate solubilization",])])
net_quant[i,34] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Extracellular phosphatases",])])
net_quant[i,35] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Edaphophilic",])])
net_quant[i,36] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Rhizophilic",])])
net_quant[i,37] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Ancestral",])])

#  number of positive edges per guild/group
Vs <- as.character(unique(c(cors2[cors2$value >0,]$Var1, cors2[cors2$value >0,]$Var2)))
net_quant[i,38] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Phospholipid turnover",])])
net_quant[i,39] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Phosphate solubilization",])])
net_quant[i,40] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Extracellular phosphatases",])])
net_quant[i,41] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Edaphophilic",])])
net_quant[i,42] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Rhizophilic",])])
net_quant[i,43] <- length(Vs[Vs%in%row.names(network_tax[network_tax$fg=="Ancestral",])])
}

#### 6. plot x stats for network quant ####
load("~/Desktop/RMBL/RMBL_gradients/sandbox/netquant_p88_pn65_pcbPhylum_amVtx.Rdata")
load("~/Desktop/RMBL/RMBL_gradients/sandbox/netquant_p88_pn6_pcbGenus_amVtx.Rdata")


net_quant$Elevation_m <- net_quant$Plot
for(i in 1:length(net_quant$Elevation_m)){
  net_quant[i,]$Elevation_m <- (md[md$Plot==net_quant[i,]$Plot,]$Elevation_m)[1]
}
net_quant$Elevation_m <- as.numeric(net_quant$Elevation_m)

library(nlme)
net_quant[is.infinite(net_quant$sign_cor_ratio)==TRUE,4] <- NA

anova(lme(sign_cor_ratio~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(A~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(mod~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(num_neg~Elevation_m + pH, random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(num_pos~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))

anova(lme(lms.V~Elevation_m + pH, random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(ecp.V~Elevation_m + pH, random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(ps.V~Elevation_m + pH, random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(edapho.V~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(rhizo.V~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(ans.V~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))

# anova(lme(edapho.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(rhizo.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(ans.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(lms.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(ecp.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(ps.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# 
# anova(lme(FB.E~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(FF.E~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
# anova(lme(BB.E~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))

anova(lme(FB.E.pos~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(FF.E.pos~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(BB.E.pos~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))

anova(lme(FB.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(FF.E.neg~Elevation_m + pH , random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))
anova(lme(BB.E.neg~Elevation_m + pH, random= ~ 1 | Gradient, data=net_quant, na.action=na.exclude))


#  melt for data plotting
net_quant.m<-melt(net_quant, id.vars = c("Plot", "Gradient", "Elevation_m", "pH", "PO4"))

Vs <- c("edapho.V", "rhizo.V", "ans.V")
nq.f<-ggplot(net_quant.m[net_quant.m$variable%in%Vs,], aes(x=Elevation_m,y=value, color = variable, linetype=variable)) + geom_point(size=3) + theme_bw() + geom_smooth(method = "lm",se=F) + 
  xlab("Elevation (m)") + ylab("Number of AM Fungal Vertices (Nodes)") + ylim(0,39)+
  scale_color_manual(values=c("#cafb6f", "#71a3ce","#9acd9b" )) + 
  scale_linetype_manual(values=c("dotted", "solid", "solid")) + 
  annotate("text", x = 3250, y = Inf, label = "Rhizophilic: F(1,13)= 26.65, p=0.0002", color = "#71a3ce", size = 5, vjust = 1.5) + 
  annotate("text", x = 3250, y = Inf, label = "Ancestral: F(1,13)= 10.78, p=0.006", color = "#9acd9b", size = 5,  vjust = 3.2) + 
  annotate("text", x = 3250, y = Inf, label = "Edaphophilic: F(1,13)= 4.36, p=0.06", color = "#cafb6f", size = 5, vjust = 4.7) + 
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12), title=element_text(size=12))

Vs <- c("ps.V", "lms.V", "ecp.V")
nq.b <- ggplot(net_quant.m[net_quant.m$variable%in%Vs,], aes(x=Elevation_m,y=value, color = variable, linetype=variable)) + geom_point(size=3) + ylim(50,318)+
  theme_bw() + geom_smooth(method = "lm",se=F) + 
  xlab("Elevation (m)") + ylab("Number of Phosphorus Cycling Bacterial Vertices (Nodes)") + 
  scale_color_manual(values=c("#ff0060", "#8f30e3","#fb40fc" )) + 
  scale_linetype_manual(values=c("dotted", "solid", "dotted")) + 
  annotate("text", x = 3250, y = Inf, label = "Phosporus Solubilization: F(1,13)= 5.23, p=0.03", color = "#8f30e3", size = 5,  vjust = 1.5) + 
  annotate("text", x = 3250, y = Inf, label = "Extracellular Phosphatases: F(1,13)= 0.14, p=0.71", color = "#fb40fc", size = 5,  vjust = 3.2) + 
  annotate("text", x = 3250, y = Inf, label = "Phospholipid Turnover: F(1,13)= 2.23, p=0.15", color = "#ff0060", size = 5,  vjust = 4.7) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12), title=element_text(size=12))

Vs <- c("num_neg", "num_pos")
nq.ne<-ggplot(net_quant.m[net_quant.m$variable%in%Vs,], aes(x=Elevation_m,y=value, color= variable, linetype= variable)) + geom_point(size=3) + theme_bw() + geom_smooth(method = "lm", se=F) + 
  xlab("Elevation (m)") + ylab("Number of Edges (Connections)") +
  scale_linetype_manual(values = c("dotted", "solid")) + 
  scale_color_manual(values=c("grey", "red")) + 
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12), title=element_text(size=12)) + 
  annotate("text", x = 3250, y = Inf, label = "Negative Edges: F(1,13)= 4.93, p=0.04", color = "red", size = 5,  vjust = 1.5) + 
  annotate("text", x = 3250, y = Inf, label = "Positive Edges: F(1,13)= 2.13, p=0.16", color = "grey", size = 5, vjust = 3.2) 

Vs <- c("A")
nq.a<-ggplot(net_quant.m[net_quant.m$variable%in%Vs,], aes(x=Elevation_m,y=value, color= variable, linetype= variable)) + geom_point(size=3) + theme_bw() + geom_smooth(method = "loess", se=F) + 
  xlab("Elevation (m)") + ylab("Assortativity") +
  scale_linetype_manual(values = c("dotted", "solid")) + 
  scale_color_manual(values=c("black")) + 
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12), title=element_text(size=12)) + 
  annotate("text", x = Inf, y = Inf, label = "F(1,15)= 3.87, p=0.04", color = "black", size = 5, hjust = 1.1, vjust = 1.5) 
summary(lm(A~poly(Elevation_est, 2, raw=TRUE), data=net_quant, na.action=na.exclude))

grid.arrange(nq.a, nq.ne , nq.f, nq.b, nrow=1)
#### make and plot network  ####

# ### A. Subset Data
# this_wk<-md[md$Elevation_est==2700,]$SampleID
# this_lab<-"2700 meters"
# # # 
# this_wk<-md[md$Elevation_est==2900,]$SampleID
# this_lab<-"2900 meters"
# # 
# this_wk<-md[md$Elevation_est==3100,]$SampleID
# this_lab<-"3100 meters"
# # # 
# this_wk<-md[md$Elevation_est==3200 |md$Elevation_est==3300, ]$SampleID
# this_lab<-"3300 meters"
# # 
# this_wk<-md[md$Elevation_est==3500,]$SampleID
# this_lab<-"3500 meters"
# # # 
# this_wk<-md[md$Elevation_est==3600 |md$Elevation_est==3700,]$SampleID
# this_lab<-"3700 meters"
this_nw<-network_tab[,colnames(network_tab)%in%this_wk]

###  B. calculate the correlation network & trim correlations based on thresholds
# net_correl <- rcorr(t(this_nw),type = "kendall")
# nc1<-(net_correl[["r"]])
# ncp<-(net_correl[["P"]])
net_correl <- cor(t(this_nw),method="kendall")
#net_correl[abs(net_correl)<p]<-0  #  same p for positive and negative cor
net_correl[p > net_correl & net_correl > p2]<-0   #  different p for positive and negative cor
net_correl[is.na(net_correl)] <- 0
nc1<-net_correl
nc1.melt <-melt(nc1)

###  C. generate network file and trim nodes
set.seed(444) # <3 <3 <3
net_work <- graph_from_adjacency_matrix(nc1,mode="upper",weighted=TRUE, diag=FALSE) #  create network file
net_work <- delete.vertices(net_work,degree(net_work)==0) #  remove nodes without edges

### D. create and add labels -- fungal families & p gene metabolic pathways
net_work_used <- data.frame(V(net_work)$name) #  extract species represented in the network
# copy over (they have to be in the same order)
net_work_used$fg <- network_tax[row.names(network_tax)%in%net_work_used$V.net_work..name,]$fg #AMF guild and PCB metabolic group
net_work_used$ft <- network_tax[row.names(network_tax)%in%net_work_used$V.net_work..name,]$ft #AMF guild and PCB Phyla
#  replace NA with unassigned
net_work_used[is.na(net_work_used$fg)==TRUE,2] <- "Unassigned"
net_work_used[is.na(net_work_used$ft)==TRUE,3] <- "Unassigned"

#  label vertex
V(net_work)$ft = net_work_used$ft #fungal families, bacterial Kingdoms
V(net_work)$fg = net_work_used$fg #fungal kingdom, functional gene group

# add color attribute
V(net_work)$color <- str_replace_all(V(net_work)$fg, colors) 
#V(net_work)$color2 <- str_replace_all(V(net_work)$Kingdom, colors2)

# ## 3d. create grouping clusters
net_work_g = net_work #  copy graph
E(net_work_g)$weight = 1 #  assign weight to edges
cluster_27<-as.data.frame(membership(cluster_louvain(net_work_g)))
V(net_work_g)$cluster <-cluster_27$x
for(i in unique(V(net_work_g)$cluster)){ #  for each Family
  GroupV = which(V(net_work_g)$cluster == i) #  pull indeces of vertexes in group i
  net_work_g = add_edges(net_work_g, #  add edges to the graph - add weight of 5 to #list of all possible combinations
                         combn(GroupV, 2),attr=list(weight=5))
}


#create a layout based on G_Grouped
LO = layout_with_kk(net_work_g)

#add sign attribute
E(net_work)$sign <- E(net_work)$weight

#  add correlation
# cors[,4]<-melt(ncp)[,3] #  add p values
# cors[is.na(cors$V4)==TRUE,]$V4 <- 1 #  turn NA p vals to 1
# cors<-cors[cors$V4 <.05,]
# E(net_work)$cor <- (cors$value)*(cors$value) #  correlations stat
# E(net_work)$ps<- cors$V4

E(net_work)$weight <- E(net_work)$weight  + 1
E(net_work)$weight[is.na(E(net_work)$weight)]<-0

#  add correlations
correlations <- E(net_work)$sign # pull out cors
signs <- ifelse(correlations > 0, "positive", "negative")# Step 2: Determine sign of each correlation
edge_colors <- ifelse(signs == "positive", "grey", "red")# Step 3: Assign colors based on sign

#  vertex sizes
sizes<-c("Edaphophilic"="5",
           "Rhizophilic"="5",
           "Ancestral"="5",
            "Unassigned"="5",
           "Phospholipid turnover"="2",
           "Extracellular phosphatases"="2",
           "Phosphate solubilization"="2",
           "Other function"="2",
         "Actinobacteria"="2", 
         "Acidobacteria"="2",
         "Proteobacteria"="2",
         "Other Bacterial Phyla"="2")
V(net_work)$size <- as.numeric(str_replace_all(V(net_work)$ft, sizes))

# plot!
par(mfrow=c(1,1), mar=c(0,0,0,0))

#n31 <- recordPlot()
plot(net_work, 
     vertex.size= V(net_work)$size, #size of nodes
     #vertex.label = V(net_work)$Family, #remove microbes names
     vertex.label = NA, #remove microbes names
     vertex.color=V(net_work)$color,
     #vertex.shape=V(net_work)$shape,
     edge.color=edge_colors, edge.width=1,
     main=NULL,
     margin=-0.07,
     layout=LO
)

legend(x = -1.75, y = .65, title_legend3,
       pch = 21, pt.bg = colors3,
       pt.cex = 2.0, cex=.75, bty = "n", ncol = 1) #adding legend


#### 8. plot everything together ####
lay <- rbind(c(1,2,3), 
             c(4,5,6,7))
             
grid.arrange(grobs=list(ggplotGrob(as.ggplot(n27)),ggplotGrob(as.ggplot(n31)),ggplotGrob(as.ggplot(n35)),ggplotGrob(nq.a), ggplotGrob(nq.ne), ggplotGrob(nq.f),  ggplotGrob(nq.b)),
             layout_matrix=lay)


#### cluster analysis ####
library(dplyr)
clustered_taxa <- merge(cluster_27,cluster_29,by="row.names",all=TRUE)
row.names(clustered_taxa)<-clustered_taxa$Row.names
clustered_taxa <-clustered_taxa[,2:3]
clustered_taxa <- merge(clustered_taxa,cluster_31,by="row.names",all=TRUE)
row.names(clustered_taxa)<-clustered_taxa$Row.names
clustered_taxa <-clustered_taxa[,2:4]
clustered_taxa <- merge(clustered_taxa,cluster_32,by="row.names",all=TRUE)
row.names(clustered_taxa)<-clustered_taxa$Row.names
clustered_taxa <-clustered_taxa[,2:5]
clustered_taxa <- merge(clustered_taxa,cluster_35,by="row.names",all=TRUE)
row.names(clustered_taxa)<-clustered_taxa$Row.names
clustered_taxa <-clustered_taxa[,2:6]
clustered_taxa <- merge(clustered_taxa,cluster_36,by="row.names",all=TRUE)
row.names(clustered_taxa)<-clustered_taxa$Row.names
clustered_taxa <-clustered_taxa[,2:7]
colnames(clustered_taxa) <- c("e27","e29","e31","e32","e35","e36")

clustered_taxa$id <- row.names(clustered_taxa)
ct.m <- melt(clustered_taxa, id.vars = "id")
ct.m$elev <- NA
ct.m[ct.m$variable=="e27",]$elev<-2700
ct.m[ct.m$variable=="e29",]$elev<-2900
ct.m[ct.m$variable=="e31",]$elev<-3100
ct.m[ct.m$variable=="e32",]$elev<-3200
ct.m[ct.m$variable=="e35",]$elev<-3500
ct.m[ct.m$variable=="e36",]$elev<-3600

ct.m[,5:12]<-NA
for(i in unique(ct.m$id)){
  ct.m[ct.m$id==i,5:12] <- network_tax[i,]
}

ggplot(ct.m, aes(x=elev, y=value, group=id)) + 
  geom_point() + geom_line() + 
  facet_wrap(facets="V12", nrow=2,ncol=3) + 
  facet_wrap(.~variable)


#### text graveyard ####
#& cors$value!=1 & cors$value !=-1
#net_quant[i,8] <- length(values)
#net_quant[i,9] <- length(unique(c(cors2$Var1, cors2$Var2)))





