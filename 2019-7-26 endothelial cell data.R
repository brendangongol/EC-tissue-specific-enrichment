########################################################################################
#### endothelial cell data from different regions of the artery GSE: "GSE3874" #########
#######################################################################################
##############################################################################
#### download data from: https://endotheliomics.shinyapps.io/endodb/ #########
##############################################################################
# Search metadata by attributes
# enter parameters: Human, heart, coronary artery
setwd("C:/Users/breng/Dropbox/Brendan Documents/Research/Research project/TIFA/TIFA EC type analysis")
library(data.table)
library(ggpubr)

#### load data ####
###################
GEOD10804 <- fread("E-GEOD-10804 data.csv")
GEOD10804m <- fread("E-GEOD-10804 metadata.csv")

GEOD21212_A <- fread("E-GEOD-21212_A data.csv")
GEOD21212_Am <- fread("E-GEOD-21212_A metadata.csv")

GEOD21212_B <- fread("E-GEOD-21212_B data.csv")
GEOD21212_Bm <- fread("E-GEOD-21212_B metadata.csv")

GEOD43475 <- fread("E-GEOD-43475 data.csv")
GEOD43475m <- fread("E-GEOD-43475 metadata.csv")

#### setup column names ####
#############################
GEOD10804m$colnames <- "Primary"
GEOD21212_Am$colnames <- "Primary"
GEOD21212_Bm$colnames <- "Primary"
GEOD43475m$colnames <- "Primary"

#### setup labels ####
######################
GEOD10804m[grepl("Fresh", GEOD10804m$SampleType),]$colnames <- "Fresh"
GEOD21212_Am[grepl("Fresh", GEOD21212_Am$SampleType),]$colnames <- "Fresh"
GEOD21212_Bm[grepl("Fresh", GEOD21212_Bm$SampleType),]$colnames <- "Fresh"
GEOD43475m[grepl("Fresh", GEOD43475m$SampleType),]$colnames <- "Fresh"

#### setup new column names ####
################################
GEOD10804m$colnames <- paste(GEOD10804m$CellType, GEOD10804m$colnames, sep = "_")
GEOD21212_Am$colnames <- paste(GEOD21212_Am$CellType, GEOD21212_Am$colnames, sep = "_")
GEOD21212_Bm$colnames <- paste(GEOD21212_Bm$CellType, GEOD21212_Bm$colnames, sep = "_")
GEOD43475m$colnames <- paste(GEOD43475m$CellType, GEOD43475m$colnames, sep = "_")

#### subset meta data columns ####
##################################
GEOD10804m <- GEOD10804m[,c(1,14)]
GEOD21212_Am <- GEOD21212_Am[,c(1,14)]
GEOD21212_Bm <- GEOD21212_Bm[,c(1,14)]
GEOD43475m <- GEOD43475m[,c(1,14)]

#### relabel colum names to the tissue types ####
##################################################
setnames(GEOD10804, c("Feature", GEOD10804m$Observation), c("Feature", 
                                                            "Corpus cavernosum endothelial cells_Primary_1", "Corpus cavernosum endothelial cells_Primary_2", "Corpus cavernosum endothelial cells_Primary_3","Corpus cavernosum endothelial cells_Primary_4", "Corpus cavernosum endothelial cells_Primary_5", 
                                                            "Umbilical vein endothelial cells_Primary_1", "Umbilical vein endothelial cells_Primary_2",    "Umbilical vein endothelial cells_Primary_3",    
                                                            "Coronary artery endothelial cells_Primary_1", "Coronary artery endothelial cells_Primary_2",   "Coronary artery endothelial cells_Primary_3",   "Coronary artery endothelial cells_Primary_4"))

setnames(GEOD21212_A, c("Feature", GEOD21212_Am$Observation), c("Feature", 
                                                                "Umbilical vein endothelial cells_Primary_1", "Umbilical vein endothelial cells_Primary_2", 
                                                                "Dermal endothelial cells_Fresh_1", "Dermal endothelial cells_Fresh_2",            
                                                                "Aorta endothelial cells_Fresh_1", "Aorta endothelial cells_Fresh_2",           
                                                                "Pulmonary artery endothelial cells_Fresh_1", "Pulmonary artery endothelial cells_Fresh_2", 
                                                                "Coronary artery endothelial cells_Fresh_1", "Coronary artery endothelial cells_Fresh_2" ))

setnames(GEOD21212_B, c("Feature", GEOD21212_Bm$Observation), c("Feature", 
                                                                "Umbilical vein endothelial cells_Primary_1", "Umbilical vein endothelial cells_Primary_2", 
                                                                "Dermal endothelial cells_Fresh_1", "Dermal endothelial cells_Fresh_2",            
                                                                "Aorta endothelial cells_Fresh_1", "Aorta endothelial cells_Fresh_2",           
                                                                "Pulmonary artery endothelial cells_Fresh_1", "Pulmonary artery endothelial cells_Fresh_2", 
                                                                "Coronary artery endothelial cells_Fresh_1", "Coronary artery endothelial cells_Fresh_2" ))

setnames(GEOD43475, c("Feature", GEOD43475m$Observation), c("Feature", 
                                                            "Hepatic artery endothelial cells_Primary_1", "Hepatic artery endothelial cells_Primary_2", "Hepatic artery endothelial cells_Primary_3",  
                                                            "Aorta endothelial cells_Primary_1", "Aorta endothelial cells_Primary_2",            
                                                            "Coronary artery endothelial cells_Primary_1", "Coronary artery endothelial cells_Primary_2", 
                                                             
                                                            "Iliac artery endothelial cells_Primary_1", "Iliac artery endothelial cells_Primary_2",    
                                                            "Iliac vein endothelial cells_Primary_3", "Iliac vein endothelial cells_Primary_4",  "Iliac vein endothelial cells_Primary_5",  
                                                            
                                                            "Pulmonary artery endothelial cells_Primary_1", "Pulmonary artery endothelial cells_Primary_2", "Pulmonary artery endothelial cells_Primary_3",
                                                            
                                                            "Pulmonary vein endothelial cells_Primary_1",   "Pulmonary vein endothelial cells_Primary_2",   
                                                            
                                                            "Umbilical artery endothelial cells_Primary_1", "Umbilical artery endothelial cells_Primary_2", "Umbilical artery endothelial cells_Primary_3",  
                                                            "Umbilical artery endothelial cells_Primary_4", "Umbilical artery endothelial cells_Primary_5", #"Umbilical artery endothelial cells_Primary_6",
                                                    
                                                            "Umbilical artery endothelial cells_Fresh_1",
                                                            "Umbilical artery endothelial cells_Fresh_2", "Umbilical artery endothelial cells_Fresh_3",   "Umbilical artery endothelial cells_Fresh_4", 
                                                             
                                                              
                                                            "Umbilical vein endothelial cells_Primary_1", "Umbilical vein endothelial cells_Primary_2", "Umbilical vein endothelial cells_Primary_3",  
                                                            "Umbilical vein endothelial cells_Primary_4", "Umbilical vein endothelial cells_Primary_5", 
                                                            
                                                            "Umbilical vein endothelial cells_Fresh_1",    
                                                            "Umbilical vein endothelial cells_Fresh_2", "Umbilical vein endothelial cells_Fresh_3", "Umbilical vein endothelial cells_Fresh_4",     
                                                            
                                                            "Hepatic vein endothelial cells_Primary_1", "Hepatic vein endothelial cells_Primary_2", "Hepatic vein endothelial cells_Primary_3"))    
#### Add dataset names to column names ####
###########################################
setnames(GEOD10804, c("Feature",colnames(GEOD10804[,2:ncol(GEOD10804)])), c("Feature",paste(colnames(GEOD10804[,2:ncol(GEOD10804)]), "GEOD10804", sep = "_")))
setnames(GEOD21212_A, c("Feature",colnames(GEOD21212_A[,2:ncol(GEOD21212_A)])), c("Feature",paste(colnames(GEOD21212_A[,2:ncol(GEOD21212_A)]), "GEOD21212_A", sep = "_")))
setnames(GEOD21212_B, c("Feature",colnames(GEOD21212_B[,2:ncol(GEOD21212_B)])), c("Feature",paste(colnames(GEOD21212_B[,2:ncol(GEOD21212_B)]), "GEOD21212_B", sep = "_")))
setnames(GEOD43475, c("Feature",colnames(GEOD43475[,2:ncol(GEOD43475)])), c("Feature",paste(colnames(GEOD43475[,2:ncol(GEOD43475)]), "GEOD43475", sep = "_")))

#### merge datasets ####
########################
mer <- merge(GEOD10804, GEOD43475, by = "Feature", all = TRUE)
mer <- merge(mer, GEOD21212_A, by = "Feature", all = TRUE)
mer <- merge(mer, GEOD21212_B, by = "Feature", all = TRUE)

# #### obtain annotation information ####
# #######################################
# library("biomaRt")
# # Values <- GEOD10804$Feature
# listMarts() # To choose BioMart database
# myMart <- useMart("ENSEMBL_MART_ENSEMBL");
# listDatasets(myMart)
# myMart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
# listAttributes(myMart)[1:100,] # Choose data types you want to download
# # myMart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast") # use if biomart sever is down.
# go <- getBM(attributes=c( "entrezgene_id", "wikigene_name", "wikigene_description"), mart=myMart)#, values = Values, filters = 'entrezgene_id') #"ensembl_gene_id", "ensembl_transcript_id", "transcription_start_site", "strand", "go_id", "definition_1006"
# write.table(go, file="GO BioMart annotations.xls", sep="\t", quote=FALSE, row.names=FALSE)

#### Annotate with Biomart annotations ####
###########################################
go <- fread("GO BioMart annotations.xls")
setnames(go, colnames(go), c("Feature", "wikigene_name", "wikigene_description"))

mer <- merge(mer, go, by = "Feature")

GEOD10804 <- merge(go, GEOD10804, by = "Feature")
GEOD21212_A <- merge(go, GEOD21212_A, by = "Feature")
GEOD21212_B <- merge(go, GEOD21212_B, by = "Feature")
GEOD43475 <- merge(go, GEOD43475, by = "Feature")

#### Quantile normalize data ####
#################################
library(preprocessCore)
boxplot(mer[,2:71])

norm_edata = normalize.quantiles(as.matrix(mer[,2:71]))

##############################
#### plot normalized data ####
##############################
boxplot(norm_edata)

#### reannotate data ####
#########################
norm_edata <- data.table(norm_edata)
setnames(norm_edata, colnames(norm_edata), colnames(mer[,2:71]))
norm_edata <- cbind(mer[,72:73], norm_edata)

#### generate graphs ####
#########################
library(ggplot2) 

names <- c("TIFA", "CASP1", "CASP4", "CASP5", "GSDMD")
gene <- norm_edata[grepl(names[1], norm_edata$wikigene_name, ignore.case = TRUE),]
gene <- gene[1,2:ncol(gene)]
gene <- melt(gene, id = 1)
title <- gene$wikigene_description[1]
gene <- gene[,2:3]
gene <- gene[!is.na(gene$value)]

p <- ggplot(data=gene, aes(x=reorder(variable, -value), y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("relative expression")+
  xlab("tissue type")+
  ggtitle(title); p
# tiff(file = paste(title, ".tiff", sep = ""), width = 3000, height = 2000, units = "px", res = 300); p; dev.off()

#########################
#### make a function ####
#########################
celltypegraph <- function(names){
gene <- norm_edata[grepl(names, norm_edata$wikigene_name, ignore.case = TRUE),]
gene <- gene[1,2:ncol(gene)]
gene <- melt(gene, id = 1)
title <- gene$wikigene_description[1]
gene <- gene[,2:3]
gene <- gene[!is.na(gene$value)]

p <- ggplot(data=gene, aes(x=reorder(variable, -value), y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("relative expression")+
  xlab("tissue type")+
  ggtitle(title); return(p)
}

celltypegraph("ATP1A1")


#################################################################################################################################
#### generalize colum names, run stats on the samples, and return the genes with the greatest coronary artery EC preference. ####
#################################################################################################################################
#### melt data ####
data_melt <- melt(norm_edata, id = c(1:2))
cols <- data_melt[!duplicated(data_melt$variable),]$variable
cols[order(cols, decreasing = TRUE)]
#### generalize column names ####
types <- c("Corpus cavernosum", "Umbilical vein", "Umbilical artery", "Coronary artery", "Hepatic artery", "Hepatic vein", 
           "Aorta", "Iliac artery", "Iliac vein", "Pulmonary artery", "Pulmonary vein", "Dermal")
data <- NULL
for(i in 1:length(types)){
  type1 <- data_melt[grepl(types[i], data_melt$variable, ignore.case = TRUE),]
  type1$variable <- as.character(type1$variable)
  type1$variable <- types[i]
  data <- rbind(data, type1)
}

#### check to make sure we have all data ####
identical(dim(data_melt), dim(data))


#### create plots ####
gene <- data[data$wikigene_name == "A1BG"]
title <- gene$wikigene_description[1]
gene <- gene[,3:4]
gene <- gene[!is.na(gene$value)]

p <- ggplot(gene, aes(x=variable, y=value))+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), binwidth = 0.05) +
  # geom_point(shape = 21, size = 3, colour = "black", fill = "#08519C")+
  scale_color_manual(values=c("#E31A1C")) +
  ggtitle(title)+ xlab("Cell type") + ylab("Expression Level") +
  theme_pubr() +
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

#

# dat <- data[data$wikigene_name == "CRABP2",]
# dat <- dat[!is.na(dat$value),]
# test <- suppressWarnings(pairwise.wilcox.test(dat$value, as.factor(dat$variable), conf.int = TRUE))
# test <- suppressWarnings(pairwise.wilcox.test(multitab2$log2_mir210_FC, as.factor(multitab2$source), conf.int = TRUE)); print("comparison between miR-210 log2(FC) and site"); test$p.value



#########################################################
#### compare coronary artery to all other cell types ####
#########################################################
datacoro <- data
datacoro[!datacoro$variable == "Coronary artery",]$variable <- "All other"
#### quantify differences ####
keys <- datacoro[!duplicated(datacoro$wikigene_name),]$wikigene_name
datastat <- NULL
for(i in 1:length(keys)){
dat <- datacoro[datacoro$wikigene_name == keys[i],]
dat <- dat[!is.na(dat$value),]
coro <- dat[dat$variable == "Coronary artery",]
other <- dat[dat$variable == "All other",]
meandiff <- mean(coro$value)- mean(other$value)
test <- suppressWarnings(wilcox.test(coro$value , other$value, conf.int = TRUE, paired = FALSE, formula = "lhs"))
dt <- data.table(keys[i], meandiff, test$p.value)
datastat <- rbind(datastat, dt)
}
setnames(datastat, colnames(datastat), c("wikigene_name", "coro_other meandiff", "p_val"))
datastat <- datastat[order(datastat$p_val, decreasing = FALSE),]
datastat
write.table(datastat, file = "./results/ranked genes.xls", sep="\t", quote=FALSE, row.names=FALSE)



dt_mer <- merge(data, datastat, by = "wikigene_name")
dt_mer <- dt_mer[dt_mer$p_val < 0.05,]
dt_mer <- dt_mer[order(dt_mer$p_val, decreasing = FALSE),]
dt_mer2 <- dt_mer
dt_mer2[!dt_mer2$variable == "Coronary artery",]$variable <- "All other"

#### single gene plots ####
DT <- dt_mer2[dt_mer2$wikigene_name == "PIR"]
DT <- DT[!is.na(DT$value),]
title <- DT$wikigene_name[1]
p <- ggplot(DT, aes(x=variable, y=value, fill = variable))+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), binwidth = 0.09) +
  # geom_point(shape = 21, size = 3, colour = "black", fill = "#08519C")+
  scale_color_manual(values=c("#E31A1C")) +
  ggtitle(title)+ xlab("Cell type") + ylab("Expression Level") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p


#### multiple gene plots ####
dtfinal <- dt_mer2[dt_mer2$p_val < 0.0005,]
dtfinal[!duplicated(dtfinal$wikigene_name),]

key <- dtfinal[!duplicated(dtfinal$wikigene_name),]$wikigene_name
DT3 <- NULL
for(a in 1:length(key)){
dt <- dtfinal[dtfinal$wikigene_name == key[a],]
dt <- dt[!is.na(dt$value),]
mea <- mean(dt[dt$variable == "All other",]$value)
  FC <- NULL
  for(i in 1:nrow(dt)){
    FC[i] <- dt$value[i]/mea
  }
dt$FC <- log2(FC)
DT3 <- rbind(DT3, dt)
}

p <- ggplot(DT3, aes(x=reorder(wikigene_name, -`coro_other meandiff`), y=FC, fill = variable))+
  geom_boxplot(position=position_dodge(0.0))+
  # geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.0), binwidth = 0.009) +
  # geom_point(shape = 21, size = 3, colour = "black", fill = "#08519C")+
  scale_color_manual(values=c("#E31A1C")) +
  ggtitle("Genes with a coronary artery enriched gene expression ")+ xlab("gene") + ylab("Expression Level") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

reorder(variable, -value)











