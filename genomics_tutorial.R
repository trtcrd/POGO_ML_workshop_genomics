###########################################################################################################

### AI/ML POGO workshop - Ostend, Belgium, 20-22th May 2019 
# Genomics - Anders Lanzén & Tristan Cordier

########################################### Setting the stage #############################################

# setting the working directory
setwd("/Users/tristan/Desktop/POGO/POGO_ML_workshop_genomics/")

# libraries we need
library('vegan')
library('matrixStats')
library('ranger')
library('irr')
library('BBI')
## for parralel CPU computation
library('doMC')
registerDoMC(cores = detectCores())

# customs functions
source("custom_functions/sml_compo.R")
source("custom_functions/plot_ml.R")
source("custom_functions/bar_plot.R")


### import our data
# the otu table
otu <- read.table("otu_table.tsv", header=TRUE, sep="\t", row.names=1)
# the metadata
met <- read.table("metadata.tsv", header=TRUE, sep="\t", dec = ".")


################################### 'Playing' and manipulating the data ###################################

### some constant variables handling
# to fecth easily the list of samples
s <- "samples_names"
# the cuttof to consider a sample
seq_depth_cutoff <- 10000

# check if empty rows or columns at the end of file
otu[(dim(otu)[1]-10):(dim(otu)[1]),1:3]  ## normally there is nothing 
otu[1:3,(dim(otu)[2]-10):(dim(otu)[2])]  ## normally Taxa and X

# extract the OTU taxonomy assignments and total abundance from the OTU table as generated from the SLIM application (Dufresne et al., 2019)
taxo <- cbind(Total_abund = rowSums(otu[,c(1:(dim(otu)[[2]]-3))]), taxon = as.character(otu[,c(dim(otu)[[2]]-2)]))
 ### and keep the otu counts only (transposed: samples as rows, variables as columns)
otu  <- data.frame(t(otu[,c(1:(dim(otu)[[2]]-3))]))

### making sure of the match between metadata and otu table
# are they there? 
table(gsub("-",".",met[,s]) %in% rownames(otu))
# and do they match exactly?
table(gsub("-",".",met[,s]) == rownames(otu))

### keep the dataset at this point for coming back to it
otu_raw <- otu

### sum of samples 
row_S <- rowSums(otu_raw)

### plot the sequencing depth
pdf("Sequencing_depth1.pdf", width=6,height=6)
plot(row_S, main="Sequencing depth - limit at 10,000 reads", ylab = "Reads counts", xlab = "sample ID", pch=16, cex=.4)
abline(h=10000, col="blue", lwd=2)
dev.off()

### plot the sequencing depth
pdf("Sequencing_depth2.pdf", width=6,height=6)
plot(log(row_S+1), main="Sequencing depth - limit at 10,000 reads", ylab = "log(Reads counts + 1)", xlab = "sample ID", yaxt="n", pch=16, cex=.4)
axis(2, at=log(c(10,100,1000,10000,100000,1000000)),labels=c(10,100,1000,10000,100000,1000000),cex.axis=0.7, las=2, ylim=log(c(10,1000000)))
abline(h=log(10000+1), col="blue", lwd=2)
dev.off()

## get samples with sequencing depth above the total average of reads per samples
otu  <- subset(otu, row_S >= seq_depth_cutoff)
met  <- subset(met, row_S >= seq_depth_cutoff)

## keep only samples for which we have metadata
otu     <- subset(otu, is.na(met$AMBI) == F)
otu_raw <- subset(otu_raw, is.na(met$AMBI) == F)
met     <- subset(met, is.na(met$AMBI) == F)

### normalize the counts (there is not "good" approach, see Weiss et al. Microbiome, 2017; Pereira et al. BMC Genomics, 2018) 
# here let's do relative abundance.. 
OTU  <- decostand(otu, method = "total")
OTUr <- otu
MET  <- met 

######### ready to go!

########################################### Real analysis begins here #############################################

############ Exploring big patterns

# ordination (NMDS)
mds <-metaMDS(OTU, distance = "bray", k = 2, binary = F) 

# explore relations between ordinations and environmental parameters
env_dat <- MET[,c("Distance_cage", "AMBI", "NSI", "ISI", "NQI1", "Shannon")]

# plot the NMDS
pdf("NMDS1.pdf", width=12,height=7)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(mds$points, col=as.numeric(MET$col_plot), pch=MET$col_plot, main = "NMDS - Bray-Curtis matrix")
leg <- unique(as.vector(MET$Locality))
legend("topright", inset=c(-0.3,0), leg, col=unique(MET$col_plot), pch=unique(MET$col_plot), box.lty=0)
legend("topright", legend = paste("stress =", round(mds$stress, 3)), box.lty=0, cex=.8)
fit <- envfit(mds, env_dat, perm = 999, display = "sites", na.rm =T)
tab <- cbind(round(fit$vectors$r,3), fit$vectors$pvals)
plot(fit, p.max = 0.05, col = "chartreuse2")
dev.off()

# Focus on pollution gradient
col_gradient <- colorRampPalette(c("blue", "green","red"), bias=1, interpolate = "linear")(600)
pdf("NMDS2.pdf", width=12,height=7)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(mds$points, col=col_gradient[as.numeric(floor(MET$AMBI*100))], pch=16, main = "NMDS - Bray-Curtis matrix")
leg <- c("low pollution", "low to moderate", "moderate to high")
legend("topright", inset=c(-0.3,0), leg, col=c("blue", "green","red"), pch=16, box.lty=0)
dev.off()

# taxonomic barplot
# high taxonomic rank
pdf("Taxo_barplot1.pdf")
bar_plot(OTUr, MET, aggreg = c("Locality"), taxo_file = taxo, title_ = "By locality", font_size_ = 0.7, tax_rank = 4)
dev.off()

# high taxonomic rank
pdf("Taxo_barplot2.pdf")
bar_plot(OTUr, MET, aggreg = c("Locality"), taxo_file = taxo, title_ = "By locality", font_size_ = 0.7, tax_rank = 7, metazoan = T)
dev.off()

############ SML 

# custom function to perform cross-validation on a fully labelled dataset
preds <- sml_compo(OTU, MET, index = "AMBI", algo = "RF", cross_val = "Locality")

# plot the independant predictions against the real values and collect stats
pdf("Predictions1.pdf", width=8,height=5)
res <- plot_ml(preds, MET, index = "AMBI", title = "Random Forest", aggreg = c("Grab", "Station", "Locality"))
dev.off()

# we can do that within a loop and export the plots
for (BI in c("AMBI", "NSI", "ISI", "NQI1", "Shannon"))
{
  ## SML on the 10 most abundant for doing it fast...
  preds <- sml_compo(OTU[,1:10], MET, index = BI, algo = "RF", cross_val = "Locality")
  ## PLOT
  plot_ml(preds, MET, index = BI, aggreg = c("Grab", "Station", "Locality"), pdf = T,
          title = paste("Predictions_", BI, ".pdf", sep=""))
}


### Variables importance - model on the full dataset 
mod <- ranger(MET[,"AMBI"] ~ ., data=OTU, mtry=floor(dim(OTU)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
imp <- tail(sort(mod$variable.importance), 50)
pdf("Variable_importance.pdf", width = 5, height=8)
p <- barplot(imp, horiz = T, xlab="Variable importance", main=paste("OTUs importance for AMBI"), las=2, cex.names = 0.6)
### look at these OTUs
tx <- taxo[names(imp),"taxon"]
taxo[names(imp),"taxon"]
## get the last rank
for (i in 1:length(tx)) tx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", fixed=TRUE)), 1)
## and add it on the plot
text(0, p, labels=tx, pos=4, cex=0.6)
dev.off()




######### WORKSHOP day one: additional scripting.

## let's make feature selection to keep only OTUs with good predictive potential and test if this improve performance. 
## the random forest algorithm had this feature built-in

## to test that (on only one farm in this exemple), we will make feature selection on a training set of 
## of all the farm but Farm_1 and predict BI values on this hold-out farm. 
TRAIN <- subset(OTU, MET$Locality != "Farm_1")
TEST <- subset(OTU, MET$Locality == "Farm_1")

MET_tr <- subset(MET, MET$Locality != "Farm_1")
MET_te <- subset(MET, MET$Locality == "Farm_1")

## let's train a model will all the features for the AMBI index
mod_all_feature_AMBI <- ranger(MET_tr[,"AMBI"] ~ ., data=TRAIN, mtry=floor(dim(TRAIN)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
## sort the OTUs by decreasing importance
imp_AMBI <- sort(mod_all_feature_AMBI$variable.importance, decreasing = T)

## same for NSI index
mod_all_feature_NSI <- ranger(MET_tr[,"NSI"] ~ ., data=TRAIN, mtry=floor(dim(TRAIN)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
## sort the OTUs by decreasing importance
imp_NSI <- sort(mod_all_feature_NSI$variable.importance, decreasing = T)

## and we will see what it the effect of the amount of feature we keep on the accuracy of predictions, from 10 to all OTUs
feat <- c(10,100,250,500,1000, ncol(TRAIN))

## collect the R2 and kappa statistics 
R2 <- array(NA, c(2, length(feat)))
colnames(R2) <- feat
rownames(R2) <- c("AMBI", "NSI")

KA <- array(NA, c(2, length(feat)))
colnames(KA) <- feat
rownames(KA) <- c("AMBI", "NSI")

for (i in feat)
{
  # get the i most important OTUs
  OTU_list_AMBI <- names(imp_AMBI[1:i])
  OTU_list_NSI  <- names(imp_NSI[1:i])
  otu_ambi <- TRAIN[,OTU_list_AMBI]
  otu_nsi  <- TRAIN[,OTU_list_NSI]
  # and keep those OTUs only in the TEST dataset
  test_ambi <- TEST[,OTU_list_AMBI]
  test_nsi  <- TEST[,OTU_list_NSI]
  
  ## SML
  for (BI in c("AMBI", "NSI"))
  {
    ## train AMBI predictive model 
    mod_AMBI <- ranger(MET_tr[,"AMBI"] ~ ., data=otu_ambi, mtry=floor(dim(otu_ambi)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
    ## train NSI predictive model 
    mod_NSI  <- ranger(MET_tr[,"NSI"] ~ ., data=otu_nsi, mtry=floor(dim(otu_nsi)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
    
    ## now make prediction on the hold-out farm
    preds_ambi <- predict(mod_AMBI, test_ambi)
    preds_nsi  <- predict(mod_NSI, test_nsi)
    
    ## get stats for each index
    res <- plot_ml(data = preds_ambi$predictions, metadata = MET_te, index = BI,  title = paste("RF_", i, "features", sep=""),
                   aggreg = c("Grab", "Station", "Locality"), taxo_group = paste("V1V2_", i, sep=""), pdf = T)
    R2[paste(BI),paste(i)] <- res$R2
    KA[paste(BI),paste(i)] <- res$KAP
  }
}


### then plot 


quartz(width = 5, height = 5)
plot(as.numeric(R2_[1,]), col="blue", cex.axis=.8,pch=1,type="o", ylim=c(min(R2_), max(R2_)), xaxt="none", ylab = "R2", xlab="Features extracted", main="R2 - random forest")
lines(as.numeric(R2_[2,]), type="o", pch=2, col="red")
lines(as.numeric(R2_[3,]), type="o", pch=3, col="green")
lines(as.numeric(R2_[4,]), type="o", pch=4, col="yellow")
lines(as.numeric(R2_[5,]), type="o", pch=5, col="black")
axis(1, at=1:ncol(R2_), cex.axis=.8,lab=gsub("F", "", colnames(R2_)))

legend("topright", rownames(R2_), cex=0.7, col=c("blue","red", "green", "yellow", "black"), pch=1:5,bty = "n")








## keeping only the best 50 OTUs
TRAIN_feat_sel <- TRAIN[,names(imp)]

predictions <- predict(mod, TEST)




##

formula <- as.formula(paste('get(nam_comp)[,index] ~ ' ,paste(dimnames(OTUtr)[[2]],collapse='+')))

set.seed(1)
assign(mod,neuralnet(formula, data = get(nam), hidden=10, err.fct="sse", algorithm = "rprop+"))

## feature selection

mod_ranger <- ranger(MET_tr[,"AMBI"] ~ ., data=TRAIN, 
                     mtry=floor(dim(TRAIN)[2]/3), num.trees = 300, 
                     importance= "impurity", write.forest = T)
## which ones to keep?
imp <- tail(sort(mod_ranger$variable.importance), 50)
## keeping only the best 50 OTUs
TRAIN_feat_sel <- TRAIN[,names(imp)]

# and retraining the model with them only
mod_features_sel <- ranger(MET_tr[,"AMBI"] ~ ., data=TRAIN_feat_sel, 
                           mtry=floor(dim(TRAIN_feat_sel)[2]/3), num.trees = 300, 
                           importance= "impurity", write.forest = T)

## now make predictions on TEST

preds_feat_sel <- predict(mod_features_sel, TEST[,names(imp)])

preds_all_features <- predict(mod_ranger, TEST)

## 
plot_ml(preds_feat_sel$predictions, MET_te, index = "AMBI", aggreg = c("Grab", "Station", "Locality"), pdf = T,
        title = paste("Predictions_feature_selected", "AMBI", ".pdf", sep=""))
dev.off()

plot_ml(preds_all_features$predictions, MET_te, index = "AMBI", aggreg = c("Grab", "Station", "Locality"), pdf = T,
        title = paste("Predictions_all_feature_selected", "AMBI", ".pdf", sep=""))



#### caret



###
mod <- randomForest(MET_tr[,"AMBI"] ~ ., data=TRAIN, 
                    mtry=floor(dim(TRAIN)[2]/3), num.trees = 300)



