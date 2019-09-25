###############################################################
##########  plot after ML
###############################################################
source("R/utils/eco_groups.R")

require(irr)

plot_ml <- function(data, metadata, xIndex, yIndex, title = NULL, aggreg = NULL, pdf = F, taxo_group = NULL) {
  
  comp <- metadata
  predictions <- data
  trainingSet <- comp[,xIndex]
  # taxo_group is used to export plots for a list of taxonomic groups
  if (is.null(taxo_group)) taxo_group <- "MLplots"
  
  ## agregate results by aggreg vector structure
  if (is.null(aggreg)) 
  {
    pred.ag <- predictions
    train.ag <- trainingSet
  }
  if (!is.null(aggreg))
  {
    l <- length(aggreg) + 1 # for indexing the right columns
    if (length(aggreg)==1) 
    {
      pred.ag <- aggregate(predictions, by=list(comp[,aggreg[1]]), FUN=mean)[,"x"]
      pred.ag_sd <- aggregate(predictions, by=list(comp[,aggreg[1]]), FUN=sd)[,"x"]
      pred.ag_sd[is.na(pred.ag_sd)] <- 0
      # know values
      train.ag <- aggregate(trainingSet, by=list(comp[,aggreg[1]]), FUN=mean)[,"x"]
      # farm name
      farm_nam <-  aggregate(predictions, by=list(comp[,aggreg[1]]), FUN=mean)[,"Group.1"]
      # col plot
      col <- aggregate(comp$col_plot, by=list(comp[,aggreg[1]]), FUN=mean)[,"x"]
    }
    if (length(aggreg)==2) 
    {
      pred.ag <- aggregate(predictions, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"x"]
      pred.ag_sd <- aggregate(predictions, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=sd)[,"x"]
      pred.ag_sd[is.na(pred.ag_sd)] <- 0
      # know values
      train.ag <- aggregate(trainingSet, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"x"]
      # farm name
      farm_nam <-  aggregate(predictions, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"Group.2"]
      # col plot
      col <- aggregate(comp$col_plot, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"x"]
    }
    if (length(aggreg)==3) 
    {
      pred.ag <- aggregate(predictions, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"x"]
      pred.ag_sd <- aggregate(predictions, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=sd)[,"x"]
      pred.ag_sd[is.na(pred.ag_sd)] <- 0
      # know values
      train.ag <- aggregate(trainingSet, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"x"]
      # farm name
      farm_nam <-  aggregate(predictions, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"Group.3"]
      # col plot
      col <- aggregate(comp$col_plot, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"x"]
    }
  }
  # for storage
  pred.cat <- c(rep(0, length(pred.ag)))
  train.cat <- c(rep(0, length(train.ag)))
  
  ### Get groups using eco_groups to-be package (#TODO)
  
  for (i in 1:length(pred.ag)){
    pred.cat[i] <- getBIGroupFromValue(pred.ag[i], bi=yIndex)
    train.cat[i] <- getBIGroupFromValue(train.ag[i], bi=xIndex)
  }

  
  ## regression
  diff <- pred.cat-train.cat 
  good_ass <- length(subset(diff, diff ==0))
  bad_ass <- length(subset(diff, diff !=0))
  

  ##### PLOT 
  
  quartz(width = 6.6, height=4)
  #par(mfrow = c(1, 2))
  mat <- rbind(c(1,2,3), c(1,2,4))
  layout(mat, widths=c(2.5,1,1))
  
  mod <- lm (pred.ag ~ train.ag)
  if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
  if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
  if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
  if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
  
  # kappa
  kap <- kappa2(cbind(pred.cat,train.cat), "squared", sort.levels=TRUE)
  if (is.na(kap$p.value) == T) kap$p.value <- 1
  if (kap$p.value >= 0.05) sigk <- "ns"
  if (kap$p.value < 0.05) sigk <- "*"
  if (kap$p.value < 0.01) sigk <- "**"
  if (kap$p.value < 0.001) sigk <- "***"
  
  plot(pred.ag ~ train.ag, xlim=c(getWeight(1,xIndex), getWeight(5,xIndex)), 
       ylim=c(getWeight(1,yIndex),getWeight(5,yIndex)), 
       pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("AMBI prediction /", title), 
       xlab = "Morphology", ylab = "Molecular")
  
  arrows(train.ag, pred.ag-pred.ag_sd, train.ag, pred.ag+pred.ag_sd, length=0.01, angle=90, code=3, col=col)
  abline(mod, col="blue")
  
  
  groupCols = c("blue","green","yellow","orange","red")
  
  if (yIndex == "Redox" | yIndex == "Redox5yAvg" | yIndex == "Redox10yAvg") {
    
    # Draw boxes around categories  = eco-groups
    for (i in c(1:5)){
      rect(-ECO_GROUPS["Redox",i+1] ,-ECO_GROUPS["Redox",i+1],
         -ECO_GROUPS["Redox",i],-ECO_GROUPS["Redox",i],
         border=groupCols[i]) 
    }
    text(700,-200, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(700,-250, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    
  } else {
    xI = getSyn(xIndex)
    yI = getSyn(yIndex)
    for (i in c(1:5)){
      rect(ECO_GROUPS[xI,i] ,ECO_GROUPS[yI,i],
           ECO_GROUPS[xI,i+1],ECO_GROUPS[yI,i+1],
           border=groupCols[i])  
    }
    # Print out R square and Kappa
    text(ECO_GROUPS[xI,6],ECO_GROUPS[yI,6]/12, 
         pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(ECO_GROUPS[xI,6],ECO_GROUPS[yI,6]/30, 
         pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
  }
  
  # second plot for legend
  par(mar=c(0,0,0,0), xpd=TRUE)
  plot(0,type="n", axes=F, xlab="", ylab="")
  legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
  
  par(mar=c(3,1,5,3))
  bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
  text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
  par(mar=c(5,1,2,3))
  bxplt <- boxplot(pred.ag - train.ag, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
  text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
  points(1, mean(pred.ag - train.ag), col="red", pch=3)
  text(x= 1.35, y= mean(pred.ag - train.ag), labels=paste(round(mean(pred.ag - train.ag), 2)), xpd=TRUE, cex=0.7, col="red")
  

  
  return(list("preds_cont" = pred.ag, "bi_values"= train.ag, "preds_cat" = NA, "labels" =NA, "R2" = paste(round(summary(mod)$adj.r.squared, 3), sig, sep=""), "KAP" = NA))
  
}
  
