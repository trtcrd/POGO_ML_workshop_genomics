###############################################################
##########  plot after ML
###############################################################



plot_ml <- function(data, metadata, index, title = NULL, aggreg = NULL, pdf = F, taxo_group = NULL) {
  
  comp <- metadata
  combi1 <- data
  combi2 <- comp[,index]
  # taxo_group is used to export plots for a list of taxonomic groups
  if (is.null(taxo_group)) taxo_group <- "MLplots"
  
  ## agregate results by aggreg vector structure
  if (is.null(aggreg)) 
  {
    combined1 <- combi1
    combined2 <- combi2
  }
  if (!is.null(aggreg))
  {
    l <- length(aggreg) + 1 # for indexing the right columns
    if (length(aggreg)==1) 
    {
      combined1 <- aggregate(combi1, by=list(comp[,aggreg[1]]), FUN=mean)[,"x"]
      combined1_sd <- aggregate(combi1, by=list(comp[,aggreg[1]]), FUN=sd)[,"x"]
      combined1_sd[is.na(combined1_sd)] <- 0
      # know values
      combined2 <- aggregate(combi2, by=list(comp[,aggreg[1]]), FUN=mean)[,"x"]
      # farm name
      farm_nam <-  aggregate(combi1, by=list(comp[,aggreg[1]]), FUN=mean)[,"Group.1"]
      # col plot
      col <- aggregate(comp$col_plot, by=list(comp[,aggreg[1]]), FUN=mean)[,"x"]
    }
    if (length(aggreg)==2) 
    {
      combined1 <- aggregate(combi1, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"x"]
      combined1_sd <- aggregate(combi1, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=sd)[,"x"]
      combined1_sd[is.na(combined1_sd)] <- 0
      # know values
      combined2 <- aggregate(combi2, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"x"]
      # farm name
      farm_nam <-  aggregate(combi1, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"Group.2"]
      # col plot
      col <- aggregate(comp$col_plot, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=mean)[,"x"]
    }
    if (length(aggreg)==3) 
    {
      combined1 <- aggregate(combi1, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"x"]
      combined1_sd <- aggregate(combi1, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=sd)[,"x"]
      combined1_sd[is.na(combined1_sd)] <- 0
      # know values
      combined2 <- aggregate(combi2, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"x"]
      # farm name
      farm_nam <-  aggregate(combi1, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"Group.3"]
      # col plot
      col <- aggregate(comp$col_plot, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=mean)[,"x"]
    }
  }
  
  
  
  # predictions
  #combined1 <- aggregate(combi1, by=list(comp$Grab, comp$Station,comp$Locality), FUN=mean)[,"x"]
  #combined1_sd <- aggregate(combi1, by=list(comp$Grab, comp$Station,comp$Locality), FUN=sd)[,"x"]
  #combined1_sd[is.na(combined1_sd)] <- 0
  # know values
  #combined2 <- aggregate(combi2, by=list(comp$Grab, comp$Station,comp$Locality), FUN=mean)[,"x"]
  # keep farm_name vector
  #farm_nam <-  aggregate(combi1, by=list(comp$Grab, comp$Station,comp$Locality), FUN=mean)[,"Group.3"]
  
  ## colors
  #col <- aggregate(comp$col_plot, by=list(comp$Grab, comp$Station,comp$Locality), FUN=mean)[,"x"]
  
  # for storage
  comb1_disc <- c(rep(0, length(combined1)))
  comb2_disc <- c(rep(0, length(combined2)))
  
  ### if AMBI
  if (index == "AMBI") 
  {
    for (i in 1:length(combined1))
    {
      if (combined1[i] < 1.2) comb1_disc[i] <- 1
      else if (combined1[i] >= 1.2 && combined1[i] < 3.3) comb1_disc[i] <- 2
      else if (combined1[i] >= 3.3 && combined1[i] < 4.3) comb1_disc[i] <- 3
      else if (combined1[i] >= 4.3 && combined1[i] < 5.5) comb1_disc[i] <- 4
      else if (combined1[i] >= 5.5 && combined1[i] <= 6) comb1_disc[i] <- 5
      else if (combined1[i] >= 6) comb1_disc[i] <- 5       # if above 6 (it can happens with MXnet for instance)
      
      if (combined2[i] < 1.2) comb2_disc[i] <- 1
      else if (combined2[i] >= 1.2 && combined2[i] < 3.3) comb2_disc[i] <- 2
      else if (combined2[i] >= 3.3 && combined2[i] < 4.3) comb2_disc[i] <- 3
      else if (combined2[i] >= 4.3 && combined2[i] < 5.5) comb2_disc[i] <- 4
      else if (combined2[i] >= 5.5 && combined2[i] <= 6) comb2_disc[i] <- 5
      else if (combined2[i] >= 6) comb2_disc[i] <- 5       # just in case here
    }
    ## regression
    diff <- comb1_disc-comb2_disc 
    good_ass <- length(subset(diff, diff ==0))
    bad_ass <- length(subset(diff, diff !=0))
  }
  
  ### if microAMBI
  if (index == "microAMBI") 
  {
    for (i in 1:length(combined1))
    {
      if (combined1[i] < 1.2) comb1_disc[i] <- 1
      else if (combined1[i] >= 1.2 && combined1[i] < 3.3) comb1_disc[i] <- 2
      else if (combined1[i] >= 3.3 && combined1[i] < 4.3) comb1_disc[i] <- 3
      else if (combined1[i] >= 4.3 && combined1[i] < 5.5) comb1_disc[i] <- 4
      else if (combined1[i] >= 5.5 && combined1[i] <= 6) comb1_disc[i] <- 5
      else if (combined1[i] > 6) comb1_disc[i] <- 5 
      
      if (combined2[i] < 1.2) comb2_disc[i] <- 1
      else if (combined2[i] >= 1.2 && combined2[i] < 3.3) comb2_disc[i] <- 2
      else if (combined2[i] >= 3.3 && combined2[i] < 4.3) comb2_disc[i] <- 3
      else if (combined2[i] >= 4.3 && combined2[i] < 5.5) comb2_disc[i] <- 4
      else if (combined2[i] >= 5.5 && combined2[i] <= 6) comb2_disc[i] <- 5
      else if (combined2[i] > 6) comb2_disc[i] <- 5 
    }
    ## regression
    diff <- comb1_disc-comb2_disc 
    good_ass <- length(subset(diff, diff ==0))
    bad_ass <- length(subset(diff, diff !=0))
  }
  
  
  ### if NSI
  if (index == "NSI") 
  {
    for (i in 1:length(combined1))
    {
      if (combined1[i] < 10) comb1_disc[i] <- 1
      else if (combined1[i] >= 10 && combined1[i] < 15) comb1_disc[i] <- 2
      else if (combined1[i] >= 15 && combined1[i] < 20) comb1_disc[i] <- 3
      else if (combined1[i] >= 20 && combined1[i] < 25) comb1_disc[i] <- 4
      else if (combined1[i] >= 25 && combined1[i] <= 31) comb1_disc[i] <- 5
      else if (combined1[i] > 31) comb1_disc[i] <- 5 
      
      if (combined2[i] < 10) comb2_disc[i] <- 1
      else if (combined2[i] >= 10 && combined2[i] < 15) comb2_disc[i] <- 2
      else if (combined2[i] >= 15 && combined2[i] < 20) comb2_disc[i] <- 3
      else if (combined2[i] >= 20 && combined2[i] < 25) comb2_disc[i] <- 4
      else if (combined2[i] >= 25 && combined2[i] <= 31) comb2_disc[i] <- 5
      else if (combined2[i] > 31) comb2_disc[i] <- 5 
    }
    ## regression
    diff <- comb1_disc-comb2_disc 
    good_ass <- length(subset(diff, diff ==0))
    bad_ass <- length(subset(diff, diff !=0))
  }
  
  ### if NQI1
  if (index == "NQI1") 
  {
    for (i in 1:length(combined1))
    {
      if (combined1[i] < 0.31) comb1_disc[i] <- 1
      else if (combined1[i] >= 0.31 && combined1[i] < 0.49) comb1_disc[i] <- 2
      else if (combined1[i] >= 0.49 && combined1[i] < 0.63) comb1_disc[i] <- 3
      else if (combined1[i] >= 0.63 && combined1[i] < 0.82) comb1_disc[i] <- 4
      else if (combined1[i] >= 0.82 && combined1[i] <= 1) comb1_disc[i] <- 5
      else if (combined1[i] > 1) comb1_disc[i] <- 5 
      
      if (combined2[i] < 0.31) comb2_disc[i] <- 1
      else if (combined2[i] >= 0.31 && combined2[i] < 0.49) comb2_disc[i] <- 2
      else if (combined2[i] >= 0.49 && combined2[i] < 0.63) comb2_disc[i] <- 3
      else if (combined2[i] >= 0.63 && combined2[i] < 0.82) comb2_disc[i] <- 4
      else if (combined2[i] >= 0.82 && combined2[i] <= 1) comb2_disc[i] <- 5
      else if (combined2[i] > 1) comb2_disc[i] <- 5 
    }
    ## regression
    diff <- comb1_disc-comb2_disc 
    good_ass <- length(subset(diff, diff ==0))
    bad_ass <- length(subset(diff, diff !=0))
  }
  
  
  ### if ISI
  if (index == "ISI") 
  {
    for (i in 1:length(combined1))
    {
      if (combined1[i] < 4.5) comb1_disc[i] <- 1
      else if (combined1[i] >= 4.5 && combined1[i] < 6.1) comb1_disc[i] <- 2
      else if (combined1[i] >= 6.1 && combined1[i] < 7.5) comb1_disc[i] <- 3
      else if (combined1[i] >= 7.5 && combined1[i] < 9.6) comb1_disc[i] <- 4
      else if (combined1[i] >= 9.6 && combined1[i] <= 13) comb1_disc[i] <- 5
      else if (combined1[i] > 13) comb1_disc[i] <- 5 
      
      if (combined2[i] < 4.5) comb2_disc[i] <- 1
      else if (combined2[i] >= 4.5 && combined2[i] < 6.1) comb2_disc[i] <- 2
      else if (combined2[i] >= 6.1 && combined2[i] < 7.5) comb2_disc[i] <- 3
      else if (combined2[i] >= 7.5 && combined2[i] < 9.6) comb2_disc[i] <- 4
      else if (combined2[i] >= 9.6 && combined2[i] <= 13) comb2_disc[i] <- 5
      else if (combined2[i] > 13) comb2_disc[i] <- 5 
    }
    ## regression
    diff <- comb1_disc-comb2_disc 
    good_ass <- length(subset(diff, diff ==0))
    bad_ass <- length(subset(diff, diff !=0))
  }
  
  if (index == "Shannon") 
  {
    for (i in 1:length(combined1))
    {
      if (combined1[i] < 0.9) comb1_disc[i] <- 1
      else if (combined1[i] >= 0.9 && combined1[i] < 1.9) comb1_disc[i] <- 2
      else if (combined1[i] >= 1.9 && combined1[i] < 3) comb1_disc[i] <- 3
      else if (combined1[i] >= 3 && combined1[i] < 4.8) comb1_disc[i] <- 4
      else if (combined1[i] >= 4.8 && combined1[i] <= 5.7) comb1_disc[i] <- 5
      else if (combined1[i] > 5.7) comb1_disc[i] <- 5 
      
      if (combined2[i] < 0.9) comb2_disc[i] <- 1
      else if (combined2[i] >= 0.9 && combined2[i] < 1.9) comb2_disc[i] <- 2
      else if (combined2[i] >= 1.9 && combined2[i] < 3) comb2_disc[i] <- 3
      else if (combined2[i] >= 3 && combined2[i] < 4.8) comb2_disc[i] <- 4
      else if (combined2[i] >= 4.8 && combined2[i] <= 5.7) comb2_disc[i] <- 5
      else if (combined2[i] > 5.7) comb2_disc[i] <- 5 
    }
    ## regression
    diff <- comb1_disc-comb2_disc 
    good_ass <- length(subset(diff, diff ==0))
    bad_ass <- length(subset(diff, diff !=0))
  }
  
  
  ## classif 
  # diff <- combined1_rf-combined2_rf 
  # good_ass <- length(subset(diff, diff ==0))
  # bad_ass <- length(subset(diff, diff !=0))

  ##### PLOT 
  # AMBI
  if (index == "AMBI") 
  {
    if(pdf)
    {
      dir.create(file.path(getwd(), taxo_group), showWarnings = F)
      pdf(file = paste(taxo_group, "/", "AMBI.pdf", sep=""), width = 6.6, height=4, useDingbats=F)
    } else quartz(width = 6.6, height=4)
    #par(mfrow = c(1, 2))
    mat <- rbind(c(1,2,3), c(1,2,4))
    layout(mat, widths=c(2.5,1,1))
    
    mod <- lm (combined1 ~ combined2)
    if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
    if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
    if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
    if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
    
    # kappa
    kap <- kappa2(cbind(comb1_disc,comb2_disc), "squared", sort.levels=TRUE) 
    if (is.na(kap$p.value) == T) kap$p.value <- 1
    if (kap$p.value >= 0.05) sigk <- "ns"
    if (kap$p.value < 0.05) sigk <- "*"
    if (kap$p.value < 0.01) sigk <- "**"
    if (kap$p.value < 0.001) sigk <- "***"
    
    #col2 <- col[as.integer(as.factor(farm_nam))]
    #col <-  comp_foram$Col_plot
    plot(combined1 ~ combined2, xlim=c(0,6), ylim=c(0,6), pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("AMBI prediction /", title), xlab = "Morphology", ylab = "Molecular")
    
    arrows(combined2, combined1-combined1_sd, combined2, combined1+combined1_sd, length=0.01, angle=90, code=3, col=col)
    abline(mod, col="blue")
    rect(0,0,1.2,1.2,border="blue")
    rect(1.2,1.2,3.3,3.3,border="green")
    rect(3.3,3.3,4.3,4.3,border="yellow")
    rect(4.3,4.3,5.5,5.5,border="orange")
    rect(5.5,5.5,6,6,border="red")
    text(6,0.5, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(6,0.2, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    
    # second plot for legend
    par(mar=c(0,0,0,0), xpd=TRUE)
    plot(0,type="n", axes=F, xlab="", ylab="")
    legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
    
    par(mar=c(3,1,5,3))
    bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
    par(mar=c(5,1,2,3))
    bxplt <- boxplot(combined1 - combined2, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
    points(1, mean(combined1 - combined2), col="red", pch=3)
    text(x= 1.35, y= mean(combined1 - combined2), labels=paste(round(mean(combined1 - combined2), 2)), xpd=TRUE, cex=0.7, col="red")
    if(pdf) dev.off()
    
  }
  
  
  # AMBI
  if (index == "microAMBI") 
  {
    if(pdf)
    {
      dir.create(file.path(getwd(), taxo_group), showWarnings = F)
      pdf(file = paste(taxo_group, "/", "microAMBI.pdf", sep=""), width = 6.6, height=4, useDingbats=F)
    } else quartz(width = 6.6, height=4)
    #par(mfrow = c(1, 2))
    mat <- rbind(c(1,2,3), c(1,2,4))
    layout(mat, widths=c(2.5,1,1))
    
    mod <- lm (combined1 ~ combined2)
    if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
    if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
    if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
    if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
    
    # kappa
    kap <- kappa2(cbind(comb1_disc,comb2_disc), "squared", sort.levels=TRUE) 
    if (is.na(kap$p.value) == T) kap$p.value <- 1
    if (kap$p.value >= 0.05) sigk <- "ns"
    if (kap$p.value < 0.05) sigk <- "*"
    if (kap$p.value < 0.01) sigk <- "**"
    if (kap$p.value < 0.001) sigk <- "***"
    
    #col2 <- col[as.integer(as.factor(farm_nam))]
    #col <-  comp_foram$Col_plot
    plot(combined1 ~ combined2, xlim=c(0,6), ylim=c(0,6), pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("microAMBI prediction /", title), xlab = "Morphology", ylab = "Molecular")
    
    arrows(combined2, combined1-combined1_sd, combined2, combined1+combined1_sd, length=0.01, angle=90, code=3, col=col)
    abline(mod, col="blue")
    rect(0,0,1.2,1.2,border="blue")
    rect(1.2,1.2,3.3,3.3,border="green")
    rect(3.3,3.3,4.3,4.3,border="yellow")
    rect(4.3,4.3,5.5,5.5,border="orange")
    rect(5.5,5.5,6,6,border="red")
    text(6,0.5, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(6,0.2, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    
    # second plot for legend
    par(mar=c(0,0,0,0), xpd=TRUE)
    plot(0,type="n", axes=F, xlab="", ylab="")
    legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
    
    par(mar=c(3,1,5,3))
    bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
    par(mar=c(5,1,2,3))
    bxplt <- boxplot(combined1 - combined2, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
    points(1, mean(combined1 - combined2), col="red", pch=3)
    text(x= 1.35, y= mean(combined1 - combined2), labels=paste(round(mean(combined1 - combined2), 2)), xpd=TRUE, cex=0.7, col="red")
    if(pdf) dev.off()
    
  }
  
  
  ##NSI
  if (index == "NSI") 
  {
    if(pdf)
    {
      dir.create(file.path(getwd(), taxo_group), showWarnings = F)
      pdf(file = paste(taxo_group, "/", "NSI.pdf", sep=""), width = 6.6, height=4, useDingbats=F)
    } else quartz(width = 6.6, height=4)
    #par(mfrow = c(1, 2))
    mat <- rbind(c(1,2,3), c(1,2,4))
    layout(mat, widths=c(2.5,1,1))
    
    mod <- lm (combined1 ~ combined2)
    if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
    if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
    if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
    if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
    
    # kappa
    kap <- kappa2(cbind(comb1_disc,comb2_disc), "squared", sort.levels=TRUE) 
    if (is.na(kap$p.value) == T) kap$p.value <- 1
    if (kap$p.value >= 0.05) sigk <- "ns"
    if (kap$p.value < 0.05) sigk <- "*"
    if (kap$p.value < 0.01) sigk <- "**"
    if (kap$p.value < 0.001) sigk <- "***"
    
    plot(combined1 ~ combined2, xlim=c(0,31), ylim=c(0,31), pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("NSI prediction /", title), xlab = "Morphology", ylab = "Molecular")
    arrows(combined2, combined1-combined1_sd, combined2, combined1+combined1_sd, length=0.01, angle=90, code=3, col=col)
    abline(mod, col="blue")
    rect(0,0,10,10,border="red")
    rect(10,10,15,15,border="orange")
    rect(15,15,20,20,border="yellow")
    rect(20,20,25,25,border="green")
    rect(25,25,31,31,border="blue")
    text(31,2.6, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(31,1, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    # second plot for legend
    par(mar=c(0,0,0,0), xpd=TRUE)
    plot(0,type="n", axes=F, xlab="", ylab="")
    legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
    
    par(mar=c(3,1,5,3))
    bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
    par(mar=c(5,1,2,3))
    bxplt <- boxplot(combined1 - combined2, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
    points(1, mean(combined1 - combined2), col="red", pch=3)
    text(x= 1.35, y= mean(combined1 - combined2), labels=paste(round(mean(combined1 - combined2), 2)), xpd=TRUE, cex=0.7, col="red")
    if(pdf) dev.off()
  }
  
  
  
  ### NQI1
  if (index == "NQI1") 
  {
    if(pdf)
    {
      dir.create(file.path(getwd(), taxo_group), showWarnings = F)
      pdf(file = paste(taxo_group, "/", "NQI1.pdf", sep=""), width = 6.6, height=4, useDingbats=F)
    } else quartz(width = 6.6, height=4)
    #par(mfrow = c(1, 2))
    mat <- rbind(c(1,2,3), c(1,2,4))
    layout(mat, widths=c(2.5,1,1))
    
    mod <- lm (combined1 ~ combined2)
    if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
    if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
    if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
    if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
    
    # kappa
    kap <- kappa2(cbind(comb1_disc,comb2_disc), "squared", sort.levels=TRUE) 
    if (is.na(kap$p.value) == T) kap$p.value <- 1
    if (kap$p.value >= 0.05) sigk <- "ns"
    if (kap$p.value < 0.05) sigk <- "*"
    if (kap$p.value < 0.01) sigk <- "**"
    if (kap$p.value < 0.001) sigk <- "***"
    
    plot(combined1 ~ combined2, xlim=c(0,1), ylim=c(0,1), pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("NQI1 prediction /", title), xlab = "Morphology", ylab = "Molecular")
    arrows(combined2, combined1-combined1_sd, combined2, combined1+combined1_sd, length=0.01, angle=90, code=3, col=col)
    abline(mod, col="blue")
    rect(0,0,.31,.31,border="red")
    rect(.31,.31,.49,.49,border="orange")
    rect(.49,.49,.63,.63,border="yellow")
    rect(.63,.63,.82,.82,border="green")
    rect(.82,.82,1,1,border="blue")
    
    text(1,.1, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(1,.05, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    
    # second plot for legend
    par(mar=c(0,0,0,0), xpd=TRUE)
    plot(0,type="n", axes=F, xlab="", ylab="")
    legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
    
    par(mar=c(3,1,5,3))
    bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
    par(mar=c(5,1,2,3))
    bxplt <- boxplot(combined1 - combined2, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
    points(1, mean(combined1 - combined2), col="red", pch=3)
    text(x= 1.35, y= mean(combined1 - combined2), labels=paste(round(mean(combined1 - combined2), 2)), xpd=TRUE, cex=0.7, col="red")
    if(pdf) dev.off()
  }
  
  
  
  #### ISI 
  if (index == "ISI") 
  {
    if(pdf)
    {
      dir.create(file.path(getwd(), taxo_group), showWarnings = F)
      pdf(file = paste(taxo_group, "/", "ISI.pdf", sep=""), width = 6.6, height=4, useDingbats=F)
    } else quartz(width = 6.6, height=4)
    #par(mfrow = c(1, 2))
    mat <- rbind(c(1,2,3), c(1,2,4))
    layout(mat, widths=c(2.5,1,1))
    
    mod <- lm (combined1 ~ combined2)
    if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
    if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
    if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
    if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
    
    # kappa
    kap <- kappa2(cbind(comb1_disc,comb2_disc), "squared", sort.levels=TRUE) 
    if (is.na(kap$p.value) == T) kap$p.value <- 1
    if (kap$p.value >= 0.05) sigk <- "ns"
    if (kap$p.value < 0.05) sigk <- "*"
    if (kap$p.value < 0.01) sigk <- "**"
    if (kap$p.value < 0.001) sigk <- "***"
    
    plot(combined1 ~ combined2, xlim=c(0,13), ylim=c(0,13), pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("ISI prediction /", title), xlab = "Morphology", ylab = "Molecular")
    arrows(combined2, combined1-combined1_sd, combined2, combined1+combined1_sd, length=0.01, angle=90, code=3, col=col)
    abline(mod, col="blue")
    rect(0,0,4.5,4.5,border="red")
    rect(4.5,4.5,6.1,6.1,border="orange")
    rect(6.1,6.1,7.5,7.5,border="yellow")
    rect(7.5,7.5,9.6,9.6,border="green")
    rect(9.6,9.6,13,13,border="blue")
    
    text(13,1.1, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(13,.4, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    
    # second plot for legend
    par(mar=c(0,0,0,0), xpd=TRUE)
    plot(0,type="n", axes=F, xlab="", ylab="")
    legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
    
    par(mar=c(3,1,5,3))
    bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
    par(mar=c(5,1,2,3))
    bxplt <- boxplot(combined1 - combined2, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
    points(1, mean(combined1 - combined2), col="red", pch=3)
    text(x= 1.35, y= mean(combined1 - combined2), labels=paste(round(mean(combined1 - combined2), 2)), xpd=TRUE, cex=0.7, col="red")
    if(pdf) dev.off()
  }
  
  #### Shannon 
  if (index == "Shannon") 
  {
    if(pdf)
    {
      dir.create(file.path(getwd(), taxo_group), showWarnings = F)
      pdf(file = paste(taxo_group, "/", "Shannon.pdf", sep=""), width = 6.6, height=4, useDingbats=F)
    } else quartz(width = 6.6, height=4)
    #par(mfrow = c(1, 2))
    mat <- rbind(c(1,2,3), c(1,2,4))
    layout(mat, widths=c(2.5,1,1))
    
    mod <- lm (combined1 ~ combined2)
    if (anova(mod)[["Pr(>F)"]][1] >= 0.05) sig <- "ns"
    if (anova(mod)[["Pr(>F)"]][1] < 0.05) sig <- "*"
    if (anova(mod)[["Pr(>F)"]][1] < 0.01) sig <- "**"
    if (anova(mod)[["Pr(>F)"]][1] < 0.001) sig <- "***"
    
    # kappa
    kap <- kappa2(cbind(comb1_disc,comb2_disc), "squared", sort.levels=TRUE) 
    if (is.na(kap$p.value) == T) kap$p.value <- 1
    if (kap$p.value >= 0.05) sigk <- "ns"
    if (kap$p.value < 0.05) sigk <- "*"
    if (kap$p.value < 0.01) sigk <- "**"
    if (kap$p.value < 0.001) sigk <- "***"
    
    plot(combined1 ~ combined2, xlim=c(0,5.7), ylim=c(0,5.7), pch =col, col= col,  xaxs="i",yaxs="i",cex=1, main = paste("Shannon prediction /", title), xlab = "Morphology", ylab = "Molecular")
    arrows(combined2, combined1-combined1_sd, combined2, combined1+combined1_sd, length=0.01, angle=90, code=3, col=col)
    abline(mod, col="blue")
    rect(0,0,0.9,0.9,border="red")
    rect(0.9,0.9,1.9,1.9,border="orange")
    rect(1.9,1.9,3,3,border="yellow")
    rect(3,3,4.8,4.8,border="green")
    rect(4.8,4.8,5.7,5.7,border="blue")
    
    text(5.7,0.5, pos=2,paste("R²=", round(summary(mod)$adj.r.squared, 3), sig, sep=""))
    text(5.7,0.2, pos=2,paste("Kappa=", round(kap$value, 3), sigk, sep=""))
    
    # second plot for legend
    par(mar=c(0,0,0,0), xpd=TRUE)
    plot(0,type="n", axes=F, xlab="", ylab="")
    legend("topleft", as.vector(unique(farm_nam)), col= unique(col), pch=unique(col), box.lty=0, inset=c(0,0.11))
    
    par(mar=c(3,1,5,3))
    bplt <- barplot(prop.table(table(diff))*100, ylim=c(0,100), xlab="Status mismatch", ylab="Percentage", cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= bplt, y= prop.table(table(diff))*100 + 5, labels=paste("n=", table(diff), sep=""), xpd=TRUE, cex=0.7)
    par(mar=c(5,1,2,3))
    bxplt <- boxplot(combined1 - combined2, cex.names = 0.7, cex.axis=0.7,cex.lab=0.5)
    text(x= 0.65, y= bxplt$stats, labels=paste(round(bxplt$stats, 2)), xpd=TRUE, cex=0.7)
    points(1, mean(combined1 - combined2), col="red", pch=3)
    text(x= 1.35, y= mean(combined1 - combined2), labels=paste(round(mean(combined1 - combined2), 2)), xpd=TRUE, cex=0.7, col="red")
    if(pdf) dev.off()
  }
  
  return(list("preds_cont" = combined1, "bi_values"= combined2, "preds_cat" = comb1_disc, "labels" = comb2_disc, "R2" = paste(round(summary(mod)$adj.r.squared, 3), sig, sep=""), "KAP" = paste(round(kap$value, 3), sigk, sep="")))
  
}
  
