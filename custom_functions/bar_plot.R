#### barplot of taxonomic rank relative abundance

bar_plot <- function(otu_table, comp_, aggreg = F, taxo_file, title_, font_size_ = 0.8, horizontal_ = T,
                     tax_rank, metazoan = F, lastrank = F , pdf = F, taxo_group = NULL, file_name = NULL, return_levels = F,
                     remove_rare = T, remove_unass = F, distinct_colors = F) {

  ### in taxo_file what is the name of the last column??
  
  otu <- otu_table
  comp <- comp_
  tax <- taxo_file
  title <- title_
  font_size <- font_size_
  horizontal <- horizontal_
  tax_rank <- tax_rank
  metazoan <- metazoan
  lastrank <- lastrank
  pdf <- pdf
  taxo_group <- taxo_group
  
  if (lastrank) { taxRank <- "last" } else { taxRank <- tax_rank }
  
  ## title
  if (lastrank) title <- paste(title, "(last rank)")
  if (!lastrank) title <- paste(title, " (rank ", tax_rank, ")", sep="")
  # 
  # otu <- OTUtpPL
  # comp <- COMP
  # tax <- TAXOtp
  # title <- paste(k,taxo_group, "-", ncol(OTUtpML),"OTUs","-",sum(OTUtpPL), "reads")
  # font_size <- 0.8
  # horizontal_ <- T
  # tax_rank <- tax_rank
  # metazoan <- F
  # lastrank <- F
  # pdf <- T
  # taxo_group <- k

  # need to remove [ ] or weird character in taxo file if any
  for (i in 1:length(tax)) tax[i] <- gsub("[^[:alnum:][:blank:]_;+?&/\\-]", "", tax[i], c) 

  # tax_rank 2 if forams..
  # tax rank 7 if metazoa 
  # tax rank 2 if bacteria
  # tax rank 4 if eukaryotes
  
  ## if metazoan keep only metazoan
  if (metazoan) 
  {
    metaz <- grep("Metazoa", tax[,"taxon"])
    otu <- otu[,metaz]
    tax <-tax[metaz,]
  }
  
  ## if remove unassigned
  if (remove_unass) 
  {
    unass <- grep("Unassigned", tax[,"taxon"], ignore.case = T)
    otu <- otu[,-unass]
    tax <-tax[-unass,]
    #tax_rank <- 7
  }
  
  # merging the data by grab 
  #otu_ <- aggregate(otu, by=list(comp, comp$Station, comp$Locality), FUN=sum)
  #otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=sum)
  
  if (aggreg==F) dat <- otu
  if (aggreg!=F)
  {
    l <- length(aggreg) + 1 # for indexing the right columns
    if (length(aggreg)==1) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]]), FUN=sum)
    if (length(aggreg)==2) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=sum)
    if (length(aggreg)==3) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=sum)
    if (length(aggreg)==4) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]], comp[,aggreg[4]]), FUN=sum)
    if (length(aggreg)==5) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]], comp[,aggreg[4]], comp[,aggreg[5]]), FUN=sum)
    # and reformatting to get single field with locality > station > grab
    dat <- otu_[,l:dim(otu_)[[2]]]
  }
  ### then sorting by station and locality 
  #dat <- otu_[with(otu_, order(Group.2, Group.3)),]
  
  # and repasting the names with all infos
  ## if reordering by station and locality ## rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[with(otu_, order(Group.2, Group.3)),c("Group.3","Group.2","Group.1")]), stringsAsFactors=FALSE))
  if (aggreg!=F)
  {
    if (length(aggreg)==1) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.1")]), stringsAsFactors=FALSE))
    if (length(aggreg)==2) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.2","Group.1")]), stringsAsFactors=FALSE))
    if (length(aggreg)==3) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.3","Group.2","Group.1")]), stringsAsFactors=FALSE))
    if (length(aggreg)==4) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.4", "Group.3","Group.2","Group.1")]), stringsAsFactors=FALSE))
    if (length(aggreg)==5) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.5", "Group.4", "Group.3","Group.2","Group.1")]), stringsAsFactors=FALSE))
  }
  # extract the "taxo_rank" assignment
  if (lastrank)
  {
    taxa <- seq(from = 1, to = dim(dat)[2])
    for (i in 1:dim(dat)[2])
    {
      if (tax[i,"taxon"] == "") tax[i,"taxon"] <- "Unassigned"
      tmp <- tail(unlist(strsplit(as.character(tax[i,"taxon"]), split=";", fixed=TRUE)), 1)
      taxa[i] <- tmp
      if (tmp == "Unassigned") taxa[i] <- "zz_Unassigned"
    }
  } else
  {
    taxa <- seq(from = 1, to = dim(dat)[2])
    for (i in 1:dim(dat)[2])
    {
      ### if the taxonomic assignment is empty -> Unassigned
      if (tax[i,"taxon"] == "") tax[i,"taxon"] <- "Unassigned" 
      ### extract the tax_rank of the assignment
      tmp <- unlist(strsplit(as.character(tax[i,"taxon"]), split=";", fixed=TRUE))[tax_rank] 
      ### if it is not empty
      if (is.na(tmp) == FALSE)
      {
        taxa[i] <- tmp
        ## zz_ to make it appear the last in the plot
        if (tmp == "Unassigned") taxa[i] <- "zz_Unassigned"
      } else {
        ### extract only the tax_rank, if lower, it goes to unassigned...
        # tmp <- tail(unlist(strsplit(as.character(tax[i,"Taxa"]), split=";", fixed=TRUE)), 1) ### previous code
        #tmp <- unlist(strsplit(as.character(tax[i,"Taxa"]), split=";", fixed=TRUE))[tax_rank]
        #if (is.na(tmp) == TRUE) taxa[i] <- "zz_Unassigned"
        #if (tmp == "Unassigned") taxa[i] <- "zz_Unassigned"
        taxa[i] <- "zz_Unassigned"
      }
    }
  }
  # creating a table for sum family
  dat_ <- as.data.frame(array(NA, c(dim(dat)[1], length(levels(as.factor(taxa))))))
  dimnames(dat_)[[1]] <- dimnames(dat)[[1]]
  ll <- levels(as.factor(taxa))
  for (i in 1:length(ll)) ll[i] <- strsplit(ll[i], split="(", fixed=TRUE)[[1]][1]    # need to remove parenthesis if any
  dimnames(dat_)[[2]] <- ll
  
  # summing the abundance by tax_rank
  for (i in ll)
  {
    tt  <- dat[,grep(i, taxa)]
    if (is.data.frame(tt)) { dat_[,i] <- rowSums(tt) } else { dat_[,i] <- tt }
  }
  
  # normalizing the otu table for % 
  dat_ <- dat_ / rowSums(dat_) *100
  dat_[dat_ == "NaN"] <- 0
  
  ## remove rare families
  if (remove_rare) 
  {
    if (ncol(dat_) > 1) dat_ <-  dat_[,colSums(dat_)>1]
    # REnormalizing the otu table for % 
    dat_ <- dat_ / rowSums(dat_) *100
    dat_[dat_ == "NaN"] <- 0
  }
  
  
  #barplot
  
  ## most distinct colors (random) 
  color <-  grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  if (distinct_colors) col <- sample(color, dim(dat_)[2])
  
  ## color ramp
  if (!distinct_colors) col <- colorRampPalette(c("black", "yellow", "blue", "darkorchid1", "red", "green", "orange", "white", "lightblue1", "lavenderblush", "olivedrab4", "seagreen1", "grey"), bias=1, interpolate = "linear")(dim(dat_)[2])

  #col <- colorRampPalette("palette")(dim(dat_)[2])
  
  if (horizontal_) 
  {
    if (pdf)
    {
      if (is.null(taxo_group) & is.null(file_name)) 
      {
        pdf(file = "taxo_group_plot.pdf")
      } else if (is.null(taxo_group) & length(file_name) > 0)
      {
        pdf(file = paste(file_name,".pdf", sep=""))
      } else {
        dir.create(file.path(getwd(), taxo_group))
        pdf(file = paste(taxo_group, "/", paste(taxo_group, "_rank_", taxRank, "_plot.pdf", sep=""), sep=""))
      }
      par(mar=c(4, 8, 2, 10), xpd=TRUE, lwd = 0.5)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=1, horiz=TRUE, cex.names = font_size, border = T, xlab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      if(remove_unass & legg[length(legg)] == "zz_Unassigned") legg[length(legg)] <- paste("Unassigned (rank ", tax_rank, ")", sep="") 
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.55,0.02), box.lty=0, border = T)
      dev.off()
    } else
    {
      quartz()
      par(mar=c(4, 8, 2, 10), xpd=TRUE, lwd = 0.5)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=1, horiz=TRUE, cex.names = font_size, border = T, xlab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      if(remove_unass & legg[length(legg)] == "zz_Unassigned") legg[length(legg)] <- paste("Unassigned (rank ", tax_rank, ")", sep="") 
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.55,0.02), box.lty=0, border = T)
    }
  } else {
    if (pdf)
    {
      if (is.null(taxo_group)) 
      {
        pdf(file = "taxo_group_plot.pdf")
      } else {
        dir.create(file.path(getwd(), taxo_group))
        pdf(file = paste(taxo_group, "/", paste(taxo_group, "_rank_", taxRank, "_plot.pdf", sep=""), sep=""))
      }
      par(mar=c(8, 4, 2, 10), xpd=TRUE, lwd = 0.5)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=3, horiz=FALSE, cex.names = font_size, border = T, ylab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      if(remove_unass & legg[length(legg)] == "zz_Unassigned") legg[length(legg)] <- paste("Unassigned (rank ", tax_rank, ")", sep="") 
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.35,0), box.lty=0, border = T)
      dev.off()
    } else
    {
      quartz()
      par(mar=c(8, 4, 2, 10), xpd=TRUE, lwd = 0.5)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=3, horiz=FALSE, cex.names = font_size, border = T, ylab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      if(remove_unass & legg[length(legg)] == "zz_Unassigned") legg[length(legg)] <- paste("Unassigned (rank ", tax_rank, ")", sep="") 
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.35,0), box.lty=0, border = T)
    }
  }
  
  if (return_levels) return(levels(as.factor(dimnames(dat_)[[2]])))
  
}
