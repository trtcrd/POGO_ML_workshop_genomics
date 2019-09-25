##############################################################
# Helper functions to get ecogroups based on their intervals #
# Anders Lanzen 2019-07-17                                   #
##############################################################

SYNONYMS = data.frame(row.names = c("PI5yAvg","PI10yAvg","Redox5yAvg","Redox10yAvg",
                               "OM5yAvg","OM10yAvg","AMBIGroup","piMetals","piHC","piOM","piOM5yAvg"),
                      bi = c("PI","PI","Redox","Redox","OM","OM","AMBI",rep("PI",4)))

# Ranges and weights for all eco-groups
# Redox values in minus e.g. +700 --> -700!
ECO_GROUPS = data.frame(row.names=c("AMBI","microgAMBI", "PI","Redox","OM","NSI","NQ1","ISI"),
                        gIMin=c(0,0,0,-700,0,0,0,0),
                        gIIMin=c(1.2,1.2,1,-500,1,10,.31,4.5),
                        gIIIMin=c(3.3,2.4,2,-300,2,15,.49,6.1),
                        gIVMin=c(4.3,3.6,3,-100,4,20,.63,7.5),
                        gVMin=c(5,4.8,4,100,8,25,.82,9.6),
                        gVMax=c(6,6,5,300,25,31,1,13),
                        gIWeight = c(0,0,0,500,0,0,0,0),
                        gIIWeight = c(1.5,1.5,1.2,300,1.5,NA,NA,NA),
                        gIIIWeight = c(3,2.2,2,0,3,NA,NA,NA),
                        gIVWeight = c(4.5,3.8,3.5,-100,6,NA,NA,NA),
                        gVWeight = c(6,6,5,-300,25,31,1,13))

## Find which BI group that a specific (predicted) value is in
getBIGroupFromValue = function(value,bi="AMBI") {
  if (!bi %in% row.names(ECO_GROUPS)) {
    if (bi %in% row.names(SYNONYMS)) {
      bi = as.character(SYNONYMS[bi,1])
    } else {
      return(NA)
    }
  }
  
  if (bi=="Redox") value = -value
  i=1
  while (i<6 & value >= ECO_GROUPS[bi,i+1]) i=i+1
  return (i)
}
  
## Find closest organism EC group based on peak
  
getECFromPeak = function(value, bi="AMBI"){
  if (!bi %in% row.names(ECO_GROUPS)) {
    if (bi %in% row.names(SYNONYMS)) {
      bi = as.character(SYNONYMS[bi,1])
    } else {
      return(NA)
    }
  }
  diffs = abs(value - ECO_GROUPS[bi,c(7:11)])
  ec = c(1:5)[diffs==min(diffs)]
  if (length(ec)>1) ec = ec[1]
  return (ec)
}

getWeight= function(group, bi="AMBI"){
  
  if (!group %in% c(1:5)) return (NA)
  
  if (!bi %in% row.names(ECO_GROUPS)) {
    if (bi %in% row.names(SYNONYMS)) {
      bi = as.character(SYNONYMS[bi,1])
    } else {
      return(NA)
    }
  }
  return (ECO_GROUPS[bi,group+6])
}

getSyn= function(bi){
  
  
  if (!bi %in% row.names(ECO_GROUPS)) {
    if (bi %in% row.names(SYNONYMS)) {
      return(as.character(rSYNONYMS[bi,1]))
    } else {
      return(NA)
    }
  }
  return (bi)
}
  