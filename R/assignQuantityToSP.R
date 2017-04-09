assignQuantityToSP<-function (cbs, dm, C=list(sps=c("SP","SP_cnv"),pms = c("PM","PM_cnv")), ambig = F) 
{
  print("Assigning copy number to SPs...")
  dm[, C$sps]=round(1000*dm[, C$sps])/1000
  sps = sort(unique(as.numeric(dm[, C$sps])))
  cols = c( colnames(cbs), paste("SP_", as.character(unique(sps)), sep = "") )
  ploidy = matrix(cbind(cbs, matrix(NaN, nrow(cbs), length(sps))), 
                  nrow = nrow(cbs), ncol = length(cols), dimnames = list(1:nrow(cbs), cols))
  
  toD = c();  allAssigned=c();
  for (k in 1:nrow(ploidy)) {
    if (mod(k, 100) == 0) {
      print(paste("Finding overlaps for CBS segment", k, "out of ", nrow(ploidy), "..."))
    }
    idx = which(dm[, "chr"] == ploidy[k, "chr"] & dm[, "startpos"] >= 
                  ploidy[k, "startpos"] & dm[, "startpos"] <= ploidy[k, "endpos"])
    if (length(idx) == 0) {
      next
    }
    for (j in idx) { ##For each mutated locus
      dmx = dm[j, ]
      for(cI in 1:length(C$sps)){ ##For each subpopulation affected by this locus
        if (is.na(C$sps[cI])) {
          next
        }
        sp = paste("SP_", dmx[C$sps[cI]], sep = "")
        thisAssigned=c(k,dmx[C$sps[cI]] ,ploidy[k,"chr"],ploidy[k,"startpos"],ploidy[k,"endpos"],dmx[C$pms[cI]]);
        allAssigned = rbind(allAssigned, thisAssigned)
        if (is.na(ploidy[k, sp]) || ploidy[k, sp] == as.double(dmx[C$pms[cI]])) {
          ploidy[k, sp] = as.double(dmx[C$pms[cI]])
        } else {
          toD = rbind(toD, thisAssigned)
        }
      }
      
    }
  }
  
  
  ######################
  ##Either remove ambiguous segments or calculate their median ploidy based on SNV ploidy 
  colName=C$pms[1];
  printErr=FALSE;
  if(length(toD)>0){
    rownames(toD)=c(1:nrow(toD))
    colnames(toD)=c("Idx","SP","chr","startpos","endpos",colName)
    colnames(allAssigned)=colnames(toD)  
    uD=as.matrix(unique(toD));
    for (i in 1:nrow(uD)){
      sp = paste("SP_", as.numeric(uD[i,"SP"]), sep = "")
      if (!ambig){
        ploidy[as.numeric(uD[i,"Idx"]),sp]=NA;
        printErr=TRUE;    
      }else{
        ii=which(as.numeric(uD[i,"SP"])==as.numeric(allAssigned[,"SP"]) & 
                   as.numeric(uD[i,"chr"])==as.numeric(allAssigned[,"chr"]) &
                   as.numeric(uD[i,"startpos"])==as.numeric(allAssigned[,"startpos"]) &
                   as.numeric(uD[i,"endpos"])==as.numeric(allAssigned[,"endpos"]));
        ploidy[as.numeric(uD[i,"Idx"]),sp]=round(median(as.numeric(allAssigned[ii,colName])));
      }
    }
  }
  
  if(printErr){
    print(paste("Ambiguous SP specific ploidies found for ",length(toD)," segment-SP pairs."));
    print("Ploidies not assigned for these segments in corresponding SPs.")
  } 
  
  print("... Done.")
  if (ambig){
    print("Warning: parameter <ambig> set to TRUE. Output includes segment-assignements where subpopulation specific ploidy is ambiguous.Recommend repeating circular binary segmentation with less stringent parameters instead, to reduce segment length and thus the prevalence of ambiguous assignements.")
  }
  
  return(ploidy)
}

