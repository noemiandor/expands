assignQuantityToSP<-function(cbs, dm, C=list(sp=c("SP","SP_cnv"), pm=c("PM","PM_cnv")), e=1, v=T){
  print("Assigning copy number to SPs...")
  #dm[, C$sp]=round(1000*dm[, C$sp])/1000
  sps = sort(unique(as.numeric(as.matrix(dm[, C$sp]))))
  cols = c( colnames(cbs), paste("SP_", as.character(unique(sps)), sep = "") )
  sp_cbs = matrix(cbind(cbs, matrix(NaN, nrow(cbs), length(sps))), 
                  nrow = nrow(cbs), ncol = length(cols), dimnames = list(1:nrow(cbs), cols))
  
  toD = c();  
  for (k in 1:nrow(sp_cbs)) {
    if (mod(k, 100) == 0) {
      .notifyUser(paste("Finding overlaps for CBS segment", k, "out of ", nrow(sp_cbs), "..."),verbose=v)
    }
    idx = which(dm[, "chr"] == sp_cbs[k, "chr"] & dm[, "startpos"] >= 
                  sp_cbs[k, "startpos"] & dm[, "startpos"] <= sp_cbs[k, "endpos"])
    if (length(idx) == 0) {
      next
    }

    for(sp in sps){
      spName= paste("SP_", sp, sep = "")
      #       iConstant=which(apply(dm[, C$sp,drop=F],1,var)==0); ##SP==SP_cnv
      #       iGoodFit=which(dm[, "%maxP"]>=1)
      idy=which(apply(dm[, C$sp,drop=F]==sp,1,any))
      idy=intersect(idx,idy); ##index of SNVs/CNVs assigned to this SP and located within segment k
      # idy=intersect(iGoodFit,idy); 
      if (length(idy) == 0) {
        next
      }
      N=sum(dm[idy,"%maxP"],na.rm=T); ##cummulative assignment confidence for these SNVs/CNVs
      
      boolI=matrix(F,length(idy),length(C$sp))
      boolI[dm[idy, C$sp,drop=F]==sp]=T
      boolI[apply(boolI,1,all),1]=F; ##take only one PM, where SP==SP_cnv
      pms=t(dm[idy,C$pm])[t(boolI)]; ##SP ploidies
      pm=sum(pms*dm[idy,"%maxP"],na.rm=T)/N; ##%maxP-weighted mean of PMs, PM_cnvs
      
      thisAssigned=c(k,sp ,sp_cbs[k,"chr"],sp_cbs[k,"startpos"],sp_cbs[k,"endpos"],pm);
      if(length(pms)>1 && var(pms,na.rm=T)>e){
        toD = rbind(toD, thisAssigned)
      }else{
        sp_cbs[k, spName]=pm
      }
    }
    
  }
  
  if(!is.null(toD) && nrow(toD)>0){
    print(paste("Ambiguous SP specific copy numbers found for ",nrow(toD)," segment-SP pairs."));
    print("Copy number not assigned for these segments in corresponding SPs.")
  } 
  
  .notifyUser("... Done.",verbose=v);
  if(e>0){
    .notifyUser("Warning: parameter <e> set to be >0. Output includes segment-assignements where subpopulation specific copy number is ambiguous. Recommend repeating circular binary segmentation with less stringent parameters instead, to reduce segment length and thus the prevalence of ambiguous assignements.",verbose=v)
  }
  
  return(sp_cbs)
}

