assignMutations<-function( dm, finalSPs, max_PM=6){
  
  addCols=c("%maxP","SP","SP_cnv","PM","PM_cnv","PM_B","PM_B_cnv");
  for (k in 1:length(addCols)){
    dm=.addColumn(dm,addCols[k],NA);
    dm[,addCols[k]]=NA; 
  }
  
  spFreq = sort(finalSPs[, "Mean Weighted",drop=F])
  precision=finalSPs[1,"precision"]
  
  success=0;
  for(k in 1:nrow(dm)){
    af=as.numeric(dm[k,"AF_Tumor"]); cnv=as.numeric(dm[k,"CN_Estimate"]); pnb=as.numeric(dm[k,"PN_B"])
    ###Get cell frequency probability solution for this locus only among previously detected SPs
    .jcall("core.utils.Common","V","setALLOWED_FREQUENCIES",.jarray(as.double(spFreq)))
    .jcall("core.utils.Common","V","setMAX_PM",as.integer(max_PM))
    o <-.jnew("core.ParallelSubpopulations",as.double(cnv), as.double(af), as.integer(pnb),as.logical(F));
    e=o$getBestFittingSubpopulation();
    bsp=e$getKey();
    
    ###Save results  
    dm[k,"%maxP"]=e$getValue();
    ##Important: `SP` is always the one carrying an SNV, may or may not carry an CNV
    dm[k,"SP"]=bsp$getCellularFrequency();
    dm[k,"PM"]=bsp$getPm();
    dm[k,"PM_B"]=bsp$getPmb();
    ##Important: `SP_cnv` is always the one carrying an CNV, if it exists. `SP_cnv` may or may not carry an SNV
    if(!is.null(bsp$getChild())){
      dm[k,"SP_cnv"]=bsp$getChild()$getCellularFrequency();
      dm[k,"PM_cnv"]=bsp$getChild()$getPm();
      dm[k,"PM_B_cnv"]=bsp$getChild()$getPmb();
    }else if(!is.null(bsp$getSibling())){
      dm[k,"SP_cnv"]=bsp$getSibling()$getCellularFrequency();
      dm[k,"PM_cnv"]=bsp$getSibling()$getPm();
      dm[k,"PM_B_cnv"]=bsp$getSibling()$getPmb();
    }else {
      dm[k,"SP_cnv"]=bsp$getParent()$getCellularFrequency();
      dm[k,"PM_cnv"]=bsp$getParent()$getPm();
      dm[k,"PM_B_cnv"]=bsp$getParent()$getPmb();
    }
    
    ###User feedback
    success=success+1;
    if (mod(k,20)==0){
      print(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ",
                  success,"/",k))
    }
  }
  
  ##Remove SPs to which no mutations were assigned
  iK=c();
  for (j in 1:nrow(finalSPs)){
    idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
    finalSPs[j,"nMutations"]=length(idx);
    if(length(idx)>0){
      iK=c(iK,j);
    }
  }
  finalSPs=finalSPs[iK,, drop=FALSE];
  
  output=list("dm"=dm,"finalSPs"=finalSPs);
  return(output);
}

