assignMutations<-function( dm, finalSPs, max_PM=6, cnvSPs=NULL, ploidy=2,verbose=T){
  
  addCols=c("%maxP","SP","SP_cnv","PM","PM_cnv","PM_B","PM_B_cnv","f");
  for (k in 1:length(addCols)){
    dm=.addColumn(dm,addCols[k],NA);
    dm[,addCols[k]]=NA; 
  }
  dm=as.data.frame(dm);
  dm$SNV_CNV_Phylo=NA
  
  spFreq = .jarray(as.double( sort(finalSPs[, "Mean Weighted",drop=F]) ))
  spcnvFreq=as.double(NULL); 
  if(!is.null(cnvSPs)){ ##Frequency spectrum for SP_cnv different from that of SP
    spcnvFreq = .jarray(as.double( sort(cnvSPs[, "Mean Weighted",drop=F]) ))
  }
  precision=finalSPs[1,"precision"]
  
  success=0;
  ##Initial mutation assignment
  for(k in 1:nrow(dm)){
    .init_rJava(spFreq,spcnvFreq,max_PM,pnb=as.numeric(dm[k,"PN_B"]),spRel=NULL)
    dm[k,]=.call_rJava(dm[k,,drop=F], ploidy=ploidy)
    
    ###User feedback
    success=success+1;
    if (mod(k,20)==0){
      .notifyUser(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ",
                  success,"/",k),verbose=verbose)
    }
  }
  
  
  ##@TODO: test
  ##Count sp pair relations among SNVs with co-occurence violation assumption
  ii=which(dm$SP!=dm$SP_cnv);   dmx=dm[ii,,drop=F]
  if(!isempty(ii)){
    print(paste("Resolving potential phylogeny conflicts among",length(ii),"loci..."))
    dmx[,c("SP","SP_cnv")]=t(apply(dmx[,c("SP","SP_cnv"),drop=F],1,sort))
    uP=unique(dmx[,c("SP","SP_cnv")]);     uPN=paste(uP[,1],uP[,2],sep="_")
    uR=c("child_parent","parent_child","siblings")
    spRel=matrix(0,length(uPN),length(uR)); rownames(spRel)=uPN; colnames(spRel)=uR
    for(i in 1:nrow(spRel)){
      eq=dmx[,c("SP","SP_cnv"),drop=F]==repmat(unlist(uP[i,]),nrow(dmx),1);
      iK=which(apply(eq ,1,all))

      fr=count(dmx$SNV_CNV_Phylo[iK])
      spRel[uPN[i],as.character(fr$x)]=fr$freq
    }
    spRel=as.data.frame(spRel); spRel$parental=apply(spRel[,c("child_parent","parent_child"),drop=F],1,sum,na.rm=T)
    .notifyUser(spRel,verbose=verbose)
    spRel=cbind(uP,spRel)

    ##Reassign mutations while accounting for paired phylogenies
    for(k in ii){
      .init_rJava(spFreq,spcnvFreq,max_PM,pnb=as.numeric(dm[k,"PN_B"]),spRel=spRel)
      dm[k,]=.call_rJava(dm[k,,drop=F], ploidy=ploidy)
    }

  }
  
  
  ##Remove SPs to which no mutations were assigned and record SNV/CNV cooccurence
  finalSPs=.addColumn(finalSPs,'snv_cnv_Co',NA);
  finalSPs=.addColumn(finalSPs,'nCNVs',NA);
  iK=c();
  for (j in 1:nrow(finalSPs)){
    idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
    finalSPs[j,"nMutations"]=length(idx);
    idx=which(dm[,"SP_cnv"]==finalSPs[j,"Mean Weighted"]);
    finalSPs[j,"nCNVs"]=length(idx);
    
    idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"] & dm[,"SP_cnv"]==finalSPs[j,"Mean Weighted"]);
    finalSPs[j,"snv_cnv_Co"]=length(idx);
  }
  finalSPs=finalSPs[finalSPs[,"nCNVs"]>0 | finalSPs[,"nMutations"]>0,, drop=FALSE];
  
  output=list("dm"=dm,"finalSPs"=finalSPs);
  return(output);
}


.init_rJava<-function(spFreq,spcnvFreq,max_PM,pnb,spRel=NULL){
  ###Get cell frequency probability solution for this locus only among previously detected SPs
  .jcall("core.utils.Common","V","setALLOWED_SP_FREQUENCIES",spFreq)
  .jcall("core.utils.Common","V","setALLOWED_SP_CNV_FREQUENCIES",spcnvFreq)
  .jcall("core.utils.Common","V","setMAX_PM",as.integer(max_PM))
  # if(pnb==1){
  #   .jcall("core.utils.Common","V","setMAX_PM",as.integer(2)); ##Limit scenarios when dealing with LOH
  # }
  if(!is.null(spRel)){
    ##Select relationships to keep:
    for(i in 1:nrow(spRel)){
      sp1=spRel$SP[i];    sp2=spRel$SP_cnv[i]
      if(spRel[i,"siblings"]>spRel[i,"parental"]){
        .jcall("core.utils.Common","V","setAllowedSiblings_SP_FREQUENCIES",as.double(sp1),as.double(sp2))
        .jcall("core.utils.Common","V","setAllowedSiblings_SP_FREQUENCIES",as.double(sp2),as.double(sp1))
      }else if(spRel[i,"siblings"]<spRel[i,"parental"]){
        .jcall("core.utils.Common","V","setAllowedSiblings_SP_FREQUENCIES"); ##This will allow only parental relations among SPs;
      }
    }
  }
}


.call_rJava<-function(dmx,ploidy){
  af=as.numeric(dmx[1,"AF_Tumor"]); cnv=as.numeric(dmx[1,"CN_Estimate"]); pnb=as.numeric(dmx[1,"PN_B"])
  o <-.jnew("core.ParallelSubpopulations",as.double(cnv), as.double(af), as.integer(pnb),as.logical(F), as.integer(ploidy));
  e=o$getBestFittingSubpopulation();
  bsp=e$getKey();
  
  ###Save results  
  dmx[1,"%maxP"]=e$getValue();
  ##Important: `SP` is always the one carrying an SNV, may or may not carry an CNV
  dmx[1,"SP"]=bsp$getCellularFrequency();
  dmx[1,"PM"]=bsp$getPm();
  dmx[1,"PM_B"]=bsp$getPmb();
  dmx[1,"f"]=dmx[1,"SP"]-bsp$getDeviation();
  ##Important: `SP_cnv` is always the one carrying an CNV, if it exists. `SP_cnv` may or may not carry an SNV
  if(!is.null(bsp$getChild())){
    dmx[1,"SP_cnv"]=bsp$getChild()$getCellularFrequency();
    dmx[1,"PM_cnv"]=bsp$getChild()$getPm();
    dmx[1,"PM_B_cnv"]=bsp$getChild()$getPmb();
    dmx[1,"SNV_CNV_Phylo"]="parent_child"
  }else if(!is.null(bsp$getSibling())){
    dmx[1,"SP_cnv"]=bsp$getSibling()$getCellularFrequency();
    dmx[1,"PM_cnv"]=bsp$getSibling()$getPm();
    dmx[1,"PM_B_cnv"]=bsp$getSibling()$getPmb();
    dmx[1,"SNV_CNV_Phylo"]="siblings"
  }else {
    dmx[1,"SP_cnv"]=bsp$getParent()$getCellularFrequency();
    dmx[1,"PM_cnv"]=bsp$getParent()$getPm();
    dmx[1,"PM_B_cnv"]=bsp$getParent()$getPmb();
    dmx[1,"SNV_CNV_Phylo"]="child_parent"
  }
  return(dmx)
}

