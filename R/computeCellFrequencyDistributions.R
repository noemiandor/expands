computeCellFrequencyDistributions<-function(dm, max_PM=6, p, min_CF=0.1, ploidy = 2, nc = 1, v = T){

  STDOUT = "cCFD.log"
  
  print("Computing cell-frequency probability distributions...")
  ##add cell-frequency that best explains CN and AF to dm
  dm=.addColumn(dm,"f",NA);
  ##compute cell-freq-probability for each mutation
  freq=t(seq(min_CF,1,by=p/10));
  
  if(nc==1){
    output=.oneThreadCellFrequencyDistributions(list(dm=dm, max_PM=max_PM, precision=p, freq=freq, min_CellFreq=min_CF, ploidy=ploidy, verbose=v))
  }else{
    ###########################
    ####Parallel processing####
    print(paste("Using ",nc," cores to calculate cell-frequency distributions..."))
    print(paste("stdout and stderr connections will be redirected to ",STDOUT))
    cl <- makeCluster(nc,outfile=STDOUT)
    ## Split mutation profile
    input=list()
    for(chr in unique(dm[,"chr"])){
      dmx=dm[dm[,"chr"]==chr,,drop=F]
      input[[as.character(chr)]]=list(dm=dmx, max_PM=max_PM, precision=p, freq=freq, min_CellFreq=min_CF, ploidy=ploidy, verbose=v);
    }
    ## Distribute jobs
    results=clusterApply(cl,input,.oneThreadCellFrequencyDistributions)
    # Gather results
    output=list(densities=c(),freq=freq,dm=c());
    for(i in 1:length(results)){
      output$densities=rbind(output$densities,results[[i]]$densities);
      output$dm=rbind(output$dm,results[[i]]$dm);
    }
    
    stopCluster(cl)
  }
  return(output)
}



.oneThreadCellFrequencyDistributions<-function(varargin){
  # library(matlab)
  # library(expands)
  dm=varargin$dm
  max_PM=varargin$max_PM
  precision=varargin$precision
  freq=varargin$freq
  min_CellFreq=varargin$min_CellFreq
  norm=T; #varargin$norm
  ploidy=varargin$ploidy
  verbose=varargin$verbose
  
  densities=matrix(matrix(NA,nrow(dm),length(freq)),nrow=nrow(dm),ncol=length(freq),dimnames=list(1:nrow(dm),freq));
  success=0; 
  errors=c();
  warnings=c();
  for (k in 1:nrow(dm)){
    output=try(cellfrequency_pdf(af=dm[k,"AF_Tumor"],cnv=dm[k,"CN_Estimate"], pnb=dm[k,"PN_B"],freq=freq, max_PM=max_PM,ploidy=ploidy),silent=TRUE);
    if(class(output)=="try-error"){
      errors=rbind(errors,output);
    }else{
      if(norm){
        output$p=output$p/sum(output$p,na.rm=T); ##under the assumption that relative rather than absolute probabilities matter for clustering
      }
      densities[k,]=output$p;
      dm[k,"f"]=output$bestF;
      success=success+1;
    } 
    
    if (mod(k,20)==0){
      .notifyUser(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ",
                  success,"/",k),verbose=verbose)
    }
  }
  failure=nrow(dm)-success;
  if (length(errors)>0){
    print(paste("Failed to find cell-frequency distribution for ",
                failure,"SNVs."));
    .notifyUser("Causes:",verbose=verbose);
    errors=unique(errors);
    for (i in 1:length(errors)){
      .notifyUser(errors[i],verbose=verbose);
    }
  }
  if (length(warnings)>0){
    warnings=unique(warnings);
    .notifyUser(paste("Additional warnings:"))
    for (i in 1:length(warnings)){
      .notifyUser(warnings[i],verbose=verbose);
    }
  }
  .notifyUser("...Done.",verbose=verbose)
  output=list(densities=densities,freq=freq,dm=dm);
  return(output);
}
