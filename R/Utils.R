.addColumn<-function(M,newCol,initVal){
  if (!any(colnames(M)==newCol)){
    if(!is.null(dim(M))){
      M=matrix(cbind(M,matrix(initVal,nrow(M),1)),nrow=nrow(M),ncol=ncol(M)+1,
               dimnames = list(rownames(M), c(colnames(M),newCol)));
    }else{
      cols=names(M);
      M=c(M,initVal);	
      names(M)=c(cols,newCol);
    }
  }
  return(M);
}





.intersect_MatlabV<-function(a,b){
  x=intersect(a,b)
  ia=match(x,a);
  ib=match(x,b);
  return(list(a=x,ia=ia,ib=ib));
}





.writeExpandsOutput<-function(X, dirF,snvF,suffix=".sps.cbs", message="Output"){
  X=as.data.frame(X)
  if(suffix==".sps.cbs"){
    X[,"LOCUS"]=paste(X$chr,":",X$startpos,"-",X$endpos,sep="")
  }else if(suffix==".sps"){
    X[,"LOCUS"]=paste(X$chr,":",X$startpos,"-",X$startpos,sep="")
  }
  output=paste(dirF, .Platform$file.sep, snvF,suffix,sep="");
  write(paste("## expands version",packageVersion("expands")), file = output, append=FALSE);
  suppressWarnings(write.table(X,file = output, append=TRUE, quote = FALSE, sep = "\t", row.names=FALSE));
  print(paste(message,"saved under",output));
}



.notifyUser<-function(message,verbose=T){
  if(verbose){
    print(message)
  }
}


.readSNVandCBS<-function(SNV,CBS,max_PM=6, min_CF=0.1,snvF=NULL, verbose){
  ##SNVs
  if (is.character(SNV) && file.exists(SNV)){
    print(paste("Running ExPANdS on: ",SNV))
    dm=read.table(SNV,sep="\t",header=TRUE,check.names = F,stringsAsFactors = F);
    if (!all(sapply(dm,is.numeric))){
      print(paste("Warning: not all columns in", SNV,"are numeric. Use only numeric data as input to ensure unexpected conversion does not occur."))
    }
    dm <- dm[ as.character(dm[,"chr"]) %in% as.character(seq(100)), ];
    dm=data.matrix(dm);
    print("Only SNVs with autosomal coordinates included.")
    if(is.null(snvF)){
      snvF=SNV;
    }
  }else if (is.matrix(SNV)){
    dm=SNV;
    if (!is.numeric(dm)) {
      print("SNV matrix has to be numeric. Likely cause: only mutations detected on autosomes accepted for ExPANdS model. Remove SNVs with allosomal and mitochondrial coordinates before you proceed.")
      return();
    }
    if (is.null(snvF)){
      snvF="out.expands";
    }
  }else{
    print("No SNVs provided. Aborting ExPANdS.");
    return();
  }
  
  ##Output
  dirF=fileparts(snvF)$pathstr;
  if (nchar(dirF)==0){
    dirF=".";
  }
  snvF=paste(fileparts(snvF)$name,fileparts(snvF)$ext,sep="");
  
  
  ##CBS
  if (is.character(CBS) && file.exists(CBS)){
    copyNumber=as.matrix(read.table(CBS,sep="\t",header=TRUE,check.names = F,stringsAsFactors = F))
    if (!all(sapply(copyNumber,is.numeric))){
      print(paste("Warning: not all columns in", CBS,"are numeric. Use only numeric data as input to ensure unexpected conversion does not occur."))
    }
  }else if (is.matrix(CBS)){
    copyNumber=CBS;
  }else{
    print("No copy number information provided. Aborting ExPANdS.")
    return();
  }
  if(!any("CN_Estimate" %in% colnames(dm))){  
    if (any(copyNumber[,"CN_Estimate"]<0) || quantile(copyNumber[,"CN_Estimate"],0.9)<1){
      print("Column <CN_Estimate> in CBS input seems to contain log-ratio entries. Absolute copy number values are required (e.g. average ~2.0 expected for predominantly diploid genomes).")
      return();
    }
    if (sum(abs(copyNumber[,"CN_Estimate"]-round(copyNumber[,"CN_Estimate"])))<1){
      print("Warning! Copy numbers have values rounded to the closest integer. Using rational positive estimates of copy numbers is recommended.")
    }
    dm=assignQuantityToMutation(dm,copyNumber,"CN_Estimate",verbose=verbose);
  }else{
    print("Using column <CN_Estimate> from <SNV> as copy number estimate. Parameter <CBS> used only for phylogeny inference.")
  }
  
  ii=which(is.na(dm[,"CN_Estimate"]));
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to unavailable copy number in that region."));
    dm=dm[-ii,];
  }
  ii=which(dm[,"CN_Estimate"]>max_PM);
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to high-level amplifications (>",max_PM, "copies) within that region. Consider increasing value of parameter max_PM to include these SNVs, provided high coverage data (> 150 fold) is available"));
    dm=dm[-ii,];
  }
  ii=which(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]<min_CF);
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to AF*CN below ", min_CF," (SNV can't be explained by an SP present in ",min_CF*100 ,"% or more of the sample)."));
    dm=dm[-ii,];
  }
  
  return(list(dm=dm, copyNumber=copyNumber, snvF=snvF))
}
