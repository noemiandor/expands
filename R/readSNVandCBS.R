.readSNVandCBS<-function(SNV,CBS,max_PM=6, min_CellFreq=0.1,snvF=NULL){
  ##SNVs
  if (is.character(SNV) && file.exists(SNV)){
    print(paste("Running ExPANdS on: ",SNV))
    dm=read.table(SNV,sep="\t",header=TRUE,stringsAsFactors = FALSE);
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
    copyNumber=as.matrix(read.table(CBS,sep="\t",header=TRUE,stringsAsFactors = FALSE))
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
      print("Column <CN_Estimate> in CBS input seems to contain log-ratio entries. Please supply absolute copy number values (e.g. average ~2.0 expected for predominantly diploid genomes).")
      return();
    }
    if (sum(abs(copyNumber[,"CN_Estimate"]-round(copyNumber[,"CN_Estimate"])))<1){
      print("Warning! Copy numbers have values rounded to the closest integer. Using rational positive estimates of copy numbers is recommended.")
    }
    dm=assignQuantityToMutation(dm,copyNumber,"CN_Estimate");
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
  ii=which(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]<min_CellFreq);
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to AF*CN below ", min_CellFreq," (SNV can't be explained by an SP present in ",min_CellFreq*100 ,"% or more of the sample)."));
    dm=dm[-ii,];
  }
  
  return(list(dm=dm, copyNumber=copyNumber, snvF=snvF))
}