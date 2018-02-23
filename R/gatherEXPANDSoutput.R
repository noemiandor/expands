gatherEXPANDSoutput<-function(outDirEXPANDS,regex=""){
  f=grep(regex,list.files(outDirEXPANDS,pattern="*.sps$",full.names = T),value=T)
  expdata=list()
  for(x in f){
    sps=read.table(x,sep="\t",header = T,check.names = F,stringsAsFactors = F);
    if(any(colnames(sps)=="LOCUS")){
      rownames(sps)=sps$LOCUS
      sps=sps[,setdiff(colnames(sps),"LOCUS")]
    }
    
    spstats=read.table(gsub(".sps$",".spstats",x),sep="\t",header = T,check.names = F,stringsAsFactors = F);
    # ##TODO --> remove this column update: --> fixed to be correct for next expands version
    # tmp=count(sps[,"SP"]); i=.intersect_MatlabV(round(tmp$x,3),round(spstats$`Mean Weighted`,3)); 
    # spstats[i$ib,"nMutations"]=tmp$freq[i$ia] 
    # ##
    
    if(file.exists(gsub(".sps$",".tree",x))){
      tree=try(phylobase::readNewick(gsub(".sps$",".tree",x)))
      treeApe=read.tree(gsub(".sps$",".tree",x))
    }else{
      tree=NULL; treeApe=NULL;
    }
    spscbs=try(read.table(gsub(".sps$",".sps.cbs",x),sep="\t",header = T,check.names = F,stringsAsFactors = F));
    if(class(spscbs)!= "try-error"){
      spscbs$seglength=1+spscbs$endpos-spscbs$startpos
      if(any(colnames(spscbs)=="LOCUS")){
        rownames(spscbs)=spscbs$LOCUS
        spscbs=spscbs[,setdiff(colnames(spscbs),"LOCUS")]
      }
    }
    sName=strsplit(fileparts(x)$name,".sps")[[1]][1]
    expdata[[sName]]=list(snv=sps,cbs=as.matrix(spscbs),spstats=as.matrix(spstats),tree=tree, treeApe=treeApe)
  }
  return(expdata)
}
