cellfrequency_pdf <-function(af,cnv,pnb,freq, max_PM=6){
  
  ###Get cell frequency probability solution for this locus
  .jcall("core.utils.Common","V","setALLOWED_FREQUENCIES",.jarray(as.double(freq)))
  .jcall("core.utils.Common","V","setMAX_PM",as.integer(max_PM))
  o <-.jnew("core.ParallelSubpopulations",as.double(cnv), as.double(af), as.integer(pnb),as.logical(T));
  cs=.jcall(o,"Ljava/util/Map;","getCellFreq2ProbabilityMap")
  
  ###Parse java object, sort and set small non-zero distance for perfect fit
  mat=.map2mat(cs)
  mat=mat[sort(mat[,"f"],index.return=T)$ix,,drop=F]
  perfectI=which(mat[,"prob"]>=.Machine$double.xmax/1E9); inperfectI=setdiff(1:nrow(mat),perfectI)
  mat[perfectI,"prob"]=max(mat[inperfectI,"prob"])
  mat[,"prob"]=mat[,"prob"]/min(mat[inperfectI,"prob"]); ##Kernel density estimation cannot handle only very large values
  
  ###Kernel density estimation
  p=density(mat[,"f"], bw = "SJ", adjust = 0.25,
            kernel = c("gaussian"),
            weights = mat[,"prob"],from=min(freq),to=max(freq))
  bestF=p$x[which.max(p$y)];
  
  output=list("p"=approx(p$x,p$y,freq)$y,"bestF"=bestF);
  return(output)
}


.map2mat<-function(cs){
  keySet<-.jrcall(cs,"keySet")
  an_iter<-.jrcall(keySet,"iterator")
  aList <- matrix(NA,.jrcall(keySet,"size"),2); colnames(aList)=c("f","prob")
  i=1
  while(.jrcall(an_iter,"hasNext")){
    e <- .jrcall(an_iter,"next");
    aList[i,]=c(e, .jrcall(cs,"get",.jnew("java.lang.Double",e)))
    i=i+1
  }
  return(aList)
}

