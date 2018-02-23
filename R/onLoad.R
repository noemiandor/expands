.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
 
  ##Make sure Java version requirement is met at runtime
  .jinit()
  jv <- .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
  if(substr(jv, 1L, 1L) == "1") {
    jvn <- as.numeric(paste0(strsplit(jv, "[.]")[[1L]][1:2], collapse = "."))
    if(jvn < 1.5) stop("Java 5.0 is needed for this package but not available")
  }
  
}
