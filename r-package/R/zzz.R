# Package init.

.onUnload <- function(libpath) {
  library.dynam.unload("misskappa", libpath)
}
