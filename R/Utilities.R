
.onUnload <- function (libpath) {
  library.dynam.unload("robin", libpath)
}