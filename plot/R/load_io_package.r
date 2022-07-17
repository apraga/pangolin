load_io_package <- function(io_type) {
  if (io_type == "hdf5") {
    if (!require("rhdf5",character.only = TRUE)) {
      stop("Package rhdf5 not found")
    }
  }
  else { }
}
