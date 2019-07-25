
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste("P values are now by default adjusted for",
                                "multiple testing using the Benjamini-Hochberg",
                                "procedure!"))
}
