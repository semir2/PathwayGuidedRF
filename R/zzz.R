
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste("Since version 0.4.0 P values are by default adjusted for",
                                "multiple testing using the Benjamini-Hochberg",
                                "procedure!"))
}
