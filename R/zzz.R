.onLoad <- function(libname, pkgname) {
    msg <- paste0(
        "==================================================================\n",
        "               Welcome to the ", pkgname, " package!\n",
        "==================================================================\n\n",
        "This package is developed by the Computational Biology Group\n",
        "of ETH Zurich and the Swiss Institute of Bioinformatics (SIB).\n\n",
        "Please cite the following paper when using this package:\n",
        "https://www.nature.com/articles/s41467-023-39400-w\n\n",
        "For any questions, please contact niko.beerenwinkel@bsse.ethz.ch\n\n"
    )
    cat(msg)
}
