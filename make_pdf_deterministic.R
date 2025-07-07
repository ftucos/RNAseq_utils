make_pdf_deterministic <- function(file) {
    # 1. strip metadata and info (including ModDate)
    system2("exiftool", c(
        "-q", "-q", # suppress normal informational messages and [minor] warnings
        "-all:all=",
        "-overwrite_original",
        file
    ))

    # 2. rewrite deterministically
    # qpdf can remove also metadata but not the ModDate from the /Info dictionary
    system2("qpdf", c(
        "--linearize", # Reorders internal PDF objects to optimize for consistent structure and fast web view
        "--deterministic-id", # uses only internal content for ID generation (no timestamp and output file name)
        file,
        "--replace-input"
    ))
}
