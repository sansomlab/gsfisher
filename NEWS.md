# gsfisher 0.2

* Added a `NEWS.md` file to track changes to the package.
* Pass `R CMD check`.
* Deprecated `readGMT` to use the established `qusage::read.gmt`.
* Added examples for usage and implicit testing of functions that do not
    rely on a working internet connection (e.g., `biomaRt`).
* Added small GMT file for testing.
* Renamed `translateGMT2mouse` to `mapENTREZhuman2mouse`.
* Renamed `fetchAnnotation` to `fetchAnnotations`.
* Split functions into multiple files by topic.
* Latest Ensembl `version` selected by `NULL`,
    consistently with `biomaRt::useEnsembl`.

# gsfisher 0.1

* Converted scripts into package.
