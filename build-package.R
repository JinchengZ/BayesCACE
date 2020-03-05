rm(list = ls())
library(devtools)
library(roxygen2)

## create vignette
# devtools::use_vignette("my-vignette")
# devtools::use_data(CODES, internal = T)
## Create testthat
## only run once the following
# devtools::use_testthat()
## run test cases of functions
# options(testthat.output_file = "test-out.xml")
devtools::load_all()
# devtools::test()
# testthat::test_dir("tests/")
# testthat::test_dir("tests/", reporter = "junit")
# check
devtools::document()
devtools::check()

# build
Sys.setenv("TAR" = "internal")
Sys.getenv("PATH")
devtools::build(manual = T)
devtools::install_local("../BayesCACE_1.0.tar.gz", dependencies = NA, upgrade = "never")
# devtools::install_github("JinchengZ/BayesCACE")

## generate the help manual.
pack <- "BayesCACE"
path <- find.package(pack)
if (file.exists(paste0(pack, ".pdf"))) {file.remove(paste0(pack, ".pdf"))}

# system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))
system(paste(file.path(R.home(""), "R"), "CMD", "Rd2pdf", shQuote(path)))


