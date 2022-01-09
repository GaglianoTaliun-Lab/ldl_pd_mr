#R/4.0.5 installation for EasyStrata on ComputeCanada
install.packages("Cairo", repo="https://RForge.net")
install.packages(c("plotrix", "data.table"))
#selected mirror  15
url <- "https://homepages.uni-regensburg.de/~wit59712/easystrata/EasyStrata_18.1.tar.gz"
pkgFile <- "EasyStrata_18.1.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)
library('EasyStrata')
EasyStrata("my_ecf_file.ecf")
EasyStrats doc: https://homepages.uni-regensburg.de/~wit59712/easystrata/EasyStrata_8.6_Commands_140615.pdf

