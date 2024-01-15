library(reticulate)
use_python("~/path/anaconda3/envs/myenv/bin/python")
use_condaenv('myenv')

install.packages("tensorflow")

tensorflow::install_tensorflow()




install.packages("tensorflow")
tensorflow::install_tensorflow(extra_packages='tensorflow-probability', version = "2.9.1")
#ensure version 2 installed
tensorflow::tf_config()
library(tensorflow)
install_tensorflow() 
install_tensorflow(method='conda')


install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.0.3.tar.gz", repo=NULL, type="source")

install.packages("devtools")
library(rlang)
library(devtools)
devtools::install_github("Irrationone/cellassign")

