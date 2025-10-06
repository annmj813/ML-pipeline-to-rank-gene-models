pkgs <- c("readr","dplyr","rsample","yardstick","tibble","yaml")
need <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(need)) install.packages(need, repos="https://cloud.r-project.org")
