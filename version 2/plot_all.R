setwd("./R") # set wd R with Code and Data inside
if (!file.exists("figures")) {
  dir.create("figures")
}
source("Code/evt_abisko.R", local = TRUE)
source("Code/pandemics.R", local = TRUE)
source("Code/Venice_png.R", local = TRUE)
source("Code/heatwaves_bw.R", local = TRUE)
source("Code/SP500.R", local = TRUE)


