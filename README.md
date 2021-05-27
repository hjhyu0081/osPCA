# os-PCA
This is a README file of the R package osPCA. In our paper, we devise the os-PCA procedure.

## Installation of package
To install our package, you may simply execute the following code.
```
# install.packages("remotes")
remotes::install_github("hjhyu0081/osPCA@main")
```

## osPCA package
For an example to perform os-PCA, we analyze rice data which is the real data example for our paper. If you want to download this rice data, put these code lines into your R console. Note that the first column is the name of each variable.
```
# install.packages("readr")
library (readr)
urlfile="https://raw.githubusercontent.com/hjhyu0081/osPCA/main/data/rice.csv"
data <- t(read.csv(url(urlfile), encoding = "UTF-8", as.is=T))
```

Main functions for this package are `con_pca`, `scatterplot_mean`, and `scoreplot`. For more detail, open osPCA.md in this directory.

