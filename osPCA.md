---
title: "os-PCA"
author: "Hyeong Jin Hyun"
date: "5/25/2021"
output: 
  html_document: 
    keep_md: true
---
# os-PCA

## Introduction

This module is to introduce a basic example to use os-PCA procedure.  
To begin with, install our package by executing this simple code.

```r
remotes::install_github("hjhyu0081/osPCA@main")
```
You may require several other packages to perform our own package. Install those packages if needed.

```r
require(osPCA)
require(dplyr)
require(quadprog)
require(ggplot2)
require(ggfittext)
require(gridExtra)
```


Next, download the rice dataset from https://github.com/hjhyu0081/osPCA and set appropriate directory to analyze this datasets. The numerical data and the variables are saved and assigned separately just for convenience.


```r
data <- scale(t(read.csv("./data/rice.csv", encoding = "UTF-8", as.is=T)[,-1]))
variables <- read.csv("./data/rice.csv")[,1]
```

If you do not want to download it into your computer, you may directly import the data from my github.

```r
library (readr)
urlfile="https://raw.githubusercontent.com/hjhyu0081/osPCA/main/data/rice.csv"
data <- scale(t(read.csv(url(urlfile), encoding = "UTF-8", as.is=T)[,-1]))
variables <- read.csv(url(urlfile))[,1]
```


Also, it is essential to set the group variables according to the group. The group variables are recommended to set by numeric variables.


```r
group<-c(rep(1,30),rep(2,30),rep(3,30),rep(4,30),rep(5,30))
```


## os-PCA using rice data

### `scatterplot_mean` and `con_pca`

It is well known that `prcomp` can be used to perform PCA. We provide a function `scatterplot_mean` that allows to draw the scatterplot of principal components according to the group, in which each x and y axes can be selected by setting the third arguments. As following, the first argument of this function should be the score object of `prcomp` function and the second argument should be assigned group variables. You can see that the PC scores are not ordered in this case.


```r
res_orgPCA <- prcomp(data, center = T, scale = T)
scatterplot_mean(res_orgPCA$x, group = group, PC = c(1,2))
```

![](osPCA_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


The function `con_pca` offers the score and basis for os-PCA, while the `prcomp` provides them for PCA. The first argument of the function `con_pca` should be numerical datasets, and the second argument should be 'constrained matrix' $A$. The third argument should be the number of components which is chosen by screeplot. Moreover, you can select the number of maximum iteration and tolerance error for the basis.


```r
res_conPCA <- con_pca(data = data, A = making_amat(group), q = 7, max.iter=100, tol = 10^{-5})
scatterplot_mean(res_conPCA$X, group = group, PC = c(1,2))
```

![](osPCA_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

The followings are screeplot code lines. As stated in our paper in Section 3.2, the screeplot function for os-PCA is pretty much similar with that of PCA. Therefore, the screeplot of the os-PCA that allow you to choose the number of PCs could be substituted with that of PCA to avoid numerically burdensome work. Here is a simple example for this work. Besides you can make your own useful screeplot function.

### Screeplot

```r
ggplot(data = data.frame(lambda = res_orgPCA$sdev^2, i = 1:length(res_orgPCA$sdev))) +
  geom_line(aes(x = i, y = lambda)) +
  geom_point(aes(x = i, y = lambda), size = 0.3)
```

![](osPCA_files/figure-html/screeplot-1.png)<!-- -->

The following codes are for Figure 3 and Figure 4 using the function `scatterplot_mean`. Figure 3 is the group means of PC scores from the original PCA of the rice data and Figure 4 is that from the os PCA.

```r
for(method in c("orgPCA","conPCA")){
  for(i in 1:4){
    for(j in (i+1):5){
      if(method == "orgPCA"){
        tempplot <- scatterplot_mean(score = res_orgPCA$x, group = group, PC = c(i,j)) + theme_bw() +
          theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "none")
        assign(paste0("g",i,j,"_",method), tempplot)
      }else{
        tempplot <- scatterplot_mean(score = get(paste0("res_",method))$X, group = group, PC = c(i,j)) + theme_bw() +
          theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "none")
        assign(paste0("g",i,j,"_",method), tempplot)
      }
    }
  }
}
```


```r
for(method in c("orgPCA","conPCA")){
  assign(paste0("graphs_",method),list())
  for(i in 1:4){
    for(j in (i+1):5){
      assign(paste0("graphs_",method), append(get(paste0("graphs_",method)), list(get(paste0("g",i,j,"_",method)))))
    }
  }
}
```

To implement the following code lines, you might have an error message to implement the following code. To resolve this, you are recommended to execute the following code. To do this, you might need to find 'osPCA.R' file in my github and download it into your apprpriate directory.


```r
source("./R/osPCA.R")
```


```r
org <- pairs(graphs_orgPCA)
```

![](osPCA_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


```r
con <- pairs(graphs_conPCA)
```

![](osPCA_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

You can see that the score of the osPCA is ordered according to the group after implementing the osPCA procedure. 

### `scoreplot`
We provide a function `scoreplot` (Figure 5). Two pairs of PC scores are presented on the plots. The gray dots display scores that are obtained from os-PCA. The black dots and red dots are from the original PCA and the rotated PC, respectively. `scoreplot` function will provide the all pairs of PC scores. The following codes will just present the 1st vs 2nd and 1st vs 3rd PC scores. 


```r
scoreplotlist<-scoreplot(data = data, group = group, q = 7)
scoreplotlist[[1]]
```

![](osPCA_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
scoreplotlist[[2]]
```

![](osPCA_files/figure-html/unnamed-chunk-13-1.png)<!-- -->
