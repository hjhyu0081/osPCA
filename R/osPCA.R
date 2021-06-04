# CheckFile by (PSO)---------------------------------------------------------------
CheckFile <-
  function(filename,
           fileext,
           filepath,
           cnt = 1,
           today = format(Sys.Date(), "%y%m%d")) {
    if (nchar(fileext) > 0) {
      fullname <- paste0(today, "_", filename, "_", cnt, ".", fileext)
      while (any(fullname %in% list.files(path = filepath))) {
        cnt <- cnt + 1
        fullname <-
          paste0(today, "_", filename, "_", cnt, ".", fileext)
      }
    } else{
      # file should be a folder
      fullname <- paste0(today, "_", filename, "_", cnt)
      while (any(fullname %in% list.files(path = filepath))) {
        cnt <- cnt + 1
        fullname <- paste0(today, "_", filename, "_", cnt)
      }
    }
    return(paste0(filepath, fullname))
  }


# A matrix making ---------------------------------------------------------
making_amat<-function(ind_grp){
  if(is.factor(ind_grp)!=TRUE) ind_grp<-as.factor(ind_grp)
  m<-length(levels(ind_grp))
  A<-matrix(0,nrow=m,ncol=length(ind_grp))
  for(i in 1:m){
    A[i,ind_grp==levels(ind_grp)[i]]<-1/table(ind_grp)[i]
    if(i>1) A[i-1,ind_grp==levels(ind_grp)[i]]<-(-1/table(ind_grp)[i])
  }
  return(A=A[-m,])
}

con_pca<-function(data, group = NULL, A = NULL, q=2, max.iter=100, tol = 10^{-5}){
  if(missing(A)) A <- making_amat(group)
  data <- scale(data);
  n<-dim(data)[1]; p<-dim(data)[2];
  res<-prcomp(data,scale=T,center=T)
  B<-t(res$rotation[,1:q]); X<-res$x[,1:q];
  iter<-1
  while(iter<max.iter){
    B.old<-B; X.old<-X;
    C<-B%*%t(data)
    X<-matrix(solve.QP(Dmat = diag(n*q), dvec=c(t(C)), Amat = -t(kronecker(diag(q),A)), bvec = rep(0,q*dim(A)[1]))$solution,
              nrow=n,ncol=q)
    res.svd<-svd(t(data)%*%X)
    B<-res.svd$v%*%t(res.svd$u)
    if(sum((B.old-B)^2)<tol) break
    iter <- iter+1
  }
  X<-data.frame(X)
  colnames(X)<-paste0("PC",1:q)
  return(list(X = X, B = B))
}

org_pca_flipped <- function(data, group, q){
  res <- prcomp(data, center = T, scale = T)
  cand <- paste0("PC",1:q)
  data_meaned <- data.frame(res$x[,1:q],group = label) %>% group_by(group) %>% summarise_at(vars(cand), funs(mean = "mean"))
  B <- t(res$rotation[,1:q]); X <- res$x[,1:q];
  X[,data_meaned[1,-1] > data_meaned[3,-1]] <- -X[,data_meaned[1,-1] > data_meaned[3,-1]]
  B[data_meaned[1,-1] > data_meaned[3,-1],] <- -B[data_meaned[1,-1] > data_meaned[3,-1],]
  return(list(X = X, B = B))
}

con_pca_flipped<-function(data, group, q=2, max.iter=100, tol = 10^{-5}){
  A <- making_amat(group)
  data <- scale(data);
  n<-dim(data)[1]; p<-dim(data)[2];
  res<-prcomp(data,scale=T,center=T)
  cand <- paste0("PC",1:q)
  data_meaned <- data.frame(res$x[,1:q],group) %>% group_by(group) %>% summarise_at(vars(cand),funs(mean="mean"))
  B<-t(res$rotation[,1:q]); X<-res$x[,1:q];
  X[,data_meaned[1,-1] > data_meaned[3,-1]] <- -X[,data_meaned[1,-1] > data_meaned[3,-1]]
  B[data_meaned[1,-1] > data_meaned[3,-1],] <- -B[data_meaned[1,-1] > data_meaned[3,-1],]

  iter<-1
  while(iter<max.iter){
    B.old<-B; X.old<-X;
    C<-B%*%t(data)
    X<-matrix(solve.QP(Dmat = diag(n*q), dvec=c(t(C)), Amat = -t(kronecker(diag(q),A)), bvec = rep(0,q*dim(A)[1]))$solution,
              nrow=n,ncol=q)
    res.svd<-svd(t(data)%*%X)
    B<-res.svd$v%*%t(res.svd$u)
    if(sum((B.old-B)^2)<tol) break
    iter <- iter+1
  }
  X<-data.frame(X)
  colnames(X)<-paste0("PC",1:q)
  return(list(X = X, B = B))
}

scatterplot_mean<-function(score, group, PC){
  data <- data.frame(score,group=as.factor(group))
  cand <- paste0("PC",PC)
  data_meaned<-data %>%
    group_by(group) %>%
    summarise_at(vars(cand), funs(mean = "mean"))
  return(ggplot(data = data, aes(x = get(cand[1]), y = get(cand[2]), col = group)) +
           geom_point(size = .5) +
           geom_point(data = data_meaned, aes (x = get(paste0(cand[1],"_mean")), y = get(paste0(cand[2],"_mean")), col = group) , size = 2) +
           geom_line(data = data_meaned, aes (x = get(paste0(cand[1],"_mean")), y = get(paste0(cand[2],"_mean"))), col = "Black") +
           labs(x = paste0(cand[1]), y = paste0(cand[2])))
}

scoreplot<-function(data, group, q, flipped = FALSE){
  if(flipped){
    res_conPCA <- con_pca_flipped(data = data, A = making_amat(group), q = q)
    res_orgPCA <- org_pca_flipped(data = data, group = label, q = q)
  }else{
    res_conPCA <- con_pca(data = data, A = making_amat(group), q = q)
    res_orgPCA <- prcomp(data, scale = T, center = T)
    res_orgPCA <- list(X = res_orgPCA$x[,1:q], B = t(res_orgPCA$rotation[,1:q]))
  }

  data_plot <- rbind(data.frame(res_conPCA$X, rotation = "os-PCA", obs = rownames(res_conPCA$X)),
                     data.frame(rbind(res_orgPCA$X, res_orgPCA$X %*% t(res_conPCA$B%*%t(res_orgPCA$B))),
                                rotation = c(rep("original", nrow(res_orgPCA$X)), rep("rotated", nrow(res_orgPCA$X))),
                                obs = rep(rownames(res_orgPCA$X),2)))

  g <- list()
  for(i in 1:(q-1)){
    for(j in (i+1):q){
      g <- append(g, list(ggplot(data = data_plot,
                                 mapping = aes(x = get(paste0("PC",i)), y = get(paste0("PC",j)),
                                               size = rotation, alpha = rotation, colour = rotation)) +
                            geom_point() +
                            scale_size_manual(values = c(2,0.5,0.5)) +
                            scale_alpha_manual(values = c(0.5,1,1)) +
                            scale_color_manual(values = c("#999999","#999999","red")) +
                            labs(title = paste0("Scatter plot of PC scores PC",i,"vs. PC",j," with q = ",q,"."),
                                 x = paste0("PC",i), y = paste0("PC",j))))
    }
  }
  return(g)
}

pairs <- function(graph){
  tmp <- ggplot_gtable(ggplot_build(graph[[1]] + theme(legend.position = "right")))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]

  for(i in 1:5){
    assign(paste0("text",i), ggplot(data = data.frame(text = paste0("PC",i), size = 1000), aes (x = text, y = 1, label = text)) + geom_fit_text(grow = T) + theme_bw() +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
    )
  }

  return(grid.arrange(text1,graph[[1]],graph[[2]],graph[[3]],graph[[4]],
                      text2,graph[[5]],graph[[6]],graph[[7]],
                      text3,graph[[8]],graph[[9]],
                      text4,graph[[10]],
                      text5,
                      legend,
                      layout_matrix = rbind(c(1,NA,NA,NA,16),
                                            c(2,6,NA,NA,NA),
                                            c(3,7,10,NA,NA),
                                            c(4,8,11,13,NA),
                                            c(5,9,12,14,15))))
}
