library(igraph)
#参数
#大小n,社区数K，实验重复次数rept
n <- 4000
K <- 2
rept <- 25
#概率矩阵P
pdata=c(3,0.5,0.5,1)
P<-matrix(pdata,K,K)
P
#度异质参量theta
theta <- 1:n
i1 <- 1:(n/4)
i2 <- (n/4+1):n
c0 <- 0.5
d0 <- 0.025*c(1,3,5,7,9)
#情形1
theta1[1:(n/4)] <- c0+(c0-d0[1])*(4*i1/n)
theta1[(n/4+1):n] <- c0+(c0-d0[1])*(4*i2/(3*n))^2
#情形2
theta2[1:(n/4)] <- c0+(c0-d0[2])*(4*i1/n)
theta2[(n/4+1):n] <- c0+(c0-d0[2])*(4*i2/(3*n))^2
#情形3
theta3[1:(n/4)] <- c0+(c0-d0[3])*(4*i1/n)
theta3[(n/4+1):n] <- c0+(c0-d0[3])*(4*i2/(3*n))^2
#情形4
theta4[1:(n/4)] <- c0+(c0-d0[4])*(4*i1/n)
theta4[(n/4+1):n] <- c0+(c0-d0[4])*(4*i2/(3*n))^2
#情形5
theta5[1:(n/4)] <- c0+(c0-d0[5])*(4*i1/n)
theta5[(n/4+1):n] <- c0+(c0-d0[5])*(4*i2/(3*n))^2

theta <- list(theta1,theta2,theta3,theta4,theta5)
#归一化
for (i in 1:5)
{
  # theta <- theta/sum(abs(theta))
  theta[[i]][1:n/4] <- theta[[i]][1:n/4]/sum(abs(theta[[i]][1:n/4]))
  theta[[i]][(n/4+1):n] <- theta[[i]][(n/4+1):n]/sum(abs(theta[[i]][(n/4+1):n]))
  theta[[i]] <- 0.8*theta[[i]]/max(theta[[i]])
  max(theta[[i]])
  min(theta[[i]])
}


#标签向量l  
l0 <- 1:n
l0[1:n/4] <- 1
l0[(n/4+1):n] <- 2


#生成主要信息矩阵omega
fomega <- function(t = 1)
{
  theta <- theta[[t]]
  omega<-matrix(0,n,n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {  
      omega[i,j] <- theta[i]*P[l0[i],l0[j]]*theta[j]
    }
  }
  return(omega)
}

#生成噪音矩阵W
fW <- function(omega)
{
  W<-matrix(0,n,n)
  for (i in 1:n)
  {
    for (j in 1:n)
    { 
      if (j>i){
        W[i,j] <- rbinom(1,1,omega[i,j])-omega[i,j]
      }
      
    }
  }
  W<-W+t(W)
  return(W)
}

#生成N（V,E)的邻接矩阵AA,并找到其最大连通分量N0(V0,E0)的邻接矩阵A0

# PAA<-omega-diag(diag(omega))
# PAA
# AA<-matrix(0,n,n)
# for (i in 1:n)
# {
#   for (j in 1:n)
#   {
#     AA[i,j] <- rbinom(1,1,PAA[i,j])
#   }
# }
# AA<-AA+W
fA <- function(omega,W)
{
  AA<-omega-diag(diag(omega))+W
  g1 <- graph_from_adjacency_matrix(AA,mode = "undirected", weighted = NULL,
                                    diag = TRUE, add.colnames = NULL, add.rownames = NA)
  V(g1)$name <- as.character(1:n)
  # plot(g1)
  
  #找最大连通分量
  if (!is.connected(g1))
  {
    cg1<-decompose.graph(g1)#cg1为连通分量列表
    cgm1<-cg1[[which.max(sapply(cg1,vcount))]]#lst[[]]提取元素；lst[]提取子列表
  }else{
    cgm1<-g1
  }
  A0<-get.adjacency(cgm1)#此时A0为dgCMatrix类下的矩阵
  A0<-as.matrix(A0)#转换为矩阵格式
  return(list(A0,g1,cgm1))
}

#取出对应连通分量的标签
fl <- function(A0,l0,g1,cgm1)
{
  n0<-dim(A0)[1]
  if (n0 != n)
  {
    vg1 <- as.vector(V(g1)$name)
    vcgm1 <- as.vector(V(cgm1)$name)
    uc <- as.integer(setdiff(vg1,vcgm1))
    l <- l0[-uc]
  }
  else{
    l <- l0
  }
  return(l)
}

ScoreError <- c()
OPCAError <- c()
for (t in 1:5)
{
  omega <- fomega(t)
  W <- fW(omega)
  lst <- fA(omega,W)
  A0 <- lst[[1]]
  g1 <- lst[[2]]
  cgm1 <- lst[[3]]
  l <- fl(A0,l0,g1,cgm1) 
  source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCORE_change.R')
  ex1_out <- SCORE(A0,K)
  ScoreErrorLst <- getError(ex1_out,l)
  ScoreError <- append(ScoreError,ScoreErrorLst$error.rate)
  source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/PCA.R')
  ex1_out <- SCORE(A0,K)
  OPCAErrorLst <- getError(ex1_out,l)
  OPCAError <- append(OPCAError,OPCAErrorLst$error.rate)
}
ScoreError
OPCAError
plot(d0,ScoreError, type='o', lwd=2, col="red",pch = 19,ylab = 'ErrorRate')
lines(d0,OPCAError,lty = 4, lwd=2, col="green")
legend(x='topleft',legend = c('SCORE','oPCA'),col = c('red','green'),lty = 1)